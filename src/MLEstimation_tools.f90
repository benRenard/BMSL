module MLEstimation_tools

!~**********************************************************************
!~* Purpose: Simple ML fitting of a dist. to some 1-d data
!~**********************************************************************
!~* Programmer: Ben Renard & Kris Kochanek, Cemagref Lyon
!~**********************************************************************
!~* Last modified: 18/08/2010
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1.
!~*		2.
!~*		3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: ML_EstimPar,GetLogLikelihood,GetUncertainty_Fisher

type::MLtype
    character(250)::distID="AintGotNoName"
    real(mrk), allocatable::InputData(:)
    integer(mik)::nData=undefIN
end type MLType
type(MLtype)::ML

Contains

subroutine GetLogLikelihood(teta,X, distID, lkh, feas, isnull,err,mess)

!^**********************************************************************
!^* Purpose: Compute Log-Likelihood
!^**********************************************************************
!^* Programmer: Ben Renard & Kris Kochanek, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 18/08/2010
!^**********************************************************************
!^* Comments: mv handled out of this sub
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.teta, par vector
!^*		2.X, data
!^*		3.distID, distribution of data
!^* OUT
!^*		1.lkh, log-likelihood
!^*		2.feas,  feasible?
!^*		3.isnull, is log-post=0?
!^*		4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetSamplePdf, GetParNumber

real(mrk), intent(in)::teta(:), X(:)
character(*), intent(in)::distID
real(mrk), intent(out), optional::lkh
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
integer(mik):: npar

err=0;mess='';feas=.true.;isnull=.false.;
if(present(lkh)) lkh=undefRN

! check size
CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
	mess='GetLogLikelihood:'//trim(mess)
	feas=.false.;return
endif
if(size(teta)/=npar) then
	feas=.false.;err=1;mess="GetLogLikelihood:FATAL:size mismatch [teta]";
	return
endif

if(present(lkh)) then
    ! get lkh
    call GetSamplePdf(DistId=DistId,x=X,par=teta,&
			    loga=.true.,pdf=lkh,feas=feas,isnull=isnull,err=err,mess=mess)
    If(err>0) then
	    mess='GetLogLikelihood:'//trim(mess)
	    feas=.false.;return
    endif
else
    return
endif

end subroutine GetLogLikelihood
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE LogLikelihoodWrapper (dataIN, dataOUT, x, feas, fx, gradFx, hessFx, err, message)

!^**********************************************************************
!^* Purpose: Wrapper to Log-Likelihood sub
!^**********************************************************************
!^* Programmer: Ben Renard & Kris Kochanek, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 18/08/2010
!^**********************************************************************
!^* Comments: mv handled out of this sub
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. dataIN, used to pass the assumed distribution and the data
!^*		2. dataOUT, unused
!^*		3. x, -> teta
!^* OUT
!^*		1. feas, feasible?
!^*		2. fx, Log-lok(teta)
!^*		3. GradFx & hessFx: unused
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type

TYPE(data_ricz_type),INTENT(in),OPTIONAL::dataIN
TYPE(data_ricz_type),INTENT(inout),OPTIONAL::dataOUT
REAL(mrk), INTENT(in) :: x(:)
LOGICAL, INTENT(out) :: feas
REAL(mrk), INTENT(out), OPTIONAL :: fx, gradFx(:), hessFx(:,:)
INTEGER(mik),INTENT(out)::err
CHARACTER(*),INTENT(out)::message
!locals
logical::isnull


CALL GetLogLikelihood(teta=x,X=ML%InputData, distID=ML%distID, lkh=fx, feas=feas, isnull=isnull,err=err,mess=message)
if(present(fx)) then
    if(isnull) fx=-HugeRe
    fx=-1*fx
endif

end SUBROUTINE LogLikelihoodWrapper
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE ML_EstimPar (X,distID,mv,OptimMethod,&
                        teta0,IterMax,&
                        ! Multistart-only optional arguments
                        nRep, HyperCube_Lo, HyperCube_Hi,&
                        MultiStartFile,&
                        OptimMethod_MultiStart,&
                        ! Optional output files
                        HistoryFile, ResultFile,&
                        teta,Hessian,err,mess)

!^**********************************************************************
!^* Purpose: ML estimation sub
!^**********************************************************************
!^* Programmer: Ben Renard & Kris Kochanek, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 18/08/2010
!^**********************************************************************
!^* Comments: mv handled here!
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.X, data
!^*		2. distID, assumed distribution
!^*		3. mv, missing value code
!^*		4. [OptimMethod], default is quasi-newton
!^*		5. [MultiStartFile], File for saving results of the MultiStart method
!^*		6. [OptimMethod_MultiStart], Optim method within MultiStart - default is quasi-newton
!^* OUT
!^*		1. teta, ML parameter
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type,data_ricz_def
use Optimisation_tools
use Distribution_tools, only: GetParNumber

real(mrk), intent(in)::X(:),mv,teta0(:)
real(mrk), intent(in), optional::HyperCube_Lo(:), HyperCube_Hi(:)
integer, intent(in)::iterMax
character(*), intent(in)::distID
integer(mik), intent(in), optional::OptimMethod,OptimMethod_MultiStart,nRep
character(*), intent(in), optional::MultiStartFile, HistoryFile, ResultFile
real(mrk), intent(out)::teta(:), Hessian(:,:)
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
integer(mik)::n0,n,nopt
logical, allocatable::mvlist(:)
real(mrk), allocatable::Xdata(:)
type(data_ricz_type)::dataIN
integer(mik)::smeth,ndim,fcalls,gcalls,hcalls,mometh, activeSet(size(teta0))
real(mrk)::msLo(size(teta0)),msHi(size(teta0)),x0(size(teta0)),fMin,&
            xLo(size(teta0)),xHi(size(teta0)),xScale(size(teta0))
character(250)::hisFile,resfile,msFile

err=0;mess='';

!optional arguments
if(present(OptimMethod)) then
    smeth=OptimMethod
else
    smeth=qndmsl_smeth
endif
if(present(OptimMethod_MultiStart)) then
    mometh=OptimMethod_MultiStart
else
    mometh=qndmsl_smeth
endif
if(present(HyperCube_Lo)) then
    msLo=HyperCube_Lo
else
    msLo=-HugeRe
endif
if(present(HyperCube_Hi)) then
    msHi=HyperCube_Hi
else
    msHi=HugeRe
endif
if(present(nRep)) then
    nopt=nRep
else
    nopt=10
endif
if(present(MultiStartFile)) then
    msFile=MultiStartFile
else
    msFile="Optim_MultiStart_Results.txt"
endif
if(present(HistoryFile)) then
    hisFile=HistoryFile
else
    hisFile="Optim_History.txt"
endif
if(present(ResultFile)) then
    resFile=ResultFile
else
    resFile="Optim_Results.txt"
endif


! handle mv
n0=size(X)
if(allocated(mvlist)) deallocate(mvlist)
allocate(mvlist(n0))
mvlist=(X/=mv)
n=count(mvlist)
if(n==0) then
	err=1;mess='ML_EstimPar:FATAL: no non-missing values!';return
endif

! populate ML
ML%nData=n
if(allocated(ML%InputData)) deallocate(ML%InputData)
allocate(ML%InputData(n))
ML%InputData=pack(X,mvlist)
ML%distID=trim(distID)

! other tunings
CALL GetParNumber(distID, nDim, err, mess)
If(err>0) then
	mess='ML_EstimPar:'//trim(mess)
	return
endif
activeSet=0
x0=teta0
xLo=-HugeRe
xHi=HugeRe
xscale=1._mrk

call optimiserWrap (smeth=smeth,nopt=nopt,evalFunc=LogLikelihoodWrapper,dataIN=dataIN,&
                          nDim=nDim,&
                          xLo=xLo,xHi=xHi,&
                          msLo=msLo,msHi=msHi,&
                          xscale=xscale,&
                          activeSet=activeSet,&
                          msFile=msFile, &
                          imeth_qn=5,gmeth=1,hmeth_qn=6,himeth_qn=1,&  ! QN options
                          mMem_lbfgs=10,&                      ! LBFGS method
                          mometh=mometh,&
                          x0=x0, &
                          xMin=teta, fMin=fMin, &
                          chkGrad=.false.,chkHess=.false.,&
                          inverseHessian=.true., hess=Hessian,&
                          iterMax=iterMax, hisFile=hisFile, &
                          iternfo=0,&
                          resFile=resFile,&
                          fcalls=fcalls,gcalls=gcalls,hcalls=hcalls,&
                          err=err, message=mess)

end SUBROUTINE ML_EstimPar
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetUncertainty_Fisher(MLestimate,Fisher,SkipUnfeas,distID,Nsim,par,err,mess)
!^**********************************************************************
!^* Purpose: Generate samples from the Gaussian-Fisher approximation
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 10/04/2012
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: 
!^**********************************************************************
!^* IN
!^*		1. MLestimate
!^*		2. Fisher, inverse Hessian from optimization
!^*		3. SkipUnfeas, logical. If TRUE, will keep trying until a feasible
!^*		   parameter is generated; if FALSE, will not test feasability 
!^*		   and may hence generate unfeasible parameters 
!^*		4. DistID, assumed distribution
!^*		5. Nsim, number of Bootstrap replicates
!^* OUT
!^*		1.par, parameter samples
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use numerix_dmsl_kit,only:normaldev
use Distribution_tools, only:GetParFeas
real(mrk), intent(in)::MLestimate(:),Fisher(:,:)
logical, intent(in)::SkipUnfeas
character(*), intent(in)::distID
integer(mik), intent(in)::Nsim
real(mrk),intent(out)::par(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik), parameter::MaxTry=100
integer(mik)::npar,i,Ntry
real(mrk)::dummy(size(MLestimate),size(MLestimate)),candid(size(MLestimate))
logical::feas

err=0;mess='';par=undefRN
! check size
npar=size(MLestimate)
if(size(Fisher,dim=1)/=npar .or. size(Fisher,dim=2)/=npar) then
    err=1;mess='GetUncertainty_Fisher:Fatal:Size mismatch[MLestimate,Fisher]';return
endif
if(size(par,dim=2)/=npar) then
    err=1;mess='GetUncertainty_Fisher:Fatal:Size mismatch[par,npar]';return
endif
if(size(par,dim=1)/=nsim) then
    err=1;mess='GetUncertainty_Fisher:Fatal:Size mismatch[par,nsim]';return
endif

do i=1,nsim
    feas=.false.;Ntry=0
    do while(feas==.false. .and. Ntry<=MaxTry)
        ! generate
        call normaldev(mean=MLestimate,covar=Fisher, do_dcmp=.true., covar_dcmp=dummy,&
                   gdev=candid, err=err, message=mess)
                   
        if(err>0) then
            mess='GetUncertainty_Fisher: '//trim(mess);return
        endif
        if(.not.SkipUnfeas) then
            feas=.true.
        else
            ! Test feasability
            call GetParFeas(DistID, candid, feas, err, mess)
            if(err>0) then
                mess='GetUncertainty_Fisher: '//trim(mess);return
            endif
        endif
        Ntry=Ntry+1
    enddo
    if(Ntry>=MaxTry) then
        mess='GetUncertainty_Fisher:WARNING:[MaxTry] exceeded - parameter likely unfeasible'
        write(*,*) trim(mess)
        err=1;return
    endif
    par(i,:)=candid
enddo
end subroutine GetUncertainty_Fisher
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end module MLEstimation_tools
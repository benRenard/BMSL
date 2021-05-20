module RegressionMCMC
!~**********************************************************************
!~* Purpose: RegressionMCMC tools
!~**********************************************************************
!~* Programmer: Xun SUN, Cemagref Lyon
!~**********************************************************************
!~* Last modified:16/12/2010
!~**********************************************************************
!~* Comments:  
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. GetLogLkh(teta,X, (CovList), distID, (linkList), lkh0, feas, isnull,err,mess)
!~*		
!~*		2. GetLogPost(teta,X, (CovList), distID, PriorList,(linkList), 
!~*		           lp, (pd), (likelihood), feas, isnull,err,mess)	
!~*	
!~*     3. RegMCMC(X, dist, PriorList, (CovList), &
!~*		          (LinkList), &
!~*		          (mv), &
!~*		          start, startStd, MinMoveRate, MaxMoveRate, DownMult, UpMult, &
!~*		          nAdapt,nCycles, OutFile, err, mess)
!~*	         
!~*     4. Bayes_Pred(File,Nburn,Nread,Nslim,& ! IN: Read properties
!~*		          Nrep, distID, & ! IN: number of reps for EACH MCMC sample + data dist
!~*		          (CovList),(newCov), Covlength, &
!~*		          modal, pred,& ! OUT: modal par estimate & replicates from the predictive
!~*		          err, mess)
!~*	  
!~*     5. Maxlkh(teta,X, (CovList), distID, (linkList), (mv), ml, (nmn), (mpar), feas, isnull,err,mess)
!~*	  
!~*     6. AICfunc(k, n, ml, out, err, mess)
!~*	  
!~*     7. BICfunc(k, n, ml, out, err, mess)
!~*	  
!~*     8. DICfunc(tetas, X, (mv), (CovList), distID, (linkList), Nburn, DIC, feas, isnull,err,mess)
!~*	  
!~*     9. DelmissingValue(X,(y), (CovList), (mv), Xn,(yn), (CovListn), (nmn), err, mess)
!~*	  
!~*     10. thetaList(teta,k, distID, (CovList), (linkList), tetaList, feas, err,mess)		     
!~**********************************************************************
use kinds_dmsl_kit ! numeric kind definitions from DMSL
use BayesianEstimation_tools, only: PriorListType, GetLogPrior 
use RegressionFunk


implicit none
Private
public :: GetLogLkh,GetLogPost,RegMCMC,Bayes_Pred,Maxlkh,BICfunc,AICfunc,DICfunc,&
                DelmissingValue,thetaList

real(mrk),allocatable:: Xdata(:)
character(100)::Xdist
type(PriorListType), allocatable:: XpriorList(:)

type, public:: CovListType
	character(100)::RegressionFunk !Model that used with the covariate
	integer(mik)::nbpar ! number of parameters used in this regression function
	real(mrk), allocatable::CovMatrix(:,:) !Covariate matrix
end type CovListType

type(CovListType), allocatable:: XCovList(:)

character(100), allocatable::LList(:)

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetLogLkh(teta,X, CovList, distID, linkList, lkh0, feas, isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute the log-likelihood with the covariates
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:31/03/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.teta, par vector
!^*		2.X, data
!^*		3.CovList, models. nb of parameters and covariates in each paramter place of distribution
!^*		   optional, no need for stationary case. but obligatory for non-stationary case
!^*		4.distID, distribution of data
!^*		5.linkList, optional, list of link function, default= 'Identity'
!^* OUT
!^*		1.lhk0, log-likelihood
!^*		2.feas,  feasible?
!^*		3.isnull, is log-post=0?
!^*		4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetPdf, GetParNumber
real(mrk), intent(in)::teta(:), x(:)
character(*), intent(in)::distID
character(*), intent(in), optional ::linkList(:)
Type(CovListType), intent(in), optional::CovList(:)
real(mrk), intent(out)::lkh0
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
character(100),allocatable::LinkL(:)
Type(CovListType), allocatable::CL(:)
integer(mik),allocatable:: nbParameter(:)
integer(mik):: npar,i,j, temp
real(mrk)::prior, lkh, fx
real(mrk), allocatable::tetaList(:,:),tetaListTemp(:,:)
logical::l

err=0;mess='';feas=.true.;isnull=.false.;lkh0=undefRN

CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
	mess='GetLogLkh:'//trim(mess)
	feas=.false.;return
endif

if(allocated(CL)) deallocate(CL)
allocate(CL(npar))

if (.not. present(CovList)) then
    if  (npar/=size(teta)) then
        err=1; mess="GetLogLkh:FATAL:missing CovList"
        return
    else

        do i=1,npar
            CL(i)%nbpar=1
            CL(i)%RegressionFunk='Identity'
            if(allocated(CL(i)%CovMatrix)) deallocate(CL(i)%CovMatrix)
            allocate(CL(i)%CovMatrix(size(X),1))
        end do
    endif
else
   !check CovList size
    if (size(CovList)/=npar) then
        err=1; mess="GetLogLkh:FATAL:Covariate Matrix size mismatch"
        return
    endif
    CL=CovList
endif

!check whether the total parameter numbers is equal to teta size
if(allocated(nbParameter)) deallocate(nbParameter)
allocate(nbParameter(npar))

do i=1,npar
    nbParameter(i)=CL(i)%nbpar
end do

if(sum(nbparameter)/=size(teta)) then
    feas=.false.;err=1;mess="GetLogLkh:FATAL:nbpar size mismatch [teta]";
    return
endif

if (allocated(LinkL)) deallocate(LinkL)
allocate(LinkL(npar))

if (.not. present(linkList)) then
    LinkL='Identity'
else
    if (size(linkList)/=npar) then
        err=1;mess="GetLogLkh:FATAL:linkList Size mismatch"
        return
    endif
    Linkl=linkList
endif

! fill in the tetaList and sum up the log-pdf of x|tetas
if (allocated(tetaList)) deallocate(tetaList)
allocate(tetaList(size(x),npar))
if (allocated(tetaListTemp)) deallocate(tetaListTemp)
allocate(tetaListTemp(size(x),npar))

lkh0=0
do i= 1, size(x)
    temp=1
    do j= 1, npar
        call ApplyRegressionFunk(RegressionFunk=CL(j)%RegressionFunk,covariate=CL(j)%CovMatrix(i,:),&
                            regressionPar=teta(temp:(temp+CL(j)%nbpar-1)),out=tetaList(i,j),err=err,mess=mess)
        if (err>0) then
        	mess='GetLogLkh:'//trim(mess)
	        feas=.false.;return
        endif
        temp=temp+CL(j)%nbpar
        
        call InverseLinkFunk(IvsLinkFunk=trim(linkl(j)),x=tetaList(i,j),fx=tetaListTemp(i,j),feas=feas,err=err,mess=mess)
        tetaList(i,j)=tetaListTemp(i,j)
        if (err>0) then
        	mess='GetLogLkh:'//trim(mess)
	        feas=.false.;return
        endif
        if(.not. feas) then
	        !err=1;mess="GetLogPost:FATAL:unfeasible "//trim(mess)
	    return
	    endif        
    end do
    
    call Getpdf(DistId=DistID,x=x(i),par=tetaList(i,:),loga=.true.,pdf=lkh,feas=feas,isnull=isnull,err=err,mess=mess)
    If(err>0) then
	    mess='GetLogLkh:'//trim(mess)
	    feas=.false.;return
    endif
    if(.not. feas) then
	    !err=1;mess="GetLogPost:FATAL:unfeasible [PriorList(i)%par]"
	    return
    endif
    
    if(isnull) return
    lkh0=lkh0+lkh
end do    
end subroutine GetLogLkh


!!!!!!!!!!!!!!!!!!!!
subroutine GetLogPost(teta,X, CovList, distID, PriorList,linkList, lp, pd, likelihood, feas, isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute unnormalized posterior with covariate
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:31/03/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.teta, par vector
!^*		2.X, data
!^*		3.CovList, models. nb of parameters and covariates in each paramter place of distribution
!^*		   optional, no need for stationary case. but obligatory for non-stationary case
!^*		4.distID, distribution of data
!^*		5.PriorList, list of priors & their par (see PriorListType above)
!^*		6.linkList, optional, list of link function, default= 'Identity'
!^* OUT
!^*		1.lp, log-post
!^*		2.pd, optional, prior value
!^*		3.likelihood, optional, likelihood value
!^*		4.feas,  feasible?
!^*		5.isnull, is log-post=0?
!^*		6.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		7.mess, error message
!^**********************************************************************

use Distribution_tools, only: GetPdf, GetParNumber
real(mrk), intent(in)::teta(:), x(:)
character(*), intent(in)::distID
character(*), intent(in), optional ::linkList(:)
Type(PriorListType), intent(in)::PriorList(:)
Type(CovListType), intent(in),optional::CovList(:)
real(mrk), intent(out)::lp
real(mrk),intent(out),optional::pd,likelihood
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
character(100),allocatable::LinkL(:)
Type(CovListType), allocatable::CL(:)
integer(mik):: npar
real(mrk)::prior, lkh0


err=0;mess='';feas=.true.;isnull=.false.;lp=undefRN

CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
	mess='GetLogPost:'//trim(mess)
	feas=.false.;return
endif

! allocate(CL)
if(allocated(CL)) deallocate(CL)
allocate(CL(npar))

!LinkList
if (allocated(LinkL)) deallocate(LinkL)
allocate(LinkL(npar))

if (.not. present(linkList)) then
    LinkL='Identity'
else
    if (size(linkList)/=npar) then
        err=1;mess="GetLogLkh:FATAL:linkList Size mismatch"
        return
    endif
    Linkl=linkList
endif

! get log-lkh
!if (present(CovList)) then
    call GetLogLkh(teta=teta, X=X, CovList=CovList, distID=distID, linkList=Linkl,&
                         lkh0=lkh0, feas=feas, isnull=isnull,err=err,mess=mess)
    If(err>0) then
	    mess='GetLogPost:'//trim(mess)
	    feas=.false.;return
    endif
    if(.not. feas) then
	    !err=1;mess="GetLogPost:FATAL:unfeasible [PriorList(i)%par]"
	    return
    endif
    if(isnull) return
!else
!    call GetLogLkh(teta=teta, X=X, distID=distID, linkList=Linkl,&
!                         lkh0=lkh0, feas=feas, isnull=isnull,err=err,mess=mess)
!    If(err>0) then
!	    mess='GetLogPost:'//trim(mess)
!	    feas=.false.;return
!    endif
!    if(.not. feas) then
!	    !err=1;mess="GetLogPost:FATAL:unfeasible [PriorList(i)%par]"
!	    return
!    endif
!    if(isnull) return
!endif
if (present(likelihood)) likelihood=lkh0  

! get log-prior
call GetLogPrior(teta,PriorList,prior, feas, isnull,err,mess)
If(err>0) then
	mess='GetLogPost:'//trim(mess)
	feas=.false.;return
endif
if(.not. feas) then
	!err=1;mess="GetLogPost:FATAL:unfeasible [PriorList(i)%par]"
	return
endif
if(isnull) return
if (present(pd)) pd=prior

lp=lkh0 +prior    
 
end subroutine GetLogPost	   


!!!!!!!!!!!!!!!!!!
subroutine RegMCMC(X, dist, PriorList, CovList, &
                 LinkList, &
                 mv, &
                 start, startStd, MinMoveRate, MaxMoveRate, DownMult, UpMult, &
                  nAdapt,nCycles, OutFile, err, mess)
                  
!^**********************************************************************
!^* Purpose: Parameter estimation using Bayes-MCMC with Covariates
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:06/01/2011
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.X, data (1-d vector)
!^*		2.distID, distribution (see distribution_tools module)
!^*		3.PriorList, list of priors & their par (see PriorListType above)
!^*		4.CovMat,optional, Covariate Matrix, with whose length= length(X)
!^*		5.CovList, models. nb of parameters and covariates in each paramter place of distribution
!^*		 for example: for a GEV distribution, there are three parameter space (location, scale, shape)
!^*		 thus the longeur of CovList=3.
!^*		6.LinkList: inverse link function 
!^*		6.mv: missing values symbol
!^*		7.nAdapt, number of iterations before adapting the jump std
!^*		8.nCycles, number of adaptations (total number of iteration is hence nAdapt*nCycles
!^*			If nCycles<1, user will be asked stop/continue
!^*		9.MinMoveRate, objective move rate (lower bound)
!^*		10.MaxMoveRate, objective move rate (upper bound)
!^*		11.DownMult, value multiplying the jump std whan move rate is too low
!^*		12.UpMult, value multiplying the jump std whan move rate is too high
!^*		13.OutFile, address of output file for MCMC samples
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.start, starting point
!^*		2.startStd, starting std for each marginal jump distribution
!^*		3.
!^*		4.
!^**********************************************************************

use MCMCStrategy_tools, only:Adaptive_Metro_OAAT
use Distribution_tools, only:  GetParNumber

real(mrk), intent(in)::X(:),MinMoveRate,MaxMoveRate,DownMult,UpMult
real(mrk), intent(in),optional::mv
Type(CovListType),intent(in),optional::CovList(:)
real(mrk), intent(inout)::start(:), startStd(:)
Character(*), intent(in)::dist
Character(*), intent(in),optional::LinkList(:)
Type(PriorListType), intent(in)::PriorList(:)
integer(mik), intent(in)::nAdapt,nCycles				
character(*), intent(in)::OutFile
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
logical,allocatable::mvlist(:)
real(mrk)::mv0
integer(mik)::n,n0,npar,nnpar,i,j
real(mrk)::fAux(3)
logical::feas, isnull
real(mrk)::fx


err=0;mess=''
feas=.true.;isnull=.false.

if (.not. present(mv)) then
    mv0=-9999
else
    mv0=mv
endif
n0=size(X)
if(allocated(mvlist)) deallocate(mvlist)
allocate(mvlist(n0))
mvlist=(X/=mv0)
n=count(mvlist)
if(n==0) then
	err=1;mess='RegMCMC:FATAL: no non-missing values!';return
endif

! handle mv

! populate global variables 
if(allocated(Xdata)) deallocate(Xdata)
allocate(Xdata(n))
Xdata=pack(X,mvlist)

Xdist=dist
npar=size(PriorList)
if(allocated(XpriorList)) deallocate(XpriorList)
allocate(XpriorList(npar))
XpriorList=PriorList

!Linklist
CALL GetParNumber(dist, nnpar, err, mess)
If(err>0) then
	mess='GetLogPost:'//trim(mess)
	feas=.false.;return
endif

if(allocated(LList)) deallocate(LList)
allocate(LList(nnpar))

if(.not. present(Linklist)) then
    LList='Identity'
else
    LList=Linklist
endif



if (.not. present(CovList)) then
    ! MCMC
    call Adaptive_Metro_OAAT(f=Post_wrapperII,x=start,&
				fx=fx,fAux=fAux,std=startStd,&
				nAdapt=nAdapt,nCycles=nCycles,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				OutFile=OutFile,err=err,mess=mess)
    If(err>0) then
	    mess='RegMCMC:'//trim(mess)
	    feas=.false.;return
    endif
    
else

!Initialize Covariate Matrix with missing values
if(allocated(XCovList)) deallocate(XCovList)
allocate(XCovList(size(CovList)))
do i=1,size(CovList)
    XCovList(i)%RegressionFunk=CovList(i)%RegressionFunk
    XCovList(i)%nbpar=CovList(i)%nbpar
    if(allocated(XCovList(i)%CovMatrix)) deallocate(XCovList(i)%CovMatrix)
    allocate( XCovList(i)%CovMatrix(n,size(CovList(i)%CovMatrix,dim=2)))
    do j=1,size(CovList(i)%CovMatrix,dim=2)
        XCovList(i)%CovMatrix(:,j)=pack(CovList(i)%CovMatrix(:,j),mvlist)
    enddo
enddo


! MCMC
call Adaptive_Metro_OAAT(f=Post_wrapper,x=start,&
				fx=fx,fAux=fAux,std=startStd,&
				nAdapt=nAdapt,nCycles=nCycles,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				OutFile=OutFile,err=err,mess=mess)
    If(err>0) then
	    mess='RegMCMC:'//trim(mess)
	    feas=.false.;return
    endif

endif

end subroutine RegMCMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Bayes_Pred(File,Nburn,Nread,Nslim,& ! IN: Read properties
								Nrep, distID, & ! IN: number of reps for EACH MCMC sample + data dist
								linkList, &
								CovList,newCov, Covlength, &
								modal, pred,& ! OUT: modal par estimate & replicates from the predictive
								err, mess)
!^**********************************************************************
!^* Purpose:  Generate deviates from the predictive distribution with covariates
!^*             and find the parameters associated to the maximum likelihood (modal)
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified: 08/03/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. File, containing MCMC samples
!^*		2. Nburn, number of lines to discard
!^*		3. Nread, number of lines to read
!^*		4. Nslim, only one line every Nslim will be used
!^*		5. Nrep, for each line, Nrep deviates will be generated
!^*(advise: keep Nrep ~ Nslim else Nrep>>Nslim will yield a huge number of replicates!)
!^*		6. distID, distribution name
!^*		7. CovList, the CovList used while doing estimation
!^*		8. NewCov, new Covariates matrix, each line contains the new covariates for the each parameter place
!^*		9. CovLength, a vector that contains the number of COVARIATES used for each parameter model
!^* OUT
!^*		1.modal, modal parameter estimate
!^*		2.pred, deviates from the predictive
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetParNumber, Generate

integer(mik), intent(in)::Nburn, Nread, Nslim, Nrep,covLength(:)
character(*),intent(in)::File, distID
character(*), intent(in), optional ::linkList(:)
Type(CovListType),intent(in)::CovList(:)
real(mrk), intent(in)::newCov(:,:)
real(mrk), intent(out)::modal(:)
real(mrk), pointer::pred(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,j, k,l, errcode, npar,n,npred, ml(1),totpar
real(mrk), allocatable::Temp(:,:), mcmc(:,:),thetaList(:),thetaListTemp(:)
logical:: feas

err=0;mess=''



open(unit=1,file=File,status='old')
read(1,*) !headers
do i=1,Nburn ! burn
	read(1,*,iostat=errcode)
	if(errcode/=0) then ! problem reading file
		err=1
		mess='Bayes_Predictive:FATAL:problem reading file during burnin - cannot proceed!'
		return
	endif
enddo

! Get number of col
CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
	mess='Bayes_Pred:'//trim(mess)
	return
endif
! check CovList size
if (size(CovList)/=npar) then
    err=1; mess="Bayes_Pred:FATAL:Covariate Matrix size mismatch"
    return
endif
!compute the total parameter number
totpar=0
do i=1, npar
    totpar=totpar+CovList(i)%nbpar
end do
!allocate Temp
if(allocated(Temp)) deallocate(Temp)
allocate(Temp(Nread,totpar+4))

do i=1,Nread ! read used lines
	k=i
	read(1,*,iostat=errcode) Temp(i,:)
	!write(*,*) Temp(i,2)
	if(errcode/=0) then ! problem reading file
		err=-1;k=k-1
		mess='Bayes_Pred:WARNING:problem reading file before [Nread] lines!'
		EXIT
	endif
enddo
close(1)

! Slim
n=k/Nslim
if(allocated(mcmc)) deallocate(mcmc)
allocate(mcmc(n,totpar+4))
mcmc=Temp(1:n:nSlim,:)

! get modal estimate
ml=Maxloc(mcmc(:,totpar+1))
if(size(modal)/=totpar) then
	err=1; mess='Bayes_Pred:SEVERE:size mismath [modal]'
	modal=undefRN
else
	modal=mcmc(ml(1),1:totpar)
endif

! check the length of covLength, that should be equal to the row number of newcov
if(size(newCov(:,1))/=size(covlength) .or. size(covlength)/=size(CovList)) then
    err=1; mess='Bayes_Pred:size mismath [newCov, covLength or CovList]'
    return;
endif

! LinkList
if(allocated(LList)) deallocate(LList)
allocate(LList(npar))

if(.not. present(Linklist)) then
    LList='Identity'
else
    LList=Linklist
endif

! get deviates from predictive
npred=n*nrep
if(associated(pred)) deallocate(pred)
allocate(pred(npred))
k=0
do i=1,n
    ! allocate thetaList and fill in it.
    if(allocated(thetaList)) deallocate(thetaList)
    allocate(thetaList(npar))
    if(allocated(thetaListTemp)) deallocate(thetaListTemp)
    allocate(thetaListTemp(npar))
    j=1
    do l=1,npar
        call ApplyRegressionFunk(RegressionFunk=CovList(l)%RegressionFunk,covariate=newcov(l,1:covLength(l)),&
                             regressionPar=mcmc(i,j:(covList(l)%nbpar+j-1)),out=thetaList(l),err=err,mess=mess)
        If(err>0) then
	        mess='Bayes_Pred:'//trim(mess)
	        return
        endif

        j=j+CovList(l)%nbpar
    
        call InverseLinkFunk(IvsLinkFunk=trim(llist(l)),x=thetaList(l),fx=thetaListTemp(l),feas=feas,err=err,mess=mess)
        thetaList(l)=thetaListTemp(l)
        if (err>0) then
        	mess='Bayes_Pred:'//trim(mess)
	        feas=.false.;return
        endif
        if(.not. feas) then
	        !err=1;mess="Bayes_Pred:FATAL:unfeasible "//trim(mess)
	    return
	    endif
    ! complete the check code. size of newcov is not good.
    end do

	do j=1,nrep
		k=k+1
		
		call Generate(DistId=distID,par=thetaList,&
					gen=pred(k),feas=feas,err=err,mess=mess)
		If(err>0) then
			mess='Bayes_Predictive: '//trim(mess);return
		endif
		if(.not. feas) then ! shouldn't happen - MCMC should NOT generate infeasible par!
			err=1;
			mess='Bayes_Predictive: unfeasible parameter!';return
		endif
	enddo
enddo


end subroutine Bayes_Pred

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Maxlkh(teta,X, CovList, distID, linkList, mv, ml, nmn, mpar, feas, isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute the maximum likelihood estimate
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:31/03/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^*		1.teta, par matrix, each row is a collection of estimated parameters
!^*		2.X, data
!^*		3.CovList, models, nb of parameters and covariates in each paramter place of distribution
!^*		   optional, no need for stationary case. but obligatory for non-stationary case
!^*		4.distID, distribution of X
!^*		5.linkList, optional, list of link function, default= 'Identity'
!^*		6.mv, optional, missing value symbol of X, default=-9999
!^* OUT
!^*		1.ml, maximum likelihood value
!^*		2.nmn, number of non-missing value in X
!^*		3.mpar, maximum likelihood estimations
!^*		3.feas,  feasible?
!^*		4.isnull, is log-post=0?
!^*		5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		6.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetPdf, GetParNumber
real(mrk), intent(in):: teta(:,:), x(:)
real(mrk), intent(in),optional::mv
character(*), intent(in)::distID
character(*), intent(in), optional ::linkList(:)
Type(CovListType), intent(in), optional::CovList(:)
real(mrk), intent(out)::ml
integer(mik),intent(out),optional::nmn
real(mrk), intent(out),optional::mpar(:)
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
logical,allocatable::mvlist(:)
real(mrk)::mv0,temp
integer(mik)::i,j,n,n0,nnpar

feas=.true.;isnull=.false.
err=0;mess=''
ml=undefRN

if (.not. present(mv)) then
    mv0=-9999._mrk
else
    mv0=mv
endif
n0=size(X)
if(allocated(mvlist)) deallocate(mvlist)
allocate(mvlist(n0))
mvlist=(X/=mv0)
n=count(mvlist)
if(n==0) then
	err=1;mess='Maxlkh:FATAL: no non-missing values!';return
endif

if (present(nmn)) nmn=n
! handle mv

! populate global variables 
if(allocated(Xdata)) deallocate(Xdata)
allocate(Xdata(n))
Xdata=pack(X,mvlist)

!Linklist
CALL GetParNumber(distID, nnpar, err, mess)
If(err>0) then
	mess='GetLogPost:'//trim(mess)
	feas=.false.;return
endif

if(allocated(LList)) deallocate(LList)
allocate(LList(nnpar))

if(.not. present(Linklist)) then
    LList='Identity'
else
    LList=Linklist
endif

if (.not. present(CovList)) then
    if (present(mpar) .and. size(mpar)/=nnpar) then
        err=1; mess='Maxlkh: output parameter size mismatch'
        return
    endif
    
    do i=1, size(teta,dim=1)
        call GetLogLkh(teta=teta(i,:),X=Xdata, distID=distID, lkh0=temp, feas=feas,isnull=isnull,err=err,mess=mess)
        If(err>0) then
	        mess='ML:'//trim(mess)
	        feas=.false.;return
        endif
        if (temp>ml) then
            ml=temp
            if (present(mpar)) mpar=teta(i,:)
        endif
    enddo 
                     
    
else
    if (present(mpar) .and. size(mpar)/=size(teta(1,:))) then
        err=1; mess='Maxlkh: output parameter size mismatch'
        return
    endif
    !Initialize Covariate Matrix with missing values
    if(allocated(XCovList)) deallocate(XCovList)
    allocate(XCovList(size(CovList)))
    do i=1,size(CovList)
       XCovList(i)%RegressionFunk=CovList(i)%RegressionFunk
       XCovList(i)%nbpar=CovList(i)%nbpar
       if(allocated(XCovList(i)%CovMatrix)) deallocate(XCovList(i)%CovMatrix)
       allocate( XCovList(i)%CovMatrix(n,size(CovList(i)%CovMatrix,dim=2)))
       do j=1,size(CovList(i)%CovMatrix,dim=2)
          XCovList(i)%CovMatrix(:,j)=pack(CovList(i)%CovMatrix(:,j),mvlist)
       enddo
    enddo


    do i=1,size(teta,dim=1)
        call GetLogLkh(teta=teta(i,:),X=Xdata, CovList=XCovList, distID=distID,linkList=LList,&
                                 lkh0=temp, feas=feas,isnull=isnull,err=err,mess=mess)
        If(err>0) then
	        mess='ML:'//trim(mess)
	        feas=.false.;return
        endif
        if (temp>ml) then
            ml=temp
            if (present(mpar)) mpar=teta(i,:)
        endif
    enddo

endif


end subroutine Maxlkh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BICfunc(k, n, ml, BIC, err, mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:  /  /2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.k, number of parameters
!^*     2.n, observation numbers
!^*     3.ml, log-maximum likelihood
!^* OUT
!^*     1.out,bic
!^*     2.err
!^*     3.mess
!^**********************************************************************
integer(mik),intent(in)::k,n
real(mrk),intent(in)::ml
real(mrk),intent(out)::BIC
integer(mik),intent(out)::err
character(*),intent(out)::mess

err=0;mess=''

if (n<=0) then
    err=1;mess='Compute BIC error, number of observatons negativel'
    return;
endif

if (k<=0) then
    err=1;mess='Compute BIC error, number of parameters negativel'
    return;
endif

BIC=k*LOG(n*1._mrk)-2*ml
end subroutine BICfunc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AICfunc(k, n, ml, AICc, err, mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:  /  /2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!!^* IN
!^*     1.k, number of parameters
!^*     2.n, observation numbers
!^*     3.ml, log-maximum likelihood
!^* OUT
!^*     1.out,bic
!^*     2.err
!^*     3.mess
!^**********************************************************************
integer(mik),intent(in)::k,n
real(mrk),intent(in)::ml
real(mrk),intent(out)::AICc

integer(mik),intent(out)::err
character(*),intent(out)::mess

err=0;mess=''

if (n<=0) then
    err=1;mess='Compute BIC error, number of observatons negativel'
    return;
endif

if (k<=0) then
    err=1;mess='Compute BIC error, number of parameters negativel'
    return;
endif

AICc=2*k-2*ml+2*k*(k+1)/(n-k-1)

end subroutine AICfunc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DICfunc(tetas, X, mv, CovList, distID, linkList,Nburn, DIC, feas, isnull,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:  /  /2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* OUT
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* INOUT
!^*     1.
!^*     2.
!^**********************************************************************
use Distribution_tools, only: GetPdf, GetParNumber
real(mrk), intent(in)::tetas(:,:), x(:)
real(mrk), intent(in),optional::mv
character(*), intent(in)::distID
character(*), intent(in), optional ::linkList(:)
Type(CovListType), intent(in), optional::CovList(:)
integer(mik),intent(in)::Nburn
real(mrk), pointer::DIC(:)
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
character(100),allocatable::LinkL(:)
Type(CovListType), allocatable::CL(:)
integer(mik),allocatable:: nbParameter(:)
integer(mik):: npar,i,j,tl
real(mrk)::missingvalue, Dtotal, D, lkh, lkhmean, pD
real(mrk),allocatable::tetatotal(:)
real(mrk),pointer::newX(:)
Type(CovListType), pointer::newCL(:)


err=0;mess='';feas=.true.;isnull=.false.

CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
	mess='DICfunc:'//trim(mess)
	feas=.false.;return
endif

if(allocated(CL)) deallocate(CL)
allocate(CL(npar))

if (.not. present(CovList)) then
    if  (npar/=size(tetas,dim=2)) then
        err=1; mess="DICfunc:FATAL:missing CovList"
        return
    else

        do i=1,npar
            CL(i)%nbpar=1
            CL(i)%RegressionFunk='Identity'
            if(allocated(CL(i)%CovMatrix)) deallocate(CL(i)%CovMatrix)
            allocate(CL(i)%CovMatrix(size(X),1))
        end do
    endif
else
   !check CovList size
    if (size(CovList)/=npar) then
        err=1; mess="DICfunc:FATAL:Covariate Matrix size mismatch"
        return
    endif
    CL=CovList
endif

!check whether the total parameter numbers is equal to teta size
if(allocated(nbParameter)) deallocate(nbParameter)
allocate(nbParameter(npar))

do i=1,npar
    nbParameter(i)=CL(i)%nbpar
end do

if(sum(nbparameter)/=size(tetas,dim=2)) then
    feas=.false.;err=1;mess="DICfunc:FATAL:nbpar size mismatch [teta]";
    return
endif

if (allocated(LinkL)) deallocate(LinkL)
allocate(LinkL(npar))

if (.not. present(linkList)) then
    LinkL='Identity'
else
    if (size(linkList)/=npar) then
        err=1;mess="DICfunc:FATAL:linkList Size mismatch"
        return
    endif
    Linkl=linkList
endif

if(present(mv)) then
    missingvalue=mv
else
    missingvalue=-9999._mrk
endif

!start compute DIC
tl=size(tetas,dim=1)
call DelmissingValue(X=x, CovList=CL, mv=missingvalue, Xn=newX, CovListn=newCL, err=err, mess=mess)
If(err>0) then
	mess='DICfunc:'//trim(mess)
	feas=.false.;return
endif

if(associated(DIC)) deallocate(DIC)
allocate(DIC(tl))
if(allocated(tetatotal)) deallocate(tetatotal)
allocate(tetatotal(size(tetas,dim=2)))

DIC(1:Nburn)=-9999._mrk
tetatotal=0._mrk
Dtotal=0._mrk

do i=(Nburn+1),tl
    call GetLogLkh(teta=tetas(i,:),X=newX, CovList=newCL, distID=distID, linkList=linkList, lkh0=lkh,& 
                    feas=feas,isnull=isnull, err=err, mess=mess)
    If(err>0) then
	    mess='DICfunc:'//trim(mess)
	    feas=.false.;return
    endif
    if(.not. feas) then
	    err=1;mess="DICfunc:FATAL:unfeasible [PriorList(i)%par]"
	    return
    endif
    if(isnull) return
    
    D=-2._mrk*lkh
    Dtotal=Dtotal+D
    tetatotal=tetatotal+tetas(i,:)
    
    call GetLogLkh(teta=(tetatotal/(i-Nburn)),X=newX, CovList=newCL, distID=distID, linkList=linkList, lkh0=lkhmean,& 
                    feas=feas,isnull=isnull, err=err, mess=mess)
    If(err>0) then
	    mess='DICfunc:'//trim(mess)
	    feas=.false.;return
    endif
    if(.not. feas) then
	    err=1;mess="DICfunc:FATAL:unfeasible [PriorList(i)%par]"
	    return
    endif
    if(isnull) return
    
    pd=Dtotal/(i-Nburn)+2._mrk*lkhmean
    DIC(i)=Dtotal/(i-Nburn)+pd
enddo

end subroutine DICfunc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DelmissingValue(X,y, CovList, mv, Xn,yn, CovListn, nmn, err, mess)
!^**********************************************************************
!^* Purpose: Delete the missing values in X
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:01/04/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.X, observation
!^*     2.CovList, optional, Covlist associated to X
!^*     3.mv, optional, missing values symbols
!^* OUT
!^*     1.Xn, result of X
!^*     2.CovListn, optional, result of CovList
!^*     3.nmn, optional, number of non-missing values
!^*     4.err
!^*     5.mess
!^**********************************************************************
real(mrk), intent(in):: x(:)
real(mrk), intent(in),optional::mv, y(:,:)
Type(CovListType), intent(in), optional::CovList(:)
real(mrk), pointer:: Xn(:)
real(mrk), pointer, optional::yn(:,:)
Type(CovListType), pointer, optional::CovListn(:)
integer(mik),intent(out),optional::nmn
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
logical,allocatable::mvlist(:)
real(mrk)::mv0
integer(mik)::i,j,n,n0


err=0;mess=''

if (.not. present(mv)) then
    mv0=-9999._mrk
else
    mv0=mv
endif
n0=size(X)
if(allocated(mvlist)) deallocate(mvlist)
allocate(mvlist(n0))
mvlist=(X/=mv0)
n=count(mvlist)
if(n==0) then
	err=1;mess='DelmissingValue:FATAL: no non-missing values!';return
endif

if (present(nmn)) nmn=n

if(associated(Xn)) deallocate(Xn)
allocate(Xn(n))
Xn=pack(X,mvlist)

!yn
if (present(yn)) then
    if(.not. present(y)) then 
        err=1;mess='DelmissingValue: missing y';return
    else
        if (size(y,dim=1)/=size(x)) then
            err=1;mess='DelmissingValue: y size mismatch with x';return
        else
            if(associated(yn)) deallocate(yn)
            allocate(yn(n,size(y,dim=2)))
            do i=1,size(y,dim=2)
                yn(:,i)=pack(y(:,i),mvlist)
            enddo
        endif
    endif
endif
       
        
!CovList
if(present(CovListn)) then
    if(.not. present(CovList)) then
        err=1;mess='DelmissingValue:FATAL: no input CovList!';return
    else
        !if (size(CovList)/=size(X)) then
            !err=1; mess='size CovList mismatch with X'
            !return
        !else
            if(associated(CovListn)) deallocate(CovListn)
            allocate(CovListn(size(CovList)))
            do i=1,size(CovList)
                CovListn(i)%RegressionFunk=CovList(i)%RegressionFunk
                CovListn(i)%nbpar=CovList(i)%nbpar
                if(allocated(CovListn(i)%CovMatrix)) deallocate(CovListn(i)%CovMatrix)
                allocate( CovListn(i)%CovMatrix(n,size(CovList(i)%CovMatrix,dim=2)))
                do j=1,size(CovList(i)%CovMatrix,dim=2)
                    CovListn(i)%CovMatrix(:,j)=pack(CovList(i)%CovMatrix(:,j),mvlist)
                enddo
            enddo
        !endif
    endif
endif
    
    
end subroutine DelmissingValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine thetaList(teta,k, distID, CovList, linkList, tetaList, feas, err,mess)
!^**********************************************************************
!^* Purpose: compute the parameters with covariates for each distribution
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:06/04/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.k, size(x)
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* OUT
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* INOUT
!^*     1.
!^*     2.
!^**********************************************************************
use Distribution_tools, only: GetParNumber
real(mrk), intent(in)::teta(:)
integer(mik),intent(in)::k
character(*), intent(in)::distID
character(*), intent(in), optional ::linkList(:)
Type(CovListType), intent(in), optional::CovList(:)
real(mrk), intent(out)::tetaList(:,:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(100),intent(out)::mess
!local
Type(CovListType), allocatable::CL(:)
integer(mik),allocatable:: nbParameter(:)
real(mrk), allocatable::tetaListTemp(:,:)
integer(mik)::i,j, temp,npar


err=0;mess='';feas=.true.

CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
	mess='thetaList:'//trim(mess)
	feas=.false.;return
endif

if(allocated(CL)) deallocate(CL)
allocate(CL(npar))

if (.not. present(CovList)) then
    if  (npar/=size(teta)) then
        err=1; mess="ThetaList:FATAL:missing CovList"
        return
    else

        do i=1,npar
            CL(i)%nbpar=1
            CL(i)%RegressionFunk='Identity'
            if(allocated(CL(i)%CovMatrix)) deallocate(CL(i)%CovMatrix)
            allocate(CL(i)%CovMatrix(k,1))
        end do
    endif
else
   !check CovList size
    if (size(CovList)/=npar) then
        err=1; mess="ThetaList:FATAL:Covariate Matrix size mismatch"
        return
    endif
    CL=CovList
endif

!check whether the total parameter numbers is equal to teta size
if(allocated(nbParameter)) deallocate(nbParameter)
allocate(nbParameter(npar))

do i=1,npar
    nbParameter(i)=CL(i)%nbpar
end do

if(sum(nbparameter)/=size(teta)) then
    feas=.false.;err=1;mess="ThetaList:FATAL:nbpar size mismatch [teta]";
    return
endif

!LinkList
if (allocated(LList)) deallocate(LList)
allocate(LList(npar))

if (.not. present(linkList)) then
    LList='Identity'
else
    if (size(linkList)/=npar) then
        err=1;mess="ThetaList:FATAL:linkList Size mismatch"
        return
    endif
    LList=linkList
endif

!check tetaList size
if (size(tetaList,dim=2)/=npar) then
    err=1;mess='thetaList:FATAL: tetaList column size mismatch!';return
endif
if (size(tetaList,dim=1)/=k) then
    err=1;mess='thetaList:FATAL: tetaList row size mismatch with size(x)!';return
endif

if (allocated(tetaListTemp)) deallocate(tetaListTemp)
allocate(tetaListTemp(k,npar))


do i= 1, k
    temp=1
    do j= 1, npar
        call ApplyRegressionFunk(RegressionFunk=CL(j)%RegressionFunk,covariate=CL(j)%CovMatrix(i,:),&
                            regressionPar=teta(temp:(temp+CL(j)%nbpar-1)),out=tetaList(i,j),err=err,mess=mess)
        if (err>0) then
        	mess='thetaList:'//trim(mess)
	        feas=.false.;return
        endif
        temp=temp+CL(j)%nbpar
        
        call InverseLinkFunk(IvsLinkFunk=trim(LList(j)),x=tetaList(i,j),fx=tetaListTemp(i,j),feas=feas,err=err,mess=mess)
        tetaList(i,j)=tetaListTemp(i,j)
        if (err>0) then
        	mess='thetaList:'//trim(mess)
	        feas=.false.;return
        endif
        if(.not. feas) then
	        !err=1;mess="thetaList:FATAL:unfeasible "//trim(mess)
	    return
	    endif        
    end do
end do

end subroutine thetaList




!!!!!!!!!
!Private!
!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!
subroutine Post_wrapper(x,feas,isnull,fx,fAux,err,mess)
!^**********************************************************************
!^* Purpose: wrapper to GetLogPost to comply with MCMC interface
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:06/01/2011
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.x, teta
!^* OUT
!^*		1. see GetLogPost
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
real(mrk),intent(in)::x(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess

if(present(fAux)) fAux=UndefRN
call GetLogPost(teta=x, X=XData, CovList=XCovList, distID=Xdist, PriorList=XpriorList, linkList=LList, lp=fx, likelihood=fAux(1),&
				feas=feas, isnull=isnull,err=err,mess=mess)
call AICfunc(k=size(x), n=size(XData), ml=fAux(1), AICc=fAux(2), err=err, mess=mess)
call BICfunc(k=size(x), n=size(XData), ml=fAux(1), BIC=fAux(3), err=err, mess=mess)
end subroutine Post_wrapper


!!!!!!!!!!!!!!!!!!!!!!!
subroutine Post_wrapperII(x,feas,isnull,fx,fAux,err,mess)
!^**********************************************************************
!^* Purpose: wrapper to GetLogPost to comply with MCMC interface
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:06/01/2011
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.x, teta
!^* OUT
!^*		1. see GetLogPost
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
real(mrk),intent(in)::x(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess

if(present(fAux)) fAux=UndefRN
call GetLogPost(teta=x, X=XData, distID=Xdist, PriorList=XpriorList, linkList=LList, lp=fx, likelihood=fAux(1), &
				feas=feas, isnull=isnull,err=err,mess=mess)
call AICfunc(k=size(x), n=size(XData), ml=fAux(1), AICc=fAux(2), err=err, mess=mess)
call BICfunc(k=size(x), n=size(XData), ml=fAux(1), BIC=fAux(3), err=err, mess=mess)

end subroutine Post_wrapperII


!!!!!!!!!!!!!!!!!!!!
!An ancient vision of getlogpost
subroutine GetLogPostII(teta,X, CovList, distID, PriorList,linkList, lp, pd, likelihood, feas, isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute unnormalized posterior with covariate
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:06/01/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.teta, par vector
!^*		2.X, data
!^*		3.CovList, models. nb of parameters and covariates in each paramter place of distribution
!^*		4.distID, distribution of data
!^*		5.PriorList, list of priors & their par (see PriorListType above)
!^* OUT
!^*		1.lp, log-post
!^*		2.pd, optional, prior value
!^*		3.likelihood, optional, likelihood value
!^*		4.feas,  feasible?
!^*		5.isnull, is log-post=0?
!^*		6.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		7.mess, error message
!^**********************************************************************

use Distribution_tools, only: GetPdf, GetParNumber
real(mrk), intent(in)::teta(:), x(:)
character(*), intent(in)::distID
character(*), intent(in), optional ::linkList(:)
Type(PriorListType), intent(in)::PriorList(:)
Type(CovListType), intent(in)::CovList(:)
real(mrk), intent(out)::lp
real(mrk),intent(out),optional::pd,likelihood
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
character(100),allocatable::LinkL(:)
integer(mik),allocatable:: nbParameter(:)
integer(mik):: npar,i,j, temp
real(mrk)::prior, lkh0, lkh, fx
real(mrk), allocatable::tetaList(:,:),tetaListTemp(:,:)
logical::l

err=0;mess='';feas=.true.;isnull=.false.;lp=undefRN

CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
	mess='GetLogPost:'//trim(mess)
	feas=.false.;return
endif


!check CovList size
if (size(CovList)/=npar) then
    err=1; mess="GetLogPost:FATAL:Covariate Matrix size mismatch"
    return
endif

!check whether the total parameter numbers is equal to teta size
if(allocated(nbParameter)) deallocate(nbParameter)
allocate(nbParameter(npar))

do i=1,npar
    nbParameter(i)=CovList(i)%nbpar
end do

if(sum(nbparameter)/=size(teta)) then
    feas=.false.;err=1;mess="GetLogPost:FATAL:nbpar size mismatch [teta]";
    return
endif

if (allocated(LinkL)) deallocate(LinkL)
allocate(LinkL(npar))

if (.not. present(linkList)) then
    LinkL='Identity'
else
    if (size(linkList)/=npar) then
        err=1;mess="GetLogPost:FATAL:linkList Size mismatch"
        return
    endif
    Linkl=linkList
endif

! get log-prior
call GetLogPrior(teta,PriorList,prior, feas, isnull,err,mess)
If(err>0) then
	mess='GetLogPost:'//trim(mess)
	feas=.false.;return
endif
if(.not. feas) then
	!err=1;mess="GetLogPost:FATAL:unfeasible [PriorList(i)%par]"
	return
endif
if(isnull) return

if (present(pd)) pd=prior

! fill in the tetaList and sum up the log-pdf of x|tetas
if (allocated(tetaList)) deallocate(tetaList)
allocate(tetaList(size(x),npar))
if (allocated(tetaListTemp)) deallocate(tetaListTemp)
allocate(tetaListTemp(size(x),npar))

lkh0=0
do i= 1, size(x)
    temp=1
    do j= 1, npar
        call ApplyRegressionFunk(RegressionFunk=CovList(j)%RegressionFunk,covariate=CovList(j)%CovMatrix(i,:),&
                            regressionPar=teta(temp:(temp+CovList(j)%nbpar-1)),out=tetaList(i,j),err=err,mess=mess)
        if (err>0) then
        	mess='GetLogPost:'//trim(mess)
	        feas=.false.;return
        endif
        temp=temp+CovList(j)%nbpar
        
        call InverseLinkFunk(IvsLinkFunk=trim(linkl(j)),x=tetaList(i,j),fx=tetaListTemp(i,j),feas=feas,err=err,mess=mess)
        tetaList(i,j)=tetaListTemp(i,j)
        if (err>0) then
        	mess='GetLogPost:'//trim(mess)
	        feas=.false.;return
        endif
        if(.not. feas) then
	        !err=1;mess="GetLogPost:FATAL:unfeasible "//trim(mess)
	    return
	    endif        
    end do
    
    call Getpdf(DistId=DistID,x=x(i),par=tetaList(i,:),loga=.true.,pdf=lkh,feas=feas,isnull=isnull,err=err,mess=mess)
    If(err>0) then
	    mess='GetLogPost:'//trim(mess)
	    feas=.false.;return
    endif
    if(.not. feas) then
	    !err=1;mess="GetLogPost:FATAL:unfeasible [PriorList(i)%par]"
	    return
    endif
    
    if(isnull) return
    lkh0=lkh0+lkh
end do    
if (present(likelihood)) likelihood=lkh0    
lp=lkh0 +prior    
 
end subroutine GetLogPostII	 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!An ancient vision of Bayes_Pred. Some errors in the following code
subroutine Bayes_Pred2(File,Nburn,Nread,Nslim,& ! IN: Read properties
								Nrep, distID, & ! IN: number of reps for EACH MCMC sample + data dist
								CovList,newCov, Covlength, &
								modal, pred,& ! OUT: modal par estimate & replicates from the predictive
								err, mess)
!^**********************************************************************
!^* Purpose:  Generate deviates from the predictive distribution with covariates
!^*             and find the parameters associated to the maximum likelihood (modal)
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified: 08/03/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. File, containing MCMC samples
!^*		2. Nburn, number of lines to discard
!^*		3. Nread, number of lines to read
!^*		4. Nslim, only one line every Nslim will be used
!^*		5. Nrep, for each line, Nrep deviates will be generated
!^*(advise: keep Nrep ~ Nslim else Nrep>>Nslim will yield a huge number of replicates!)
!^*		6. distID, distribution name
!^*		7. CovList, the CovList used while doing estimation
!^*		8. NewCov, new Covariates matrix, each line contains the new covariates for the each parameter place
!^*		9. CovLength, a vector that contains the number of COVARIATES used for each parameter model
!^* OUT
!^*		1.modal, modal parameter estimate
!^*		2.pred, deviates from the predictive
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetParNumber, Generate

integer(mik), intent(in)::Nburn, Nread, Nslim, Nrep,covLength(:)
character(*),intent(in)::File, distID
Type(CovListType),intent(in)::CovList(:)
real(mrk)::newCov(:,:)
real(mrk), intent(out)::modal(:)
real(mrk), pointer::pred(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,j, k, errcode, npar,n,npred, ml(1),totpar
real(mrk), allocatable::Temp(:,:), mcmc(:,:),thetaList(:)
logical:: feas

err=0;mess=''


open(unit=1,file=File,status='old')
read(1,*) !headers
do i=1,Nburn ! burn
	read(1,*,iostat=errcode)
	if(errcode/=0) then ! problem reading file
		err=1
		mess='Bayes_Predictive:FATAL:problem reading file during burnin - cannot proceed!'
		return
	endif
enddo

! Get number of col
CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
	mess='Bayes_Pred:'//trim(mess)
	return
endif
! check CovList size
if (size(CovList)/=npar) then
    err=1; mess="Bayes_Pred:FATAL:Covariate Matrix size mismatch"
    return
endif
!compute the total parameter number
totpar=0
do i=1, npar
    totpar=totpar+CovList(i)%nbpar
end do
!allocate Temp
if(allocated(Temp)) deallocate(Temp)
allocate(Temp(Nread,totpar+4))

do i=1,Nread ! read used lines
	k=i
	read(1,*,iostat=errcode) Temp(i,:)
	!write(*,*) Temp(i,2)
	if(errcode/=0) then ! problem reading file
		err=-1;k=k-1
		mess='Bayes_Pred:WARNING:problem reading file before [Nread] lines!'
		EXIT
	endif
enddo
close(1)

! Slim
n=k/Nslim
if(allocated(mcmc)) deallocate(mcmc)
allocate(mcmc(n,totpar+4))
mcmc=Temp(1:n:nSlim,:)

! get modal estimate
ml=Maxloc(mcmc(:,totpar+1))
if(size(modal)/=totpar) then
	err=1; mess='Bayes_Pred:SEVERE:size mismath [modal]'
	modal=undefRN
else
	modal=mcmc(ml(1),1:totpar)
endif

! check the length of covLength, that should be equal to the row number of newcov
if(size(newCov(:,1))/=size(covlength) .or. size(covlength)/=size(CovList)) then
    err=1; mess='Bayes_Pred:size mismath [newCov, covLength or CovList]'
    return;
endif


! allocate thetaList and fill in it.
if(allocated(thetaList)) deallocate(thetaList)
allocate(thetaList(npar))
j=1
do i=1,npar
    call ApplyRegressionFunk(RegressionFunk=CovList(i)%RegressionFunk,covariate=newcov(i,1:covLength(i)),&
                             regressionPar=mcmc(i,j:(covList(i)%nbpar+j-1)),out=thetaList(i),err=err,mess=mess)
    If(err>0) then
	    mess='Bayes_Pred:'//trim(mess)
	    return
    endif

    j=CovList(i)%nbpar+j
    
! complete the check code. size of newcov is not good.
end do

! get deviates from predictive
npred=n*nrep
if(associated(pred)) deallocate(pred)
allocate(pred(npred))
k=0
do i=1,n
	do j=1,nrep
		k=k+1
		
		call Generate(DistId=distID,par=thetaList,&
					gen=pred(k),feas=feas,err=err,mess=mess)
		If(err>0) then
			mess='Bayes_Predictive: '//trim(mess);return
		endif
		if(.not. feas) then ! shouldn't happen - MCMC should NOT generate infeasible par!
			err=1;
			mess='Bayes_Predictive: unfeasible parameter!';return
		endif
	enddo
enddo


end subroutine Bayes_Pred2






end module RegressionMCMC
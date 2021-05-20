module MomentEstimation_tools

!~**********************************************************************
!~* Purpose: Moment estimators & derivatives (L-moment) 
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 04/04/2012
!~**********************************************************************
!~* Comments: 
!~* 
!~* 
!~* 
!~* 
!~* 
!~* 
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. 
!~*		2. 
!~*		3. 
!~*		4. 
!~**********************************************************************
use kinds_dmsl_kit ! numeric kind definitions from DMSL
use EmpiricalStats_tools, only:GetEmpiricalStats,Get4LMoments
implicit none
Private
public :: GetMomentEstimate,GetMOMEstimate,GetUncertainty_Bootstrap,&
          GetUncertainty_ParametricBootstrap

contains

subroutine GetMomentEstimate(DistID,X,Mtype,par,err,mess)
!^**********************************************************************
!^* Purpose: Moment Estimation, with optional 'Mtype' for switching 
!^*          between standard MOM and L-moment
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/04/2012
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: 
!^**********************************************************************
!^* IN
!^*		1. DistID, assumes distribution
!^*		2. X, sample used to estimate
!^*		3. [Mtype], 'MOM'(default),'LMOM'
!^* OUT
!^*		1.par, parameter estimate
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
Use Distribution_tools
character(*), intent(in)::DistID
character(*), intent(in),optional::Mtype
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
character(250)::MT

err=0;mess='';par=UndefRN
! handle optional MT
if(present(Mtype)) then
    MT=Mtype
else
    MT='MOM'
endif
select case(trim(MT))
case('MOM')
    call GetMOMEstimate(DistID,X,par,err,mess)
    if(err>0) then
        mess='GetMomentEstimate:'//trim(mess);return
    endif
case('LMOM')
    call GetLMOMEstimate(DistID,X,par,err,mess)
    if(err>0) then
        mess='GetMomentEstimate:'//trim(mess);return
    endif
case default
    err=1;mess='GetMomentEstimate:Fatal:Unavailable [M-type]';return
end select
end subroutine GetMomentEstimate
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetMOMEstimate(DistID,X,par,err,mess)
!^**********************************************************************
!^* Purpose: standard Moment Estimation
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/04/2012
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: 
!^**********************************************************************
!^* IN
!^*		1. DistID, assumes distribution
!^*		2. X, sample used to estimate
!^* OUT
!^*		1.par, parameter estimate
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
Use Distribution_tools
character(*), intent(in)::DistID
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::npar
real(mrk)::par_dummy(size(par))

err=0;mess='';par=UndefRN
! check size
call GetParNumber(DistID, npar, err, mess)
if(err>0) then
    mess='GetMOMEstimate:'//trim(mess);return
endif
if(size(par)/=npar) then
    err=1;mess='GetMOMEstimate:Fatal:Size mismatch [par,DistID]';return
endif
select case(DistID)
case(Gaussian)
    call Gaussian_MOM(X,par,err,mess)
    if(err>0) then
        mess='GetMOMEstimate:'//trim(mess);return
    endif
case(Gaussian_reparameterized)
    call Gaussian_MOM(X,par_dummy,err,mess)
    if(err>0) then
        mess='GetMOMEstimate:'//trim(mess);return
    endif
    par(1)=par_dummy(1)
    if(par_dummy(1)==0._mrk) then
        err=1;mess='GetMOMEstimate:Fatal:zero-mean';return
    else
        par(2)=par_dummy(2)/par_dummy(1)
    endif
case(Gumbel)
    call Gumbel_MOM(X,par,err,mess)
    if(err>0) then
        mess='GetMOMEstimate:'//trim(mess);return
    endif
case(Gumbel_reparameterized)
    call Gumbel_MOM(X,par_dummy,err,mess)
    if(err>0) then
        mess='GetMOMEstimate:'//trim(mess);return
    endif
    par(1)=par_dummy(1)
    if(par_dummy(1)==0._mrk) then
        err=1;mess='GetMOMEstimate:Fatal:zero-location';return
    else
        par(2)=par_dummy(2)/par_dummy(1)
    endif
case(PearsonIII)
    call PearsonIII_MOM(X,par,err,mess)
    if(err>0) then
        mess='GetMOMEstimate:'//trim(mess);return
    endif
case(PearsonIII_reparameterized)
    call PearsonIII_MOM(X,par_dummy,err,mess)
    if(err>0) then
        mess='GetMOMEstimate:'//trim(mess);return
    endif
    par(1)=par_dummy(1)
    if(par_dummy(1)==0._mrk) then
        err=1;mess='GetMOMEstimate:Fatal:zero-location';return
    else
        par(2)=par_dummy(2)/par_dummy(1)
    endif
    par(3)=par_dummy(3)
case(GEV)
    call GEV_MOM(X,par,err,mess)
    if(err>0) then
        mess='GetMOMEstimate:'//trim(mess);return
    endif
case(GEV_reparameterized)
    call GEV_MOM(X,par_dummy,err,mess)
    if(err>0) then
        mess='GetMOMEstimate:'//trim(mess);return
    endif
    par(1)=par_dummy(1)
    if(par_dummy(1)==0._mrk) then
        err=1;mess='GetMOMEstimate:Fatal:zero-location';return
    else
        par(2)=par_dummy(2)/par_dummy(1)
    endif
    par(3)=par_dummy(3)
case default
    err=1;mess='GetMOMEstimate:Fatal:Unavailable [DistID]';return
end select
end subroutine GetMOMEstimate
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetLMOMEstimate(DistID,X,par,err,mess)
!^**********************************************************************
!^* Purpose: L-Moment Estimation
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/04/2012
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: 
!^**********************************************************************
!^* IN
!^*		1. DistID, assumes distribution
!^*		2. X, sample used to estimate
!^* OUT
!^*		1.par, parameter estimate
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
Use Distribution_tools
character(*), intent(in)::DistID
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::npar
real(mrk)::par_dummy(size(par))

err=0;mess='';par=UndefRN
! check size
call GetParNumber(DistID, npar, err, mess)
if(err>0) then
    mess='GetLMOMEstimate:'//trim(mess);return
endif
if(size(par)/=npar) then
    err=1;mess='GetLMOMEstimate:Fatal:Size mismatch [par,DistID]';return
endif
select case(DistID)
case(Gaussian)
    call Gaussian_LMOM(X,par,err,mess)
    if(err>0) then
        mess='GetLMOMEstimate:'//trim(mess);return
    endif
case(Gaussian_reparameterized)
    call Gaussian_LMOM(X,par_dummy,err,mess)
    if(err>0) then
        mess='GetLMOMEstimate:'//trim(mess);return
    endif
    par(1)=par_dummy(1)
    if(par_dummy(1)==0._mrk) then
        err=1;mess='GetLMOMEstimate:Fatal:zero-mean';return
    else
        par(2)=par_dummy(2)/par_dummy(1)
    endif
case(Gumbel)
    call Gumbel_LMOM(X,par,err,mess)
    if(err>0) then
        mess='GetLMOMEstimate:'//trim(mess);return
    endif
case(Gumbel_reparameterized)
    call Gumbel_LMOM(X,par_dummy,err,mess)
    if(err>0) then
        mess='GetLMOMEstimate:'//trim(mess);return
    endif
    par(1)=par_dummy(1)
    if(par_dummy(1)==0._mrk) then
        err=1;mess='GetLMOMEstimate:Fatal:zero-location';return
    else
        par(2)=par_dummy(2)/par_dummy(1)
    endif
case(PearsonIII)
    call PearsonIII_LMOM(X,par,err,mess)
    if(err>0) then
        mess='GetLMOMEstimate:'//trim(mess);return
    endif
case(PearsonIII_reparameterized)
    call PearsonIII_LMOM(X,par_dummy,err,mess)
    if(err>0) then
        mess='GetLMOMEstimate:'//trim(mess);return
    endif
    par(1)=par_dummy(1)
    if(par_dummy(1)==0._mrk) then
        err=1;mess='GetLMOMEstimate:Fatal:zero-location';return
    else
        par(2)=par_dummy(2)/par_dummy(1)
    endif
    par(3)=par_dummy(3)
case(GEV)
    call GEV_LMOM(X,par,err,mess)
    if(err>0) then
        mess='GetLMOMEstimate:'//trim(mess);return
    endif
case(GEV_reparameterized)
    call GEV_LMOM(X,par_dummy,err,mess)
    if(err>0) then
        mess='GetLMOMEstimate:'//trim(mess);return
    endif
    par(1)=par_dummy(1)
    if(par_dummy(1)==0._mrk) then
        err=1;mess='GetLMOMEstimate:Fatal:zero-location';return
    else
        par(2)=par_dummy(2)/par_dummy(1)
    endif
    par(3)=par_dummy(3)
case default
    err=1;mess='GetLMOMEstimate:Fatal:Unavailable [DistID]';return
end select
end subroutine GetLMOMEstimate
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetUncertainty_Bootstrap(DistID,X,Mtype,Nsim,par,err,mess)
!^**********************************************************************
!^* Purpose: Generate (standard) Bootstrap samples of parameter estimates
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 06/04/2012
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: 
!^**********************************************************************
!^* IN
!^*		1. DistID, assumed distribution
!^*		2. X, sample used to estimate
!^*		3. [Mtype], 'MOM'(default),'LMOM'
!^*		4. Nsim, number of Bootstrap replicates
!^* OUT
!^*		1.par, parameter estimate
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
character(*), intent(in)::DistID
character(*), intent(in),optional::Mtype
real(mrk), intent(in)::X(:)
integer(mik), intent(in)::Nsim
real(mrk), intent(out)::par(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::i
real(mrk)::Boot(size(X))

err=0;mess='';par=UndefRN
! check size
if(size(par,dim=1)/=Nsim) then
    err=1;mess='GetUncertainty_Bootstrap:Fatal:Size mismatch [par,Nsim]';return
endif
do i=1,Nsim
    !generate bootstrap sample
    call CreateBootstrapSample(X=X,Boot=Boot,err=err,mess=mess)
    if(err>0) then
        mess="GetUncertainty_Bootstrap: "//trim(mess);return
    endif
    ! Estimate
    call GetMomentEstimate(DistID=DistID,X=Boot,Mtype=Mtype,par=par(i,:),err=err,mess=mess)
    if(err>0) then
        mess="GetUncertainty_Bootstrap: "//trim(mess);return
    endif
enddo
end subroutine GetUncertainty_Bootstrap
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetUncertainty_ParametricBootstrap(DistID,X,Mtype,Nsim,par0,par,err,mess)
!^**********************************************************************
!^* Purpose: Generate Parametric Bootstrap samples of parameter estimates
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 06/04/2012
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: 
!^**********************************************************************
!^* IN
!^*		1. DistID, assumed distribution
!^*		2. X, sample used to estimate
!^*		3. [Mtype], 'MOM'(default),'LMOM' 
!^*		4. Nsim, number of Bootstrap replicates
!^*		5. [par0], parameter estimate. Samples will be generated from  
!^*		           the distribution with parameters par0 (if present)  
!^*		           If not-present, par0 is re-estimated within this sub 
!^* OUT
!^*		1.par, parameter estimate
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use Distribution_tools
character(*), intent(in)::DistID
character(*), intent(in),optional::Mtype
real(mrk), intent(in)::X(:)
real(mrk), intent(in),optional::par0(:)
integer(mik), intent(in)::Nsim
real(mrk), intent(out)::par(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::i
real(mrk)::Boot(size(X)),parZero(size(par,dim=2))
logical:: feas

err=0;mess='';par=UndefRN
! check size
if(size(par,dim=1)/=Nsim) then
    err=1;mess='GetUncertainty_ParametricBootstrap:Fatal:Size mismatch [par,Nsim]';return
endif
if(present(par0)) then
    if(size(par0)/=size(parZero)) then
        err=1;mess='GetUncertainty_ParametricBootstrap:Fatal:Size mismatch [par0,parZero]';return
    endif
    parZero=par0
else
    call GetMomentEstimate(DistID=DistID,X=X,Mtype=Mtype,par=parZero,err=err,mess=mess)
    if(err>0) then
        mess="GetUncertainty_ParametricBootstrap:"//trim(mess);return
    endif
endif
do i=1,Nsim
    !generate bootstrap sample
    call GenerateSample(DistId=DistId,par=parZero,gen=Boot,feas=feas,err=err,mess=mess)
    if(err>0) then
        mess="GetUncertainty_ParametricBootstrap:"//trim(mess);return
    endif
    if (.not.feas) then
        err=1;mess='GetUncertainty_ParametricBootstrap:Fatal:unfeasible [parZero]';return
    endif
    ! Estimate
    call GetMomentEstimate(DistID=DistID,X=Boot,Mtype=Mtype,par=par(i,:),err=err,mess=mess)
    if(err>0) then
        mess="GetUncertainty_ParametricBootstrap:"//trim(mess);return
    endif
enddo
end subroutine GetUncertainty_ParametricBootstrap
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! private
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine CreateBootstrapSample(X,Boot,err,mess)
use numerix_dmsl_kit,only:uniran
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::Boot(Size(X))
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::n,Indx(size(X)),i

err=0;mess=''
n=size(X)
do i=1,n
    call uniran(Indx(i),1,n)
enddo
Boot=X(Indx)
end subroutine CreateBootstrapSample
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine Gaussian_MOM(X,par,err,mess)
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::mean,std

err=0;mess='';par=UndefRN
call GetEmpiricalStats(x=X,mean=mean,std=std,&
                       err=err,mess=mess)
if(err>0) then
    mess='Gaussian_MOM:'//trim(mess);return
endif
par(1)=mean
par(2)=std
end subroutine Gaussian_MOM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine Gaussian_LMOM(X,par,err,mess)
use utilities_dmsl_kit, only:pi
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::LMOM(4)

err=0;mess='';par=UndefRN
call Get4LMoments(X,LMOM,err,mess)
if(err>0) then
    mess='Gaussian_LMOM:'//trim(mess);return
endif
par(1)=LMOM(1)
par(2)=LMOM(2)*sqrt(pi)
end subroutine Gaussian_LMOM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine Gumbel_MOM(X,par,err,mess)
use utilities_dmsl_kit, only:pi,euler
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::mean,std

err=0;mess='';par=UndefRN
call GetEmpiricalStats(x=X,mean=mean,std=std,&
                       err=err,mess=mess)
if(err>0) then
    mess='Gumbel_MOM:'//trim(mess);return
endif
par(2)=(sqrt(6._mrk)/pi)*std
par(1)=mean-par(2)*euler
end subroutine Gumbel_MOM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine Gumbel_LMOM(X,par,err,mess)
use utilities_dmsl_kit, only:euler
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::LMOM(4)

err=0;mess='';par=UndefRN
call Get4LMoments(X,LMOM,err,mess)
if(err>0) then
    mess='Gumbel_LMOM:'//trim(mess);return
endif
par(2)=LMOM(2)/log(2._mrk)
par(1)=LMOM(1)-euler*par(2)
end subroutine Gumbel_LMOM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine PearsonIII_MOM(X,par,err,mess)
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::mean,std,skew

err=0;mess='';par=UndefRN
call GetEmpiricalStats(x=X,mean=mean,std=std,skewness=skew,&
                       err=err,mess=mess)
if(err>0) then
    mess='PearsonIII_MOM:'//trim(mess);return
endif
if(skew==0._mrk) then
    err=1;mess='PearsonIII_MOM:Fatal:Skew=0, MOM estimator undefined';return
endif
par(3)=4._mrk/skew**2
par(2)=sign(1._mrk,skew)*(std/sqrt(par(3)))
par(1)=mean-par(2)*par(3)
end subroutine PearsonIII_MOM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine PearsonIII_LMOM(X,par,err,mess)
use utilities_dmsl_kit, only:pi,gammaln
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::LMOM(4),z,a,gamma,sigma,mu

err=0;mess='';par=UndefRN
call Get4LMoments(X,LMOM,err,mess)
if(err>0) then
    mess='PearsonIII_LMOM:'//trim(mess);return
endif
if(abs(LMOM(3))<1._mrk/3._mrk) then
    z=3._mrk*pi*LMOM(3)**2
    a=(1._mrk+0.2906_mrk*z)/(z+0.1882_mrk*z**2+0.0442_mrk*z**3)
else
    z=1._mrk-abs(LMOM(3))
    a=(0.36067_mrk*z-0.59567_mrk*z**2+0.25361_mrk*z**3)/&
    (1._mrk-2.78861_mrk*z+2.56096_mrk*z**2-0.77045_mrk*z**3)
endif
gamma=2._mrk*sign(1._mrk,LMOM(3))/sqrt(a)
sigma=LMOM(2)*sqrt(pi*a)*exp(gammaln(a)-gammaln(0.5_mrk+a))
mu=LMOM(1)
par(1)=mu-2._mrk*(sigma/gamma)
par(2)=0.5*sigma*gamma
par(3)=4._mrk/gamma**2
end subroutine PearsonIII_LMOM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GEV_LMOM(X,par,err,mess)
use utilities_dmsl_kit, only:gammaln
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::LMOM(4),c

err=0;mess='';par=UndefRN
call Get4LMoments(X,LMOM,err,mess)
if(err>0) then
    mess='GEV_LMOM:'//trim(mess);return
endif
c=2._mrk/(3._mrk+LMOM(3))-log(2._mrk)/log(3._mrk)
par(3)=7.8590_mrk*c+2.9554_mrk*c**2
par(2)=(LMOM(2)*par(3))/((1._mrk-2._mrk**(-1._mrk*par(3)))*exp(gammaln(1._mrk+par(3))))
par(1)=LMOM(1)+(par(2)/par(3))*(1._mrk-exp(gammaln(1._mrk+par(3))))
end subroutine GEV_LMOM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GEV_MOM(X,par,err,mess)
use numerix_dmsl_kit, only:zbrent
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit, only:gammaln
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::mean,std,skew,fzero,g1,g2
integer(mik)::fcalls
type(data_ricz_type)::dataIN
real(mrk), parameter::SmallShape=-1._mrk/3._mrk+0.000001_mrk,&
                      LargeShape=4._mrk,&
                      TolX=0.001_mrk,TolF=0.001_mrk,&
                      xscale=0.1_mrk,fscale=1._mrk
integer(mik), parameter::itmax=10000

err=0;mess='';par=UndefRN
call GetEmpiricalStats(x=X,mean=mean,std=std,skewness=skew,&
                       err=err,mess=mess)
if(err>0) then
    mess='GEV_MOM:'//trim(mess);return
endif
dataIN%rs1=skew
! Solve Equation for Shape Estimate
call  zbrent(evalFunc=GEV_MOM_ShapeRoot,dataIN=dataIN,&
             x1=SmallShape,x2=LargeShape,&
             tolX=tolX,tolF=tolF,&
             xscale=xscale,fscale=fscale,&
             itmax=itmax,&
             xroot=par(3),froot=fzero,fcalls=fcalls,&
             err=err,message=mess)
if(err>0) then
    mess='GEV_MOM:'//trim(mess);return
endif
g1=exp(gammaln(1._mrk+par(3)))
g2=exp(gammaln(1._mrk+2._mrk*par(3)))
par(2)=(std*abs(par(3)))/sqrt(g2-g1**2)
par(1)=mean-par(2)/par(3)+(par(2)/par(3))*g1
end subroutine GEV_MOM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GEV_MOM_ShapeRoot(dataIN,dataOUT,x,feas,fx,dfdxV,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  use utilities_dmsl_kit, only:gammaln
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,dfdxV(:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  !locals
  real(mrk)::skew,g1,g2,g3,Tskew
  real(mrk), parameter::TinyShape=0.001_mrk,GumbelSkew=1.139547_mrk
 
  message='';err=0;feas=.true.
  if(present(dfdxV)) dfdxV=UndefRN
  if(x<=-1._mrk/3._mrk) then
      feas=.false.;return
  endif
  skew=dataIN%rs1 ! empirical skew
  ! compute theoretical skewness
  g1=exp(gammaln(1._mrk+x))
  g2=exp(gammaln(1._mrk+2._mrk*x))
  g3=exp(gammaln(1._mrk+3._mrk*x))
  if(abs(x)<=TinyShape) then ! switch to Gumbel Formula
    Tskew=GumbelSkew
  else
    Tskew=-1._mrk*sign(1._mrk,x)*(g3-3._mrk*g1*g2+2._mrk*g1**3)/((g2-g1**2)**1.5_mrk)
  endif
  fx=Tskew-skew 
end subroutine GEV_MOM_ShapeRoot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end module MomentEstimation_tools
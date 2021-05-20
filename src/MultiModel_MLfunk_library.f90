module MultiModel_MLfunk_library

!~**********************************************************************
!~* Purpose: Catalogue of Max-Lkh functions (for COMPLEX)
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 27/07/2015
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
use types_dmsl_kit, only:data_ricz_type
use Optimisation_tools
implicit none
Private
public :: MLfunk,&
          COMPLEX_GetPattern

!==================!
! TYPES DEFINITION !
!==================!

type,public:: BasicDataType
    integer(mik)::Nt=undefIN,Nd=undefIN
    real(mrk), allocatable::Z(:,:) ! Nt*Nd, data
    real(mrk), allocatable::Zaux(:,:) ! used to store any other auxiliary variable
end type BasicDataType

type,public:: SingleMdataType
    integer(mik)::Nrep=undefIN
    type(BasicDataType), allocatable::rep(:)
    logical, allocatable::DirtyRep(:) ! where true, rep is affected by observation errors
end type SingleMdataType

type,public:: InferType ! inference properties
    character(250)::ID="aintgotnoname"
    integer(mik)::Npar=undefIN
end type InferType

!==================!
! GLOBAL VARIABLES !
!==================!
Character(250), parameter, PUBLIC:: & ! Models' catalogue
            MM_IID='MM_IID',& ! IID data, distribution is given in option%cs1
            MM_MGaussian='MM_MultivariateGaussian',& ! IID multivariate gaussian data
            MM_COMPLEX='MM_COMPLEX'   ! COMPLEX model

type(SingleMdataType),public::GX
real(mrk),allocatable,public::GSigmaInv(:,:)
real(mrk),public::GSigmaDet ! log-determinant
type(InferType),public::Ginfer
Type(data_ricz_type),public::Goption

!============!
! PARAMETERS !
!============!
real(mrk),parameter::scalecoeff=1._mrk,scaledefault=1._mrk
! Optimizer options
integer(mik),parameter::Optim_meth=qndmsl_smeth,Optim_nstart=1,Optim_activeSet=0,Optim_iterMax=100000
real(mrk), parameter::Optim_msLo=-HugeRe,Optim_msHi=HugeRe
character(250),parameter::Optim_msFile="Optim_MultiStart_Results.txt",&
                          Optim_hisFile="Optim_History.txt",&
                          Optim_resFile="Optim_Results.txt"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============!
! PUBLIC SUBS !
!=============!
subroutine MLfunk(par,eps,fmax,n,err,mess)
!^**********************************************************************
!^* Purpose: ML parameter estimation
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^* OUT
!^*		1.par (size Npar) Max-lkh estimates
!^*		2.eps (size Npar) Max-lkh estimates of obs errors epsilon
!^*		3.fmax, Max-lkh value
!^*		4.n, total number of individuals for computing Lkh
!^*		5.err, error code
!^*		6.mess, error message
!^**********************************************************************
real(mrk),intent(out)::par(:),eps(:),fmax
integer(mik),intent(out)::n,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='MLfunk'
real(mrk)::par0(size(par)),parscale(size(par)),parlow(size(par)),parhigh(size(par))
real(mrk)::eps0(size(eps)),fmin
logical::feas,isnull
real(mrk),allocatable::x0(:),xLo(:),xHi(:),xScale(:),xMin(:),Hessian(:,:),mslo(:),mshi(:)
integer(mik),allocatable::activeset(:)
integer(mik)::npar,np,fcalls,gcalls,hcalls

err=0;mess='';n=sum(GX%Rep(:)%Nt)
par=undefRN;eps=undefRN;fmax=undefRN
npar=size(par)

! Get rough estimation
call RoughEstimate(GX,Goption,par0,parscale,parlow,parhigh,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

! SPEED-UP: stop here for cases where rough estimation is the ML estimator
select case(trim(Ginfer%ID))
case(MM_IID)
    if( (.not.any(GX%DirtyRep)) .and. (Goption%cs1=="Gaussian")) then ! no need to optimize!
        par=par0
        call LLfunk_iid(par,eps,fmax,feas,isnull,err,mess)
        return
    endif
case(MM_MGaussian)
    if( (.not.any(GX%DirtyRep)) ) then ! no need to optimize!
        par=par0
        call LLfunk_mgauss(par,eps,fmax,feas,isnull,err,mess)
        return
    endif
end select

! OPTIMIZATION
if(any(GX%DirtyRep)) then ! need to concatenate par and eps
    np=2*npar
    allocate(x0(np),xLo(np),xHi(np),xScale(np),xMin(np),Hessian(np,np),mslo(np),mshi(np),activeset(np))
    x0(1:npar)=par0;x0((npar+1):)=0._mrk;
    xscale(1:npar)=parscale;xscale((npar+1):)=parscale;
    xLo(1:npar)=parlow;xLo((npar+1):)=-HugeRe;
    xHi(1:npar)=parhigh;xHi((npar+1):)=HugeRe;
else ! par only
    np=npar
    allocate(x0(np),xLo(np),xHi(np),xScale(np),xMin(np),Hessian(np,np),mslo(np),mshi(np),activeset(np))
    x0=par0;xscale=parscale;xLo=parlow;xHi=parhigh
endif
mslo=Optim_msLo;mshi=Optim_msHi;
activeset=Optim_activeSet
call optimiserWrap (smeth=Optim_meth,& ! optim method
                    nopt=Optim_nstart,& ! multi-start
                    evalFunc=wrap,& ! objective function
                    nDim=np,&
                    xLo=xLo,xHi=xHi,&
                    msLo=mslo,msHi=mshi,&
                    xscale=xscale,&
                    activeSet=activeset,&
                    msFile=Optim_msFile, &
                    imeth_qn=5,gmeth=1,hmeth_qn=6,himeth_qn=1,&  ! QN options
                    mMem_lbfgs=10,&                      ! LBFGS method
                    mometh=Optim_meth,&
                    x0=x0, &
                    xMin=xMin, fMin=fMin, &
                    chkGrad=.false.,chkHess=.false.,&
                    inverseHessian=.false., hess=Hessian,&
                    iterMax=Optim_iterMax, hisFile=Optim_hisFile, &
                    iternfo=0,&
                    resFile=Optim_resFile,&
                    fcalls=fcalls,gcalls=gcalls,hcalls=hcalls,&
                    err=err, message=mess)
if (err/=0) then;write(*,*) 'WARNING: Optimization failed';endif
fmax=-1._mrk*fmin
if(any(GX%DirtyRep)) then ! need to separate par and eps
    par=xMin(1:npar)
    eps=xMin((npar+1):)
else
    par=xMin
endif
deallocate(x0,xLo,xHi,xScale,xMin,Hessian,mslo,mshi,activeset)

end subroutine MLfunk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============!
! GENERAL SUBS !
!==============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LLfunk(par,eps,LL,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Log-likelihood function
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.par (size Npar) parameters
!^*		2.eps (size Npar) obs errors epsilon
!^* OUT
!^*		1. LL, log-likelihood
!^*		2. feas, feasability of parameters 
!^*		3. isNull, is likelihood=0?
!^*		4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5. mess, error message
!^**********************************************************************
real(mrk),intent(in)::par(:),eps(:)
real(mrk),intent(out)::LL
logical,intent(out)::feas,isnull
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='LLfunk'

select case(trim(Ginfer%ID))
case(MM_IID)
    call LLfunk_iid(par,eps,LL,feas,isnull,err,mess)
case(MM_MGaussian)
    call LLfunk_mgauss(par,eps,LL,feas,isnull,err,mess)
case(MM_COMPLEX)
    call LLfunk_COMPLEX(par,eps,LL,feas,isnull,err,mess)
case default
    err=1;mess=trim(procname)//":Fatal:unknown infer%ID"
end select

end subroutine LLfunk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RoughEstimate(X,option,par,parscale,parlow,parhigh,err,mess)
!^**********************************************************************
!^* Purpose: rough estimate
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.X, single model object
!^*		2.option, options
!^* OUT
!^*		1.par (size Npar) rough estimates
!^*		1.parscale (size Npar) scale
!^*		1.parlow (size Npar) lower bound
!^*		1.parhigh (size Npar) higher bound
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
use Distribution_tools, only:GetRoughEstimate
type(SingleMdataType), intent(in)::X
Type(data_ricz_type), intent(in)::option
real(mrk),intent(out)::par(:),parscale(:),parlow(:),parhigh(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='RoughEstimate'

select case(trim(Ginfer%ID))
case(MM_IID)
    call RoughEstimate_iid(X,option,par,parscale,parlow,parhigh,err,mess)
case(MM_MGaussian)
    call RoughEstimate_mgauss(X,option,par,parscale,parlow,parhigh,err,mess)
case(MM_COMPLEX)
    call RoughEstimate_COMPLEX(X,option,par,parscale,parlow,parhigh,err,mess)
case default
    err=1;mess=trim(procname)//":Fatal:unknown infer%ID"
end select

end subroutine RoughEstimate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wrap(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
! interface to make optimizer happy
use kinds_dmsl_kit
use types_dmsl_kit,only:data_ricz_type
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x(:)
logical(mlk),intent(out)::feas
real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
logical::isnull
integer(mik)::nx
real(mrk),allocatable::par(:),eps(:)

err=0;message='';feas=.true.
nx=size(x)
if(.not.present(fx)) return
if(any(GX%DirtyRep)) then
  allocate(par(nx/2),eps(nx/2))
  par=x(1:nx/2);eps=x((nx/2+1):)
else
  allocate(par(nx),eps(nx))
  par=x;eps=0._mrk
endif
call LLfunk(par=par,eps=eps,LL=fx,feas=feas,isnull=isnull,err=err,mess=message)
if(err>0) return
if(present(fx)) then
  if(isnull .or. (.not.feas)) then
    fx=-HugeRe
  endif
  fx=-1*fx
endif
if(present(gradFx)) gradFx=undefRN
if(present(hessFx)) hessFx=undefRN
endsubroutine wrap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================!
! INFER-SPECIFIC subs !
!=====================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LLfunk_iid(par,eps,LL,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Log-likelihood for iid data 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.par (size Npar) parameters
!^*		2.eps (size Npar) obs errors epsilon
!^* OUT
!^*		1. LL, log-likelihood
!^*		2. feas, feasability of parameters 
!^*		3. isNull, is likelihood=0?
!^*		4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5. mess, error message
!^**********************************************************************
use Distribution_tools, only:GetSamplePdf
real(mrk),intent(in)::par(:),eps(:)
real(mrk),intent(out)::LL
logical,intent(out)::feas,isnull
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='LLfunk_iidGauss'
integer(mik)::i
real(mrk)::pdf,param(size(par))

err=0;mess='';feas=.true.;isnull=.false.;LL=0._mrk
do i=1,GX%Nrep
    if(GX%DirtyRep(i)) then ! add observation errors
        param=par+eps
    else ! keep original par values
        param=par
    endif
    call GetSamplePdf(DistId=Goption%cs1,x=GX%rep(i)%Z(:,1),par=param,&
             loga=.true.,pdf=pdf,feas=feas,isnull=isnull,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if((.not.feas) .or. isnull) return
    LL=LL+pdf
enddo
! Add hierarchical component if needed
if(any(GX%DirtyRep)) then 
    LL=LL+mnorm0_logp(eps,GSigmaInv,GSigmaDet)
endif

end subroutine LLfunk_iid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LLfunk_mgauss(par,eps,LL,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Log-likelihood for iid multivariate gaussian data 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 13/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.par (size Npar) parameters
!^*		2.eps (size Npar) obs errors epsilon
!^* OUT
!^*		1. LL, log-likelihood
!^*		2. feas, feasability of parameters 
!^*		3. isNull, is likelihood=0?
!^*		4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5. mess, error message
!^**********************************************************************
use linalg_dmsl_kit, only:choles_invrt
real(mrk),intent(in)::par(:),eps(:)
real(mrk),intent(out)::LL
logical,intent(out)::feas,isnull
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='LLfunk_mgauss'
integer(mik)::i,j,k,Nd,compt
real(mrk)::pdf,SigmaDet,param(size(par)),mu(size(GX%rep(1)%Z,dim=2)),&
           sigma(size(GX%rep(1)%Z,dim=2),size(GX%rep(1)%Z,dim=2)),&
           sigmainv(size(GX%rep(1)%Z,dim=2),size(GX%rep(1)%Z,dim=2))
logical::posdef

err=0;mess='';feas=.true.;isnull=.false.;LL=0._mrk
Nd=size(GX%rep(1)%Z,dim=2)
do i=1,GX%Nrep
    if(GX%DirtyRep(i)) then ! add observation errors
        param=par+eps
    else ! keep original par values
        param=par
    endif
    ! unpack param
    mu=param(1:Nd)
    do j=1,Nd
        sigma(j,j)=param(Nd+j)
    enddo
    compt=2*Nd
    do j=2,Nd
        do k=1,j-1
            compt=compt+1
            sigma(j,k)=param(compt)
            sigma(k,j)=sigma(j,k)
        enddo
    enddo
    ! invert sigma
    posdef=.true.
    call choles_invrt(a=Sigma,ainv=SigmaInv,posDefinite=posdef,logDet=SigmaDet,err=err,message=mess)
    if(.not.posdef) then
        feas=.false.;return
    endif
    if(err/=0) then
        mess=trim(procname)//':'//trim(mess);feas=.false.;return
    endif
    do j=1,size(GX%rep(i)%Z,dim=1)
        pdf=mnorm0_logp(GX%rep(i)%Z(j,:)-mu,SigmaInv,SigmaDet)
        LL=LL+pdf
    enddo
enddo
! Add hierarchical component if needed
if(any(GX%DirtyRep)) then 
    pdf=mnorm0_logp(eps,GSigmaInv,GSigmaDet)
    LL=LL+pdf
endif

end subroutine LLfunk_mgauss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LLfunk_COMPLEX(par,eps,LL,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Log-likelihood for COMPLEX model 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 13/08/2015
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.par (size Npar) parameters
!^*		2.eps (size Npar) obs errors epsilon
!^* OUT
!^*		1. LL, log-likelihood
!^*		2. feas, feasability of parameters 
!^*		3. isNull, is likelihood=0?
!^*		4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5. mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:twoPi
use Distribution_tools, only:GetSamplePdf,GetPdf

real(mrk),intent(in)::par(:),eps(:)
real(mrk),intent(out)::LL
logical,intent(out)::feas,isnull
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk),parameter::ro_epsilon=0.005_mrk
real(mrk),allocatable::mu(:),sigma(:),res(:),innov(:)
real(mrk)::ro,pdf,maxtime,param(size(par))
integer(mik)::i,j,k,Nt,mu_Nseason,mu_Nleg,sig_Nseason,sig_Nleg,Nd

! Initialise
err=0;mess='';isNull=.false.;feas=.true.;LL=0._mrk
! Make sense of options
maxtime=Goption%rs1
mu_Nseason=Goption%is1
mu_Nleg=Goption%is2
sig_Nseason=Goption%is3
sig_Nleg=Goption%is4
Nd=size(GX%rep(1)%Z,dim=2)
! compute likelihoods
do i=1,GX%Nrep
    if(GX%DirtyRep(i)) then ! add observation errors
        param=par+eps
    else ! keep original par values
        param=par
    endif
    Nt=GX%rep(i)%Nt
    k=0
    allocate(mu(Nt),sigma(Nt),res(Nt),innov(Nt-1))
    do j=1,Nd
        ! get mean
        call COMPLEX_GetPattern(k=k,par=param,time=GX%rep(i)%Zaux(:,1),maxtime=maxtime,&
                                Nseason=mu_Nseason,Nleg=mu_Nleg,pattern=mu,&
                                feas=feas,err=err,mess=mess)
        if(.not.feas) return
        ! get Standard deviation
        call COMPLEX_GetPattern(k=k,par=param,time=GX%rep(i)%Zaux(:,1),maxtime=maxtime,&
                                Nseason=sig_Nseason,Nleg=sig_Nleg,pattern=sigma,&
                                feas=feas,err=err,mess=mess)
        if(.not.feas) return
        sigma=exp(sigma)
        ! Get residuals and innovations
        ! compute log-lkh
        res=(GX%rep(i)%Z(:,j)-mu)/sigma
        ! -----------------
        ! Evin et al's approach - only works for AR(1)
        k=k+1;ro=param(k)
        if(abs(ro)>=1._mrk-ro_epsilon) then;feas=.false.;return;endif
        pdf = -1._mrk*sum(log(sigma))& ! Jacobian
        -0.5_mrk*log(twoPi)-0.5_mrk*res(1)**2& ! N(0,1) marginal for first step
        -0.5_mrk*real(Nt-1,mrk)*(log(twoPi)+log(1._mrk-ro**2))-0.5_mrk*sum( ((res(2:Nt)-ro*res(1:Nt-1))**2)/(1._mrk-ro**2) ) ! conditional: N(ro*res(t-1),sqrt(1-ro**2)
        LL=LL+pdf
    enddo
    deallocate(mu,sigma,res,innov)
enddo
! Add hierarchical component if needed
if(any(GX%DirtyRep)) then 
    pdf=mnorm0_logp(eps,GSigmaInv,GSigmaDet)
    LL=LL+pdf
endif
end subroutine LLfunk_COMPLEX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RoughEstimate_iid(X,option,par,parscale,parlow,parhigh,err,mess)
!^**********************************************************************
!^* Purpose: rough estimate for iid data 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.X, single model object
!^*		2.option, options
!^* OUT
!^*		1.par (size Npar) rough estimates
!^*		1.parscale (size Npar) scale
!^*		1.parlow (size Npar) lower bound
!^*		1.parhigh (size Npar) higher bound
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
use Distribution_tools, only:GetRoughEstimate
type(SingleMdataType), intent(in)::X
Type(data_ricz_type), intent(in)::option
real(mrk),intent(out)::par(:),parscale(:),parlow(:),parhigh(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='RoughEstimate_iid'
integer(mik)::i,k
real(mrk)::foo(sum(X%rep(:)%Nt))

err=0;mess='';
par=undefRN;parscale=undefRN;parlow=undefRN;parhigh=undefRN
! pack all replications
k=0
do i=1,X%nrep
    foo((k+1):(k+X%rep(i)%Nt))=X%rep(i)%Z(:,1)
    k=k+X%rep(i)%Nt    
enddo
! Get rough estimate
call GetRoughEstimate(DistID=option%cs1,X=foo,par=par,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
parscale=scalecoeff*abs(par)
where(parscale==0._mrk) parscale=scaledefault
parlow=-HugeRe
parhigh=HugeRe
end subroutine RoughEstimate_iid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RoughEstimate_mgauss(X,option,par,parscale,parlow,parhigh,err,mess)
!^**********************************************************************
!^* Purpose: rough estimate for iid multivariate gaussian data 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 13/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.X, single model object
!^*		2.option, options
!^* OUT
!^*		1.par (size Npar) rough estimates
!^*		1.parscale (size Npar) scale
!^*		1.parlow (size Npar) lower bound
!^*		1.parhigh (size Npar) higher bound
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
use numerix_dmsl_kit, only:getmeanvar
type(SingleMdataType), intent(in)::X
Type(data_ricz_type), intent(in)::option
real(mrk),intent(out)::par(:),parscale(:),parlow(:),parhigh(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='RoughEstimate_mgauss'
integer(mik)::i,j,k,Nd
real(mrk)::foo(sum(X%rep(:)%Nt),size(X%rep(1)%Z,dim=2)),&
           mu(size(X%rep(1)%Z,dim=2)),&
           sigma(size(X%rep(1)%Z,dim=2),size(X%rep(1)%Z,dim=2))

err=0;mess='';
par=undefRN;parscale=undefRN;parlow=undefRN;parhigh=undefRN
Nd=size(X%rep(1)%Z,dim=2)
! pack all replications
k=0
do i=1,X%nrep
    foo((k+1):(k+X%rep(i)%Nt),:)=X%rep(i)%Z
    k=k+X%rep(i)%Nt    
enddo
! Get mean and covar
call getmeanvar(x=foo,package=2,mean=mu,covar=sigma,method="c",err=err,message=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
! pack into par
par(1:Nd)=mu
do i=1,Nd
    par(Nd+i)=sigma(i,i)
enddo
k=2*Nd
do i=2,Nd
    do j=1,i-1
        k=k+1
        par(k)=sigma(i,j)
    enddo
enddo
parscale=scalecoeff*abs(par)
where(parscale==0._mrk) parscale=scaledefault
parlow=-HugeRe
parhigh=HugeRe
end subroutine RoughEstimate_mgauss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RoughEstimate_COMPLEX(X,option,par,parscale,parlow,parhigh,err,mess)
!^**********************************************************************
!^* Purpose: rough estimate for COMPLEX model 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 13/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.X, single model object
!^*		2.option, options
!^* OUT
!^*		1.par (size Npar) rough estimates
!^*		1.parscale (size Npar) scale
!^*		1.parlow (size Npar) lower bound
!^*		1.parhigh (size Npar) higher bound
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
use linalg_dmsl_kit, only:choles_invrt
use EmpiricalStats_tools, only:GetEmpiricalStats

type(SingleMdataType), intent(in)::X
Type(data_ricz_type), intent(in)::option
real(mrk),intent(out)::par(:),parscale(:),parlow(:),parhigh(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='RoughEstimate_COMPLEX'
integer(mik)::i,j,k,Nt,mu_Nseason,mu_Nleg,sig_Nseason,sig_Nleg,Nd
real(mrk),allocatable::time(:),B(:,:),Binv(:,:),res(:),stdres(:),foo(:)
real(mrk)::moy,sdev,maxtime
logical,allocatable::mask(:)
logical::pd

err=0;mess='';
par=undefRN;parscale=undefRN;parlow=undefRN;parhigh=undefRN
! Make sense of options
maxtime=Goption%rs1
mu_Nseason=Goption%is1
mu_Nleg=Goption%is2
sig_Nseason=Goption%is3
sig_Nleg=Goption%is4
Nt=X%rep(1)%Nt
Nd=size(X%rep(1)%Z,dim=2)
k=0
allocate(mask(Nt),time(Nt),B(Nt,mu_Nseason+mu_Nleg),Binv(mu_Nseason+mu_Nleg,mu_Nseason+mu_Nleg),res(Nt),stdres(Nt))
time=X%rep(1)%Zaux(:,1)
do j=1,Nd
    moy=sum(X%rep(1)%Z(:,j))/real(Nt,mrk)
    sdev=sqrt(sum( (X%rep(1)%Z(:,j)-moy)**2 )/real(Nt,mrk))
    ! Start by estimating S & L for the mean using the standard OLS estimate
    B=0._mrk
    do i=1,mu_Nleg
        B(:,i)=LegendrePoly(2._mrk*(time/maxtime)-1._mrk,i)
    enddo
    do i=1,mu_Nseason
        mask=(time-floor(time))<real(i,mrk)/real(mu_Nseason,mrk) .and. (time-floor(time))>=real(i-1,mrk)/real(mu_Nseason,mrk)
        where (mask) B(:,i+mu_Nleg)=1._mrk
    enddo
    call choles_invrt(a=matmul(transpose(B),B),ainv=Binv,posDefinite=pd,err=err,message=mess)
    if(err/=0 .or. .not.pd)then
        err=1;mess=trim(procname)//':'//'BtB not invertible';return
    endif
    par( (k+1):(k+mu_Nleg+mu_Nseason) )=matmul(matmul(Binv,transpose(B)),X%rep(1)%Z(:,j))
    ! Get the residuals and estimate seasonnal sdevs
    res=X%rep(1)%Z(:,j)-matmul(B,par( (k+1):(k+mu_Nleg+mu_Nseason) ))
    k=k+mu_Nleg+mu_Nseason
    if(sig_Nleg>0) then
        par( (k+1):(k+sig_Nleg) )=0._mrk
        k=k+sig_Nleg
    endif
    do i=1,sig_Nseason
        mask=(time-floor(time))<real(i,mrk)/real(sig_Nseason,mrk) .and. (time-floor(time))>=real(i-1,mrk)/real(sig_Nseason,mrk)
        if(allocated(foo)) deallocate(foo);allocate(foo(count(mask)))
        foo=pack(res,mask)
        call GetEmpiricalStats(x=foo, std=par(k+1),err=err,mess=mess)
        where(mask) stdres=res/par(k+1)
        par(k+1)=log(par(k+1))
        k=k+1
    enddo
    ! Get the lag-1 coeff of standardized residuals
    moy=sum(stdres)/real(Nt,mrk)
    par(k+1)=dot_product(stdres(1:(Nt-1))-moy,stdres(2:Nt)-moy)/dot_product(stdres-moy,stdres-moy)
    k=k+1
enddo
! parLo and parHi
parscale=scalecoeff*abs(par)
where(parscale==0._mrk) parscale=scaledefault
parlow=-HugeRe
parhigh=HugeRe
deallocate(mask,time,B,Binv,res,stdres)
end subroutine RoughEstimate_COMPLEX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============!
! MISC. SUBS !
!============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function mnorm0_logp(x,SigmaInv,logDet)
! log-pdf of a multinormal dist with mean zero
use utilities_dmsl_kit,only:half,lntwopi,quadform
real(mrk), intent(in)::x(:),SigmaInv(:,:),logDet
real(mrk)::mnorm0_logp
integer(mik)::n
n=size(x)
mnorm0_logp=-half*(quadform(x,SigmaInv)+real(n,mrk)*lntwopi+logDet)
end function mnorm0_logp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elemental function LegendrePoly(x,k)
! Legendre polynomials, useful for describing non-linear trends
real(mrk),intent(in)::x
integer(mik),intent(in)::k
real(mrk):: LegendrePoly
select case (k)
case(0)
    LegendrePoly=1._mrk
case(1)
    LegendrePoly=x
case(2)
    LegendrePoly=0.5_mrk*(3._mrk*x**2-1._mrk)
case(3)
    LegendrePoly=0.5_mrk*(5._mrk*x**3-3._mrk*x)
case(4)
    LegendrePoly=0.125_mrk*(35._mrk*x**4-30._mrk*x**2+3._mrk)
case(5)
    LegendrePoly=0.125_mrk*(63._mrk*x**5-70._mrk*x**3+15._mrk*x)
case default
    LegendrePoly=undefRN
end select
end function LegendrePoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function SeasonnalSignal(x,Nseason,par)
! linear interpolation between seasonal values stored in par
use utilities_dmsl_kit, only:quickLinInterp
integer(mik),intent(in)::Nseason
real(mrk),intent(in)::x,par(:)
real(mrk):: SeasonnalSignal
integer(mik)::i
SeasonnalSignal=quickLinInterp(x=x-floor(x),xx=(/(real(i,mrk)/real(Nseason,mrk),i=0,Nseason)/),&
                               yy=(/par,par(1)/), guessI=1+floor(Nseason*x) )
end function SeasonnalSignal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine COMPLEX_GetPattern(k,par,time,maxtime,Nseason,Nleg,pattern,feas,err,mess)
! Get temporal pattern using trend & seasonnality models
integer(mik), intent(inout)::k
integer(mik),intent(in)::Nseason,Nleg
real(mrk),intent(in)::par(:),time(:),maxtime
real(mrk),intent(out)::pattern(:)
logical,intent(out)::feas
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::i,t,Nt
real(mrk)::Leg(size(time),Nleg)

err=0;mess='';feas=.true.
Nt=size(time)
! Legendre polynomials
if(Nleg>0) then
    do i=1,Nleg
        Leg(:,i)=LegendrePoly(2._mrk*(time/maxtime)-1._mrk,i)
    enddo
    pattern=matmul(Leg,par( (k+1):(k+Nleg) ))
else
    pattern=0._mrk
endif
k=k+Nleg
! Seasonnality cycle
if(Nseason>0) then
    do t=1,Nt
        pattern(t)=pattern(t)+SeasonnalSignal(time(t),Nseason,par( (k+1):(k+Nseason) ))
    enddo
    k=k+Nseason
endif
end subroutine COMPLEX_GetPattern
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module MultiModel_MLfunk_library
module MLfunk_catalogue

!~**********************************************************************
!~* Purpose: Catalogue of Max-Lkh functions (for COMPLEX)
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 04/08/2014
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
public :: MLfunk_PCs,PCs_GetPattern

type,public:: BasicDataType
    integer(mik)::Nt=undefIN,Nd=undefIN
    real(mrk), allocatable::Z(:,:) ! Nt*Nd, data
    real(mrk), allocatable::Zaux(:,:) ! used to store any other auxiliary variable
end type BasicDataType

type,public:: SingleMdataType
    integer(mik)::Nrep=undefIN
    type(BasicDataType), allocatable::rep(:)
    logical, allocatable::DirtyRep(:) ! if true, will add JPV's "observation errors"
end type SingleMdataType
type(SingleMdataType),public::Mdata

type,public:: MoptionsType
    integer(mik)::npar=undefIN
    integer(mik)::N0=undefIN,N1=undefIN,N2=undefIN,N3=undefIN,&
                  N4=undefIN,N5=undefIN,N6=undefIN,N7=undefIN
    real(mrk)::R0=undefRN
    real(mrk),allocatable::V0(:)
    real(mrk),allocatable::M0(:,:)
end type MoptionsType
Type(MoptionsType),public::Mopt

real(mrk), public::DirtySdev
real(mrk),parameter::ro_epsilon=0.05_mrk
logical, parameter::DoSmartStart=.true., DoQuickEstimation=.false.

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine MLfunk_PCs(par,fmax,n,err,mess)
!^**********************************************************************
!^* Purpose: independent stat. models for several PCs of a PCA 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/08/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^* OUT
!^*		1.par (size Npar) Max-lkh estimates
!^*		2.fmax, Max-lkh value
!^*		3.n, total number of individuals for computing Lkh
!^*		4.err, error code
!^*		5.mess, error message
!^**********************************************************************
use kinds_dmsl_kit
use Optimisation_tools
real(mrk),intent(out)::par(:),fmax
integer(mik),intent(out)::n,err
character(*),intent(out)::mess
! locals
real(mrk),allocatable::xlo(:),xhi(:),mslo(:),mshi(:),&
           x0(:),xscale(:),teta(:),Hessian(:,:)
integer(mik),allocatable::activeSet(:)     
real(mrk)::fmin           
integer(mik)::meth,fcalls,gcalls,hcalls,Npar
character(250)::hisFile,resfile,msFile
logical::feas, isnull

err=0;mess='';n=sum(Mdata%Rep(:)%Nt)
! handle dirtyfication
if(any(Mdata%DirtyRep)) then
    Npar=2*Mopt%Npar
else
    Npar=Mopt%Npar
endif
allocate(xlo(Npar),xhi(Npar),mslo(Npar),mshi(Npar),activeSet(Npar),&
           x0(Npar),xscale(Npar),teta(Npar),Hessian(Npar,Npar))
! Optimize!
activeSet=0
call PCs_GetStartingPoint(any(Mdata%DirtyRep),x0,xscale,xLo,XHi)
if(DoQuickEstimation) then
    par=x0
    call LLfunk_PCs(x=par,feas=feas,isnull=isnull,fx=fmax,err=err,mess=mess)
    if (err/=0 .or. (.not.feas) .or. isnull) then;write(*,*) 'WARNING: Impossible estimate';endif
else
    msLo=-HugeRe
    msHi=HugeRe
    msFile="Optim_MultiStart_Results.txt"
    hisFile="Optim_History.txt"
    resFile="Optim_Results.txt"
    meth=qndmsl_smeth ! lbfgs_smeth !
    call optimiserWrap (smeth=meth,& ! optim method
                        nopt=1,& ! multi-start
                        evalFunc=ObjFunk_PCs,& ! objective function
                        nDim=Npar,&
                        xLo=xLo,xHi=xHi,&
                        msLo=msLo,msHi=msHi,&
                        xscale=xscale,&
                        activeSet=activeSet,&
                        msFile=msFile, &
                        imeth_qn=5,gmeth=1,hmeth_qn=6,himeth_qn=1,&  ! QN options
                        mMem_lbfgs=10,&                      ! LBFGS method
                        mometh=meth,&
                        x0=x0, &
                        xMin=teta, fMin=fMin, &
                        chkGrad=.false.,chkHess=.false.,&
                        inverseHessian=.true., hess=Hessian,&
                        iterMax=100000, hisFile=hisFile, &
                        iternfo=0,&
                        resFile=resFile,&
                        fcalls=fcalls,gcalls=gcalls,hcalls=hcalls,&
                        err=err, message=mess)
    if (err/=0) then;write(*,*) 'WARNING: Optimization failed';endif
    fmax=-1._mrk*fmin
    par=teta
endif
deallocate(xlo,xhi,mslo,mshi,activeSet,&
           x0,xscale,teta,Hessian)
end subroutine MLfunk_PCs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Private subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LLfunk_PCs(x,feas,isnull,fx,fAux,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 01/08/2014
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
!^*		1. feas, feasability of parameter set x
!^*		2. isNull, is fx=0?
!^*		3. fx, f(x)
!^*		4. fAux, unused
!^*		5. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		6. mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:twoPi
use DIstribution_tools, only:GetSamplePdf,GetPdf

real(mrk),intent(in)::x(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk),allocatable::mu(:),sigma(:),res(:),innov(:),par(:)
real(mrk)::ro,fx0
integer(mik)::i,k,j,Nt,Nlag,Npar

! Initialise
err=0;mess='';isNull=.false.;feas=.true.;fx=0._mrk
if(present(fAux)) fAux=UndefRN
Npar=Mopt%Npar
Nlag=Mopt%N7
do j=1,Mdata%Nrep
    if(allocated(par)) deallocate(par);allocate(par(Npar))
    if(Mdata%DirtyRep(j)) then
        par=x(1:Npar)+x( (Npar+1):) !x(1:Npar)*(1._mrk+x( (Npar+1):)) ! Dirtyfication!
    else
        par=x(1:Npar)
    endif
    Nt=Mdata%rep(j)%Nt
    allocate(mu(Nt),sigma(Nt),res(Nt),innov(Nt-Nlag))
    ! get mean
    k=0
    call PCs_GetPattern(k=k,DoWhat=(/Mopt%N1,Mopt%N2/),x=par,time=Mdata%rep(j)%Zaux(:,1),&
                        maxtime=Mopt%R0,Nseason=Mopt%N5,pattern=mu,feas=feas,err=err,mess=mess)
    if(.not.feas) return
    ! get Standard deviation
    call PCs_GetPattern(k=k,DoWhat=(/Mopt%N3,Mopt%N4/),x=par,time=Mdata%rep(j)%Zaux(:,1),&
                        maxtime=Mopt%R0,Nseason=Mopt%N5,pattern=sigma,feas=feas,err=err,mess=mess)
    if(.not.feas) return
    sigma=exp(sigma)
    ! Get residuals and innovations
    ! compute log-lkh
    res=(Mdata%rep(j)%Z(:,Mopt%N0)-mu)/sigma
    if(Nlag>0) then
        if(Mopt%N6==0) then
        ! -----------------
            ! Evin et al's approach - only works for AR(1)
            ro=par(Npar)
            if(abs(ro)>=1._mrk-ro_epsilon) then;feas=.false.;return;endif
            fx0 = -1._mrk*sum(log(sigma))& ! Jacobian
            -0.5_mrk*log(twoPi)-0.5_mrk*res(1)**2& ! N(0,1) marginal for first step
            -0.5_mrk*real(Nt-1,mrk)*(log(twoPi)+log(1._mrk-ro**2))-0.5_mrk*sum( ((res(2:Nt)-ro*res(1:Nt-1))**2)/(1._mrk-ro**2) ) ! conditional: N(ro*res(t-1),sqrt(1-ro**2)
        else
        ! -----------------
            ! Attempt at AR(p) - buggy...
            innov=res( (Nlag+1):Nt)
            do i=1,Nlag
                if(abs(par(Npar-Nlag+i))>=1._mrk-ro_epsilon) then;feas=.false.;return;endif
                innov=innov-par(Npar-Nlag+i)*res( (Nlag+1-i):(Nt-i) )
            enddo
            fx0=-1._mrk*sum(log(sigma((Nlag+1):Nt)))& ! Jacobian
            -0.5_mrk*real(Nt-Nlag,mrk)*log(twoPi)-0.5_mrk*sum( innov**2 )
        endif
    else
        fx0 = -1._mrk*sum(log(sigma))& ! Jacobian
        -0.5_mrk*real(Nt-1,mrk)*(log(twoPi))-0.5_mrk*sum( res**2 ) ! conditional: N(ro*res(t-1),sqrt(1-ro**2)
    endif
    fx=fx+fx0
    deallocate(mu,sigma,res,innov,par)
enddo
if(any(Mdata%DirtyRep)) then ! add log-prior on obs error
    !call GetSamplePdf(DistId='Gaussian',x=x( (Npar+1):),par=(/0._mrk,DirtySdev/),&
    !                  loga=.true.,pdf=fx0,feas=feas,isnull=isnull,err=err,mess=mess)
    do i=1,Npar
        call GetPdf(DistId='Gaussian',x=x(Npar+i),par=(/0._mrk,DirtySdev*Mopt%V0(i)/),&
                          loga=.true.,pdf=fx0,feas=feas,isnull=isnull,err=err,mess=mess)
        if(err/=0) then;mess='LLfunk_PCs:'//trim(mess);return;endif
        if(.not.feas) return
        fx=fx+fx0
    enddo
endif
end subroutine LLfunk_PCs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ObjFunk_PCs(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
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
  
  err=0;message='';feas=.true.
  if(.not.present(fx)) return
  call LLfunk_PCs(x=x,feas=feas,isnull=isnull,fx=fx,err=err,mess=message)
  if(err>0) return
  if(present(fx)) then
    if(isnull) fx=-HugeRe
    fx=-1*fx
endif
endsubroutine ObjFunk_PCs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PCs_GetPattern(k,DoWhat,x,time,maxtime,Nseason,pattern,feas,err,mess)
! Get temporal pattern using trend & seasonnality models
use utilities_dmsl_kit,only:twoPi
integer(mik), intent(inout)::k
integer(mik),intent(in)::DoWhat(:),Nseason
real(mrk),intent(in)::x(:),time(:),maxtime
real(mrk),intent(out)::pattern(:)
logical,intent(out)::feas
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk),parameter::eps=0.1_mrk
integer(mik)::i,t,Nt
real(mrk)::Leg(size(time),DoWhat(1))
real(mrk)::ones(size(time)),a0,a1,a2

err=0;mess='';feas=.true.
Nt=size(time)
! Legendre polynomials
if(DoWhat(1)>0) then
    do i=1,DoWhat(1);Leg(:,i)=LegendrePoly(2._mrk*(time/maxtime)-1._mrk,i);enddo
    pattern=matmul(Leg,x( (k+1):(k+DoWhat(1)) ))
else
    pattern=0._mrk
endif
k=k+DoWhat(1)
! Annual cycle
if(DoWhat(2)>0) then
!    a0=x(k+1);a1=x(k+2);a2=1._mrk
!    if(a0<0._mrk .or. a1<-0.5_mrk*a2 .or. a1>0.5_mrk*a2) then;feas=.false.;return;endif
!!    if(a0<0._mrk .or. a1<0._mrk .or. a1>a2) then;feas=.false.;return;endif
!    pattern=pattern+a0*cos((twoPi/a2)*(time-a1))
!    k=k+2
    do t=1,Nt
        pattern(t)=pattern(t)+SeasonnalSignal(time(t),Nseason,x( (k+1):(k+Nseason) ))
    enddo
    k=k+Nseason
endif
!! sub-Annual cycle
!if(DoWhat(3)>0) then
!    a0=x(k+1);a1=x(k+2);a2=x(k+3)
!    if(a0<0._mrk .or. a1<-0.5_mrk*a2 .or. a1>0.5_mrk*a2 .or. a2>=1._mrk-eps) then;feas=.false.;return;endif
!!    if(a0<0._mrk .or. a1<0_mrk .or. a1>a2 .or. a2>=1._mrk-eps) then;feas=.false.;return;endif
!    pattern=pattern+a0*cos((twoPi/a2)*(time-a1))
!    k=k+3
!endif
end subroutine PCs_GetPattern
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PCs_GetStartingPoint(AnyDirty,x0,xscale,xLo,xHi)
! quite crappy, may certainly be improved...
use linalg_dmsl_kit, only:choles_invrt
use EmpiricalStats_tools, only:GetEmpiricalStats
real(mrk),intent(out)::x0(:),xscale(:),xLo(:),xHi(:)
logical(mrk), intent(in)::AnyDirty
!locals
integer(mik)::k,i,err
logical,allocatable::mask(:)
integer(mik)::Nt,Nleg,Nseason,Nlag,Npar,WhichComp
real(mrk), allocatable::time(:), B(:,:),Binv(:,:),res(:),stdres(:),foo(:)
real(mrk)::moy,sdev,maxtime
logical::pd
character(250)::mess

k=0;xscale=1._mrk;x0=0._mrk
Nt=Mdata%rep(1)%Nt
WhichComp=Mopt%N0
Nleg=Mopt%N1
Nseason=Mopt%N2
Npar=Mopt%Npar
Nlag=Mopt%N7
maxtime=Mopt%R0
allocate(mask(Nt),time(Nt),B(Nt,Nseason+Nleg),Binv(Nseason+Nleg,Nseason+Nleg),res(Nt),stdres(Nt))
time=Mdata%rep(1)%Zaux(:,1)
moy=sum(Mdata%rep(1)%Z(:,WhichComp))/real(Nt,mrk)
sdev=sqrt(sum( (Mdata%rep(1)%Z(:,WhichComp)-moy)**2 )/real(Nt,mrk))
! x0 and xscale
if(DoSmartStart) then
    ! Start by estimating S & L for the mean using the standard OLS estimate
    B=0._mrk
    do i=1,Nleg
        B(:,i)=LegendrePoly(2._mrk*(time/maxtime)-1._mrk,i)
    enddo
    do i=1,Nseason
        mask=(time-floor(time))<real(i,mrk)/real(Nseason,mrk) .and. (time-floor(time))>=real(i-1,mrk)/real(Nseason,mrk)
        where (mask) B(:,i+Nleg)=1._mrk
    enddo
    call choles_invrt(a=matmul(transpose(B),B),ainv=Binv,posDefinite=pd,err=err,message=mess)
    if(err/=0 .or. .not.pd) return
    x0( (k+1):(k+Nleg+Nseason) )=matmul(matmul(Binv,transpose(B)),Mdata%rep(1)%Z(:,WhichComp))
    ! Get the residuals and estimate seasonnal sdevs
    res=Mdata%rep(1)%Z(:,WhichComp)-matmul(B,x0( (k+1):(k+Nleg+Nseason) ))
    k=k+Nleg+Nseason
    do i=1,Nseason
        mask=(time-floor(time))<real(i,mrk)/real(Nseason,mrk) .and. (time-floor(time))>=real(i-1,mrk)/real(Nseason,mrk)
        if(allocated(foo)) deallocate(foo);allocate(foo(count(mask)))
        foo=pack(res,mask)
        call GetEmpiricalStats(x=foo, std=x0(k+1),err=err,mess=mess)
        where(mask) stdres=res/x0(k+1)
        x0(k+1)=log(x0(k+1))
        k=k+1
    enddo
    ! Get the lag-1 coeff of standardized residuals
    if(Nlag>0) then
        moy=sum(stdres)/real(Nt,mrk)
        x0(k+1)=dot_product(stdres(1:(Nt-1))-moy,stdres(2:Nt)-moy)/dot_product(stdres-moy,stdres-moy)
        if(Nlag>1) x0((k+2):)=0._mrk
    endif
    xscale(1:Npar)=10._mrk*Mopt%V0
else
    if(Nleg>0) then;x0((k+1):(k+Nleg))=0._mrk;k=k+Nleg;endif
    do i=1,Nseason
        mask=(time-floor(time))<real(i,mrk)/real(Nseason,mrk) .and. (time-floor(time))>=real(i-1,mrk)/real(Nseason,mrk)
        x0(k+1)=sum(Mdata%rep(1)%Z(:,WhichComp),mask)/real(count(mask),mrk);
        xscale(k+1)=sdev;k=k+1
    enddo
    do i=1,Nseason
        mask=(time-floor(time))<real(i,mrk)/real(Nseason,mrk) .and. (time-floor(time))>=real(i-1,mrk)/real(Nseason,mrk)
        x0(k+1)=sqrt(sum( (Mdata%rep(1)%Z(:,WhichComp)-x0(k+1-Nseason))**2,mask )/real(count(mask),mrk))
        x0(k+1)=log(x0(k+1));
        xscale(k+1)=abs(x0(k+1))
        k=k+1
    enddo
    if(Nlag>0) then;x0( (Npar-Nlag+1):Npar)=0.1_mrk;xscale( (Npar-Nlag+1):Npar)=0.1_mrk;endif
endif
deallocate(mask,time,B,Binv,res,stdres)
if(AnyDirty) then
    x0((Npar+1):)=0._mrk
    xscale((Npar+1):)=DirtySdev*Mopt%V0
endif

! xLo and xHi
xLo=-HugeRe;xHi=HugeRe
if(Nlag>0) then;xLo( (Npar-Nlag+1):Npar)=-1._mrk+ro_epsilon;xHi( (Npar-Nlag+1):Npar)=1._mrk-ro_epsilon;endif
if(AnyDirty) then
    xLo((Npar+1):)=-HugeRe
    xHi((Npar+1):)=HugeRe
endif
end subroutine PCs_GetStartingPoint
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module MLfunk_catalogue
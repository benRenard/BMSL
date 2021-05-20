module MultiModel_tools

!~**********************************************************************
!~* Purpose: tools for multi-model analysis (weighting etc)
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
use MultiModel_MLfunk_library

implicit none
Private
public :: GetGCMWeight

type,public:: MultiMdataType
    integer(mik)::Nm=undefIN ! Note: model Nm+1 for observation
    type(SingleMdataType), allocatable::model(:) ! Warning: size Nm+1!!! 
end type MultiMdataType

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetGCMWeight(X,Sigma,infer,option,AICweight,BICweight,relweight,&
                        H,aic,bic,par,eps,fmax,pk,N,err,mess)
!^**********************************************************************
!^* Purpose: Get weights according the method developped in COMPLEX
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
!^*		1.X, multimodel data object
!^*		2.Sigma, observation errors covariance matrix
!^*		3.infer, inference setup
!^*		4.option, options passed to MLfunk
!^* OUT
!^*		1.AICweight (size Nm+1) AIC weights
!^*		2.BICweight (size Nm+1) BIC weights
!^*		3.relweight (size Nm) relative weights (H0 excluded)
!^*		4.H (size Nm+1) likelihood associated with Hk
!^*		5.aic (size Nm+1)
!^*		6.bic (size Nm+1)
!^*		7.par (size Npar*Nm+1) Max-lkh estimates, par(:,Nm+1) is for observations
!^*		8.eps (size Npar*Nm+1) Max-lkh estimates of epsilon
!^*		9.fmax (size Nm+1) Max-lkh value, fmax(Nm+1) is for observations
!^*		10.pk (size Nm+1) total number of parameters for each hypothesis
!^*		11.N (size Nm+1) total number of samples for each hypothesis
!^*		12.err, error code
!^*		13.mess, error message
!^**********************************************************************
use linalg_dmsl_kit,only:choles_invrt
use types_dmsl_kit, only:data_ricz_type
use utilities_dmsl_kit,only:number_string
type(MultiMdataType), intent(in)::X
real(mrk), intent(in)::Sigma(:,:)
type(InferType), intent(in)::infer
Type(data_ricz_type), intent(in)::option
real(mrk), intent(out)::AICweight(X%Nm+1),BICweight(X%Nm+1),relweight(X%Nm),&
                        AIC(X%Nm+1),BIC(X%Nm+1),&
                        fmax(X%Nm+1),&
                        par(infer%Npar,X%Nm+1),eps(infer%Npar,X%Nm),H(X%Nm+1)
integer(mik), intent(out)::pk(X%Nm+1),N(X%Nm+1),err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='GetGCMWeight'
integer(mik)::Nm,Npar,m,Nhat(X%Nm+1),Nt,Nd,Nt2,Nd2,r,i
real(mrk)::w(X%Nm+1),rw(X%Nm),Hini(X%Nm)
real(mrk)::dummypar(infer%Npar),dummyeps(infer%Npar)
logical::posdef

err=0;mess='';
AICweight=undefRN;BICweight=undefRN;relweight=undefRN;
aic=undefRN;bic=undefRN;fmax=undefRN
par=undefRN;H=undefRN
! Pilot global variables
Goption=option
Ginfer=infer
Nm=X%Nm
Npar=infer%Npar
if(allocated(GSigmaInv)) deallocate(GSigmaInv);allocate(GSigmaInv(Npar,Npar))
posdef=.true.
call choles_invrt(a=Sigma,ainv=GSigmaInv,posDefinite=posdef,logDet=GSigmaDet,err=err,message=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not.posdef) then;err=1;mess=trim(procname)//':Fatal:Sigma is not positive definite';return;endif

do m=Nm+1,1,-1
    ! Get Mak-lkh estimates for each individual model
    !write(*,*) '*************************************************'
    !write(*,*) 'Processing model:',m,'/',Nm,'+1'
    ! Pilot global object GX describing data 
    call PilotGX(X=X,m=m,addobs=.false.)
    ! perform estimation
    call MLfunk(par=par(:,m),eps=dummyeps,fmax=fmax(m),n=Nhat(m),err=err,mess=mess)
    if(err>0) then
        mess=trim(procname)//':'//trim(mess)
        write(*,*) 'WARNING: optim failed with message:'
        write(*,*) trim(mess)
        !return
    endif
    ! Hypothesis Hk: model k is correct, re-estimate with observations
    if(m<=Nm) then
        call PilotGX(X=X,m=m,addobs=.true.)
        ! perform estimation
        call MLfunk(par=dummypar,eps=eps(:,m),fmax=Hini(m),n=N(m),err=err,mess=mess)
        if(err>0) then
            mess=trim(procname)//':'//trim(mess)
            write(*,*) 'WARNING: optim failed with message:'
            write(*,*) trim(mess)
            !return
        endif
    endif
enddo
! Correct H & M
do m=1,Nm
    H(m)=Hini(m)+sum(fmax)-fmax(m)-fmax(Nm+1)
    N(m)=N(m)+sum(Nhat)-Nhat(m)-Nhat(Nm+1)
enddo
! H0 ypothesis: all models are wrong
H(Nm+1)=sum(fmax);N(Nm+1)=sum(Nhat)
! Compute criteria
pk=Npar*(Nm+1);pk(Nm+1)=Npar*(Nm+2)
aic=-2._mrk*H+2._mrk*pk
bic=-2._mrk*H+pk*log(1._mrk*N)
w=aic-minval(aic);AICweight=exp(-0.5_mrk*w)/sum(exp(-0.5_mrk*w))
w=bic-minval(bic);BICweight=exp(-0.5_mrk*w)/sum(exp(-0.5_mrk*w))
rw=bic(1:Nm)-minval(bic(1:Nm));relweight=exp(-0.5_mrk*rw)/sum(exp(-0.5_mrk*rw))

end subroutine GetGCMWeight
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!==============!
! PRIVATE SUBS !
!==============!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PilotGX(X,m,addobs)
!^**********************************************************************
!^* Purpose: Pilot global object GX
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
!^*		1.X, multimodel data object
!^*		2.m, model index
!^*		3.addobs, Add obs as the last replication?
!^* OUT
!^**********************************************************************
type(MultiMdataType), intent(in)::X
integer(mik), intent(in)::m
logical, intent(in)::addobs
! locals
integer(mik)::upto,r,Nt,Nd,Nt2,Nd2,Nm

Nm=X%Nm
if(addobs) then
    GX%nrep=X%model(m)%Nrep+1
else
    GX%nrep=X%model(m)%Nrep
endif
if(allocated(GX%rep)) deallocate(GX%rep);allocate(GX%rep(GX%nrep))
if(allocated(GX%DirtyRep)) deallocate(GX%DirtyRep);allocate(GX%DirtyRep(GX%nrep))
if(addobs) then
    GX%DirtyRep=.false.
    upto=GX%nrep-1
else
    GX%DirtyRep=m>X%Nm
    upto=GX%nrep
endif
do r=1,upto
    call PilotGX_rep(X,m,r)
enddo

if(addobs) then
    r=GX%nrep
    Nt=size(X%model(Nm+1)%rep(1)%Z,dim=1)
    Nd=size(X%model(Nm+1)%rep(1)%Z,dim=2)
    GX%rep(r)%Nt=Nt;GX%rep(r)%Nd=Nd
    if(allocated(GX%rep(r)%Z)) deallocate(GX%rep(r)%Z)
    allocate(GX%rep(r)%Z(Nt,Nd))
    GX%rep(r)%Z=X%Model(Nm+1)%rep(1)%Z
    Nt2=size(X%model(Nm+1)%rep(1)%Zaux,dim=1)
    Nd2=size(X%model(Nm+1)%rep(1)%Zaux,dim=2)
    if(allocated(GX%rep(r)%Zaux)) deallocate(GX%rep(r)%Zaux)
    allocate(GX%rep(r)%Zaux(Nt2,Nd2))
    GX%rep(r)%Zaux=X%Model(Nm+1)%rep(1)%Zaux
    GX%DirtyRep(r)=.true.
endif
end subroutine PilotGX
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PilotGX_rep(X,m,r)
!^**********************************************************************
!^* Purpose: Pilot a single rep in global object GX
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
!^*		1.X, multimodel data object
!^*		2.m, model index
!^*		3.r, rep index
!^* OUT
!^**********************************************************************
type(MultiMdataType), intent(in)::X
integer(mik), intent(in)::m,r
!locals
integer(mik)::Nt,Nd,Nt2,Nd2

Nt=size(X%model(m)%rep(r)%Z,dim=1)
Nd=size(X%model(m)%rep(r)%Z,dim=2)
GX%rep(r)%Nt=Nt;GX%rep(r)%Nd=Nd
if(allocated(GX%rep(r)%Z)) deallocate(GX%rep(r)%Z)
allocate(GX%rep(r)%Z(Nt,Nd))
GX%rep(r)%Z=X%Model(m)%rep(r)%Z
Nt2=size(X%model(m)%rep(r)%Zaux,dim=1)
Nd2=size(X%model(m)%rep(r)%Zaux,dim=2)
if(allocated(GX%rep(r)%Zaux)) deallocate(GX%rep(r)%Zaux)
allocate(GX%rep(r)%Zaux(Nt2,Nd2))
GX%rep(r)%Zaux=X%Model(m)%rep(r)%Zaux

end subroutine PilotGX_rep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module MultiModel_tools
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
use MLfunk_catalogue

implicit none
Private
public :: GetGCMWeight

type,public:: MultiMdataType
    integer(mik)::Nm=undefIN ! Note: model Nm+1 for observation
    type(SingleMdataType), allocatable::model(:) ! Warning: size Nm+1!!! 
end type MultiMdataType

logical, parameter::SaveEpsilon=.true.
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetGCMWeight(X,Moptions,ObsStdError,MLfunk,weight,relweight,H,bic,par,fmax,err,mess)
!^**********************************************************************
!^* Purpose: Get weights according the method developped in COMPLEX
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
!^*		1.X, multimodel data object
!^*		2.MLfunk, Max-lkh function
!^* OUT
!^*		1.weight (size Nm+1) weights
!^*		2.relweight (size Nm) relative weights (H0 excluded)
!^*		3.H (size Nm+1) likelihood associated with Hk
!^*		4.bic (size Nm+1)
!^*		5.par (size Npar*Nm+1) Max-lkh estimates, par(:,Nm+1) is for observations
!^*		6.fmax (size Nm+1) Max-lkh value, fmax(Nm+1) is for observations
!^*		7.err, error code
!^*		8.mess, error message
!^**********************************************************************
use MLfunk_catalogue
use utilities_dmsl_kit,only:number_string
type(MultiMdataType), intent(in)::X
Type(MoptionsType), intent(in)::Moptions
real(mrk), intent(in)::ObsStdError
real(mrk), intent(out)::weight(X%Nm+1),relweight(X%Nm),bic(X%Nm+1),&
                        fmax(X%Nm+1),par(Moptions%Npar,X%Nm+1),H(X%Nm+1)
integer(mik), intent(out)::err
character(*),intent(out)::mess
interface
    subroutine MLfunk(par,fmax,n,err,mess)
        use kinds_dmsl_kit
        use MLfunk_catalogue
        real(mrk),intent(out)::par(:),fmax
        integer(mik),intent(out)::n,err
        character(*),intent(out)::mess
    end subroutine MLfunk
end interface
! locals
character(250),parameter::procname='GetGCMWeight'
integer(mik)::Nm,Npar,m,Nhat(X%Nm+1),N(X%Nm+1),Nt,Nd,Nt2,Nd2,r,i
real(mrk)::pk(X%Nm+1),w(X%Nm+1),rw(X%Nm),Hini(X%Nm)
real(mrk), allocatable::dummypar(:)

err=0;mess='';weight=undefRN;relweight=undefRN;bic=undefRN;fmax=undefRN
par=undefRN;H=undefRN
Mopt=Moptions
Nm=X%Nm;Npar=Mopt%Npar
DirtySdev=ObsStdError

do m=Nm+1,1,-1
    ! Get Mak-lkh estimates for each individual model
    write(*,*) '*************************************************'
    write(*,*) 'Processing model:',m,'/',Nm,'+1'
    ! Pilot global object Mdata
    Mdata%nrep=X%model(m)%Nrep
    if(allocated(Mdata%rep)) deallocate(Mdata%rep);allocate(Mdata%rep(Mdata%nrep))
    if(allocated(Mdata%DirtyRep)) deallocate(Mdata%DirtyRep);allocate(Mdata%DirtyRep(Mdata%nrep))
    if(m<=Nm) then
        Mdata%DirtyRep=.false.
        if(allocated(dummypar)) deallocate(dummypar);allocate(dummypar(Moptions%Npar))
    else
        Mdata%DirtyRep=.true.
        if(allocated(dummypar)) deallocate(dummypar);allocate(dummypar(2*Moptions%Npar))
    endif
    do r=1,Mdata%nrep
        Nt=size(X%model(m)%rep(r)%Z,dim=1)
        Nd=size(X%model(m)%rep(r)%Z,dim=2)
        Mdata%rep(r)%Nt=Nt;Mdata%rep(r)%Nd=Nd
        if(allocated(Mdata%rep(r)%Z)) deallocate(Mdata%rep(r)%Z)
        allocate(Mdata%rep(r)%Z(Nt,Nd))
        Mdata%rep(r)%Z=X%Model(m)%rep(r)%Z
        Nt2=size(X%model(m)%rep(r)%Zaux,dim=1)
        Nd2=size(X%model(m)%rep(r)%Zaux,dim=2)
        if(allocated(Mdata%rep(r)%Zaux)) deallocate(Mdata%rep(r)%Zaux)
        allocate(Mdata%rep(r)%Zaux(Nt2,Nd2))
        Mdata%rep(r)%Zaux=X%Model(m)%rep(r)%Zaux
    enddo
    ! perform estimation
    call MLfunk(par=dummypar,fmax=fmax(m),n=Nhat(m),err=err,mess=mess)
    par(:,m)=dummypar(1:Moptions%Npar)
    if(err>0) then
        mess=trim(procname)//':'//trim(mess)
        write(*,*) 'WARNING: optim failed with message:'
        write(*,*) trim(mess)
        !return
    endif
    ! Hypothesis Hk: model k is correct, re-estimate with observations
    if(allocated(dummypar)) deallocate(dummypar);allocate(dummypar(2*Moptions%Npar))
    if(m<=Nm) then
        Mdata%nrep=X%model(m)%Nrep+1
        if(allocated(Mdata%rep)) deallocate(Mdata%rep);allocate(Mdata%rep(Mdata%nrep))
        if(allocated(Mdata%DirtyRep)) deallocate(Mdata%DirtyRep);allocate(Mdata%DirtyRep(Mdata%nrep))
        Mdata%DirtyRep=.false.
        do r=1,Mdata%nrep-1
            Nt=size(X%model(m)%rep(r)%Z,dim=1)
            Nd=size(X%model(m)%rep(r)%Z,dim=2)
            Mdata%rep(r)%Nt=Nt;Mdata%rep(r)%Nd=Nd
            if(allocated(Mdata%rep(r)%Z)) deallocate(Mdata%rep(r)%Z)
            allocate(Mdata%rep(r)%Z(Nt,Nd))
            Mdata%rep(r)%Z=X%Model(m)%rep(r)%Z
            Nt2=size(X%model(m)%rep(r)%Zaux,dim=1)
            Nd2=size(X%model(m)%rep(r)%Zaux,dim=2)
            if(allocated(Mdata%rep(r)%Zaux)) deallocate(Mdata%rep(r)%Zaux)
            allocate(Mdata%rep(r)%Zaux(Nt2,Nd2))
            Mdata%rep(r)%Zaux=X%Model(m)%rep(r)%Zaux
        enddo
        ! Add obs as just another rep
        r=Mdata%nrep
        Nt=size(X%model(Nm+1)%rep(1)%Z,dim=1)
        Nd=size(X%model(Nm+1)%rep(1)%Z,dim=2)
        Mdata%rep(r)%Nt=Nt;Mdata%rep(r)%Nd=Nd
        if(allocated(Mdata%rep(r)%Z)) deallocate(Mdata%rep(r)%Z)
        allocate(Mdata%rep(r)%Z(Nt,Nd))
        Mdata%rep(r)%Z=X%Model(Nm+1)%rep(1)%Z
        Nt2=size(X%model(Nm+1)%rep(1)%Zaux,dim=1)
        Nd2=size(X%model(Nm+1)%rep(1)%Zaux,dim=2)
        if(allocated(Mdata%rep(r)%Zaux)) deallocate(Mdata%rep(r)%Zaux)
        allocate(Mdata%rep(r)%Zaux(Nt2,Nd2))
        Mdata%rep(r)%Zaux=X%Model(Nm+1)%rep(1)%Zaux
        Mdata%DirtyRep(r)=.true.
        ! perform estimation
        call MLfunk(par=dummypar,fmax=Hini(m),n=N(m),err=err,mess=mess)
        if(err>0) then
            mess=trim(procname)//':'//trim(mess)
            write(*,*) 'WARNING: optim failed with message:'
            write(*,*) trim(mess)
            !return
        endif
        if(SaveEpsilon) then
            open(unit=666,file='Epsilon_M'//trim(number_string(m))//'.txt',status='replace')
            do i=1,Moptions%Npar
                write(666,'(2e20.10)') dummypar(Moptions%Npar+i),DirtySdev*Mopt%V0(i)
            enddo
            close(666)
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
bic=-2._mrk*H+pk*log(1._mrk*N)
w=bic-minval(bic);weight=exp(-0.5_mrk*w)/sum(exp(-0.5_mrk*w))
rw=bic(1:Nm)-minval(bic(1:Nm));relweight=exp(-0.5_mrk*rw)/sum(exp(-0.5_mrk*rw))

end subroutine GetGCMWeight
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module MultiModel_tools
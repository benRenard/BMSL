module MCMCStrategy_tools

!~**********************************************************************
!~* Purpose: MCMC strategies
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref
!~**********************************************************************
!~* Last modified:
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1.
!~*    2.
!~*    3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use MCMC_tools


implicit none
Private
public :: Adaptive_Metro_OAAT, Adaptive_Metro

Contains


subroutine Adaptive_Metro(f,x,fx,fAux,covar,scale,nAdapt,nCycles,&
                MinMoveRate,MaxMoveRate,DownMult,UpMult,&
                OutFile,headers,err,mess)

!^**********************************************************************
!^* Purpose: Metropolis with adaptive scale factor
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:09/11/2008
!^**********************************************************************
!^* Comments: Only the scale factor is adapted - var is choleskized
!^* out of this sub
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.
!^*    2.
!^*    3.
!^* OUT
!^*    1.
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^* INOUT
!^*    1.
!^*    2.
!^*    3.
!^**********************************************************************
use utilities_dmsl_kit, only:number_string
!Use RFortran

real(mrk),intent(inout)::x(:),covar(:,:),scale
integer(mik),intent(in)::nAdapt,nCycles
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in),optional::OutFile,headers(:)
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
logical::feas,isnull,again
character(100), allocatable::head(:), fAuxHead(:)
logical::move
real(mrk)::moverate
integer(mik)::n,i,p,k,compt
real(mrk)::postlist(nAdapt)
character(250)::fmt

interface
    subroutine f(x,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface

!Handle Output file headers
n=size(x,dim=1)
if(present(OutFile)) then
    open(unit=1,file=OutFile,status='replace')
    if(present(headers)) then ! trust user - no size check at the moment
        k=size(headers)
        fmt='(' // trim(number_string(k)) // '(A13," "))'
        write(1,trim(fmt)) headers
    else
        if(allocated(head)) deallocate(head)
        allocate(head(n))
        do i=1,n
            head(i)='par'//trim(number_string(i))
        enddo
        if(present(fAux)) then
            p=size(fAux,dim=1)
            if(allocated(fAuxHead)) deallocate(fAuxHead)
            allocate(fAuxHead(p))
            do i=1,p
                fAuxhead(i)='fAux'//trim(number_string(i))
            enddo
            fmt='(' // trim(number_string(n+1+p)) // '(A13," "))'
            write(1,trim(fmt)) head(:),'Objectivefunction',fAuxHead(:)
            deallocate(fAuxHead)
        else
            fmt='(' // trim(number_string(n+1)) // '(A13," "))'
            write(1,trim(fmt)) head(:),'Objectivefunction'
        endif
        deallocate(head)
    endif
endif

! Check starting value
call f(x=x,feas=feas,isnull=isnull,fx=fx,fAux=fAux,err=err,mess=mess)
If(err>0) then
    mess='Adaptive_Metro'//trim(mess)
    if(present(OutFile)) close(1)
    return
endif
if((.not.feas).or.isnull) then
    err=1;mess='Adaptive_Metro:Fatal:Unfeasible starting point'
    if(present(OutFile)) close(1)
    return
endif
if(present(OutFile)) then
    if(present(fAux)) then
        fmt='(' // trim(number_string(n+1+p)) // 'e14.6)'
        write (1,trim(fmt)) x,fx,fAux
    else
        fmt='(' // trim(number_string(n+1)) // 'e14.6)'
        write (1,trim(fmt)) x,fx
    endif
endif

!ok = Rinit();ok=Reval('X11()')

again=.true.;compt=0
do while(again)
    !Start sampling
    moverate=0._mrk
    do i=1,nAdapt
        call Metro_GaussJump(f,x,fx,fAux,(scale**2)*covar,move,err,mess)
        If(err>0) then
            mess='Adaptive_Metro'//trim(mess)
            if(present(OutFile)) close(1)
            return
        endif
        postlist(i)=fx
        !Write2File
        if(present(OutFile)) then
            if(present(fAux)) then
                fmt='(' // trim(number_string(n+1+p)) // 'e14.6)'
                write (1,trim(fmt)) x,fx,fAux
            else
                fmt='(' // trim(number_string(n+1)) // 'e14.6)'
                write (1,trim(fmt)) x,fx
            endif
        endif
        !Update Move Rates
        if(move) moverate=moverate+1._mrk/nAdapt
    enddo
!    ok=Reval('layout(c(1, 2))')
!    ok=Rput('mr', moverate)
!    ok=Rput('compt', compt+1)
!    ok=Reval('mrl[compt]<-mr')
!    ok=Reval('matplot(mrl,type="l", xlab="Parameter", ylab="Move Rate", main="Move Rates", ylim=c(0, 1))')
!    ok=Rput('post', postlist)
!    ok=Reval('matplot(post,type="l", xlab="Iteration", ylab="Log-posterior", main="Log-posterior")')

    ! Adapt Jump standard deviation
    if(moverate<=MinMoveRate) scale=scale*DownMult
    if(moverate>=MaxMoveRate) scale=scale*UpMult
    compt=compt+1
    if(nCycles>=1) then
        if(compt>=nCycles) again=.false.
    else
        write(*,*) 'Again?'
        read(*, *) again
    endif
enddo

if(present(OutFile)) close(1)
!ok = Rclose()  ! End R session

end subroutine Adaptive_Metro

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Adaptive_Metro_OAAT(f,x,fx,fAux,std,nAdapt,nCycles,&
                MinMoveRate,MaxMoveRate,DownMult,UpMult,&
                OutFile,headers,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.
!^*    2.
!^*    3.
!^* OUT
!^*    1.
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^* INOUT
!^*    1.
!^*    2.
!^*    3.
!^**********************************************************************
use utilities_dmsl_kit, only:number_string
!Use RFortran

real(mrk),intent(inout)::x(:),std(:)
integer(mik),intent(in)::nAdapt,nCycles
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in),optional::OutFile,headers(:)
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
logical::feas,isnull,again
character(100), allocatable::head(:), fAuxHead(:)
logical,allocatable::move(:)
real(mrk),allocatable::moverate(:)
integer(mik)::n,i,p,k,compt
real(mrk)::postlist(nAdapt)
character(250)::fmt

interface
    subroutine f(x,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface

!Handle Output file headers
n=size(x,dim=1)
if(present(fAux)) then;p=size(fAux,dim=1);else;p=0;endif

if(present(OutFile)) then
    open(unit=1,file=OutFile,status='replace')
    if(present(headers)) then ! trust user - no size check at the moment
        k=size(headers)
        fmt='(' // trim(number_string(k)) // '(A13," "))'
        write(1,trim(fmt)) headers
    else
        if(allocated(head)) deallocate(head)
        allocate(head(n))
        do i=1,n
            head(i)='par'//trim(number_string(i))
        enddo
        if(present(fAux)) then
            if(allocated(fAuxHead)) deallocate(fAuxHead)
            allocate(fAuxHead(p))
            do i=1,p
                fAuxhead(i)='fAux'//trim(number_string(i))
            enddo
            fmt='(' // trim(number_string(n+1+p)) // '(A13," "))'
            write(1,trim(fmt)) head(:),'Objectivefunction',fAuxHead(:)
            deallocate(fAuxHead)
        else
            fmt='(' // trim(number_string(n+1)) // '(A13," "))'
            write(1,trim(fmt)) head(:),'Objectivefunction'
        endif
        deallocate(head)
    endif
endif

! Check starting value
call f(x=x,feas=feas,isnull=isnull,fx=fx,fAux=fAux,err=err,mess=mess)
If(err>0) then
    mess='Adaptive_Metro_OAAT'//trim(mess)
    if(present(OutFile)) close(1)
    return
endif
if((.not.feas).or.isnull) then
    err=1;mess='Adaptive_Metro_OAAT:Fatal:Unfeasible starting point'
    if(present(OutFile)) close(1)
    return
endif
if(present(OutFile)) then
    if(present(fAux)) then
        fmt='(' // trim(number_string(n+1+p)) // 'e14.6)'
        write (1,trim(fmt)) x,fx,fAux
    else
        fmt='(' // trim(number_string(n+1)) // 'e14.6)'
        write (1,trim(fmt)) x,fx
    endif
endif
! Allocate MoveRate stuff
if(allocated(move)) deallocate(move)
allocate(move(n))
if(allocated(moverate)) deallocate(moverate)
allocate(moverate(n))

!ok = Rinit();ok=Reval('X11()')

again=.true.;compt=0;!ok=Reval("mrl<-c()")
do while(again)
    !Start sampling
    moverate=0._mrk
    do i=1,nAdapt
        call Metro_OAAT_GaussJump(f,x,fx,fAux,std,move,err,mess)
        If(err>0) then
            mess='Adaptive_Metro_OAAT'//trim(mess)
            if(present(OutFile)) close(1)
            return
        endif
        postlist(i)=fx
        !Write2File
        if(present(OutFile)) then
            if(present(fAux)) then
                fmt='(' // trim(number_string(n+1+p)) // 'e14.6)'
                write (1,trim(fmt)) x,fx,fAux
            else
                fmt='(' // trim(number_string(n+1)) // 'e14.6)'
                write (1,trim(fmt)) x,fx
            endif
        endif
        !Update Move Rates
        where(move) moverate=moverate+1._mrk/nAdapt
    enddo
!    ok=Reval('layout(c(1, 2))')
!    ok=Rput('mr', moverate)
!    ok=Reval('matplot(mr,type="l", xlab="Parameter", ylab="Move Rate", main="Move Rates", ylim=c(0, 1))')
!    ok=Rput('post', postlist)
!    ok=Reval('matplot(post,type="l", xlab="Iteration", ylab="Log-posterior", main="Log-posterior")')

    ! Adapt Jump standard deviation
    where(moverate<=MinMoveRate) std=std*DownMult
    where(moverate>=MaxMoveRate) std=std*UpMult
    compt=compt+1
    if(nCycles>=1) then
        if(compt>=nCycles) again=.false.
    else
        write(*,*) 'Again?'
        read(*, *) again
    endif
enddo

if(present(OutFile)) close(1)
!ok = Rclose()  ! End R session

end subroutine Adaptive_Metro_OAAT


end module MCMCStrategy_tools

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
                OutFile,MonitorFile,headers,err,mess)

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
use utilities_dmsl_kit, only:number_string,getSpareUnit

real(mrk),intent(inout)::x(:),covar(:,:),scale
integer(mik),intent(in)::nAdapt,nCycles
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in),optional::OutFile,MonitorFile,headers(:)
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
logical::feas,isnull,again
character(100), allocatable::head(:), fAuxHead(:)
logical::move
real(mrk)::moverate
integer(mik)::n,i,p,k,compt,untOut,untMonitor
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
    call getSpareUnit(untOut,err,mess)
    If(err>0) then;mess='Adaptive_Metro:'//trim(mess);return;endif
    open(unit=untOut,file=OutFile,status='replace')
    if(present(headers)) then ! trust user - no size check at the moment
        k=size(headers)
        fmt='(' // trim(number_string(k)) // '(A13," "))'
        write(untOut,trim(fmt)) headers
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
            write(untOut,trim(fmt)) head(:),'Objectivefunction',fAuxHead(:)
            deallocate(fAuxHead)
        else
            fmt='(' // trim(number_string(n+1)) // '(A13," "))'
            write(untOut,trim(fmt)) head(:),'Objectivefunction'
        endif
        deallocate(head)
    endif
endif

!Handle Monitoring file
if(present(MonitorFile)) then
    call getSpareUnit(untMonitor,err,mess)
    If(err>0) then
        mess='Adaptive_Metro:'//trim(mess)
        if(present(OutFile)) close(untOut)
        return
    endif
    open(unit=untMonitor,file=MonitorFile,status='replace')
endif

! Check starting value
call f(x=x,feas=feas,isnull=isnull,fx=fx,fAux=fAux,err=err,mess=mess)
If(err>0) then
    mess='Adaptive_Metro:'//trim(mess)
    if(present(OutFile)) close(untOut)
    if(present(MonitorFile)) close(untMonitor)
    return
endif
if((.not.feas).or.isnull) then
    err=1;mess='Adaptive_Metro:Fatal:Unfeasible starting point'
    if(present(OutFile)) close(untOut)
    if(present(MonitorFile)) close(untMonitor)
    return
endif
if(present(OutFile)) then
    if(present(fAux)) then
        fmt='(' // trim(number_string(n+1+p)) // 'e14.6)'
        write (untOut,trim(fmt)) x,fx,fAux
    else
        fmt='(' // trim(number_string(n+1)) // 'e14.6)'
        write (untOut,trim(fmt)) x,fx
    endif
endif

again=.true.;compt=0
if(present(MonitorFile)) call writeMonitor(untMonitor,nAdapt*compt,nAdapt*nCycles)
do while(again)
    !Start sampling
    moverate=0._mrk
    do i=1,nAdapt
        call Metro_GaussJump(f,x,fx,fAux,(scale**2)*covar,move,err,mess)
        If(err>0) then
            mess='Adaptive_Metro:'//trim(mess)
            if(present(OutFile)) close(untOut)
            if(present(MonitorFile)) close(untMonitor)
            return
        endif
        postlist(i)=fx
        !Write2File
        if(present(OutFile)) then
            if(present(fAux)) then
                fmt='(' // trim(number_string(n+1+p)) // 'e14.6)'
                write (untOut,trim(fmt)) x,fx,fAux
            else
                fmt='(' // trim(number_string(n+1)) // 'e14.6)'
                write (untOut,trim(fmt)) x,fx
            endif
        endif
        !Update Move Rates
        if(move) moverate=moverate+1._mrk/nAdapt
    enddo

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
    if(present(MonitorFile)) call writeMonitor(untMonitor,nAdapt*compt,nAdapt*nCycles)
enddo

if(present(OutFile)) close(untOut)
if(present(MonitorFile)) close(untMonitor)

end subroutine Adaptive_Metro

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Adaptive_Metro_OAAT(f,x,fx,fAux,std,nAdapt,nCycles,&
                MinMoveRate,MaxMoveRate,DownMult,UpMult,&
                OutFile,MonitorFile,headers,err,mess)

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
use utilities_dmsl_kit, only:number_string,getSpareUnit

real(mrk),intent(inout)::x(:),std(:)
integer(mik),intent(in)::nAdapt,nCycles
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in),optional::OutFile,MonitorFile,headers(:)
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
logical::feas,isnull,again
character(100), allocatable::head(:), fAuxHead(:)
logical,allocatable::move(:)
real(mrk),allocatable::moverate(:)
integer(mik)::n,i,p,k,compt,untOut,untMonitor
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
    call getSpareUnit(untOut,err,mess)
    If(err>0) then;mess='Adaptive_Metro_OAAT:'//trim(mess);return;endif
    open(unit=untOut,file=OutFile,status='replace')
    if(present(headers)) then ! trust user - no size check at the moment
        k=size(headers)
        fmt='(' // trim(number_string(k)) // '(A13," "))'
        write(untOut,trim(fmt)) headers
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
            write(untOut,trim(fmt)) head(:),'Objectivefunction',fAuxHead(:)
            deallocate(fAuxHead)
        else
            fmt='(' // trim(number_string(n+1)) // '(A13," "))'
            write(untOut,trim(fmt)) head(:),'Objectivefunction'
        endif
        deallocate(head)
    endif
endif

!Handle Monitoring file
if(present(MonitorFile)) then
    call getSpareUnit(untMonitor,err,mess)
    If(err>0) then
        mess='Adaptive_Metro:'//trim(mess)
        if(present(OutFile)) close(untOut)
        return
    endif
    open(unit=untMonitor,file=MonitorFile,status='replace')
endif

! Check starting value
call f(x=x,feas=feas,isnull=isnull,fx=fx,fAux=fAux,err=err,mess=mess)
If(err>0) then
    mess='Adaptive_Metro_OAAT:'//trim(mess)
    if(present(OutFile)) close(untOut)
    if(present(MonitorFile)) close(untMonitor)
    return
endif
if((.not.feas).or.isnull) then
    err=1;mess='Adaptive_Metro_OAAT:Fatal:Unfeasible starting point'
    if(present(OutFile)) close(untOut)
    if(present(MonitorFile)) close(untMonitor)
    return
endif
if(present(OutFile)) then
    if(present(fAux)) then
        fmt='(' // trim(number_string(n+1+p)) // 'e14.6)'
        write (untOut,trim(fmt)) x,fx,fAux
    else
        fmt='(' // trim(number_string(n+1)) // 'e14.6)'
        write (untOut,trim(fmt)) x,fx
    endif
endif
! Allocate MoveRate stuff
if(allocated(move)) deallocate(move)
allocate(move(n))
if(allocated(moverate)) deallocate(moverate)
allocate(moverate(n))

again=.true.;compt=0
if(present(MonitorFile)) call writeMonitor(untMonitor,nAdapt*compt,nAdapt*nCycles)
do while(again)
    !Start sampling
    moverate=0._mrk
    do i=1,nAdapt
        call Metro_OAAT_GaussJump(f,x,fx,fAux,std,move,err,mess)
        If(err>0) then
            mess='Adaptive_Metro_OAAT:'//trim(mess)
            if(present(OutFile)) close(untOut)
            if(present(MonitorFile)) close(untMonitor)
            return
        endif
        postlist(i)=fx
        !Write2File
        if(present(OutFile)) then
            if(present(fAux)) then
                fmt='(' // trim(number_string(n+1+p)) // 'e14.6)'
                write (untOut,trim(fmt)) x,fx,fAux
            else
                fmt='(' // trim(number_string(n+1)) // 'e14.6)'
                write (untOut,trim(fmt)) x,fx
            endif
        endif
        !Update Move Rates
        where(move) moverate=moverate+1._mrk/nAdapt
    enddo

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
    if(present(MonitorFile)) call writeMonitor(untMonitor,nAdapt*compt,nAdapt*nCycles)
enddo

if(present(OutFile)) close(untOut)
if(present(MonitorFile)) close(untMonitor)

end subroutine Adaptive_Metro_OAAT

!*********!
! PRIVATE !
!*********!

subroutine writeMonitor(unt,i,n)
use utilities_dmsl_kit, only:number_string
integer(mik),intent(in)::unt,i,n
rewind(unt)
write (unt,'(A)') trim(number_string(i))//'/'//trim(number_string(n))
end subroutine


end module MCMCStrategy_tools

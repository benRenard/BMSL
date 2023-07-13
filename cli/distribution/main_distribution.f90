program distribution

use kinds_dmsl_kit
use utilities_dmsl_kit,only:countSubstringInString,number_string
use Distribution_tools

implicit none

!-----------------------
! Constants
character(len_stdStrD),parameter::version="0.1.0 July 2023"
character(250),parameter::fmt_numeric='e24.15',fmt_string='A24'
integer(mik),parameter::nsim_def=100,nx_def=100
character(len_vLongStr),parameter::resFile_def='distribution_result.txt'
real(mrk),parameter:: low_def=0.001_mrk,high_def=0.999_mrk
!-----------------------
! Distribution specification
character(len_stdStrL)::distID,act
real(mrk),allocatable::par(:)
integer(mik)::npar,nsim
!-----------------------
! Results
real(mrk),allocatable::res(:),xgrid(:)
real(mrk):: low,high
integer(mik)::nx
character(len_vLongStr)::resFile
!-----------------------
! Misc.
integer(mik)::i,err,narg,k
character(len_vLongStr)::mess,arg
logical::feas,isnull
!-----------------------

! init
distID='';npar=undefIN;act=''
nsim=nsim_def
resFile=resFile_def
nx=undefIN

! Interpret command line arguments
i=1
narg=command_argument_count()
do while (i<=narg)
     call get_command_argument(i,arg)
     select case (arg)
     case ('-name', '--name')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            DistID=trim(arg)
            i=i+1
        else
            call consoleMessage(-1,'-name requires the name of the distribution as a string')
        endif
     case ('-par', '--parameters')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            npar=countSubstringInString(arg,',')+1
            allocate(par(npar)); par=undefRN
            read(arg,*,iostat=err) par
            if(err/=0) then
                call consoleMessage(-1,'-par requires the list of parameters (comma-separated list of reals)')
            endif
            i=i+1
        else
            call consoleMessage(-1,'-par requires the list of parameters (comma-separated list of reals)')
        endif
     case ('-act', '--action')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            act=trim(arg)
            i=i+1
        else
            call consoleMessage(-1,'-act requires the action to be performed as a string in: d,p,q,r')
        endif
     case ('-v', '--version')
        write(*,*) 'version: ', trim(version)
        STOP
     case ('-h', '--help')
        call printHelp()
        STOP
     case ('-n', '--nsim')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            read(arg,*,iostat=err) nsim
            if(err/=0) then
                call consoleMessage(-1,'-n requires the number of simulated values as a positive integer')
            endif
            i=i+1
        else
            call consoleMessage(-1,'-n requires the number of simulated values as a positive integer')
        endif
     case ('-rf', '--result')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            resFile=trim(arg)
            i=i+1
        else
            call consoleMessage(-1,'-rf requires the name of the results file as a path to a file')
        endif
     case ('-x', '--xgrid')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            k=countSubstringInString(arg,',')
            if(k/=2) then
                call consoleMessage(-1,'-x should define a sequence with 3 comma-separated numbers as follows: low,high,nvalues ')
            endif
            read(arg,*,iostat=err) low,high,nx
            if(err/=0) then
                call consoleMessage(-1,'-x should define a sequence with 3 comma-separated numbers as follows: low,high,nvalues ')
            endif
            i=i+1
        else
            call consoleMessage(-1,'-x should define a sequence with 3 comma-separated numbers as follows: low,high,nvalues ')
        endif
     case default
        write(*,*) 'Unrecognized command-line option: ', trim(arg)
        call printHelp()
        call fatalExit
     end select
  end do

! Check compulsory arguments are present
if(act=='' .or. distID=='' .or. npar==undefIN) then
    call consoleMessage(-1,'one of the compulsory arguments -name (name of the distribution), &
                           &-par (parameters) or -act (action to perform) is missing.')
endif

! Check number of parameters and feasibility
call GetParNumber(DistID, i, err, mess)
if(err>0) then; call consoleMessage(-1,trim(mess));endif
if(i/=npar) then
    call consoleMessage(-1,'npar='//trim(number_string(npar))//' is not compatible with distribution '//trim(distID))
endif
call GetParFeas(DistID,par,feas,err,mess)
if(err>0) then; call consoleMessage(-1,trim(mess));endif
if(.not.feas) then
    call consoleMessage(-1,'Unfeasible parameter vector')
endif

! Handle xgrid
if(nx==undefIN) then ! xgrid is unspecified, use defaults
    nx=nx_def
    allocate(xgrid(nx))
    if(act=='d' .or. act=='p') then ! grid of values is needed
        call GetQuantile(DistId=DistId,p=low_def,par=par,q=low,feas=feas,err=err,mess=mess)
        if(err>0) then; call consoleMessage(-1,trim(mess));endif
        call GetQuantile(DistId=DistId,p=high_def,par=par,q=high,feas=feas,err=err,mess=mess)
        if(err>0) then; call consoleMessage(-1,trim(mess));endif
        call getGrid(low,high,nx,xgrid)
    else if (act=='q') then ! grid of probabilities is needed
        call getGrid(low_def,high_def,nx,xgrid)
    else ! no grid needed - do nothing
    endif
else ! use user-provided xgrid
    allocate(xgrid(nx))
    call getGrid(low,high,nx,xgrid)
endif

! Do as required
select case (trim(act))
    case ('d','p','q')
        allocate(res(nx))
        select case (trim(act))
            case('d')
                do i=1,nx
                    call GetPdf(DistId=DistId,x=xgrid(i),par=par,loga=.false.,pdf=res(i),feas=feas,isnull=isnull,err=err,mess=mess)
                    if(err>0) then; call consoleMessage(-1,trim(mess));endif
                enddo
            case('p')
                do i=1,nx
                    call GetCdf(DistId=DistId,x=xgrid(i),par=par,cdf=res(i),feas=feas,err=err,mess=mess)
                    if(err>0) then; call consoleMessage(-1,trim(mess));endif
                enddo
            case('q')
                do i=1,nx
                    call GetQuantile(DistId=DistId,p=xgrid(i),par=par,q=res(i),feas=feas,err=err,mess=mess)
                    if(err>0) then; call consoleMessage(-1,trim(mess));endif
                enddo
        end select
    case ('r')
        allocate(res(nsim))
        call GenerateSample(DistId=DistId,par=par,gen=res,feas=feas,err=err,mess=mess)
        if(err>0) then; call consoleMessage(-1,trim(mess));endif
    case default
        call consoleMessage(-1,'Unrecognized action: '//trim(act)//'. See distribution --help.')
        call fatalExit
end select

! Write result to file and stop
open(unit=1,file=resFile,status='REPLACE')
if(act=='r') then
    do i=1,size(res)
        write(1,'('//trim(fmt_numeric)//')') res(i)
    enddo
else
    do i=1,size(res)
        write(1,'('//trim(fmt_numeric)//','//trim(fmt_numeric)//')') xgrid(i),res(i)
    enddo
endif
close(1)
call consoleMessage(1,resFile)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine printHelp()
!^**********************************************************************
!^* Purpose: print help on command line arguments in console
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified:13/07/2023
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************

    write(*,'(a)') 'usage: distribution -name XXX -par XXX -a XXX [OPTIONS]'
    write(*,'(a)') '  -name, --name:..............name of the distribution, character string'
    write(*,'(a)') '  -par, --parameters:.........parameters, comma-separated list of reals'
    write(*,'(a)') '  -act, --action:.............action, character in d/p/q/r indicating the action to perform'
    write(*,'(a)') '                              d=density (pdf), p=non-exceddance probability (cdf), '
    write(*,'(a)') '                              q=quantile, r=random realizations'
    write(*,'(a)') 'available options:'
    write(*,'(a)') '  -v, --version:..............print version information and exit'
    write(*,'(a)') '  -h, --help:.................print help and exit'
    write(*,'(a)') '  -n, --nsim XXX:.............number of realizations when -act is r. Default 1000'
    write(*,'(a)') '  -x, --xgrid XXX:............computation grid when -act is d,p or q, in the form: low,high,nvalues'
    write(*,'(a)') '  -rf, --result XXX:..........path to results file. Default distribution_result.txt'
end subroutine printHelp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine consoleMessage(id,mess)
!^**********************************************************************
!^* Purpose: console messages printed during execution
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:13/07/2023
!^**********************************************************************
!^* Comments: negative id's are for fatal error messages
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. id, message id
!^*    2. mess, message that will be parsed into the printed message
!^*             - typically an error message issued by some subroutine
!^**********************************************************************

integer(mik),intent(in)::id
character(*), intent(in)::mess

select case(id)
case(1)
    write(*,*) '*********************************'
    write(*,*) 'Done!'
    write(*,*) 'Results written to file = ', trim(mess)
    write(*,*) '*********************************'
    write(*,*) ''
case(2)
    ! available
case(3)
    ! available
case (-1) ! Fatal Error - general
    write(*,*) ""
    write(*,*) "A FATAL ERROR has occured:"
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call fatalExit
case(-8) ! Fatal Error - writing to a file
    write(*,*) ""
    write(*,*) "A FATAL ERROR has occured"
    write(*,*) "while writting to the file:"
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call fatalExit
case default
    write(*,*) ""
    write(*,*) "a FATAL ERROR has occured"
    write(*,*) "with unknown error ID..."
    write(*,*) "Execution will stop"
    write(*,*) ""
    call fatalExit
end select
end subroutine consoleMessage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fatalExit
! action taken on fatal error
!read(*,*)
STOP
end subroutine fatalExit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getGrid(low,high,nx,xgrid)
real(mrk),intent(in)::low,high
integer(mik), intent(in)::nx
real(mrk),intent(out)::xgrid(nx)
! locals
integer(mik)::i

if(nx>1) then
    forall(i=1:nx) xgrid(i)=low+(i-1)*(high-low)/(nx-1)
else
    xgrid=low
endif

end subroutine getGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program distribution


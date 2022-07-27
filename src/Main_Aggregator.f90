program Main_Aggregator

! use what?
use kinds_dmsl_kit
use BaRatin_tools
use Dates_tools
use TimeSeries_tools
use aggregation_tools
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! constants
character(250),parameter::  Config_file="Config_BaRatinAggregator.txt"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variables read from config file
character(250):: Tfile_IN,Tfile_OUT,Sfile_IN,Sfile_OUT,scale,funk
integer(mik)::Tfile_nhead,Sfile_nhead,nT,Tformat_IN,Tformat_OUT,nspag,scaleMult,option
integer(mik), dimension(6)::FirstDate, LastDate
logical::safeMV,safeInterpol
real(mrk)::mvcode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! locals
integer(mik)::i,err,qcOption,ccOption
character(250)::mess
type(OneDTimeSeriesType)::TSin,TSout
type(DateType),allocatable::dates(:)
type(DateType)::d0,dn
real(mrk),allocatable::spagIN(:,:),spagOUT(:,:)
real(mrk)::tstart,tend
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Welcome user and read Workspace
call BaRatin_ConsoleMessage(101,'')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read configuration file
call BaRatinConfig_Read_Aggregation(trim(Config_file),&
                                   Tfile_IN,Tfile_nhead,nT,Tformat_IN,Tfile_OUT,Tformat_OUT,&
                                   Sfile_IN,Sfile_nhead,nspag,mvcode,Sfile_OUT,&
                                   FirstDate,LastDate,scale,scaleMult,&
                                   funk,safeMV,safeInterpol,option,&
                                   err,mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif

! Read input files: time
TSin%n=nT;allocate(TSin%ts(nT));TSin%nameSeries="Input"
if(allocated(dates)) deallocate(dates);allocate(dates(nT))
call ReadTimeFile(file=trim(Tfile_IN),fType=Tformat_IN,nhead=Tfile_nhead,dates=dates,err=err,mess=mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
TSin%ts(:)%date=dates
deallocate(dates)

! Read input files: spaghettis
write(*,*) 'Reading spaghetti file...:'
call CPU_TIME(tstart)
if(allocated(spagIN)) deallocate(spagIN);allocate(spagIN(nT,nspag))
call ReadSpagFile(file=trim(Sfile_IN),nhead=Sfile_nhead,nT=nT,nspag=nspag,spag=spagIN,err=err,mess=mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
call CPU_TIME(tend)
write(*,*) 'Done! in:',tend-tstart, ' seconds'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop on each spaghetti and aggregate
write(*,*) 'Aggregating spaghettis...:'
call CPU_TIME(tstart)
if(safeMV) then; qcOption=INFECT_qc;else; qcOption=IGNORE_qc;endif
if(safeInterpol) then; ccOption=CONSERV_cc;else; ccOption=LIBERAL_cc;endif
d0%year=FirstDate(1);d0%month=FirstDate(2);d0%day=FirstDate(3);d0%hour=FirstDate(4);d0%minute=FirstDate(5);d0%second=FirstDate(6);
dn%year=LastDate(1);dn%month=LastDate(2);dn%day=LastDate(3);dn%hour=LastDate(4);dn%minute=LastDate(5);dn%second=LastDate(6);
TSin%mv=mvcode
do i=1,nspag
    write(*,*) 'Processing spaghetti:',i,'/',nspag
    TSin%ts(:)%q=spagIN(:,i)
    call TSAggregate(operateur=funk,qcOption=qcOption,ccOption=ccOption,&
                        TSin=TSin,scale=scale,scaleMult=scaleMult,&
                        DateListOption=option,FirstDate=d0,LastDate=dn,&
                        TSout=TSout,err=err,message=mess)
    if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
    if(.not.allocated(spagOUT)) then ! first pass - allocate spagOUT and record TSout%dates
        allocate(spagOUT(TSout%n,nspag));spagOUT=undefRN;
        if(allocated(dates)) deallocate(dates);allocate(dates(TSout%n));dates=TSout%ts(:)%date
    endif
    spagOUT(:,i)=TSout%ts(:)%q
    deallocate(TSout%ts)
enddo
call CPU_TIME(tend)
write(*,*) 'Done! in:',tend-tstart, ' seconds'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write to file
write(*,*) 'Writing aggregated spaghettis to file...:'
call CPU_TIME(tstart)
call WriteTimeFile(file=trim(Tfile_OUT),fType=Tformat_OUT,dates=dates,err=err,mess=mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
call WriteSpagFile(file=trim(Sfile_OUT),spag=spagOUT,err=err,mess=mess)
if(err>0) then; call BaRatin_ConsoleMessage(-1,trim(mess));endif
call CPU_TIME(tend)
write(*,*) 'Done! in:',tend-tstart, ' seconds'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Good bye!
call BaRatin_ConsoleMessage(999, '')
read(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadTimeFile(file,fType,nhead,dates,err,mess)
! Note: not integrated into module BaRatin_tools to avoid introducing dates into BaRatin
use utilities_dmsl_kit,only:getSpareUnit,ReplacedString
character(*), intent(in)::file
integer(mik),intent(in)::fType,nhead
type(DateType), intent(out)::dates(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
character(250), parameter::procname="ReadTimeFile"
integer(mik)::i,n,unt
character(250)::foo

err=0;mess=''
n=size(dates)
call getspareunit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
do i=1,nhead
    read(unt,*,iostat=err) ! skip
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo
do i=1,n    
    read(unt,'(A)',iostat=err) foo
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
    foo=ReplacedString(foo,"'","")
    foo=ReplacedString(foo,'"','')
    call UnformatDate(trim(foo),fType,dates(i),err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo
close(unt)
end subroutine ReadTimeFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadSpagFile(file,nhead,nT,nspag,spag,err,mess)
! Note: not integrated into module BaRatin_tools because still experimental: may need to load the file column-y-column...
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik),intent(in)::nhead,nT,nspag
real(mrk), intent(out)::spag(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
character(250), parameter::procname="ReadSpagFile"
integer(mik)::i,n,unt
real(mrk):: foo

err=0;mess=''
call getspareunit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
do i=1,nhead
    read(unt,*,iostat=err) ! skip
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo
do i=1,nT    
    read(unt,*,iostat=err) foo,spag(i,:)
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo
close(unt)
end subroutine ReadSpagFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteTimeFile(file,fType,dates,err,mess)
! Note: not integrated into module BaRatin_tools to avoid introducing dates into BaRatin
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik),intent(in)::fType
type(DateType), intent(in)::dates(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
character(250), parameter::procname="WriteTimeFile"
integer(mik)::i,n,unt
character(250)::foo

err=0;mess=''
n=size(dates)
call getspareunit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='replace', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
do i=1,n    
    call FormatDate(dates(i),fType,foo,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    write(unt,'(A)',iostat=err) foo
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Read,trim(file));endif
enddo
close(unt)
end subroutine WriteTimeFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteSpagFile(file,spag,err,mess)
! Note: not integrated into module BaRatin_tools because still experimental: may need to load the file column-y-column...
use utilities_dmsl_kit,only:getSpareUnit,number_string
character(*), intent(in)::file
real(mrk), intent(in)::spag(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
character(250), parameter::procname="WriteSpagFile"
integer(mik)::i,p,n,unt
character(250)::fmt,head(size(spag,dim=2))

err=0;mess=''
n=size(spag,dim=1)
p=size(spag,dim=2)
fmt="("//trim(number_string(p))//trim(BaRatin_realFmt)//")"
call getspareunit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='replace', iostat=err)
if(err>0) then;call BaRatin_ConsoleMessage(messID_Open,trim(file));endif
do i=1,p;head(i)="Q_"//trim(number_string(i));enddo
write(unt,'(<p>A15)') head
do i=1,n    
    write(unt,trim(fmt),iostat=err) spag(i,:)
    if(err>0) then;call BaRatin_ConsoleMessage(messID_Write,trim(file));endif
enddo
close(unt)
end subroutine WriteSpagFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program Main_Aggregator

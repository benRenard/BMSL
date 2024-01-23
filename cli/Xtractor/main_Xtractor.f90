program Xtractor

use kinds_dmsl_kit
use DataRW_tools, only: DTSRead, WriteSeparatedFile
use TimeSeries_tools
use HydroIndices_tools, only:AnnualIndices_daily
use Dates_tools

implicit none

! character(250), parameter::exepath="..\..\..\Core\Xtractor\" ! Debug/release
! character(250), parameter::exepath="..\..\Core\Xtractor\" ! Core
character(250), parameter::exepath="" ! Core
character(250), parameter::ConfigFile="Xtractor_Config.txt"
character(250), parameter::sep=";" 
integer(mik),parameter::err_ok=0,err_openConfig=1,err_readConfig=2,err_unknownFilter=3,&
    err_unknownXtractor=4,err_readInput=5,err_unknownInputFormat=6,err_runFilter=7,&
    err_runXtractor=8,err_DorV=9,err_writeOutput=10,err_unknownOutputFormat=11

type(OneDTimeSeriesType)::y0,yIn,QcodeIn
integer(mik)::err, i, n, dStart, mStart, dEnd, mEnd
character(250):: InFile,OutFile, InFormat, OutFormat, mess, varname, &
                  Xwhat, head(3), filter
type(DateType)::YearStart, YearEnd
type(DateType), pointer :: youtDate(:),youtDate3(:,:,:)
real(mrk), pointer::yout(:), pmv(:), yout2(:,:)
integer(mik), pointer::YearList(:)
real(mrk), allocatable::mat(:,:),filterPar(:)
real(mrk)::Xpar, Xpar3
character(250)::Xpar2

! Open and read config file
open(unit=1, file=trim(exepath)//trim(ConfigFile),status="OLD",iostat=err)
if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem opening config file';call exit(err_openConfig);endif
read(1,'(A250)', iostat=err) InFile
if(err/=0) then; write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
read(1,'(A250)', iostat=err) InFormat
if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
read(1,'(A250)', iostat=err) OutFile
if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
read(1,'(A250)', iostat=err) OutFormat
if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
read(1,'(A250)', iostat=err) filter
if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
select case(filter)
  case("none")
    if(allocated(filterPar)) deallocate(filterPar)
    allocate(filterPar(0))
    read(1,*, iostat=err) 
    if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
  case("MovingAverage", "MovingMin", "MovingMax", "SPA")
    if(allocated(filterPar)) deallocate(filterPar)
    allocate(filterPar(1))
    read(1,*, iostat=err) filterPar
    if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
  case("BFS")
    if(allocated(filterPar)) deallocate(filterPar)
    allocate(filterPar(2))
    read(1,*, iostat=err) filterPar
    if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
  case default
    write(*,*) 'Xtractor:FATAL:Unknown Filter';call exit(err_unknownFilter);
end select
read(1,'(A250)', iostat=err) Xwhat
if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
select case(Xwhat)
  case("Amax", "AmaxDate", "Amean", "Amin", "AminDate", "Asum")
    read(1,*, iostat=err)
    if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
  case("Aquantile", "DroughtDuration", "DroughtVolumeDeficit", "HighFlowDuration", "HighFlowVolume")
    read(1,*, iostat=err) XPar
    if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
  case("DroughtCenterOfMass", "HighFlowCenterOfMass")
    read(1,*, iostat=err) XPar, Xpar3, Xpar2
    if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
  case default
    write(*,*) 'Xtractor:FATAL:Unknown Xtraction';call exit(err_unknownXtractor);
end select
read(1,*, iostat=err) dStart, mStart
if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
read(1,*, iostat=err) dEnd, mEnd
if(err/=0) then;write(*,*) 'Xtractor:FATAL:Problem reading config file';call exit(err_readConfig);endif
read(1,'(A250)', iostat=err) varname
if(err/=0) varname=Xwhat
close(1)

! Read series in input format
select case(InFormat)
  case("bay_dts")
    call DTSRead(file=InFile,y=y0,Qcode=QcodeIn,&
                     err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem reading Input file';call exit(err_readInput);endif
  case default
    write(*,*) 'Xtractor:FATAL:Unknown input format';call exit(err_unknownInputFormat);
end select

! Apply Filter
select case(filter)
  case("none")
    yin=y0
  case("MovingAverage")
    call TSFiltering_daily(filter=MovingAverage,param=filterPar,&
                           xin=y0,xout=yin,err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem applying filter';call exit(err_runFilter);endif
  case("MovingMin")
    call TSFiltering_daily(filter=MovingMin,param=filterPar,&
                           xin=y0,xout=yin,err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem applying filter';call exit(err_runFilter);endif
  case("MovingMax")
    call TSFiltering_daily(filter=MovingMax,param=filterPar,&
                           xin=y0,xout=yin,err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem applying filter';call exit(err_runFilter);endif
  case("SPA")
    call TSFiltering_daily(filter=SPA,param=filterPar,&
                           xin=y0,xout=yin,err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem applying filter';call exit(err_runFilter);endif
  case("BFS")
    call TSFiltering_daily(filter=BFS,param=filterPar,&
                           xin=y0,xout=yin,err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem applying filter';call exit(err_runFilter);endif
  case default
    write(*,*) 'Xtractor:FATAL:Unknown Filter';call exit(err_unknownFilter);
end select


!Extraction
YearStart%day=dStart
YearStart%month=mStart
YearEnd%day=dEnd
YearEnd%month=mEnd
YearStart%year=1
YearEnd%year=1

select case(trim(Xwhat))
  
  case("Amax")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  AMax=yout,& !Annual maxima
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif

  case("AmaxDate")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  AMaxDate=youtDate,& !Date of Annual maxima
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    n=size(youtDate)
    if(associated(yout)) deallocate(yout)
    allocate(yout(n))
    do i=1,n
		  if (youtDate(i)%year==undefIN) then
			  yout(i)=yin%mv
	    else
			  yout(i)= Date2Days(youtDate(i)-Datetype(YearList(i),YearStart%month,YearStart%day,0,0,0))
		  endif
    enddo

  case("Amean")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  Modul=yout,& !Annual mean
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
  
  case("Amin")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  AMin=yout,& !Annual minima
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif

  case("AminDate")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  AMinDate=youtDate,& !Date of Annual minima
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    n=size(youtDate)
    if(associated(yout)) deallocate(yout)
    allocate(yout(n))
    do i=1,n
		  if (youtDate(i)%year==undefIN) then
			  yout(i)=yin%mv
	    else
			  yout(i)= Date2Days(youtDate(i)-Datetype(YearList(i),YearStart%month,YearStart%day,0,0,0))
		  endif
    enddo

  case("Aquantile")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  pval=(/Xpar/),Aquant=yout2,& !Annual quantile
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    if(associated(yout)) deallocate(yout)
    allocate(yout(size(yout2,dim=1)))
    yout=yout2(:,1)

  case("Asum")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  Total=yout,& !Annual sum
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
  
  case("DroughtDuration")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  lowT=(/Xpar/),lowD=yout2,& !Drought duration
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    if(associated(yout)) deallocate(yout)
    allocate(yout(size(yout2,dim=1)))
    yout=yout2(:,1)
  
  case("DroughtVolumeDeficit")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  lowT=(/Xpar/),lowV=yout2,& !Drought Volume Deficit
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    if(associated(yout)) deallocate(yout)
    allocate(yout(size(yout2,dim=1)))
    yout=yout2(:,1)
  
  case("DroughtCenterOfMass")
    select case (Xpar2)
    case("D")
      call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  lowT=(/Xpar/),lowP=(/Xpar3/), lowComD=youtDate3,& !LowComD
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
      if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    case("V")
      call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  lowT=(/Xpar/),lowP=(/Xpar3/), lowComV=youtDate3,& !LowComD
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
      if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    case default
      write(*,*) 'Xtractor:FATAL:Unknown parameter for DroughtCenterOfMass [par2 be D or V]';call exit(err_DorV);
    end select
    n=size(youtDate3, dim=1)
    if(associated(yout)) deallocate(yout)
    allocate(yout(n))
    do i=1,n
		  if (youtDate3(i,1,1)%year==undefIN) then
			  yout(i)=yin%mv
	    else
			  yout(i)= Date2Days(youtDate3(i,1,1)-Datetype(YearList(i),YearStart%month,YearStart%day,0,0,0))
		  endif
    enddo

  case("HighFlowDuration")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  hiT=(/Xpar/),hiD=yout2,& !High Flow duration
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    if(associated(yout)) deallocate(yout)
    allocate(yout(size(yout2,dim=1)))
    yout=yout2(:,1)
  
  case("HighFlowVolume")
    call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  hiT=(/Xpar/),hiV=yout2,& !High flow Volume 
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
    if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    if(associated(yout)) deallocate(yout)
    allocate(yout(size(yout2,dim=1)))
    yout=yout2(:,1)
  
  case("HighFlowCenterOfMass")
    select case (Xpar2)
    case("D")
      call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  hiT=(/Xpar/),hiP=(/Xpar3/), hiComD=youtDate3,& !LowComD
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
      if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    case("V")
      call AnnualIndices_daily(xin=yin,&
				  YearStart=YearStart,YearEnd=YearEnd, & !define season of interest (IN)
				  hiT=(/Xpar/),hiP=(/Xpar3/), hiComV=youtDate3,& !LowComD
				  Pmv=pmv,& !percentage of missing values (% of season length)
				  YearList=YearList,& !convention: year from 10/1998 to 03/1999 is 1998
				  err=err,mess=mess)
      if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem during Xtraction';call exit(err_runXtractor);endif
    case default
      write(*,*) 'Xtractor:FATAL:Unknown parameter for DroughtCenterOfMass [par2 be D or V]';call exit(err_DorV);
    end select
    n=size(youtDate3, dim=1)
    if(associated(yout)) deallocate(yout)
    allocate(yout(n))
    do i=1,n
		  if (youtDate3(i,1,1)%year==undefIN) then
			  yout(i)=yin%mv
	    else
			  yout(i)= Date2Days(youtDate3(i,1,1)-Datetype(YearList(i),YearStart%month,YearStart%day,0,0,0))
		  endif
    enddo

  case default
    write(*,*) 'Xtractor:FATAL:Unknown Xtraction';call exit(err_unknownXtractor);
    call exit(-1)

end select

n=size(YearList)
if(allocated(mat)) deallocate(mat)
allocate(mat(n,3))
mat(:,1)=YearList
mat(:,2)=yout
mat(:,3)=pmv
head=(/'year',trim(varname),'%missing'/)
!! Write series in output format
select case(OutFormat)
  case("bay_var") 
  call WriteSeparatedFile(file=OutFile,sep=sep,y=mat,headers=head,err=err,mess=mess)
  if(err>0) then;write(*,*) 'Xtractor:FATAL:Problem Writing Output file';call exit(err_writeOutput);endif
  case default
    write(*,*) 'Xtractor:FATAL:Unknown Output format';call exit(err_unknownOutputFormat);
end select

write(*,*) 'Xtractor ran successfully'
call exit(err_ok)
end program Xtractor
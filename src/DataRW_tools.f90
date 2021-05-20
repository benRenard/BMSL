module DataRW_tools

!~**********************************************************************
!~* Purpose: Read/Write various data format encountered here and there...
!~* After reading, put the data in a variable of type TimeSeries - 
!~* module TimeSeries_tools is therefore needed 
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~*             Antoine Bard, Cemagref Lyon.
!~**********************************************************************
!~* Last modified: 18/06/2013
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1.ReadSQR: SQR format (MF)
!~*    2.ReadMF: Daily rainfall (MF)
!~*    3.ReadDTG_R: Daily Rainfall (EDF-DTG)
!~*    4.ReadDTG_Q: Daily Runoff (EDF-DTG)
!~*    5.ReadHYDROII: Daily Runoff (HYDRO II)
!~*    6.ExtraRead: Daily Data, either P or Q (EXTRAFLO format)
!~*    7.ExtraWrite: Daily Data, either P or Q (EXTRAFLO format)
!~*    8.DTSRead: Daily Data
!~*    9.DTSWrite: Daily Data
!~*    10.DatWrite: Data matrix
!~*    11.DatRead: Data matrix
!~*    12.HodgkinsWrite: csv-like format requested by Glenn Hodgkins
!~*    13.ReadRHN: read in RHN format as defined by Paul Whitfield
!~*    14.ReadMalekiFormat: read in Maleki's format
!~*    14.ReadXceedances_GlennFormat
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use TimeSeries_tools

implicit none
Private
public ::ReadSQR,ExtraWrite,ReadDTG_R, ReadDTG_Q, ReadMF, ExtraRead, ReadHYDROII,ReadDLY,&
         DTSWrite, DTSRead, DatWrite,DatRead,HodgkinsWrite,AdaptAlpRead,ReadOz,ReadRHN,ReadRHN4,ReadMalekiFormat,&
         ReadXceedances_GlennFormat,LoadFromRepository,ReadSeparatedFile,WriteSeparatedFile

type, public:: MetadataType
    integer(mik) :: X=-99,Y=-99,Z=-99
    character(250)::name='AintGotNoName'
    character(25)::Zc='Zunknown'
    real(mrk)::Atopo=undefRN,Ahydro=undefRN ! Catchment areas
end type MetadataType

interface WriteSeparatedFile
  module procedure WriteSeparatedFile_r,WriteSeparatedFile_i
endinterface WriteSeparatedFile

contains

subroutine ReadSQR(file,y,Qcode,metadata,err,mess)

!^**********************************************************************
!^* Purpose: read data in SQR formar (MeteoFrance)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 23/02/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^*    2.
!^*    3.
!^* OUT
!^*    1.y, 1D-time series of daily data
!^*    2.Qcode, 1D-time series of quality code
!^*    3.metadata
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,Qcode
type(MetadataType),intent(out)::metadata
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
!!!!!
real(mrk),parameter::SQR_mv=-999._mrk
integer(mik), parameter:: nHeaders=11, nSkip=4
!!!!!
integer(mik):: j, errcode, compt, nb, crapI, date, qual, qual2
real(mrk)::crapR, val
character(80)::string
character(59)::commune, lieudit

err=0;mess='';
100 format(I8,I9,I5,I6,I6,I6,I6,F6.2,F7.2,F7.2,F6.2,F9.1,I5,F9.1,I5,I5,F7.1,I4,I3)

open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='ReadSQR:FATAL:error opening file';
    return;
endif

! First read just 2 count number of days
do j=1, nHeaders !headers
    read(1,*)
enddo
compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! read metadata
rewind(1)
do j=1,nskip
    read(1,*)
enddo
read(1,'(A80)') string
read(string(21:79),'(A59)') commune
read(1,'(A80)') string
read(string(21:79),'(A59)') lieudit
y%NameSeries=trim(trim(commune)//'*'//trim(lieudit))
Qcode%NameSeries=trim(trim(commune)//'*'//trim(lieudit)//trim('_Qcode'))
metadata%name=y%NameSeries
read(1,'(A80)') string
read(string(21:79),'(I12)') nb
metadata%Z=nb
read(1,'(A80)') string
read(string(21:79),'(I12)') nb
metadata%X=100*nb
read(1,'(A80)') string
read(string(21:79),'(I12)') nb
metadata%Y=100*nb

! read data
rewind(1)
do j=1,nHeaders !headers
    read(1,*)
enddo
do j=1,y%n
    read(1,100,iostat=errcode) date,crapI, crapI, crapI, crapI, crapI, crapI, &
        crapR, crapR, crapR, crapR, crapR, crapI, crapR, crapI, crapI,val,qual,qual2
    if(errcode/=0) exit
    ! date
    y%ts(y%n-j+1)%date%year=int(date/10000)
    y%ts(y%n-j+1)%date%month=int((date-y%ts(y%n-j+1)%date%year*10000)/100)
    y%ts(y%n-j+1)%date%day=int(date-y%ts(y%n-j+1)%date%year*10000-y%ts(y%n-j+1)%date%month*100)
    ! Data
    if(val==SQR_mv) then ! missing value
        y%ts(y%n-j+1)%q=y%mv
    else
        y%ts(y%n-j+1)%q=val
    endif
    Qcode%ts(y%n-j+1)%q=10._mrk*qual2+qual
enddo
Qcode%ts(:)%date=y%ts(:)%date
close(1)

end subroutine ReadSQR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadDTG_R(file,y,Qcode,metadata,err,mess)

!^**********************************************************************
!^* Purpose: read data in EDF-DTG format (Rainfall, Gottardi dataset)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 23/02/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^* OUT
!^*    1.y, 1D-time series of daily data
!^*    2.Qcode, 1D-time series of quality code
!^*    3.metadata
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,Qcode
type(MetadataType),intent(out)::metadata
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
!!!!!
real(mrk),parameter::DTG_mv=-9999.00_mrk
!!!!!
integer(mik):: j, errcode, compt, nheaders
real(mrk)::val
character(80)::string
character(12)::DateString

err=0;mess='';

100 format(A12,F9.2)
open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='ReadDTG_R:FATAL:error opening file';
    return;
endif

! First find nb of headers lines
compt=0;string='';
do while(string(1:1)/='.')
    read(1,*) string
    compt=compt+1
enddo
nHeaders=compt;
compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! metadata - missing in DTG format
y%NameSeries='R'
Qcode%NameSeries='R_Qcode'

! read data
rewind(1)
do j=1,nHeaders !headers
    read(1,*)
enddo
do j=1,y%n
    read(1,100,iostat=errcode) DateString,val
    if(errcode/=0) exit

    ! date
    read(DateString(2:3),'(I2)') y%ts(j)%date%day
    read(DateString(5:6),'(I2)') y%ts(j)%date%month
    read(DateString(8:11),'(I4)') y%ts(j)%date%year

    ! Data
    if(val==DTG_mv) then ! missing value
        y%ts(j)%q=y%mv
    else
        y%ts(j)%q=val
    endif
    Qcode%ts(j)%q=Qcode%mv
enddo
Qcode%ts(:)%date=y%ts(:)%date
close(1)

end subroutine ReadDTG_R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadDTG_Q(file,y,Qcode,metadata,err,mess)

!^**********************************************************************
!^* Purpose: read data in EDF-DTG format (Runoff)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 11/05/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^* OUT
!^*    1.y, 1D-time series of daily data
!^*    2.Qcode, 1D-time series of quality code
!^*    3.metadata
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,Qcode
type(MetadataType),intent(out)::metadata
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
!!!!!
real(mrk),parameter::DTG_mv=-9999.00_mrk
integer(mik), parameter:: nHeaders=5
!!!!!
integer(mik):: j, errcode, compt
real(mrk)::val
character(80)::string
character(12)::DateString

err=0;mess='';

100 format(A12,F9.3)
open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='ReadDTG_R:FATAL:error opening file';
    return;
endif

! First find nb of headers lines
compt=0;string='';
do j=1,nHeaders
    read(1,*)
enddo

compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! metadata - missing in DTG format
y%NameSeries='Q'
Qcode%NameSeries='Q_Qcode'

! read data
rewind(1)
do j=1,nHeaders !headers
    read(1,*)
enddo
do j=1,y%n
    read(1,100,iostat=errcode) DateString,val
    if(errcode/=0) exit

    ! date
    read(DateString(2:3),'(I2)') y%ts(j)%date%day
    read(DateString(5:6),'(I2)') y%ts(j)%date%month
    read(DateString(8:11),'(I4)') y%ts(j)%date%year

    ! Data
    if(val==DTG_mv) then ! missing value
        y%ts(j)%q=y%mv
    else
        y%ts(j)%q=val
    endif
    Qcode%ts(j)%q=Qcode%mv
enddo
Qcode%ts(:)%date=y%ts(:)%date
close(1)

end subroutine ReadDTG_Q

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadMF(file,y,Qcode,metadata,err,mess)

!^**********************************************************************
!^* Purpose: read PJ data in Meteo France format
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 23/02/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^* OUT
!^*    1.y, 1D-time series of daily data
!^*    2.Qcode, 1D-time series of quality code
!^*    3.metadata
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,Qcode
type(MetadataType),intent(out)::metadata
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
!!!!!
real(mrk),parameter::MF_mv=-9999.00_mrk
!!!!!
integer(mik):: j, errcode, compt, date, np
real(mrk)::val
character(250)::string


err=0;mess='';

100 format(I8,A250)
open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='ReadMF:FATAL:error opening file';
    return;
endif

! First find nb of days
compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! metadata - missing in MF format
y%NameSeries='R'
Qcode%NameSeries='R_Qcode'

! read data
rewind(1)
do j=1,y%n
    read(1,100,iostat=errcode) date,string
    if(errcode/=0) exit

    ! date
    y%ts(j)%date%year=int(date/10000)
    y%ts(j)%date%month=int((date-y%ts(j)%date%year*10000)/100)
    y%ts(j)%date%day=int(date-y%ts(j)%date%year*10000-y%ts(j)%date%month*100)

    ! Data
    np=index(string(2:), ';')
    read(string(2:np),*) val
    if(val==MF_mv) then ! missing value
        y%ts(j)%q=y%mv
    else
        y%ts(j)%q=val
    endif
    Qcode%ts(j)%q=Qcode%mv
enddo
Qcode%ts(:)%date=y%ts(:)%date
close(1)


end subroutine ReadMF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadHYDROII(file,y,Qcode,metadata,err,mess)

!^**********************************************************************
!^* Purpose: read daily runoff data in HYDROII format
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 15/04/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^* OUT
!^*    1.y, 1D-time series of daily data
!^*    2.Qcode, 1D-time series of quality code
!^*    3.metadata
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,Qcode
type(MetadataType),intent(out)::metadata
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik), parameter::nHeaders=3
character(1), parameter::sep=';'
integer(mik)::errcode, compt,j, np1, np2
character::code*3, line*80, deb*80, bigline*250, stationCode*250, stationName*250,valid*1
real(mrk)::alt

err=0;mess=''
open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='ReadHYDROII:FATAL:error opening file';
    return
endif

! First find nb of days
do j=1, nHeaders !headers
    read(1,*)
enddo
compt=0;errcode=0;code='QJO'
do while(errcode==0 .and. code=='QJO') !just count number of days
    compt=compt+1
    read(1,'(A3)',iostat=errcode) code
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! Read Metadata
rewind(1)
read(1,*);read(1,*)
read(1,'(A250)') bigline
! Station Code
np1=index(bigline, sep);np2=np1+index(bigline(np1+1:), sep)
if(np2==np1+1) then
    stationCode='AintGotNoCode'
else
    read(bigline(np1+1:np2-1),'(A)') stationCode
endif
bigline=bigline(np2+1:)
! Station Name
np1=index(bigline, sep)
if(np1==1) then
    stationName='AintGotNoName'
else
    read(bigline(1:np1-1),'(A)') stationName
endif
bigline=bigline(np1+1:)
metadata%name=trim(StationCode)//'_'//trim(stationName)
!Altitude
np1=index(bigline, sep)
if(np1>1) then
    read(bigline(1:np1-1),'(A)') line
    read(line, '(F20.3)') alt
    metadata%Z=floor(alt)
    metadata%Zc=trim(line)
endif
bigline=bigline(np1+1:)
!Catchment area (hydro)
np1=index(bigline, sep)
if(np1>1) then
    read(bigline(1:np1-1),'(A)') line
    read(line, '(F20.3)') metadata%Ahydro
endif
bigline=bigline(np1+1:)
!Catchment area (topo)
np1=index(bigline, sep)
if(np1>1) then
    read(bigline(1:np1-1),'(A)') line
    read(line, '(F20.3)') metadata%Atopo
endif
bigline=bigline(np1+1:)

! read data
do j=1,y%n
    READ(1,'(A80)') line
    !check for daily mean discharge code= QJO
    read(line(1:3), '(A3)') code
    if (code/='QJO') then
        err=1;mess='ReadHYDROII:Fatal:QJOexpected'
        return
    end if
    read(line(14:17), '(I4)') y%ts(j)%date%year
    read(line(18:19), '(I2)') y%ts(j)%date%month
    read(line(20:21), '(I2)') y%ts(j)%date%day
    line=line(23:80)
    np1=index(line, sep)
    if(np1<=1) then
        y%ts(j)%q=y%mv
        Qcode%ts(j)%q=-1
    else
        read(line(1:np1-1), '(A)') deb
        read(deb, '(F20.3)') y%ts(j)%q
        y%ts(j)%q=0.001_mrk*y%ts(j)%q !in  m3.s-1
        read(line( (np1+3):(np1+3) ), '(A1)') valid
        if(valid=='0' .or. valid=='I' .or. valid=='S') then
            y%ts(j)%q=y%mv
            Qcode%ts(j)%q=-1
        else if(valid=='5') then
            Qcode%ts(j)%q=1
        else
            Qcode%ts(j)%q=0
        endif
    end if
enddo

!Qcode
Qcode%ts(:)%date=y%ts(:)%date
!Qcode%ts(:)%q=y%mv
close(1)

end subroutine ReadHYDROII

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadDLY(file,ReadWhat,y,Qcode,metadata,err,mess)

!^**********************************************************************
!^* Purpose: read daily format from GHCN dataset (NOAA)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon, coded @Columbia
!^**********************************************************************
!^* Last modified: 04/02/2013
!^**********************************************************************
!^* Comments: Qcode not handled properly yet - conservative approach 
!^* forcing a mv if any flag is nontrivial. But output Qcode is junk 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Proper handling of flags through Qcode
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^*    2.ReadWhat, what variable? default 'PRCP'
!^*    (see ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt)
!^* OUT
!^*    1.y, 1D-time series of daily data
!^*    2.Qcode, 1D-time series of quality code
!^*    3.metadata
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use Dates_tools, only:DaysInMonth
use TimeSeries_tools, only:FillTheGap_daily
character(*), intent(in)::file
character(*), intent(in), optional::ReadWhat
type(OneDTimeSeriesType), intent(out)::y,Qcode
type(MetadataType),intent(out)::metadata
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
type(OneDTimeSeriesType)::y0
integer(mik), parameter::nHeaders=0,dlyMVflag=-9999
integer(mik)::errcode, compt,i,j,year,month,nDays,nLines,val
character(269)::line
character(250)::RW
character(1)::mflag,qflag

10 format(A269)
11 format(A1)
220 format(I2)
240 format(I4)
250 format(I5)
err=0;mess=''
if(present(ReadWhat)) then;RW=ReadWhat;else;RW='PRCP';endif
open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode;mess='ReadDLY:FATAL:error opening file';return
endif

! First find nb of days
do j=1, nHeaders !headers
    read(1,*)
enddo
compt=0;errcode=0;nLines=0
do while(errcode==0) !just count number of days
    read(1,10,iostat=errcode) line
    if(errcode/=0) exit
    nLines=nLines+1
    if(line(18:21)==trim(RW)) then
        read(line(12:15),240) year
        read(line(16:17),220) month
        compt=compt+DaysInMonth(month,year)
    endif
enddo

! start populating TS type
y0%n=compt;
if(associated(y0%ts)) nullify(y0%ts)
allocate(y0%ts(compt))
if(compt==0) return

! read data
rewind(1)
do j=1, nHeaders !headers
    read(1,*)
enddo
compt=0;
do j=1,nLines
    read(1,10,iostat=errcode) line
    if(line(18:21)==trim(RW)) then ! process monthly line
        read(line(12:15),240) year
        read(line(16:17),220) month
        nDays=DaysInMonth(month,year)
        do i=1,nDays
            compt=compt+1
            read(line( (21+(8*(i-1))+1) : (21+(8*(i-1))+5) ),250) val
            read(line( (21+(8*(i-1))+6) : (21+(8*(i-1))+6) ),11) mflag
            read(line( (21+(8*(i-1))+7) : (21+(8*(i-1))+7) ),11) qflag
            y0%ts(compt)%date%year=year
            y0%ts(compt)%date%month=month
            y0%ts(compt)%date%day=i
            if(val==dlyMVflag .or. mflag/=' ' .or. qflag/=' ') then
                y0%ts(compt)%q=y%mv
            else
                if(trim(RW)=='PRCP') then
                    y0%ts(compt)%q=real(val,mrk)/10._mrk
                else
                    y0%ts(compt)%q=real(val,mrk)
                endif
            endif
        enddo
    endif
enddo
close(1)

! fill gaps
Call FillTheGap_daily(xin=y0,xout=y,err=err,mess=mess)
if(err/=0) then;mess='ReadDLY:'//trim(mess);return;endif
!Qcode
Qcode%n=y%n
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(Qcode%ts(y%n))
Qcode%ts(:)%date=y%ts(:)%date
Qcode%ts(:)%q=y%mv

end subroutine ReadDLY

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadOz(file,y,Qcode,err,mess)

!^**********************************************************************
!^* Purpose: read data in OZ format (Rainfall, SEQ dataset)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon, coded @Columbia
!^**********************************************************************
!^* Last modified: 18/01/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^* OUT
!^*    1.y, 1D-time series of daily data
!^*    2.Qcode, 1D-time series of quality code
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,Qcode
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
!!!!!
real(mrk),parameter::Oz_mv=99999.9_mrk
!!!!!
integer(mik):: j, errcode, compt
real(mrk)::val
character(12)::DateString

err=0;mess='';

100 format(A8,F8.1)
open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='ReadOz:FATAL:error opening file';
    return;
endif

! First find nb of days
read(1,*) ! headers
compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! read data
rewind(1)
read(1,*)!headers
do j=1,y%n
    read(1,100,iostat=errcode) DateString,val
    if(errcode/=0) exit
    ! date
    read(DateString(7:8),'(I2)') y%ts(j)%date%day
    read(DateString(5:6),'(I2)') y%ts(j)%date%month
    read(DateString(1:4),'(I4)') y%ts(j)%date%year

    ! Data
    if(val==Oz_mv) then ! missing value
        y%ts(j)%q=y%mv
    else
        y%ts(j)%q=val
    endif
    Qcode%ts(j)%q=Qcode%mv ! No quality code in this format
enddo
Qcode%ts(:)%date=y%ts(:)%date
close(1)

end subroutine ReadOz

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ExtraWrite(y,Qcode,file,type)

!^**********************************************************************
!^* Purpose: write data into EXTRAFLO format
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 23/02/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.y, data time series
!^*    2.Qcode, quality code time series
!^*    3.file, data file
!^*    4.type, 'R' or 'Q' for rainfall/runoff (default rainfall)
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(in)::y,QCode
character(1), optional::type
!locals
integer(mik)::i
character(1)::typ

if(present(type)) then
    typ=type
else
    typ='R'
endif

100 format(I8,I8,I8,F10.3,I8) !data
200 format(A8,A8,A8,A10,A8) !headers

open(unit=1, file=file,status='REPLACE')
!headers
if(typ=='R') then
    write(1,200) 'Year','Month','Day','R (mm)','QCode'
else
    write(1,200) 'Year','Month','Day','Q (m3/s)','QCode'
endif
do i=1,y%n
    write(1,100) y%ts(i)%date%year,y%ts(i)%date%month,&
                y%ts(i)%date%day, y%ts(i)%q, int(anint(Qcode%ts(i)%q))
enddo
close(1)

end subroutine ExtraWrite

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ExtraRead(file,y,Qcode,err,mess)

!^**********************************************************************
!^* Purpose: Read from EXTRAFLO format
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 24/02/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, destination file
!^* OUT
!^*    1.y, data time series
!^*    2.Qcode, quality code time series
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,QCode
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i, errcode,compt

err=0;mess=''
100 format(I8,I8,I8,F10.3,F8.0) !data

open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='ExtraRead:FATAL:error opening file';
    return;
endif

! First read just 2 count number of days
read(1,*) !headers
compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! read data
rewind(1)
read(1,*) !headers
do i=1,y%n
    read(1,100) y%ts(i)%date%year,y%ts(i)%date%month,&
                y%ts(i)%date%day, y%ts(i)%q, Qcode%ts(i)%q
enddo
Qcode%ts(:)%date=y%ts(:)%date
close(1)

end subroutine ExtraRead

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DTSWrite(y,Qcode,file,varname,err,mess)

!^**********************************************************************
!^* Purpose: write data into DTS format (Daily Time Series)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 18/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.y, data time series
!^*    2.Qcode, quality code time series
!^*    3.file, data file
!^*    4.[varname] (default 'value[]')
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(in)::y,QCode
character(*), intent(in), optional::varname
character(*), intent(out)::mess
integer(mik), intent(out)::err
!locals
integer(mik)::i
character(250)::vname

err=0;mess=''
if(present(varname)) then
    vname=varname
else
    vname='value[]'
endif

100 format(I8,I8,I8,F15.3,I8) !data
200 format(A8,A8,A8,A15,A8) !headers

open(unit=1, file=file,status='REPLACE',iostat=err)
if (err/=0) then
    mess='DTSWrite:FATAL:error opening file';
    return;
endif

!headers
write(1,200) 'Year','Month','Day',trim(vname),'QCode'
!data
do i=1,y%n
    write(1,100) y%ts(i)%date%year,y%ts(i)%date%month,&
                y%ts(i)%date%day, y%ts(i)%q, int(anint(Qcode%ts(i)%q))
enddo
close(1)

end subroutine DTSWrite

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine oldDTSRead(file,y,Qcode,err,mess)

!^**********************************************************************
!^* Purpose: Read from DTS format (Daily Time Series)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 18/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, destination file
!^* OUT
!^*    1.y, data time series
!^*    2.Qcode, quality code time series
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,QCode
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik), parameter::nhead=4
character(1),parameter::sep=";"
integer(mik)::i, errcode,compt

err=0;mess=''
100 format(I8,I8,I8,F15.3,F8.0) !data
open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='DTSRead:FATAL:error opening file';
    return;
endif

! First read just 2 count number of days
do i=1,nhead;read(1,*);enddo !headers
compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! read data
rewind(1)
do i=1,nhead;read(1,*);enddo !headers
do i=1,y%n
    read(1,100) y%ts(i)%date%year,y%ts(i)%date%month ,&
                y%ts(i)%date%day, y%ts(i)%q, Qcode%ts(i)%q
enddo
Qcode%ts(:)%date=y%ts(:)%date
close(1)

end subroutine oldDTSRead

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DTSRead(file,y,Qcode,err,mess)

!^**********************************************************************
!^* Purpose: Read from DTS format (Daily Time Series)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 18/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, destination file
!^* OUT
!^*    1.y, data time series
!^*    2.Qcode, quality code time series
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,QCode
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik), parameter::nhead=4,ncol=5
character(1),parameter::sep=";"
character(250), parameter::procname="DTSRead"
real(mrk),pointer::mat(:,:)

err=0;mess=''

call ReadSeparatedFile(file=file,sep=sep,nhead=nhead,ncol=ncol,y=mat,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

! start populating TS type
y%n=size(mat,1);Qcode%n=y%n
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(y%n),Qcode%ts(Qcode%n))

y%ts(:)%date%year=int(mat(:,1),mik)
y%ts(:)%date%month=int(mat(:,2),mik)
y%ts(:)%date%day=int(mat(:,3),mik)
y%ts(:)%q=real(mat(:,4),mrk)
Qcode%ts(:)%q=real(mat(:,5),mrk)

end subroutine DTSRead

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadSeparatedFile(file,sep,nhead,ncol,y,err,mess)
!^**********************************************************************
!^* Purpose: Read data in a separator format
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 12/01/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, data file
!^*    2.sep, separator
!^*    3.nhead, number of header rows
!^*    4.ncol, number of columns
!^* OUT
!^*    1.y, result matrix
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use kinds_dmsl_kit
use utilities_dmsl_kit, only:getspareunit,getNumItemsInFile
character(*),intent(in)::file,sep
integer(mik), intent(in)::nhead,ncol
real(mrk), pointer::y(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::Nrow,unt,i,j,first,last,indx
character(250)::jchar
character(15000)::buffer

!Init
err=0;mess='';
! Read Station file
call getspareunit(unt,err,mess)
if(err/=0) then;mess='ReadSeparatedFile:'//trim(mess);return;endif
open(unit=unt,file=trim(file),status='old')
call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=nhead,nitems=Nrow,postPos=0,jchar=jchar,err=err,message=mess)
if(err/=0) then;mess='ReadSeparatedFile:'//trim(mess);close(unt);return;endif
! allocate y
if(associated(y)) nullify(y);allocate(y(Nrow,Ncol))
! headers
do i=1,nhead;read(unt,*);enddo
! read data
do i=1,Nrow
    read(unt,'(A15000)') buffer
    do j=1,ncol
        indx=index(buffer,sep)
        if(indx==0) then
            read(buffer,*) y(i,j)
        else
            first=1;last=indx-1
            read(buffer(first:last),*) y(i,j)
            buffer=buffer(indx+1:)
        endif
    enddo
enddo
close(unt)
end subroutine ReadSeparatedFile
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine WriteSeparatedFile_r(file,sep,y,headers,append,err,mess)
!^**********************************************************************
!^* Purpose: Write data in a separator format
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 12/01/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, data file
!^*    2.sep, separator
!^*    3.y, result matrix
!^*    4.[headers]
!^*    5.[append], default false - if true append at end of existing file
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use kinds_dmsl_kit
use utilities_dmsl_kit, only:getspareunit
character(*),intent(in)::file,sep
character(*),intent(in),optional::headers(:)
logical,intent(in),optional::append
real(mrk), intent(in)::y(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
logical,parameter::append_def=.false.
integer(mik)::i,j,unt,nrow,ncol
logical::app

!Init
err=0;mess='';
if(present(append)) then;app=append;else;app=append_def;endif
call getspareunit(unt,err,mess)
if(err/=0) then;mess='WriteSeparatedFile_r:'//trim(mess);return;endif
if(app) then
    open(unit=unt,file=trim(file),position='append')
else
    open(unit=unt,file=trim(file),status='replace')
endif
nrow=size(y,1);ncol=size(y,2)
! headers
if(present(headers)) then
    do j=1,size(headers)
        write(unt,'(2A)',advance='NO') trim(headers(j)),trim(sep)
    enddo
    write(unt,'(A)',advance='YES') ""
endif
! read data
do i=1,Nrow
    do j=1,ncol
        write(unt,'(e14.6,A)',advance='NO') y(i,j),trim(sep)
    enddo
    write(unt,'(A)',advance='YES') ""
enddo
close(unt)
end subroutine WriteSeparatedFile_r
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WriteSeparatedFile_i(file,sep,y,headers,append,err,mess)
!^**********************************************************************
!^* Purpose: Write data in a separator format, integers
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 21/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, data file
!^*    2.sep, separator
!^*    3.y, result matrix
!^*    4.[headers]
!^*    5.[append], default false - if true append at end of existing file
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use kinds_dmsl_kit
use utilities_dmsl_kit, only:getspareunit
character(*),intent(in)::file,sep
character(*),intent(in),optional::headers(:)
logical,intent(in),optional::append
integer(mik), intent(in)::y(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
logical,parameter::append_def=.false.
integer(mik)::i,j,unt,nrow,ncol
logical::app

!Init
err=0;mess='';
if(present(append)) then;app=append;else;app=append_def;endif
call getspareunit(unt,err,mess)
if(err/=0) then;mess='WriteSeparatedFile_i:'//trim(mess);return;endif
if(app) then
    open(unit=unt,file=trim(file),position='append')
else
    open(unit=unt,file=trim(file),status='replace')
endif
nrow=size(y,1);ncol=size(y,2)
! headers
if(present(headers)) then
    do j=1,size(headers)
        write(unt,'(2A)',advance='NO') trim(headers(j)),trim(sep)
    enddo
    write(unt,'(A)',advance='YES') ""
endif
! read data
do i=1,Nrow
    do j=1,ncol
        write(unt,'(I0,A)',advance='NO') y(i,j),trim(sep)
    enddo
    write(unt,'(A)',advance='YES') ""
enddo
close(unt)
end subroutine WriteSeparatedFile_i
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DatWrite(y,file,headers,err,mess)

!^**********************************************************************
!^* Purpose: write data into Dat format (Data matrix)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 19/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.y, data matrix
!^*    2.file, data file
!^*    3.[headers] (default 'var1, var2 etc...')
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use  utilities_dmsl_kit,only: number_string

character(*), intent(in)::file
real(mrk), intent(in)::y(:,:)
character(*), intent(in), optional::headers(:)
character(*), intent(out)::mess
integer(mik), intent(out)::err
!locals
integer(mik)::i,n,p
character(250), allocatable::head(:)
character(250)::i_str,fmt

err=0;mess=''
n=size(y,dim=1)
p=size(y,dim=2)

!handle headers
if(present(headers)) then
  if(size(headers)/=p) then
    err=1;mess='DatWrite:FATAL:incorrect size(headers)';return
  endif
  if(allocated(head)) deallocate(head)
  allocate(head(p))
  head=headers
else
  if(allocated(head)) deallocate(head)
  allocate(head(p))
    do i=1,p
    write(i_str,'(i0)') i
    head(i)='var_'//trim(i_str)
  enddo
endif

open(unit=1, file=file,status='REPLACE',iostat=err)
if (err/=0) then
    mess='DatWrite:FATAL:error opening file';
    return;
endif

!headers
fmt='('//trim(number_string(p))//'(A14," "))'
write(1,trim(fmt)) head
!data
fmt='('//trim(number_string(p))//'E15.6)'
do i=1,n
    write(1,trim(fmt)) y(i,:)
enddo
close(1)
deallocate(head)

end subroutine DatWrite

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DatRead(file,ncol,y,headers,err,mess)

!^**********************************************************************
!^* Purpose: Read data into Dat format (Data matrix)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 31/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    2.file, data file
!^* OUT
!^*    1.y, data matrix
!^*    2.headers
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************

character(*), intent(in)::file
integer(mik), intent(in)::ncol
real(mrk), pointer::y(:,:)
character(*), intent(out)::headers(ncol)
character(*), intent(out)::mess
integer(mik), intent(out)::err
!locals
integer(mik)::i,n,compt,errcode

err=0;mess=''

open(unit=1, file=file,status='OLD',iostat=err)
if (err/=0) then
    mess='DatRead:FATAL:error opening file';
    return;
endif

! First read just 2 count number of days
read(1,*) !headers
compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
n=compt-1
rewind(1)

if(associated(y)) nullify(y)
allocate(y(n,ncol))
!headers
read(1,*) headers
!data
do i=1,n
    read(1,*) y(i,:)
enddo
close(1)

end subroutine DatRead

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine HodgkinsWrite(file,y,ID,err,mess)

!^**********************************************************************
!^* Purpose: csv-like format requested by Glenn Hodkins
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/06/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.y, data matrix
!^*    2.file, data file
!^*    3.ID
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use Dates_tools, only:FormatDate

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(in)::y
character(*), intent(in)::ID
character(*), intent(out)::mess
integer(mik), intent(out)::err
!locals
character(250),parameter::header='Agency,ID,Day,Month,Year,Flow,Flag'
integer(mik)::i,flag

err=0;mess=''

!100 format(A,',',A,',',A8,',',F15.3,',',I1) !data
100 format(A,',',A,',',I2.2,',',I2.2,',',I4.4,',',F15.3,',',I1) !data

open(unit=1, file=file,status='REPLACE',iostat=err)
if (err/=0) then
    mess='HodgkinsWrite:FATAL:error opening file';
    return;
endif

!headers
write(1,'(A150)') header
!data
do i=1,y%n
    ! call FormatDate(y%ts(i)%date,5,fDate,err,mess)
    !if (err/=0) then
    !    mess='HodgkinsWrite:FATAL:error converting date';
    !    return;
    !endif
    if(y%ts(i)%q/=y%mv) then
        flag=1
    else
        flag=0
    endif        
    !write(1,100) 'HYDRO',trim(ID),trim(fdate),y%ts(i)%q,flag
    write(1,100) 'HYDRO',trim(ID),y%ts(i)%date%day,y%ts(i)%date%month,&
                y%ts(i)%date%year,y%ts(i)%q,flag
enddo
close(1)

end subroutine HodgkinsWrite

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine AdaptAlpRead(y,Qcode,file,err,mess)

!^**********************************************************************
!^* Purpose: Reader for AdaptAlp format defined by Antoine Bard
!^**********************************************************************
!^* Programmer: Antoine Bard, Irstea Lyon - slightly modified by Ben Renard
!^**********************************************************************
!^* Last modified: Bring to this module on 9/10/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, data file
!^* OUT
!^*    1.y, time series
!^*    1.Qcode, quality codes - nor useful for AdaptAlp format, kept for consistency
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,Qcode
integer(mik),intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik):: i,k,n,errcode
character line*80,deb*50

err=0;mess='';n=0;errcode=0

OPEN(1,status='OLD',FILE=trim(file), IOSTAT=errcode, ERR=99)
read(1,'(a80)')
read(1,'(a80)')
read(1,'(a80)')
!get data dimention 
do while((errcode==0))
read(1,*,iostat=errcode)
n=n+1
end do
CLOSE(1)
! Define attributes of series
if(associated(y%ts)) nullify(y%ts)
allocate(y%ts(n-1))
y%n=n-1
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(Qcode%ts(n-1))
Qcode%n=n-1

OPEN(1,status='OLD',FILE=trim(file), IOSTAT=errcode, ERR=98)
do i=1,3;read(1,*);end do
errcode=0;k=1;
do while ((errcode==0).and.(k<=y%n))
    READ(1,'(A80)',iostat=errcode,ERR=97) line
    read(line(1:4),'(I4)') y%ts(k)%date%year
    read(line(6:7),'(I2)') y%ts(k)%date%month
    read(line(9:10),'(I2)') y%ts(k)%date%day
    deb=line(12:80)
    if (len_trim(deb)==0)then
        y%ts(k)%q=y%mv
    else
        Read(deb, '(F20.3)') y%ts(k)%q
        if (y%ts(k)%q < 0.)then
            y%ts(k)%q=y%mv
        endif
    endif
    k=k+1
end do
CLOSE(1)
Qcode%ts(:)%date=y%ts(:)%date
Qcode%ts(:)%q=Qcode%mv
return
97 mess='AdaptAlp error reading data'
err=1;return
98 mess='AdaptAlp error reading station file - 2nd opening'
err=1;return
99 mess='AdaptAlp error reading station file - 1st opening '
err=1;return

end subroutine AdaptAlpRead

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadRHN(file,y,Qcode,err,mess)

!^**********************************************************************
!^* Purpose: Read from RHNRead format (Daily Q, as defined by Paul Whitfield)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 18/06/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, source file
!^* OUT
!^*    1.y, data time series
!^*    2.Qcode, quality code time series
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,QCode
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i, errcode,compt
character(250)::crap1, crap2, crap3, crap4

err=0;mess=''

open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='RHNRead:FATAL:error opening file';
    return;
endif

! First read just 2 count number of days
read(1,*) !headers
compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! read data
rewind(1)
read(1,*) !headers
do i=1,y%n
    read(1,*) crap1,crap2,crap3,y%ts(i)%date%year,y%ts(i)%date%month,&
                y%ts(i)%date%day, y%ts(i)%q, crap4
    if(y%ts(i)%q<0._mrk) y%ts(i)%q=y%mv
enddo
Qcode%ts(:)%date=y%ts(:)%date
close(1)

end subroutine ReadRHN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ReadRHN4(file,y,Qcode,err,mess)

!^**********************************************************************
!^* Purpose: Read from RHN4 format (Daily Q, as defined by Paul Whitfield)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 20/04/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, source file
!^* OUT
!^*    1.y, data time series
!^*    2.Qcode, quality code time series
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************

character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,QCode
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i, errcode,compt
character(250)::crap1, crap2, crap3, crap4, crap5, crap6,foo1,foo2,foo3,foo4

err=0;mess=''

open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='RHNRead:FATAL:error opening file';
    return;
endif

! First read just 2 count number of days
read(1,*) !headers
compt=0;errcode=0
do while(errcode==0) !just count number of days
    compt=compt+1
    read(1,*,iostat=errcode)
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))

! read data
rewind(1)
read(1,*) !headers
do i=1,y%n
    read(1,*) crap1,crap2,crap3,crap4,foo1,foo2,&
                foo3,crap5,foo4,crap6
    if(trim(foo4)/='NA') then
        read(foo4,*) y%ts(i)%q
        if(y%ts(i)%q<0._mrk) y%ts(i)%q=y%mv
    else
        y%ts(i)%q=y%mv
    endif
    read(foo1,*,iostat=err) y%ts(i)%date%year
    if(err/=0) then;mess='ReadRHN4: read error';return;endif
    read(foo2,*,iostat=err) y%ts(i)%date%month
    if(err/=0) then;mess='ReadRHN4: read error';return;endif
    read(foo3,*,iostat=err) y%ts(i)%date%day
    if(err/=0) then;mess='ReadRHN4: read error';return;endif
enddo
Qcode%ts(:)%date=y%ts(:)%date
close(1)

end subroutine ReadRHN4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadMalekiFormat(file,y,Qcode,err,mess)

!^**********************************************************************
!^* Purpose: read daily runoff/rainfall data in Maleki's format
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 09/08/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^* OUT
!^*    1.y, 1D-time series of daily data
!^*    2.Qcode, 1D-time series of quality code
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
character(*), intent(in)::file
type(OneDTimeSeriesType), intent(out)::y,Qcode
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik), parameter::nHeaders=0

integer(mik)::errcode, compt,j, np0,np1, np2
character(250):: line,cdate,cvalue
logical::isempty

err=0;mess=''
open(unit=1,file=file, status='OLD',iostat=errcode)
! handle OPEN error
if (errcode/=0) then
    err=errcode
    mess='ReadMalekiFormat:FATAL:error opening file';
    return
endif

! First find nb of days
do j=1, nHeaders !headers
    read(1,*)
enddo
compt=0;errcode=0;isempty=.false.
do while(errcode==0 .and. (.not. isempty)) !just count number of days
    compt=compt+1
    read(1,'(A)',iostat=errcode) line
    isempty=(trim(line)=='') .or. (trim(line)==char(9))
enddo
compt=compt-1

! start populating TS type
y%n=compt;Qcode%n=compt;y%mv=-9999._mrk
if(associated(y%ts)) nullify(y%ts)
if(associated(Qcode%ts)) nullify(Qcode%ts)
allocate(y%ts(compt),Qcode%ts(compt))
y%ts(:)%q=UndefRN

rewind(1)
! read data
do j=1,y%n
    READ(1,'(A)') line
    np0=index(line, char(9))
    cdate=line(1:(np0-1))
    cvalue=line((np0+1):)
    ! work out date
    np1=index(cdate, '.')
    read(cdate(1:(np1-1)),*) y%ts(j)%date%day
    np2=index(cdate(np1+1:), '.')
    read(cdate((np1+1):(np1+np2-1)),*) y%ts(j)%date%month
    read(cdate((np1+np2+1):),*) y%ts(j)%date%year
    
    if( (trim(cvalue)=='') .or. (trim(cvalue)==char(9))) then
        y%ts(j)%q=y%mv
    else
        read(cvalue,*) y%ts(j)%q
    endif
            
    
    !    read(line((np1+1):), '(A)') cvalue
    
!
!    !check for daily mean discharge code= QJO
!    read(line(1:3), '(A3)') code
!    if (code/='QJO') then
!        err=1;mess='ReadHYDROII:Fatal:QJOexpected'
!        return
!    end if
!    read(line(14:17), '(I4)') y%ts(j)%date%year
!    read(line(18:19), '(I2)') y%ts(j)%date%month
!    read(line(20:21), '(I2)') y%ts(j)%date%day
!    line=line(23:80)
!    np1=index(line, sep)
!    if(np1<=1) then
!        y%ts(j)%q=y%mv
!    else
!        read(line(1:np1-1), '(A)') deb
!        read(deb, '(F20.3)') y%ts(j)%q
!        y%ts(j)%q=0.001_mrk*y%ts(j)%q !in  m3.s-1
!    end if
enddo

!Qcode
Qcode%ts(:)%date=y%ts(:)%date
Qcode%ts(:)%q=y%mv
close(1)

end subroutine ReadMalekiFormat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadXceedances_GlennFormat(file,M,k,n,StationList,YearList,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 28/08/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^* OUT
!^*    1.M, 0/1 matrix (missing = -99)
!^*    2.k, time series of number of successes
!^*    3.n, time series of number of trials
!^*    4.StationList, 
!^*    5.YearList, 
!^*    6.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    7.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:countSubstringInString,replacedString,getNumItemsInFile
character(*), intent(in)::file
integer(mik),pointer::M(:,:),k(:),n(:),YearList(:)
character(*),pointer::StationList(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
character(1), parameter::sep=','
integer(mik), parameter::nskip=1,mvcode=-99
integer(mik)::i,j,Nstation,Nt
character(20000)::line
character(250)::craps
character(50),allocatable::foo(:)

err=0;mess='';
! Get number of stations
open(unit=1,file=trim(file),status='old')
read(1,'(A)') line
Nstation=countSubstringInString(trim(line),sep)
! read station names
if(associated(StationList)) nullify(StationList);allocate(StationList(Nstation))
read(line,*) craps,StationList
do i=1,Nstation
    StationList(i)='Floods.'//replacedString(StationList(i),'-','.')//'.csv'
enddo
! Get number of years
call getNumItemsInFile(unt=1,preRwnd=.true.,nskip=nskip,nitems=Nt,postPos=nskip,jchar=craps,err=err,message=mess)
if(err>0) then
    mess="ReadXceedances_GlennFormat: "//trim(mess);return
endif
! Read data
if(associated(M)) nullify(M);allocate(M(Nt,Nstation))
if(allocated(foo)) deallocate(foo);allocate(foo(Nstation))
if(associated(k)) nullify(k);allocate(k(Nt))
if(associated(n)) nullify(n);allocate(n(Nt))
if(associated(YearList)) nullify(YearList);allocate(YearList(Nt))

do i=1,Nt
    read(1,*) YearList(i),foo
    do j=1,Nstation
        select case(foo(j))
        case('0');M(i,j)=0
        case('1');M(i,j)=1
        case('NA');M(i,j)=mvcode
        case default
            err=1;mess='ReadXceedances_GlennFormat:fatal:unrecognized item';return
        end select        
    enddo
    n(i)=count(M(i,:)>=0)
    k(i)=count(M(i,:)==1)
enddo
close(1)

end subroutine ReadXceedances_GlennFormat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoadFromRepository(Repo,StationListFile,Extension,ReadSub,y,Qcode,err,mess)

!^**********************************************************************
!^* Purpose: Load datafiles from a repository
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 16/11/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.Repo, full path 2 repository
!^*    2.StationListFile, file containing the name of the files to read in the repository
!^*    1.Extension, extension to add to all files
!^*    1.ReadSub, reading subroutine - see interface below
!^* OUT
!^*    1.y, vector of time series containing data [dim: number of files in StationListFile] 
!^*    2.Qcode, vector of quality code time series
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:getspareunit,getNumItemsInFile
character(*), intent(in)::Repo,StationListFile,Extension
type(OneDTimeSeriesType), pointer::y(:),Qcode(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
interface
    subroutine ReadSub(file,y,Qcode,err,mess)
        use kinds_dmsl_kit
        use TimeSeries_tools, only:OneDTimeSeriesType
        character(*),intent(in)::file
        type(OneDTimeSeriesType), intent(out)::y,Qcode
        integer(mik), intent(out)::err
        character(*),intent(out)::mess
    end subroutine ReadSub
end interface
!locals
integer(mik)::unt,Nx,i
character(250)::jchar
character(250),allocatable::StationList(:)

!Init
err=0;mess='';
! Read Station file
call getspareunit(unt,err,mess)
if(err/=0) then;mess='LoadFromRepository:'//trim(mess);return;endif
open(unit=unt,file=trim(StationListFile),status='old')
call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=0,nitems=Nx,postPos=0,jchar=jchar,err=err,message=mess)
if(err/=0) then;mess='LoadFromRepository:'//trim(mess);return;endif
! Allocate stuff & read StationList
if(allocated(StationList)) deallocate(StationList);allocate(StationList(Nx))
if(associated(y)) nullify(y);allocate(y(Nx))
if(associated(Qcode)) nullify(Qcode);allocate(Qcode(Nx))
do i=1,Nx; read(unt,'(A)') StationList(i);enddo
close(unt)
! Start Reading repository
do i=1,Nx
    call ReadSub(file=trim(Repo)//trim(StationList(i))//trim(Extension),y=y(i),Qcode=Qcode(i),err=err,mess=mess)
    if(err/=0) then;mess='LoadFromRepository:'//trim(mess);return;endif
    y(i)%NameSeries=trim(StationList(i))
    Qcode(i)%NameSeries=trim(StationList(i))
enddo

end subroutine LoadFromRepository


end module DataRW_tools

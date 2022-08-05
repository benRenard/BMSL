module TimeSeries_tools

!~**********************************************************************
!~* Purpose: Define the "TimeSeries" Type and implement associated tools
!~**********************************************************************
!~* Programmer: Ben Renard & Antoine Bard, Cemagref Lyon
!~**********************************************************************
!~* Last modified:28/07/2011
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1. TSFiltering_daily: various filters (eg moving average, BFS, SPA, etc.)
!~*    2. FillTheGap_daily: fill-incomplete daily series, by replacing non-present dates by the date+a missing value code.
!~*    3. Hourly2Daily: transform an hourly TS into a daily one.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use Dates_tools, only:DateType
implicit none
Private
public :: & ! Aggregation
            TSAggregate_preprocess, TSAggregate,Hourly2Daily,Daily2Monthly,&
            GetMonthlyStat, GetSeasonStat,GetMonthlyInterannualStat,RemoveMonthlyClimato,&
            ! filtering
            TSFiltering_daily, FillTheGap_daily,&
            ! misc. utilities
            GetValueFromDate,SetValueToDate,ExtractOneDTimeSeries,GatherMultiDTimeSeries,&
            GatherMultiDTimeSeries_Matrix

integer(mik), parameter, PUBLIC:: & ! Defines names of the filters
            MovingAverage=1, MovingMin=2, MovingMax=3,&
            BFS=4,& !BaseFlow Separation, see Tallaksen and Van Lanen 2004
            SPA=5 !Secant Peak Algorithm, see Tallaksen and Van Lanen 2004

! defaults
integer(mik), parameter::qc_def=QC_PRESENT,cc_def=QC_PRESENT
real(mrk), parameter:: mv_def=-99._mrk

! Definition of 1-D and Multi-D time series types

! Intermediate types definition TS1D and TSMultiD
type, public:: TS1D
    type(DateType)::date ! t
    real(mrk) :: q=undefRN ! q(t) - 1D
    integer(mik)::qc=qc_def ! quality code (t)
    integer(mik)::cc=cc_def ! continuity code (t)
end type TS1D

type, public:: TSMultiD
    type(DateType) :: date ! t
    real(mrk), pointer :: q(:) => NULL() ! q(t) - Multidimensional
    integer(mik), pointer :: qc(:) => NULL() ! quality code(t) - Multidimensional
    integer(mik), pointer :: cc(:) => NULL() ! continuity code(t) - Multidimensional
end type TSMultiD

type, public :: OneDTimeSeriesType ! 1-Dimensional
    character(len=250)::NameSeries
    type(TS1D), pointer::ts(:) => NULL() ! time series (date(i),q(i)), i=1...n
    integer(mik) :: n=undefIN ! time series length
    real(mrk) :: mv=mv_def ! missing values code
end type OneDTimeSeriesType

type, public :: TimeSeriesType ! Most general type, multidimensional
    character(len=250)::NameSeries
    type(TSMultiD), pointer::ts(:) => NULL() ! time series (date(i),q(i,k)), i=1...n,k=1...Nvariables
    integer(mik) :: n=UndefIN ! time series length
    real(mrk) :: mv=mv_def ! missing values code
end type TimeSeriesType

! Overloads
interface TSAggregate
  module procedure TSAggregate_engine,TSAggregate_ez
endinterface TSAggregate

!!!!!!!!!!
Contains
!!!!!!!!!!

!-----------------------------
! Aggregation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TSAggregate_preprocess(TSin,scale,scaleMult,&
                        DateListOption,FirstDate,LastDate,&
                        x,y,qc,cc,grid,err,message)

!^**********************************************************************
!^* Purpose: Pre-processing before aggregation; in particular, returns
!^*          x and grid in absolute time format, and quality/continuity
!^*          codes and y in array format.
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:11/05/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. TSin, input time series
!^*    2. scale of the grid (one of YMDhms)
!^*    3. [scaleMult] (e.g. for 15 minutes step, scaleMult = 15, scale = 'm')
!^*    4. [DateListOption]
!^*    5. [FirstDate] beginning of the grid; is absent set to first date in TSin
!^*    6. [LastDate] end of the grid; is absent set to last date in TSin
!^* OUT
!^*    1. x, dates of the time series in absolute time
!^*    2. y, values of the time series
!^*    3. qc, quality codes
!^*    4. cc, continuity codes
!^*    5. grid, target grid on which aggregation will be performed (in absolute time)
!^*    6.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    7.mess, error message
!^**********************************************************************
use Dates_tools, only:DateList,AbsoluteTime
use Aggregation_tools, only:CC_CONT

type(OneDTimeSeriesType), intent(in)::TSin
character(*), intent(in)::scale
integer(mik), intent(in),optional::scaleMult
integer(mik),intent(in),optional::DateListOption
type(DateType),intent(in),optional::FirstDate,LastDate
real(mrk),pointer::x(:),y(:),grid(:)
integer(mik), pointer::qc(:),cc(:)
integer(mik), intent(out)::err
character(100),intent(out)::message
! locals
character(*),parameter::procnam='TSAggregate_preprocess'
integer(mik)::i,n,ngrid
type(DateType)::d0,dn
type(DateType), pointer::list(:)

err=EXIT_SUCCESS; message=procnam//'/ok';
n=TSin%n
if(present(FirstDate)) then; d0=FirstDate; else; d0=TSin%ts(1)%date; endif
if(present(LastDate)) then; dn=LastDate; else; dn=TSin%ts(n)%date; endif
! create grid
call DateList(d0,dn,scale,scaleMult,list,err,message)
if(err/=EXIT_SUCCESS)then;message='f-'//procnam//'&'//trim(message);return;endif
ngrid=size(list)
! fill in output arrays
allocate(x(n),y(n),qc(n),cc(n),grid(ngrid))
x=AbsoluteTime(TSin%ts(:)%date,DateListOption,d0)
y=TSin%ts(:)%q
do i=1,n
    qc(i)=TSin%ts(i)%qc
    if(y(i)==TSin%mv) qc(i)=QC_MISSVAL
    cc(i)=TSin%ts(i)%cc
enddo
grid=AbsoluteTime(list,DateListOption,d0)
deallocate(list)
end subroutine TSAggregate_preprocess
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TSAggregate_engine(operateur,opID,TSin,scale,scaleMult,&
                       aux,DateListOption,FirstDate,LastDate,&
                       TSout,tx,tgrid,err,message)

!^**********************************************************************
!^* Purpose: Generic aggregator of an input time series using a regular grid
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:10/05/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: handle vectorial operators
!^**********************************************************************
!^* IN
!^*    1. operateur, interface to the operator to be applied
!^*    2. opID, properties of the operator
!^*    3. TSin, input time series
!^*    4. scale of the grid (one of YMDhms)
!^*    5. [scaleMult] (e.g. for 15 minutes step, scaleMult = 15, scale = 'm')
!^*    6. [aux] auxiliaries of the operator
!^*    7. [DateListOption]
!^*    8. [FirstDate] beginning of the grid; is absent set to first date in TSin
!^*    9. [LastDate] end of the grid; is absent set to last date in TSin
!^* OUT
!^*    1. TSout, output aggregated time series
!^*    2. [tx], dates of the time series in absolute time
!^*    3. [tgrid], target grid on which aggregation will be performed (in absolute time)
!^*    4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5. message, error message
!^**********************************************************************
use Dates_tools, only:DateList,InXDays
use Aggregation_tools, only:aggregator,CC_CONT

interface
  pure subroutine operateur(opID,y,x,xmin,xmax,qc,cc,aux,&
                            yOut,xOut,qcOut,ccOut,err,message)
  use kinds_dmsl_kit
  implicit none
  integer(mik),intent(in)::opID(:)
  real(mrk),intent(in)::y(:),x(:),xmin,xmax
  integer(mik),intent(in),optional::qc(:),cc(:)
  real(mrk),intent(in),optional::aux(:)
  real(mrk),intent(out)::yout(:),xout
  integer(mik),intent(out),optional::qcOut(:),ccOut(:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine operateur
endinterface
integer(mik),intent(in)::opID(:)
type(OneDTimeSeriesType), intent(in)::TSin
character(*), intent(in)::scale
integer(mik), intent(in),optional::scaleMult
real(mrk),intent(in),optional::aux(:)
integer(mik),intent(in),optional::DateListOption
type(DateType), intent(in), optional::FirstDate,LastDate
type(OneDTimeSeriesType), intent(out)::TSout
real(mrk),pointer,optional::tgrid(:),tx(:)
integer(mik), intent(out)::err
character(*),intent(out)::message
! locals
integer(mik), parameter::opDim=1 ! hard-coded dimensionality of the operator - need to think how to handle it properly
character(*),parameter::procnam='TSAggregate_engine'
integer(mik)::i,ngrid
real(mrk), pointer::x(:),y(:),grid(:)
real(mrk), allocatable::xOut(:),yOut(:,:)
integer(mik), pointer::qc(:),cc(:)
integer(mik), allocatable::qcOut(:,:),ccOut(:,:)
type(DateType)::d0

err=EXIT_SUCCESS; message=procnam//'/ok';
if(present(FirstDate)) then; d0=FirstDate; else; d0=TSin%ts(1)%date; endif
! pre-processor
call TSAggregate_preprocess(TSin,scale,scaleMult,DateListOption,FirstDate,LastDate,&
                       x,y,qc,cc,grid,err,message)
if(err/=EXIT_SUCCESS)then;message='f-'//procnam//'&'//trim(message);return;endif
! return continuous-time tx and tgrid if requested
if(present(tx)) then;allocate(tx(TSin%n));tx=x;endif;
if(present(tgrid)) then;allocate(tgrid(size(grid)));tgrid=grid;endif;
! Aggregate
ngrid=size(grid)
allocate(xOut(ngrid-1),yOut(ngrid-1,opDim),qcOut(ngrid-1,opDim),ccOut(ngrid-1,opDim))
call aggregator(operateur,opID,y,x,grid,qc,cc,aux,& !IN
                yOut,xOut,qcOut,ccOut,err,message) !OUT
if(err/=EXIT_SUCCESS)then;message='f-'//procnam//'&'//trim(message);return;endif
! pack back to TSout
TSout%n=ngrid-1
TSout%NameSeries=trim(TSin%NameSeries)//'_aggregated'
TSout%mv=TSin%mv
allocate(TSout%ts(TSout%n))
do i=1,TSout%n
    TSout%ts(i)%date=InXDays(d0,xOut(i))
    TSout%ts(i)%q=yOut(i,opDim)
    TSout%ts(i)%qc=qcOut(i,opDim)
    if(qcOut(i,opDim)/=QC_PRESENT) TSout%ts(i)%q=TSout%mv
    TSout%ts(i)%cc=ccOut(i,opDim)
enddo
deallocate(x,y,grid,xOut,yOut,qc,cc,qcOut,ccOut)
end subroutine TSAggregate_engine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine TSAggregate_ez(operateur,qcOption,ccOption,TSin,scale,scaleMult,&
                        aux,DateListOption,FirstDate,LastDate,&
                        TSout,tx,tgrid,err,message)

!^**********************************************************************
!^* Purpose: easy version with default options
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon, suggested by the Big-DK, UoA
!^**********************************************************************
!^* Last modified:11/05/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: handle vectorial operators
!^**********************************************************************
!^* IN
!^*    1. operateur, as a string!
!^*    2. [qcOption] - see Aggregation_tools
!^*    3. [ccOption] - see Aggregation_tools
!^*    4. TSin, input time series
!^*    5. scale of the grid (one of YMDhms)
!^*    6. [scaleMult] (e.g. for 15 minutes step, scaleMult = 15, scale = 'm')
!^*    7. [aux] auxiliaries of the operator
!^*    8. [DateListOption]
!^*    9. [FirstDate] beginning of the grid; is absent set to first date in TSin
!^*    10. [LastDate] end of the grid; is absent set to last date in TSin
!^* OUT
!^*    1. TSout, output aggregated time series
!^*    2. [tx], dates of the time series in absolute time
!^*    3. [tgrid], target grid on which aggregation will be performed (in absolute time)
!^*    4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5. message, error message
!^**********************************************************************
use Aggregation_tools, only:basic_op,LINEAR_interpol,INFECT_qc,CONSERV_cc,AVE_funk,&
                        SUM_funk,MIN_funk,MAX_funk,ix_interpol,ix_funk,ix_qc,ix_cc
use utilities_dmsl_kit, only:changeCase

character(*),intent(in)::operateur
integer(mik),intent(in),optional ::qcOption,ccOption
type(OneDTimeSeriesType), intent(in)::TSin
character(*), intent(in)::scale
integer(mik), intent(in),optional::scaleMult
real(mrk),intent(in),optional::aux(:)
integer(mik),intent(in),optional::DateListOption
type(DateType), intent(in), optional::FirstDate,LastDate
type(OneDTimeSeriesType), intent(out)::TSout
real(mrk),pointer,optional::tgrid(:),tx(:)
integer(mik), intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam='TSAggregate_ez'
integer(mik), parameter::interpol=LINEAR_interpol
integer(mik), parameter::qc_def=INFECT_qc, cc_def=CONSERV_cc
integer(mik)::opID(4),qc,cc

err=EXIT_SUCCESS; message=procnam//'/ok';
if(present(qcOption)) then; qc=qcOption; else; qc=qc_def;endif
if(present(ccOption)) then; cc=ccOption; else; cc=cc_def;endif
opID(ix_interpol)=interpol;opID(ix_qc)=qc;opID(ix_cc)=cc;
selectcase(changeCase(operateur,'l'))
case('sum','total','integral')
    opID(ix_funk)=SUM_funk
case('ave','average','mean')
    opID(ix_funk)=AVE_funk
case('min','mini','minimum')
    opID(ix_funk)=MIN_funk
case('max','maxi','maximum')
    opID(ix_funk)=MAX_funk
case default
    err=EXIT_FAILURE;message='f-'//procnam//'/unknown[operateur='//trim(operateur)//']';return
endselect

call TSAggregate_engine(basic_op,opID,TSin,scale,scaleMult,aux,DateListOption,FirstDate,LastDate,&
                       TSout,tx,tgrid,err,message)
if(err/=EXIT_SUCCESS)then;message='f-'//procnam//'&'//trim(message);return;endif

! pre-processor
end subroutine TSAggregate_ez
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Hourly2Daily(xin,DayStart,h2dFunk,xout,pmv,err,mess)

!^**********************************************************************
!^* Purpose: transform an hourly TS into a daily one, by summing OR
!^*          averaging OR MAXing OR MINing hourly values
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 28/07/2011
!^**********************************************************************
!^* Comments: By convention, the day containing hours from
!^*           (day j, hour DayStart) to (day j+1, hour DayStart-1)
!^*           is day j in the daily time series.
!^*           Be careful with missing values: this sub will compute the
!^*           daily mean/sum/max/min even if there is a single hourly
!^*           value within the day! It's your responsability to decide
!^*           what to do with missing value, using pmv
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin, hourly time series
!^*    2. <DayStart>, integer, first hour of the day; default 0
!^*    3. <h2dFunk>, character, which function to transform from hourly
!^*       to daily? default 'mean', available 'sum', 'max', 'min'
!^* OUT
!^*    1.xout, daily time series
!^*    2.pmv, daily time series containing the percentage of missing
!^*           values within each day
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use Dates_tools, only:DaysBetween, IsInSeason, yesterday,tomorrow,OneHourAgo

type(OneDTimeSeriesType), intent(in)::xin
integer(mik), optional, intent(in)::DayStart
character(*), optional, intent(in)::h2dFunk
type(OneDTimeSeriesType), intent(out)::xout,pmv
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::DS,nH,nD,i,ndata,compt,up
character(250)::h2d
type(DateType)::currentday
logical, allocatable::mask1(:),mask2(:)
type(OneDTimeSeriesType)::QJ

err=0;mess=''

! Handle default values for optional inputs !
if(present(DayStart)) then
    DS=DayStart;
else
    DS=0
endif
if(present(h2dFunk)) then
    h2d=h2dFunk;
else
    h2d='mean'
endif

! Compute number of days and start populating xout & pmv
nH=size(xin%ts)
nD=ceiling(DaysBetween(xin%ts(1)%date,xin%ts(nH)%date))+2
xout%n=nD
pmv%n=nD
if(associated(xout%ts)) nullify(xout%ts)
allocate(xout%ts(nD))
if(associated(pmv%ts)) nullify(pmv%ts)
allocate(pmv%ts(nD))
xout%mv=xin%mv
pmv%mv=xin%mv
xout%NameSeries=trim(xin%NameSeries)//'_daily'//trim(h2D)
pmv%NameSeries=trim(xin%NameSeries)//'_daily'//trim(h2D)//'_pmv'


currentday=yesterday(xin%ts(1)%date)
compt=1
do i=1,nD
    pmv%ts(i)%date=DateType(currentday%year,currentday%month,currentday%day,0,0,0)
    xout%ts(i)%date=DateType(currentday%year,currentday%month,currentday%day,0,0,0)
    up=min(compt+23,nH)
    ! extract non-missing hourly values within current day
    if(allocated(mask1)) deallocate(mask1)
    allocate(mask1(up-compt+1))
    if(allocated(mask2)) deallocate(mask2)
    allocate(mask2(up-compt+1))
    mask1=IsInSeason(date=xin%ts(compt:up)%date,&
                    Sstart=DateType(currentday%year,currentday%month,currentday%day,DS,0,0),&
                    Send=OneHourAgo(tomorrow(DateType(currentday%year,currentday%month,currentday%day,DS,0,0))),&
                    year=currentday%year)
    mask2=(xin%ts(compt:up)%q/=xin%mv)
    ndata=count(mask1.and.mask2)

    if(ndata==0) then ! no data this day
        pmv%ts(i)%q=100._mrk
        xout%ts(i)%q=xout%mv
    else ! At least 1 hourly value - proceed to computation
        pmv%ts(i)%q=100._mrk*(1._mrk-real(ndata,mrk)/24._mrk)
!        if (compt+23>nH) then ! reach end of hourly series
!            if(associated(QJ%ts)) nullify(QJ%ts)
!            ndata=count(mask1 .and. mask2 .and. (/(compt+j,j=0,23)/)<=nH)
!            allocate(QJ%ts(ndata))
!            QJ%ts=pack(xin%ts(compt:(compt+23)),mask1 .and. mask2 .and. (/(compt+j,j=0,23)/)<=nH) !QJ= hourly data for the current day
!        else
            if(associated(QJ%ts)) nullify(QJ%ts)
            allocate(QJ%ts(ndata))
            QJ%ts=pack(xin%ts(compt:up),mask1.and.mask2) !QJ= hourly data for the current day
!        endif
        selectcase(h2d)
        case('mean') !default
            xout%ts(i)%q=sum(QJ%ts(:)%q)/ndata
        case('sum')
            xout%ts(i)%q=sum(QJ%ts(:)%q)
        case('min')
            xout%ts(i)%q=minval(QJ%ts(:)%q)
        case('max')
            xout%ts(i)%q=maxval(QJ%ts(:)%q)
        case default
            err=1;mess='Hourly2Daily:FATAL:Unavailable [h2dFunk]';return
        end select
    endif
    currentday=tomorrow(currentday)
    compt=compt+count(mask1)
enddo

end subroutine Hourly2Daily
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Daily2Monthly(xin,d2mFunk,xout,pmv,err,mess)

!^**********************************************************************
!^* Purpose: transform an daily TS into a monthly one, by summing OR
!^*          averaging OR MAXing OR MINing daily values
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 28/07/2011
!^**********************************************************************
!^* Comments: Be careful with missing values: this sub will compute the
!^*           monthly mean/sum/max/min even if there is a single daily
!^*           value within the month! It's your responsability to decide
!^*           what to do with missing value, using pmv
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin, hourly time series
!^*    2. <d2mFunk>, character, which function to transform from daily
!^*       to monthly? default 'mean', available 'sum', 'max', 'min'
!^* OUT
!^*    1.xout, daily time series
!^*    2.pmv, daily time series containing the percentage of missing
!^*           values within each day
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use Dates_tools, only:DaysBetween, DaysInMonth

type(OneDTimeSeriesType), intent(in)::xin
character(*), optional, intent(in)::d2mFunk
type(OneDTimeSeriesType), intent(out)::xout,pmv
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::nM,nD,i,ndata,compt,up
character(250)::d2m
integer(mik)::currentmonth, currentyear
logical, allocatable::mask1(:),mask2(:)
type(OneDTimeSeriesType)::QM

err=0;mess=''

! Handle default values for optional inputs !
if(present(d2mFunk)) then
    d2m=d2mFunk;
else
    d2m='mean'
endif

! Compute number of days and start populating xout & pmv
nD=size(xin%ts)
nM=ceiling(real(ceiling(DaysBetween(xin%ts(1)%date,xin%ts(nD)%date)),mrk)/30.5_mrk)+2
xout%n=nM
pmv%n=nM
if(associated(xout%ts)) nullify(xout%ts)
allocate(xout%ts(nM))
if(associated(pmv%ts)) nullify(pmv%ts)
allocate(pmv%ts(nM))
xout%mv=xin%mv
pmv%mv=xin%mv
xout%NameSeries=trim(xin%NameSeries)//'_monthly'//trim(d2m)
pmv%NameSeries=trim(xin%NameSeries)//'_monthly'//trim(d2m)//'_pmv'


currentmonth=xin%ts(1)%date%month
currentyear=xin%ts(1)%date%year
compt=1
do i=1,nM
    pmv%ts(i)%date=DateType(currentyear,currentmonth,1,0,0,0)
    xout%ts(i)%date=DateType(currentyear,currentmonth,1,0,0,0)
    up=min(compt+31,nD)
    ! extract non-missing daily values within current month
    if(allocated(mask1)) deallocate(mask1)
    allocate(mask1(up-compt+1))
    if(allocated(mask2)) deallocate(mask2)
    allocate(mask2(up-compt+1))
    mask1= (xin%ts(compt:up)%date%year==currentyear) .and. (xin%ts(compt:up)%date%month==currentmonth)
    !IsInSeason(date=xin%ts(compt:up)%date,&
    !                Sstart=DateType(currentday%year,currentday%month,currentday%day,DS,0,0),&
    !                Send=OneHourAgo(tomorrow(DateType(currentday%year,currentday%month,currentday%day,DS,0,0))),&
    !                year=currentday%year)
    mask2=(xin%ts(compt:up)%q/=xin%mv)
    ndata=count(mask1.and.mask2)

    if(ndata==0) then ! no data this month
        pmv%ts(i)%q=100._mrk
        xout%ts(i)%q=xout%mv
    else ! At least 1 daily value - proceed to computation
        pmv%ts(i)%q=100._mrk*( 1._mrk-real(ndata,mrk)/real(DaysInMonth(currentmonth,currentyear),mrk) )
        if(associated(QM%ts)) nullify(QM%ts)
            allocate(QM%ts(ndata))
            QM%ts=pack(xin%ts(compt:up),mask1.and.mask2) !QM= daily data for the current month
        selectcase(d2m)
        case('mean') !default
            xout%ts(i)%q=sum(QM%ts(:)%q)/ndata
        case('sum')
            xout%ts(i)%q=sum(QM%ts(:)%q)
        case('min')
            xout%ts(i)%q=minval(QM%ts(:)%q)
        case('max')
            xout%ts(i)%q=maxval(QM%ts(:)%q)
        case default
            err=1;mess='Daily2monthly:FATAL:Unavailable [d2mFunk]';return
        end select
    endif
    if(currentmonth==12) then
        currentmonth=1
        currentyear=currentyear+1
    else
        currentmonth=currentmonth+1
    endif
    compt=compt+count(mask1)
enddo

end subroutine Daily2Monthly
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetMonthlyClimato(X,month,dates,mu,std,err,mess)
!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:12/11/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. X, big matrix
!^*    2. month
!^*    3. dates
!^* OUT
!^*    1. mu, monthly means
!^*    2. std, monthly std
!^*    3. err
!^*    4. mess
!^**********************************************************************
use Dates_tools
use numerix_dmsl_kit, only:getmeanvar
type(DateType),intent(in)::dates(:)
real(mrk), intent(in)::X(:,:)
integer(mik), intent(in)::month
real(mrk),intent(out)::mu(size(X,dim=2)),std(size(X,dim=2))
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::i,n
logical::mask(size(dates))
real(mrk), allocatable::foo(:)

err=0;mess=''
mask=dates(:)%month==month
n=count(mask)
do i=1,size(X,dim=2)
    if(allocated(foo)) deallocate(foo);allocate(foo(n))
    foo=pack(X(:,i),mask)
    call getmeanvar(x=foo,mean=mu(i),var=std(i),method='f',err=err,message=mess)
    if(err/=0) then;mess='GetMonthlyClimato:'//trim(mess);return;endif
    std(i)=sqrt(std(i))
enddo
end subroutine GetMonthlyClimato
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RemoveMonthlyClimato(X,dates,err,mess)
!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:12/11/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. dates
!^* OUT
!^*    1. err
!^*    2. mess
!^* INOUT
!^*    1. X, big matrix
!^**********************************************************************
use Dates_tools
type(DateType),intent(in)::dates(:)
real(mrk), intent(inout)::X(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::i,month
real(mrk)::mu(12,size(X,dim=2)),std(12,size(X,dim=2))

err=0;mess=''

do month=1,12
    call GetMonthlyClimato(X,month,dates,mu(month,:),std(month,:),err,mess)
    if(err/=0) then;mess='RemoveMonthlyClimato:'//trim(mess);return;endif
enddo
do i=1,size(X,dim=1)
    X(i,:)=(X(i,:)-mu(dates(i)%month,:))/std(dates(i)%month,:)
enddo
end subroutine RemoveMonthlyClimato
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetMonthlyStat(XIN,datesIN,stat,XOUT,datesOUT,err,mess)
!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:12/11/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. XIN, big daily matrix
!^*    2. datesIN
!^*    3. stat, apply what? eg 'min' or 'max'
!^* OUT
!^*    1. XOUT, big monthly matrix
!^*    2. datesOUT
!^*    3. err
!^*    4. mess
!^**********************************************************************
use Dates_tools
type(DateType),intent(in)::datesIN(:)
real(mrk), intent(in)::XIN(:,:)
character(*), intent(in)::stat
real(mrk), pointer::XOUT(:,:)
type(DateType),pointer::datesOUT(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik), parameter::ndaymin=27
integer(mik)::i,j,n,ncol,p,k,m
logical::mask(size(datesIN))
real(mrk), allocatable::foo(:)

err=0;mess='';
n=size(datesIN);ncol=size(XIN,dim=2)
! first loop to get number of complete months
k=0
do i=datesIN(1)%year,datesIN(n)%year
    do j=1,12
        mask=(datesIN(:)%year==i).and.(datesIN(:)%month==j)
        p=count(mask)
        if(p>=ndaymin) k=k+1
    enddo
enddo
if(associated(XOUT)) nullify(XOUT);allocate(XOUT(k,ncol))
if(associated(datesOUT)) nullify(datesOUT);allocate(datesOUT(k))
k=0
do i=datesIN(1)%year,datesIN(n)%year
    write(*,*) i
    do j=1,12
        mask=(datesIN(:)%year==i).and.(datesIN(:)%month==j)
        p=count(mask)
        if(p>=ndaymin) then
            k=k+1
            datesOUT(k)=DateType(i,j,1,0,0,0)
            do m=1,ncol
                if(allocated(foo)) deallocate(foo);allocate(foo(p))
                foo=pack(XIN(:,m),mask)
                select case(stat)
                case('min')
                    XOUT(k,m)=minval(foo)
                case('max')
                    XOUT(k,m)=maxval(foo)
                case('mean')
                    XOUT(k,m)=sum(foo)/real(p,mrk)
                case default
                    err=1;mess='GetMonthlyStat:Fatal:unknown stat';return
                endselect
            enddo
        endif
    enddo
enddo

end subroutine GetMonthlyStat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetSeasonStat(XIN,datesIN,stat,SeasonStart,SeasonEnd,IsMonthly,XOUT,datesOUT,err,mess)
!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:12/11/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. XIN, big daily matrix
!^*    2. datesIN
!^*    3. stat, apply what? eg 'min' or 'max'
!^*    4-5. SeasonStart,SeasonEnd
!^*    6. [IsMontly]
!^* OUT
!^*    1. XOUT, big seasonal matrix
!^*    2. [DatesOUT], dates for each line (only year is relevant)
!^*    3. err
!^*    4. mess
!^**********************************************************************
use Dates_tools
type(DateType),intent(in)::datesIN(:),SeasonStart,SeasonEnd
real(mrk), intent(in)::XIN(:,:)
character(*), intent(in)::stat
logical, intent(in), optional::IsMonthly
type(DateType),pointer, optional::datesOUT(:)
real(mrk), pointer::XOUT(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
real(mrk), parameter::pdaymin=0.9_mrk
integer(mik)::i,n,ncol,p,k,m
real(mrk)::SeasonDuration
logical::mask(size(datesIN)),IsM
real(mrk), allocatable::foo(:)

err=0;mess='';
if(present(IsMonthly)) then
    IsM=IsMonthly
else
    IsM=.false.
endif
n=size(datesIN);ncol=size(XIN,dim=2)
SeasonDuration=Date2Days(SeasonEnd.Yminus.SeasonStart)
if(IsM) SeasonDuration=SeasonDuration/( 365.25_mrk/real(12,mrk))
! first loop to get number of complete seasons
k=0
do i=datesIN(1)%year,datesIN(n)%year
    mask=IsInSeason(datesIN,SeasonStart,SeasonEnd,i)
    p=count(mask)
    if(real(p,mrk)>=pdaymin*real(SeasonDuration,mrk)) k=k+1
enddo
if(associated(XOUT)) nullify(XOUT);allocate(XOUT(k,ncol))
if(present(DatesOUT)) then
    if(associated(DatesOUT)) nullify(DatesOUT);allocate(DatesOUT(k))
endif
k=0
do i=datesIN(1)%year,datesIN(n)%year
    mask=IsInSeason(datesIN,SeasonStart,SeasonEnd,i)
    p=count(mask)
    if(p>=pdaymin*SeasonDuration) then
        k=k+1
        if(present(DatesOUT)) DatesOUT(k)=DateType(i,1,1,0,0,0)
        do m=1,ncol
            if(allocated(foo)) deallocate(foo);allocate(foo(p))
            foo=pack(XIN(:,m),mask)
            select case(stat)
            case('min')
                XOUT(k,m)=minval(foo)
            case('max')
                XOUT(k,m)=maxval(foo)
            case('mean')
                XOUT(k,m)=sum(foo)/real(p,mrk)
            case default
                err=1;mess='GetSeasonStat:Fatal:unknown stat';return
            endselect
        enddo
    endif
enddo

end subroutine GetSeasonStat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetMonthlyInterannualStat(y,stat,out,err,mess)
!^**********************************************************************
!^* Purpose:
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
!^*    1. y, daily time series
!^*    2. stat, apply what? eg 'mean' or 'bfi'
!^* OUT
!^*    1. out, vector of size 12 cntaining the interanual stat
!^*    2. err
!^*    3. mess
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalStats
type(OneDTimeSeriesType), intent(in)::y
character(*), intent(in)::stat
real(mrk), intent(out)::out(12)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::month,n,p
logical, allocatable::mask(:)
real(mrk), allocatable::z(:)
! BEGIN WhereAmax only
integer(mik)::FirstYear, LastYear, year,mloc(1),MaxCount(12),monthofmax
integer(mik), parameter::nDayMin=365/2
! END WhereAmax only

! Init
out=y%mv;err=0;mess=''
n=size(y%ts);
if(allocated(mask)) deallocate(mask);allocate(mask(n))

! Optional prelimininaries
if(trim(stat)=='WhereAmax') then !extract annual maxima
    MaxCount=0
    FirstYear=y%ts(1)%date%year; LastYear=y%ts(y%n)%date%year;
    do year=FirstYear,LastYear
        mask=(y%ts(:)%date%year==year) .and.(y%ts(:)%q/=y%mv)
        if(count(mask)>nDayMin) then
            mloc=maxloc(y%ts(:)%q,mask)
            monthofmax=y%ts(mloc(1))%date%month
            MaxCount(monthofmax)=MaxCount(monthofmax)+1
        endif
    enddo
endif

! Apply stat
do month=1,12
    mask=(y%ts(:)%date%month==month).and.(y%ts(:)%q/=y%mv)
    p=count(mask)
    select case(trim(stat))
    case('mean')
        if(p==0) cycle
        out(month)=sum(y%ts(:)%q,mask)/real(p,mrk)
    case('sdev')
        if(p==0) cycle
        if(allocated(z)) deallocate(z);allocate(z(p))
        z=pack(y%ts(:)%q,mask)
        call GetEmpiricalStats(x=z,std=out(month),err=err,mess=mess)
        if(err>0) then;mess='GetMonthlyInterannualStat:FATAL:'//trim(mess);return;endif
    case('cv')
        if(p==0) cycle
        if(allocated(z)) deallocate(z);allocate(z(p))
        z=pack(y%ts(:)%q,mask)
        call GetEmpiricalStats(x=z,cv=out(month),err=err,mess=mess)
        if(err>0) then;mess='GetMonthlyInterannualStat:FATAL:'//trim(mess);return;endif
    case('WhereAmax')
        out(month)=1._mrk*MaxCount(month)
    case default
        err=1;mess='GetMonthlyInterannualStat:FATAL:Unknown[stat]';return
    end select
enddo

end subroutine GetMonthlyInterannualStat

!-----------------------------
! Filtering

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine TSFiltering_daily(filter,param,xin,xout,err,mess)

!^**********************************************************************
!^* Purpose: applies different kind of filtering to a daily time series
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:11/05/2009
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
integer(mik), intent(in)::filter
real(mrk), intent(in), optional::param(:)
type(OneDTimeSeriesType), intent(in)::xin
type(OneDTimeSeriesType), intent(out)::xout
integer(mik), intent(out)::err
character(*),intent(out)::mess

err=0;mess='';
if(associated(xout%ts)) nullify(xout%ts)

select case (filter)
case(MovingAverage)
    if (.not. present(param)) then
        err=1;mess='TSFiltering_daily:MovingAverage:Fatal:param is not optional for this filter';
        return
    else
        if (size(param)/=1) then
            err=1;mess='TSFiltering_daily:MovingAverage:Fatal:length(param) should be 1';
            return
        endif
        call MovingMean(xin=xin,d=param(1),xout=xout,err=err,mess=mess)
    endif
case(MovingMin)
    if (.not. present(param)) then
        err=1;mess='TSFiltering_daily:MovingMin:Fatal:param is not optional for this filter';
        return
    else
        if (size(param)/=1) then
            err=1;mess='TSFiltering_daily:MovingMin:Fatal:length(param) should be 1';
            return
        endif
        call MovingMinimum(xin=xin,d=param(1),xout=xout,err=err,mess=mess)
    endif
case(MovingMax)
    if (.not. present(param)) then
        err=1;mess='TSFiltering_daily:MovingMax:Fatal:param is not optional for this filter';
        return
    else
        if (size(param)/=1) then
            err=1;mess='TSFiltering_daily:MovingMax:Fatal:length(param) should be 1';
            return
        endif
        call MovingMaximum(xin=xin,d=param(1),xout=xout,err=err,mess=mess)
    endif
case(BFS)
    if (.not. present(param)) then
        call BaseFlowSeparation(xin=xin,xout=xout,err=err,mess=mess)
    else
        if (size(param)/=2) then
            err=1;mess='TSFiltering_daily:BFS:Fatal:length(param) should be 2';
            return
        endif
        call BaseFlowSeparation(xin=xin,blocLength=param(1),turningPoint=param(2),&
                                xout=xout,err=err,mess=mess)
    endif
case(SPA)
    if (.not. present(param)) then
        err=1;mess='TSFiltering_daily:SPA:Fatal:param is not optional for this filter';
        return
    else
        if (size(param)/=1) then
            err=1;mess='TSFiltering_daily:SPA:Fatal:length(param) should be 1';
            return
        endif
        call SequentPeakAlgorithm(xin=xin,Q0=param(1),xout=xout,err=err,mess=mess)
    endif
case default
    err=1;mess='TSFiltering_daily:Fatal:Unknown Filter';return
end select
end subroutine TSFiltering_daily
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine FillTheGap_daily(xin,xout,err,mess)

!^**********************************************************************
!^* Purpose: replace missing days in a time series by missing value
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:23/02/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.xin, time series with missing days
!^* OUT
!^*    1.xout, no missing days but missing values instead
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^* INOUT
!^*    1.
!^*    2.
!^*    3.
!^**********************************************************************

use Dates_tools, only: DaysBetween,DailyList, operator(==)

type(OneDTimeSeriesType), intent(in)::xin
type(OneDTimeSeriesType), intent(out)::xout
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::n,i,k

err=0;mess=''
! Populate TS type
n=nint(DaysBetween(xin%ts(1)%date,xin%ts(xin%n)%date))+1
xout%n=n
if(associated(xout%ts)) nullify(xout%ts)
allocate(xout%ts(n))
xout%mv=xin%mv;xout%NameSeries=trim(xin%NameSeries)//'_Filled'

!create list of dates for xout
xout%ts(:)%date=DailyList(xin%ts(1)%date,xin%ts(xin%n)%date)
! Fill
k=1
do i=1,n
    if(xout%ts(i)%date==xin%ts(k)%date) then
        xout%ts(i)%q=xin%ts(k)%q
        k=k+1
    else
        xout%ts(i)%q=xout%mv;
    endif
enddo

end subroutine FillTheGap_daily

!-----------------------------
! Misc. utilities

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetValueFromDate(xin,date,val,DateFound,err,mess)

!^**********************************************************************
!^* Purpose: Extract a value from a time series given a date
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon/Columbia
!^**********************************************************************
!^* Last modified: 05/09/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin, time series
!^*    2. date
!^* OUT
!^*    1.val, value if date found, UndefRN otherwise
!^*    2.DateFound, was the date present in the time series?
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use Dates_tools
type(OneDTimeSeriesType), intent(in)::xin
type(DateType), intent(in)::date
real(mrk), intent(out)::val
logical, intent(out)::DateFound
integer(mik), intent(out)::err
character(*), intent(out)::mess
! local
logical, allocatable::mask(:)
integer(mik)::i,n
integer(mik)::foo

err=0;mess='';val=undefRN;DateFound=.false.
if(allocated(mask)) deallocate(mask)
allocate(mask(xin%n))
do i=1,xin%n
    mask(i)=(xin%ts(i)%date==date)
    if(mask(i)) foo=i
enddo
n=count(mask)
if(n==0) then
    return
elseif (n==1) then
    val=xin%ts(foo)%q
    DateFound=.true.
    return
else
    DateFound=.true.
    err=1;mess='GetValueFromDate:several values found:not handled yet'
    return
endif
end subroutine GetValueFromDate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SetValueToDate(xin,date,val,DateFound,err,mess)

!^**********************************************************************
!^* Purpose: set a value at a given date in a time series
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 20/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. date
!^*    2.val, value
!^* OUT
!^*    2.DateFound, was the date present in the time series?
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^* INOUT
!^*    1. xin, time series
!^**********************************************************************
use Dates_tools
type(OneDTimeSeriesType), intent(inout)::xin
type(DateType), intent(in)::date
real(mrk), intent(in)::val
logical, intent(out)::DateFound
integer(mik), intent(out)::err
character(*), intent(out)::mess
! local
logical, allocatable::mask(:)
integer(mik)::i,n
integer(mik)::foo

err=0;mess='';DateFound=.false.
if(allocated(mask)) deallocate(mask)
allocate(mask(xin%n))
do i=1,xin%n
    mask(i)=(xin%ts(i)%date==date)
    if(mask(i)) foo=i
enddo
n=count(mask)
if(n==0) then
    return
elseif (n==1) then
    xin%ts(foo)%q=val
    DateFound=.true.
    return
else
    DateFound=.true.
    err=1;mess='SetValueToDate:several values found:not handled yet'
    return
endif
end subroutine SetValueToDate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ExtractOneDTimeSeries(TS,i,TSi,err,mess)
!^**********************************************************************
!^* Purpose: Extract the ith component of multivariate TS
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon/Columbia
!^**********************************************************************
!^* Last modified: 11/01/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. TS, time series
!^*    2. i,component to be extracted
!^* OUT
!^*    1.TSi, resultulting One-D TS
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
type(TimeSeriesType), intent(in)::TS
integer(mik), intent(in)::i
type(OneDTimeSeriesType), intent(out)::TSi
integer(mik), intent(out)::err
character(*), intent(out)::mess
! local
integer(mik)::t

err=0;mess=''
TSi%n=TS%n;TSi%mv=TS%mv
if(associated(TSi%ts)) nullify(TSi%ts);allocate(TSi%ts(TSi%n))
do t=1,TS%n
    if(i>size(TS%ts(t)%q)) then; mess='ExtractOneDTimeSeries: [i] too large';err=1;return;endif
    TSi%ts(t)%date=TS%ts(t)%date
    TSi%ts(t)%q=TS%ts(t)%q(i)
enddo

end subroutine ExtractOneDTimeSeries
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GatherMultiDTimeSeries(TS1D,TS,err,mess)
!^**********************************************************************
!^* Purpose: Assemble a multivariate TS from a vector of 1D TS
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 22/11/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. TS1D, vector of 1D time series
!^* OUT
!^*    1.TS, resultulting Multi-D TS
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use Dates_tools,only:DaysBetween,FindExtremalDates,DailyList
type(OneDTimeSeriesType), intent(in)::TS1D(:)
type(TimeSeriesType), intent(out)::TS
integer(mik), intent(out)::err
character(*), intent(out)::mess
! local
integer(mik)::t,s,Nx,from
type(datetype)::oldest,newest,oldcrap,newcrap,DateList(size(TS1D))

err=0;mess=''
Nx=size(TS1D)
do s=1,Nx;DateList(s)=TS1D(s)%ts(1)%date;enddo
call FindExtremalDates(DateList,oldest,newcrap,err,mess)
if(err/=0) then;mess='GatherMultiDTimeSeries:'//trim(mess);return;endif
do s=1,Nx;DateList(s)=TS1D(s)%ts(TS1D(s)%n)%date;enddo
call FindExtremalDates(DateList,oldcrap,newest,err,mess)
if(err/=0) then;mess='GatherMultiDTimeSeries:'//trim(mess);return;endif
TS%n=nint(DaysBetween(oldest,newest))+1
if(associated(TS%ts)) nullify(TS%ts);allocate(TS%ts(TS%n))
TS%ts(:)%date=DailyList(oldest,newest)
do t=1,TS%n
    if(associated(TS%ts(t)%q)) nullify(TS%ts(t)%q);allocate(TS%ts(t)%q(Nx))
    TS%ts(t)%q=TS%mv
enddo

do s=1,Nx
    from=nint(DaysBetween(oldest,TS1D(s)%ts(1)%date))
    do t=1,TS1D(s)%n
        TS%ts(from+t)%q(s)=TS1D(s)%ts(t)%q
    enddo
enddo

end subroutine GatherMultiDTimeSeries
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GatherMultiDTimeSeries_Matrix(TS1D,RemoveMV,TS,dates,err,mess)
!^**********************************************************************
!^* Purpose: Assemble a matrix TS from a vector of 1D TS
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 22/11/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. TS1D, vector of 1D time series
!^*    2. [RemoveMV], default .false.
!^* OUT
!^*    1.TS, resultulting matrix
!^*    2.dates
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
type(OneDTimeSeriesType), intent(in)::TS1D(:)
logical, intent(in), optional::RemoveMV
real(mrk), pointer::TS(:,:)
type(datetype), pointer::dates(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! local
logical::RmV,keepme
type(TimeSeriesType)::MultiD_TS
integer(mik)::n,Nx,t,k
real(mrk), allocatable::dummy(:,:)
type(datetype), allocatable::DummyDates(:)

err=0;mess=''
if(present(RemoveMV)) then
    RMV=RemoveMV
else
    RMV=.false.
endif
call GatherMultiDTimeSeries(TS1D,MultiD_TS,err,mess)
if(err/=0) then;mess='GatherMultiDTimeSeries_Matrix:'//trim(mess);return;endif
n=MultiD_TS%n
Nx=size(TS1D)
k=0
if(allocated(dummy)) deallocate(dummy);allocate(dummy(n,Nx))
if(allocated(DummyDates)) deallocate(DummyDates);allocate(DummyDates(n))
do t=1,n
    if(RMV) then
        keepme=all(MultiD_TS%ts(t)%q/=MultiD_TS%mv)
    else
        keepme=.true.
    endif
    if(keepme) then
        k=k+1
        dummy(k,:)=MultiD_TS%ts(t)%q
        DummyDates(k)=MultiD_TS%ts(t)%date
    endif
enddo
if(associated(TS)) nullify(TS);allocate(TS(k,Nx))
if(associated(dates)) nullify(dates);allocate(dates(k))
TS=dummy(1:k,:);dates=dummydates(1:k)
end subroutine GatherMultiDTimeSeries_Matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!
! PRIVATE !
!!!!!!!!!!!

pure subroutine MovingMean(xin,d,xout,err,mess)
!^**********************************************************************
!^* Purpose: Moving average on d previous time steps
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:13/05/2009
!^**********************************************************************
!^* Comments: 1/ doesnt work for variable time step series
!^* 2/ d is floored to closer integer
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin, raw time series
!^*    2. d, moving window length
!^* OUT
!^*    1.xout, filtered series
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
type(OneDTimeSeriesType), intent(in)::xin
real(mrk), intent(in)::d
type(OneDTimeSeriesType), intent(out)::xout
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::k,i
real(mrk), allocatable::window(:)
err=0;mess=''

if (xin%n==undefIN) then
    err=1;mess='MovingMean:Fatal:xin%n undefined';return
endif

if (associated(xout%ts)) deallocate(xout%ts)
allocate(xout%ts(xin%n))
xout%n=xin%n
xout%mv=xin%mv
xout%NameSeries=trim(xin%Nameseries)//'_MovAv'
xout%ts(:)%date=xin%ts(:)%date

k=int(d+0.99999_mrk) ! trick to avoid flooring at [d-1] due to machine precision
if (k<1) then
    err=2;mess='MovingMean:Fatal:Window length cannot be smaller than 1';return
endif

if (k>1) xout%ts(1:k-1)%q=xout%mv
do i=k,xin%n
    allocate(window(k))
    window=xin%ts(i-k+1:i)%q
    if (any(window(:)==xin%mv))then
        xout%ts(i)%q=xin%mv
    else
        xout%ts(i)%q=sum(window)/(1._mrk*k)
    endif
    deallocate(window)
enddo
end subroutine MovingMean
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine MovingMinimum(xin,d,xout,err,mess)
!^**********************************************************************
!^* Purpose: Moving minimum on d previous time steps
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:18/05/2009
!^**********************************************************************
!^* Comments: 1/ doesnt work for variable time step series
!^* 2/ d is floored to closest integer
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin, raw time series
!^*    2. d, moving window length
!^* OUT
!^*    1.xout, filtered series
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
type(OneDTimeSeriesType), intent(in)::xin
real(mrk), intent(in)::d
type(OneDTimeSeriesType), intent(out)::xout
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::k,i
real(mrk), allocatable::window(:)
err=0;mess=''

if (xin%n==undefIN) then
    err=1;mess='MovingMinimum:Fatal:xin%n undefined';return
endif

if (associated(xout%ts)) deallocate(xout%ts)
allocate(xout%ts(xin%n))
xout%n=xin%n
xout%mv=xin%mv
xout%NameSeries=trim(xin%Nameseries)//'_MovMin'
xout%ts(:)%date=xin%ts(:)%date

k=int(d+0.99999_mrk) ! trick to avoid flooring at [d-1] due to machine precision
if (k<1) then
    err=2;mess='MovingMinimum:Fatal:Window length cannot be smaller than 1';return
endif

if (k>1) xout%ts(1:k-1)%q=xout%mv
do i=k,xin%n
    allocate(window(k))
    window=xin%ts(i-k+1:i)%q
    if (any(window(:)==xin%mv))then
        xout%ts(i)%q=xin%mv
    else
        xout%ts(i)%q=minval(window)
    endif
    deallocate(window)
enddo
end subroutine MovingMinimum
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine MovingMaximum(xin,d,xout,err,mess)
!^**********************************************************************
!^* Purpose: Moving maximum on d previous time steps
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:18/05/2009
!^**********************************************************************
!^* Comments: 1/ doesnt work for variable time step series
!^* 2/ d is floored to closest integer
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin, raw time series
!^*    2. d, moving window length
!^* OUT
!^*    1.xout, filtered series
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
type(OneDTimeSeriesType), intent(in)::xin
real(mrk), intent(in)::d
type(OneDTimeSeriesType), intent(out)::xout
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::k,i
real(mrk), allocatable::window(:)
err=0;mess=''

if (xin%n==undefIN) then
    err=1;mess='MovingMaximum:Fatal:xin%n undefined';return
endif

if (associated(xout%ts)) deallocate(xout%ts)
allocate(xout%ts(xin%n))
xout%n=xin%n
xout%mv=xin%mv
xout%NameSeries=trim(xin%Nameseries)//'_MovMax'
xout%ts(:)%date=xin%ts(:)%date

k=int(d+0.99999_mrk) ! trick to avoid flooring at [d-1] due to machine precision
if (k<1) then
    err=2;mess='MovingMaximum:Fatal:Window length cannot be smaller than 1';return
endif

if (k>1) xout%ts(1:k-1)%q=xout%mv
do i=k,xin%n
    allocate(window(k))
    window=xin%ts(i-k+1:i)%q
    if (any(window(:)==xin%mv))then
        xout%ts(i)%q=xin%mv
    else
        xout%ts(i)%q=maxval(window)
    endif
    deallocate(window)
enddo
end subroutine MovingMaximum
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine BaseFlowSeparation(xin,blocLength,turningPoint,xout,err,mess)
!^**********************************************************************
!^* Purpose: BaseFlow Separation Algorithm
!^* cf Tallaksen and Van Lanen 2004
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:18/05/2009
!^**********************************************************************
!^* Comments: 1/ Default values for optional blocLength and turningPoint
!^* are 5 and 0.9 respectively, cf Tallaksen and Van Lanen 2004
!^* 2/ blocLength is floored to closest integer
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin, raw time series
!^*    2. <blocLength>
!^*    3. <turningPoint>
!^* OUT
!^*    1.xout, filtered series
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
type(OneDTimeSeriesType), intent(in)::xin
real(mrk), intent(in), optional::blocLength, turningPoint
type(OneDTimeSeriesType), intent(out)::xout
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::k,i,j,n,ntp, bl, indx(1)
real(mrk), allocatable::b1(:),b2(:),b3(:)
integer(mik), allocatable::tplist(:)
real(mrk)::bloc, tp, mini, a, b

err=0;mess=''

if (xin%n==undefIN) then
    err=1;mess='BFS:Fatal:xin%n undefined';return
endif

n=xin%n
if (associated(xout%ts)) deallocate(xout%ts)
allocate(xout%ts(n))
xout%n=n
xout%mv=xin%mv
xout%NameSeries=trim(xin%Nameseries)//'_BFS'
xout%ts(:)%date=xin%ts(:)%date

!Handle optional arguments
if(present(blocLength)) then
    bloc=blocLength
else
    bloc=5._mrk
endif
if(present(turningPoint)) then
    tp=turningPoint
else
    tp=0.9_mrk
endif

k=int(bloc+0.99999_mrk) ! trick to avoid flooring at [d-1] due to machine precision
if (k<1) then
    err=2;mess='BFS:Fatal:Bloc length cannot be smaller than 1';return
endif

xout%ts(:)%q=xout%mv
xout%ts(1)%q=xin%ts(1)%q;xout%ts(n)%q=xin%ts(n)%q
if (allocated(tplist)) deallocate(tplist)
allocate(tplist(n))
if (allocated(b1)) deallocate(b1)
allocate(b1(k))
if (allocated(b2)) deallocate(b2)
allocate(b2(k))
if (allocated(b3)) deallocate(b3)
allocate(b3(k))
ntp=1;tplist(1)=1
bl=floor((1.*n)/(1.*k))

do i=1, bl-2
    b1=xin%ts(k*(i-1)+1:k*i)%q
    b2=xin%ts(k*(i)+1:k*(i+1))%q
    b3=xin%ts(k*(i+1)+1:k*(i+2))%q

    if(any(b1==xin%mv).or.any(b2==xin%mv).or.any(b3==xin%mv)) then
        xout%ts(k*i+3)%q=xout%mv
        CYCLE
    endif
    if(tp*minval(b2)<=minval(b1)) then
        if (tp*minval(b2)<=minval(b3)) then
            indx=minloc(b2)
            mini=minval(b2)
            xout%ts(k*i+indx(1))%q=mini
            ntp=ntp+1
            tplist(ntp)=k*i+indx(1)
        endif
    endif
enddo

ntp=ntp+1;tplist(ntp)=n
! linear interpolation between turning points
do i=1, ntp-1
    if(tplist(i+1)-tplist(i)>1) then
        if(any(xin%ts(tplist(i):tplist(i+1))%q==xin%mv)) CYCLE
        a=(xout%ts(tplist(i+1))%q-xout%ts(tplist(i))%q)/(1.0_mrk*(tplist(i+1)-tplist(i)))
        b=(xout%ts(tplist(i))%q*tplist(i+1)-xout%ts(tplist(i+1))%q*tplist(i))/(1.0_mrk*(tplist(i+1)-tplist(i)))
        do j=tplist(i), tplist(i+1)
            xout%ts(j)%q=a*j+b
        enddo
    endif
enddo
end subroutine BaseFlowSeparation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine SequentPeakAlgorithm(xin, Q0, xout, err, mess)
!^**********************************************************************
!^* Purpose: sequent peak algorithm compute storage volume from daily discharge and needed rate Q0
!^**********************************************************************
!^* Programmer: Antoine Bard
!^**********************************************************************
!^* Last modified: 19/05/09
!^**********************************************************************
!^* Comments: works
!^* Ben's Comments: tricky handling of mv...
!^* Here, storage is forced back to zero every time a mv is encountered
!^* Consequently, it is strongly advised to discard years containing at
!^* least one mv during the analysis.
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin : dailydischarge series
!^*    2. Q0 : needed rate
!^* OUT
!^*    1. xout : dailystorage
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************

type(OneDtimeseriestype), intent(in)::xin
type(OneDtimeseriestype), intent(out)::xout
real(mrk), intent(in)::Q0
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::n, i
real(mrk)::previousS

err=0;mess=''

if (xin%n==undefIN) then
    err=1;mess='SPA:Fatal:xin%n undefined';return
endif
n=xin%n

if (associated(xout%ts)) deallocate(xout%ts)
allocate(xout%ts(n))
xout%n=n
xout%mv=xin%mv
xout%NameSeries=trim(xin%Nameseries)//'_SPA'
xout%ts(:)%date=xin%ts(:)%date

!null initial storage
if (xin%ts(1)%q==xin%mv) then ! MV
    xout%ts(1)%q=xout%mv;previousS=0._mrk
else
    If (Q0-xin%ts(1)%q >0._mrk) then
        xout%ts(1)%q=Q0-xin%ts(1)%q
    else
        xout%ts(1)%q=0._mrk
    end if
    previousS=xout%ts(1)%q
endif

do i=2,n
    If (xin%ts(i)%q==xin%mv) then
        xout%ts(i)%q=xout%mv;previousS=0._mrk
    else
        If (previousS+Q0-xin%ts(i)%q >0._mrk) then
            xout%ts(i)%q=previousS+Q0-xin%ts(i)%q
        else
            xout%ts(i)%q=0._mrk
        end if
        previousS=xout%ts(i)%q
    endif
end do
end subroutine SequentPeakAlgorithm
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module TimeSeries_tools

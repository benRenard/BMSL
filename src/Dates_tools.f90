module Dates_tools

!~**********************************************************************
!~* Purpose: Various tools for handling dates & durations
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified: 25/07/2011
!~**********************************************************************
!~* Comments: only works for year>0 for the time being
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*  . See "public" list below - names should be explicit enough...
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none

Private
public :: &
        ! basic operators
        operator(==), operator(/=), operator(<=), operator(>=),&
        operator(<), operator(>), operator(-), operator(+), & ! Overload of common operators
        operator(.Yeq.), operator(.Yneq.), operator(.Yle.), operator(.Yge.),&
        operator(.Ylt.), operator(.Ygt.), operator(.Yminus.),& ! Within-year operators
        ! logical funk
        IsDateValid,IsInSeason,&
        ! leap years utilities
        IsLeapYear, CountLeapYears, DaysInMonth,& 
        ! approximate date-duration conversions
        Days2Date, Date2Days,&
        ! absolute time conversions
        AbsoluteTime,DaysBetween,DaysFromRef,&
        ! stepping
        InOneYear, OneYearAgo, InNYears, & ! monthly stepping
        InOneMonth, OneMonthAgo, InNMonths, & ! monthly stepping
        Tomorrow,Yesterday,InNDays,InXDays,& ! daily stepping
        InOneHour,OneHourAgo,InNhours,& ! hourly stepping
        InOneMinute,OneMinuteAgo,InNMinutes,& ! minutely stepping
        InOneSecond,OneSecondAgo,InNSeconds,& ! hourly stepping
        ! lists
        DateList,DailyList,HourlyList,&
        ! Formatting
        FormatDate,UnformatDate,&
        ! Misc.
        FindExtremalDates,GetCommonDates
        
type, public :: DateType ! Used for dates AND for durations
    integer(mik) :: year = undefIN ! no negative year, ie only A.D. years - 0 allowed for durations
    integer(mik) :: month = 1
    integer(mik) :: day = 1
    integer(mik) :: hour = 0 ! note: 24:00 not allowed (00:00 instead)
    integer(mik) :: minute = 0
    integer(mik) :: second = 0
end type DateType

! Interfaces for overloading ==,/=,<, etc...
! I. Date comparison
Interface operator (==)
    module procedure Simultaneous
end interface
Interface operator (/=)
    module procedure Distinct
end interface
Interface operator (<=)
    module procedure Before
end interface
Interface operator (<)
    module procedure StrictlyBefore
end interface
Interface operator (>)
    module procedure StrictlyAfter
end interface
Interface operator (>=)
    module procedure After
end interface
! II. Within-year dates comparison
Interface operator (.Yeq.)
    module procedure IsIdenticalInYear
end interface
Interface operator (.Yneq.)
    module procedure IsDistinctInYear
end interface
Interface operator (.Yle.)
    module procedure IsEarlierInYear
end interface
Interface operator (.Ylt.)
    module procedure IsStrictlyEarlierInYear
end interface
Interface operator (.Yge.)
    module procedure IsLaterInYear
end interface
Interface operator (.Ygt.)
    module procedure IsStrictlyLaterInYear
end interface
! III. Dates algebra
Interface operator (-) ! See Comment in sub MinusOverload
    module procedure MinusOverload
end interface
Interface operator (.Yminus.) ! See Comment in sub MinusOverload
    module procedure DifferenceInYear
end interface
Interface operator (+)
    module procedure PlusOverload
end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------
! Basic operators

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function StrictlyBefore(date1,date2)

!#**********************************************************************
!#* Purpose: is date1 strictly before date2?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:04/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#*  3.
!#* OUT
!#*  1.logical StrictlyBefore
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical StrictlyBefore

if( (.not. IsDateValid(date1)) .or. (.not. IsDateValid(date2)) ) then
    StrictlyBefore=.false.
    return
endif

if(date1%year<date2%year) then
    StrictlyBefore=.true.;return
elseif(date1%year>date2%year) then
    StrictlyBefore=.false.;return
endif

if(date1%month<date2%month) then
    StrictlyBefore=.true.;return
elseif(date1%month>date2%month) then
    StrictlyBefore=.false.;return
endif

if(date1%day<date2%day) then
    StrictlyBefore=.true.;return
elseif(date1%day>date2%day) then
    StrictlyBefore=.false.;return
endif

if(date1%hour<date2%hour) then
    StrictlyBefore=.true.;return
elseif(date1%hour>date2%hour) then
    StrictlyBefore=.false.;return
endif

if(date1%minute<date2%minute) then
    StrictlyBefore=.true.;return
elseif(date1%minute>date2%minute) then
    StrictlyBefore=.false.;return
endif

if(date1%second<date2%second) then
    StrictlyBefore=.true.;return
elseif(date1%second>date2%second) then
    StrictlyBefore=.false.;return
endif

if(date1%second==date2%second) StrictlyBefore=.false.

end function StrictlyBefore

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function Simultaneous(date1,date2)

!#**********************************************************************
!#* Purpose: are dates equal?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified: 04/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.Simultaneous
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical Simultaneous

if( (.not. IsDateValid(date1)) .or. (.not. IsDateValid(date2)) ) then
    Simultaneous=.false.
    return
endif

if ( (date1%year/=date2%year) .or. (date1%month/=date2%month) .or. &
     (date1%day/=date2%day) .or. (date1%hour/=date2%hour) .or. &
     (date1%minute/=date2%minute) .or. (date1%second/=date2%second) ) then
    Simultaneous=.false.
else
    Simultaneous=.true.
endif

end function Simultaneous

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function Distinct(date1,date2)

!#**********************************************************************
!#* Purpose: are dates distinct?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified: 04/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.Distinct
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical Distinct

if( (.not. IsDateValid(date1)) .or. (.not. IsDateValid(date2)) ) then
    Distinct=.false.
    return
endif

if ( Simultaneous(date1,date2) ) then
    Distinct=.false.
else
    Distinct=.true.
endif

end function Distinct

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function Before(date1,date2)

!#**********************************************************************
!#* Purpose: is date1 before (or simultaneous to) date2?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:04/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#*  3.
!#* OUT
!#*  1.logical Before
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical Before

if( (.not. IsDateValid(date1)) .or. (.not. IsDateValid(date2)) ) then
    Before=.false.
    return
endif

if( StrictlyBefore(date1,date2) .or. Simultaneous(date1,date2) ) then
    Before=.true.
else
    Before=.false.
endif

end function Before

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function After(date1,date2)

!#**********************************************************************
!#* Purpose: is date1 After (or simultaneous to) date2?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:04/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#*  3.
!#* OUT
!#*  1.logical After
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical After

if( (.not. IsDateValid(date1)) .or. (.not. IsDateValid(date2)) ) then
    After=.false.
    return
endif

if( .not. StrictlyBefore(date1,date2) ) then
    After=.true.
else
    After=.false.
endif

end function After

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function StrictlyAfter(date1,date2)

!#**********************************************************************
!#* Purpose: is date1 Strictly After date2?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:04/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.logical StrictlyAfter
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical StrictlyAfter

if( (.not. IsDateValid(date1)) .or. (.not. IsDateValid(date2)) ) then
    StrictlyAfter=.false.
    return
endif

if( .not. Before(date1,date2) ) then
    StrictlyAfter=.true.
else
    StrictlyAfter=.false.
endif

end function StrictlyAfter

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function MinusOverload(d1,d2)

!#**********************************************************************
!#* Purpose: overload d1-d2 where d1 and d2 are dates
!#* returns the duration between d1 and d2 in a date format
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments: d1-d2=d2-d1 - order doesnt matter
!#* Assumes constant-length year of 365.25 days. Consequently,
!#* 01/01/2003-01/01/2002 is NOT 1year0month0days0hour0min0sec
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List: Alternative version with variable-length year
!#**********************************************************************
!#* IN
!#*  1.d1
!#*  2.d2
!#* OUT
!#*  1.MinusOverload
!#**********************************************************************
type(dateType),intent(in)::d1,d2
type(dateType)::MinusOverload
!locals
real(mrk)::x

x=DaysBetween(d1,d2)
MinusOverload=Days2Date(x)

end function MinusOverload

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function PlusOverload(d1,d2)

!#**********************************************************************
!#* Purpose: overload d1+d2 where d1 and d2 are durations in date format
!#* returns the total duration d1+d2 in a date format
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.d1
!#*  2.d2
!#* OUT
!#*  1.PlusOverload
!#**********************************************************************
type(dateType),intent(in)::d1,d2
type(dateType)::PlusOverload

PlusOverload=Days2Date(Date2Days(d1)+Date2Days(d2))

end function PlusOverload

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function IsEarlierInYear(date1,date2)

!#**********************************************************************
!#* Purpose: test whether date1 <= date2 WITHIN a year
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.IsEarlierInYear
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical:: IsEarlierInYear
!locals
type(DateType)::dummyD1,dummyD2

dummyD1=date1;dummyD2=date2
dummyD1%year=2000;dummyD2%year=2000;
IsEarlierInYear=dummyD1<=dummyD2

end function IsEarlierInYear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function IsStrictlyEarlierInYear(date1,date2)

!#**********************************************************************
!#* Purpose: test whether date1 < date2 WITHIN a year
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.IsStrictlyEarlierInYear
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical:: IsStrictlyEarlierInYear
!locals
type(DateType)::dummyD1,dummyD2

dummyD1=date1;dummyD2=date2
dummyD1%year=2000;dummyD2%year=2000;
IsStrictlyEarlierInYear=dummyD1<dummyD2

end function IsStrictlyEarlierInYear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function IsStrictlyLaterInYear(date1,date2)

!#**********************************************************************
!#* Purpose: test whether date1 > date2 WITHIN a year
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.IsStrictlyLaterInYear
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical:: IsStrictlyLaterInYear
!locals
type(DateType)::dummyD1,dummyD2

dummyD1=date1;dummyD2=date2
dummyD1%year=2000;dummyD2%year=2000
IsStrictlyLaterInYear=dummyD1>dummyD2

end function IsStrictlyLaterInYear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function IsLaterInYear(date1,date2)

!#**********************************************************************
!#* Purpose: test whether date1 >= date2 WITHIN a year
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.IsLaterInYear
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical:: IsLaterInYear
!locals
type(DateType)::dummyD1,dummyD2

dummyD1=date1;dummyD2=date2
dummyD1%year=2000;dummyD2%year=2000
IsLaterInYear=dummyD1>=dummyD2

end function IsLaterInYear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function IsIdenticalInYear(date1,date2)

!#**********************************************************************
!#* Purpose: test whether date1 == date2 WITHIN a year
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.IsIdenticalInYear
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical:: IsIdenticalInYear
!locals
type(DateType)::dummyD1,dummyD2

dummyD1=date1;dummyD2=date2
dummyD1%year=2000;dummyD2%year=2000
IsIdenticalInYear=(dummyD1==dummyD2)

end function IsIdenticalInYear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function IsDistinctInYear(date1,date2)

!#**********************************************************************
!#* Purpose: test whether date1 /= date2 WITHIN a year
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.IsDistinctInYear
!#**********************************************************************
type(DateType), intent(in)::date1,date2
logical:: IsDistinctInYear
!locals
type(DateType)::dummyD1,dummyD2

dummyD1=date1;dummyD2=date2
dummyD1%year=2000;dummyD2%year=2000
IsDistinctInYear=(dummyD1/=dummyD2)

end function IsDistinctInYear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function DifferenceInYear(date1,date2)

!#**********************************************************************
!#* Purpose: compute date1-date2 WITHIN a year;
!#* eg. December - March = 9
!#* March - December = 3
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon, & Xun Sun
!#**********************************************************************
!#* Last modified:22/11/2011
!#**********************************************************************
!#* Comments: Count based on a leap year, ie 366 days/year
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.DifferenceInYear
!#**********************************************************************
type(DateType), intent(in)::date1,date2
type(DateType):: DifferenceInYear
!locals
type(DateType)::dummyD1,dummyD2

dummyD1=date1;dummyD2=date2
dummyD1%year=2;dummyD2%year=2
if(dummyD1<dummyD2) then !March - December (ie season from Dec to Mar)
    dummyD2%year=1
else !December - March (ie season from Mar 2 Dec)
    dummyD2%year=2
endif

DifferenceInYear=dummyD1-dummyD2

end function DifferenceInYear


!-----------------------------
! Logical funks

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function IsDateValid(date)

!#**********************************************************************
!#* Purpose: Is date a valid date?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:04/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1. logical IsDateValid
!#**********************************************************************
type(DateType),intent(in)::date
logical IsDateValid

IsDateValid=.true.

! year >=0
if (date%year<0) then
    IsDateValid=.false.
    return
endif

! Month between 1 and 12
if ((date%month>12).or.(date%month<1)) then
    IsDateValid=.false.
    return
endif

! day - slightly more tricky...
if (date%day<1) then
    IsDateValid=.false.
    return
endif
select case (date%month)
case(1,3,5,7,8,10,12)
    if (date%day>31) then
        IsDateValid=.false.
        return
    endif
case(4,6,9,11)
    if (date%day>30) then
        IsDateValid=.false.
        return
    endif
case(2)
    if(IsLeapYear(date%year)) then
        if (date%day>29) then
            IsDateValid=.false.
            return
        endif
    else
        if (date%day>28) then
            IsDateValid=.false.
            return
        endif
    endif
case default ! should never happen
    IsDateValid=.false.
    return
end select

!hour in (0;23)
if ((date%hour>23).or.(date%hour<0)) then
    IsDateValid=.false.
    return
endif

!minute in (0;59)
if ((date%minute>59).or.(date%minute<0)) then
    IsDateValid=.false.
    return
endif

!second in (0;59)
if ((date%second>59).or.(date%second<0)) then
    IsDateValid=.false.
    return
endif

end function IsDateValid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function IsInSeason(date,Sstart,Send,year)

!#**********************************************************************
!#* Purpose: test whether date is included in the season define by
!#* Sstart & Send
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments: Sstart might be after Send (eg october to march).
!#* year is an optional argument: if not present, the function tests
!#* whether date is in the season, whatever the year
!#* if year provided, the function tests whether date is in the season
!#* of a given year
!#* Convention: season from march to october 1960 is hydrological year 1960
!#* season from october 1960 to march 1961 is hydrological year 1960
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.Sstart
!#*  3.Send
!#*  4.<year>
!#* OUT
!#*  1.IsInSeason
!#**********************************************************************
type(DateType), intent(in)::date,Sstart,Send
integer,intent(in),optional::year
logical:: IsInSeason
!locals
type(DateType)::dummyS,dummyE


if( (.not. IsDateValid(Sstart)) .or. (.not. IsDateValid(Send)) .or. (.not. IsDateValid(date))) then
    IsInSeason=.false.
    return
endif

dummyS=Sstart;dummyE=SEnd

if(Sstart .Yle. Send) then !march to october
    if(present(year)) then
        dummyS%year=year;dummyE%year=year
        IsInSeason= (date>=dummyS) .and. (date<=dummyE)
    else
        IsInSeason= (date .Yge. dummyS) .and. (date .Yle. dummyE)
    endif
else ! october to march
    if(present(year)) then
        dummyS%year=year;dummyE%year=year+1
        IsInSeason= (date>=dummyS) .and. (date<=dummyE)
    else
        IsInSeason= ((date .Yge. dummyS) .and. (date .Yle. DateType(1,12,31,23,59,59))) &
            .or. ((date .Yle. dummyE) .and. (date .Yge. DateType(1,1,1,0,0,0)))

    endif
endif

end function IsInSeason

!-----------------------------
! Leap years utilities

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function IsLeapYear(year)

!#**********************************************************************
!#* Purpose:
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.year
!#* OUT
!#*  1.IsLeapYear
!#**********************************************************************
integer(mik),intent(in)::year
logical IsLeapYear

if(year==0) then
    IsLeapYear=.false.
    return
endif
if ((mod(year,4)==0 .and. mod(year,100)/=0) .or. (mod(year,400)==0) ) then
    IsLeapYear=.true.
else
    IsLeapYear=.false.
endif

end function IsLeapYear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function CountLeapYears(year1,year2)

!#**********************************************************************
!#* Purpose: Count leap years between year1 and year2 (both included)
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:04/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.year1
!#*  2.year2
!#* OUT
!#*  1.CountLeapYears
!#**********************************************************************
integer(mik),intent(in)::year1,year2
integer(mik) CountLeapYears
!locals
integer(mik),allocatable:: ylist(:)
integer(mik)::n,i

n=abs(year1-year2)+1
if(allocated(ylist)) deallocate(ylist)
allocate(ylist(n))
if(year1<=year2) then
    ylist=(/(year1+i-1,i=1,n)/)
else
    ylist=(/(year2+i-1,i=1,n)/)
endif
CountLeapYears=count(IsLeapYear(ylist))
deallocate(ylist)
end function CountLeapYears

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function DaysInMonth(month,year)

!#**********************************************************************
!#* Purpose: How many days in this month of this year?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:28/07/2011
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.month
!#*  2.year
!#* OUT
!#*  1.Number of days
!#**********************************************************************
integer(mik),intent(in)::month,year
integer(mik) DaysInMonth

select case (month)
case(1,3,5,7,8,10,12)
    DaysInMonth=31
case(4,6,9,11)
    DaysInMonth=30
case(2)
    if(IsLeapYear(year)) then
        DaysInMonth=29
    else
        DaysInMonth=28
    endif
end select

end function DaysInMonth

!--------------------------------------
! approximate date-duration conversions

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function Date2Days(duration)

!#**********************************************************************
!#* Purpose: convert a duration from date format 2 days
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments: Since date is considered here as a duration, the following
!#* approximations are made: 1year=365.25 days and 1 month=30.4375 days
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.duration
!#* OUT
!#*  1.Date2Days
!#**********************************************************************
type(DateType), intent(in)::duration
real(mrk)::Date2Days

Date2Days=365.25_mrk*(duration%year)+ &
          30.4375_mrk*(duration%month)+ &
          1._mrk*(duration%day)+ &
         (1._mrk/24._mrk)*duration%hour+ &
         (1._mrk/1440._mrk)*duration%minute+ &
         (1._mrk/86400._mrk)*duration%second

end function Date2Days

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function Days2Date(d)

!#**********************************************************************
!#* Purpose: convert a duration from days 2 date format
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments: Since d is considered here as a duration, the following
!#* approximations are made: 1year=365.25 days and 1 month=30.4375 days.
!#* Likely roundoff error since date format does not go below seconds
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.d
!#* OUT
!#*  1.Days2Date
!#**********************************************************************
real(mrk), intent(in)::d
type(DateType)::Days2Date
!locals
real(mrk)::rem

Days2Date%year=int(d/365.25_mrk);rem=mod(d,365.25_mrk)
Days2Date%month=int(rem/30.4375_mrk);rem=mod(rem,30.4375_mrk)
Days2Date%day=int(rem/1._mrk);rem=mod(rem,1._mrk)
Days2Date%hour=int(rem/(1._mrk/24._mrk));rem=mod(rem,1._mrk/24._mrk)
Days2Date%minute=int(rem/(1._mrk/1440_mrk));rem=mod(rem,1._mrk/1440._mrk)
Days2Date%second=int(rem/(1._mrk/86400._mrk))

end function Days2Date

!--------------------------------------
! Absolute time conversions

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function AbsoluteTime(date,option,FirstDate)

!#**********************************************************************
!#* Purpose: transform an array of dates into absolute time
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified:10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.[option] needed to balance speed vs. numerical accuracy
!#*            +1 = forward computation, growing inaccuracy
!#*            -1 = backward computation, exact at the end of the series but inaccurate when moving back
!#*             0 = first part (~4/5) forward, last part backward (based on empirical observations)
!#*           666 = (default) full treatment, by far the most accurate, but slower (~10 times)
!#*  3.[FirstDate] where is t=0? if absent set to first element of array "date"
!#* OUT
!#*  1.AbsoluteTime
!#**********************************************************************
type(dateType), intent(in)::date(:)
integer(mik), intent(in), optional::option
type(dateType), intent(in),optional::FirstDate
real(mrk)::AbsoluteTime(size(date))
!locals
integer(mik), parameter::option_def=666
real(mrk), parameter::frac=0.8_mrk ! fraction of the series to be computed forward
integer(mik)::i,n,opt,lim
real(mrk)::dummyfor(size(date)),dummybak(size(date))
type(DateType)::d0

AbsoluteTime=undefRN
if(present(option)) then; opt=option;else;opt=option_def;endif
if(present(FirstDate)) then; d0=FirstDate;else;d0=date(1);endif
n=size(date);dummyfor=0._mrk;dummybak=0._mrk

if(opt==1) then ! forward
    dummyfor(1)=DaysFromRef(date(1),d0)
    do i=2,n
        dummyfor(i)=dummyfor(i-1)+DaysFromRef(date(i),date(i-1))
    enddo
endif

if(opt==-1) then ! backward
    dummybak(n)=DaysFromRef(date(n),d0)
    do i=n-1,1,-1
        dummybak(i)=dummybak(i+1)-DaysFromRef(date(i+1),date(i))
    enddo
endif

if(opt==0) then ! mixed
    lim=floor(real(n,mrk)*frac)
    dummyfor(1)=DaysFromRef(date(1),d0)
    do i=2,lim
        dummyfor(i)=dummyfor(i-1)+DaysFromRef(date(i),date(i-1))
    enddo
    dummybak(n)=DaysFromRef(date(n),d0)
    do i=n-1,lim+1,-1
        dummybak(i)=dummybak(i+1)-DaysFromRef(date(i+1),date(i))
    enddo
endif

select case(opt)
case(-1);AbsoluteTime=dummybak    
case(0);AbsoluteTime=dummybak+dummyfor
case(1);AbsoluteTime=dummyfor
case(666);AbsoluteTime=DaysFromRef(date,d0)
case default;AbsoluteTime=undefRN
endselect

end function AbsoluteTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function DaysBetween(date1,date2)

!#**********************************************************************
!#* Purpose: Duration in days between date1 and date2
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:05/02/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.DaysBetween
!#**********************************************************************
type(DateType), intent(in)::date1,date2
real(mrk)::DaysBetween

!DaysBetween=abs(DaysFromZero(date1)-DaysFromZero(date2))
DaysBetween=abs(DaysFromRef(date1,date2))

end function DaysBetween

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function DaysFromRef(date,ref)

!#**********************************************************************
!#* Purpose: compute the duration (in days) between date and ref
!#*          if(date<ref), duration is negative 
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified:10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.ref
!#* OUT
!#*  1.DaysFromRef
!#**********************************************************************
type(dateType), intent(in)::date,ref
real(mrk)::DaysFromRef
!locals
integer(mik)::nleap,n
type(dateType)::d0,d
real(mrk)::k,z,z0,zend

if(date>=ref) then
    d0=ref;d=date;k=1._mrk
else
    d0=date;d=ref;k=-1._mrk
endif

DaysFromRef=0._mrk
z0=DaysFromJanuaryFirst(d0)
if(isLeapYear(d0%year)) then;zend=366._mrk;else;zend=365._mrk;endif
z=DaysFromJanuaryFirst(d)
if(d0%year==d%year) then
    DaysFromRef=k*(z-z0);return
else if(d0%year==d%year-1) then
    DaysFromRef=k*(zend-z0+z)
else
    n=d%year-d0%year-1
    nleap=CountLeapYears(d0%year+1,d%year-1)
    DaysFromRef=k*(zend-z0+z+real(nleap,mrk)*366._mrk+real(n-nleap,mrk)*365._mrk)
endif

end function DaysFromRef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function DaysFromZero(date)

!#**********************************************************************
!#* Purpose: compute the duration (in days) from Christ's birth
!#* 01/01/01 00:00:00
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:04/02
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.DaysFromZero
!#**********************************************************************
type(dateType), intent(in)::date
real(mrk)::DaysFromZero
!locals
integer(mik)::nLeap,n,DinM(12)

DaysFromZero=0._mrk
n=date%year-1
nLeap=CountLeapYears(1,date%year-1)
DaysFromZero=DaysFromZero+nLeap*366._mrk+(n-nLeap)*365._mrk

if(IsLeapYear(date%year)) then
    DinM=(/31,29,31,30,31,30,31,31,30,31,30,31/)
else
    DinM=(/31,28,31,30,31,30,31,31,30,31,30,31/)
endif
if(date%month>1) DaysFromZero=DaysFromZero+1._mrk*sum(DinM(1:date%month-1))
DaysFromZero=DaysFromZero+1._mrk*(date%day-1)+ &
             (1._mrk/24._mrk)*date%hour+ &
             (1._mrk/1440._mrk)*date%minute+ &
             (1._mrk/86400._mrk)*date%second

end function DaysFromZero

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function DaysFromJanuaryFirst(date)
!#**********************************************************************
!#* Purpose: compute the duration (in days) from 01/01 00:00:00 of the
!#* current year
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified:10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.DaysFromJanuaryFirst
!#**********************************************************************
type(dateType), intent(in)::date
real(mrk)::DaysFromJanuaryFirst
!locals
integer(mik)::i

DaysFromJanuaryFirst=0._mrk
do i=1,date%month-1
    DaysFromJanuaryFirst=DaysFromJanuaryFirst+real(DaysInMonth(i,date%year),mrk)
enddo
DaysFromJanuaryFirst=DaysFromJanuaryFirst+real(date%day-1,mrk)+ &
             (1._mrk/24._mrk)*real(date%hour,mrk)+ &
             (1._mrk/1440._mrk)*real(date%minute,mrk)+ &
             (1._mrk/86400._mrk)*real(date%second,mrk)
end function DaysFromJanuaryFirst

!--------------------------------------
! stepping

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function InOneYear(date,option)

!#**********************************************************************
!#* Purpose: what will be the date in 1 year?
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified:10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.[option]; if 0 (default), does not check that the result is a feasible date
!#*              if 1, look for the closest feasible date (in the past)
!#* OUT
!#*  1.date InOneYear
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in),optional::option
type(DateType):: InOneYear
! locals
integer(mik), parameter::option_default=0
integer(mik)::opt
type(DateType)::dummy

if(present(option)) then; opt=option;else;opt=option_default;endif
dummy=date
dummy%year=date%year+1;
if(opt==0) then
    InOneYear=dummy
else
    InOneYear=ClosestFeasibleDate(dummy)
endif
end function InOneYear

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function OneYearAgo(date,option)

!#**********************************************************************
!#* Purpose: what was the date 1 year ago?
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified:10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.[option]; if 0 (default), does not check that the result is a feasible date
!#*              if 1, look for the closest feasible date (in the past)
!#* OUT
!#*  1.date OneYearAgo
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in),optional::option
type(DateType):: OneYearAgo
! locals
integer(mik), parameter::option_default=0
integer(mik)::opt
type(DateType)::dummy

if(present(option)) then; opt=option;else;opt=option_default;endif
dummy=date
dummy%year=date%year-1;
if(opt==0) then
    OneYearAgo=dummy
else
    OneYearAgo=ClosestFeasibleDate(dummy)
endif
end function OneYearAgo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function InNYears(date,N,option)
!#**********************************************************************
!#* Purpose: what will be the date in N years? (same day & hour)
!#* (or what was the date N months ago if N<0)
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified:10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.N, >0 or <0
!#*  3.[option]; if 0 (default), does not check that the result is a feasible date
!#*              if 1, look for the closest feasible date (in the past)
!#* OUT
!#*  1.InNYears
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in)::N
integer(mik), intent(in),optional::option
type(DateType):: InNYears
!locals
integer(mik), parameter::option_default=0
integer(mik)::opt
type(DateType)::dummy

if(present(option)) then; opt=option;else;opt=option_default;endif
dummy=date
dummy%year=dummy%year+N
if(opt==0) then
    InNYears=dummy
else
    InNYears=ClosestFeasibleDate(dummy)
endif
end function InNYears

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function InOneMonth(date,option)

!#**********************************************************************
!#* Purpose: what will be the date in 1 month?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.[option]; if 0 (default), does not check that the result is a feasible date
!#*              if 1, look for the closest feasible date (in the past)
!#* OUT
!#*  1.date InOneMonth
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in),optional::option
type(DateType):: InOneMonth
! locals
integer(mik), parameter::option_default=0
integer(mik)::opt
type(DateType)::dummy

if(present(option)) then; opt=option;else;opt=option_default;endif
dummy=date
if(date%month==12) then
    dummy%month=1;dummy%year=date%year+1
else
    dummy%month=date%month+1;
endif
if(opt==0) then
    InOneMonth=dummy
else
    InOneMonth=ClosestFeasibleDate(dummy)
endif
end function InOneMonth

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function OneMonthAgo(date,option)

!#**********************************************************************
!#* Purpose: what was the date 1 month ago?
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified:10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.[option]; if 0 (default), does not check that the result is a feasible date
!#*              if 1, look for the closest feasible date (in the past)
!#* OUT
!#*  1.date OneMonthAgo
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in),optional::option
type(DateType):: OneMonthAgo
! locals
integer(mik), parameter::option_default=0
integer(mik)::opt
type(DateType)::dummy

if(present(option)) then; opt=option;else;opt=option_default;endif
dummy=date
if(date%month==1) then
    dummy%month=12;dummy%year=date%year-1
else
    dummy%month=date%month-1;
endif
if(opt==0) then
    OneMonthAgo=dummy
else
    OneMonthAgo=ClosestFeasibleDate(dummy)
endif
end function OneMonthAgo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function InNMonths(date,N,option)
!#**********************************************************************
!#* Purpose: what will be the date in N months? (same day & hour)
!#* (or what was the date N months ago if N<0)
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified:10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.N, >0 or <0
!#*  3.[option]; if 0 (default), does not check that the result is a feasible date
!#*              if 1, look for the closest feasible date (in the past)
!#* OUT
!#*  1.InNMonths
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in)::N
integer(mik), intent(in),optional::option
type(DateType):: InNMonths
!locals
integer(mik), parameter::option_default=0
integer(mik)::opt,i
type(DateType)::dummy

if(present(option)) then; opt=option;else;opt=option_default;endif
dummy=date
if(N>0) then
    do i=1,N
        if( (dummy%month+1)<=12) then
            dummy%month=dummy%month+1
        else
            dummy%month=1
            dummy%year=dummy%year+1
        endif
    enddo
elseif(N<0) then
    do i=1,abs(N)
        if( (dummy%month-1)>=1) then
            dummy%month=dummy%month-1
        else
            dummy%month=12
            dummy%year=dummy%year-1
        endif
    enddo
endif
if(opt==0) then
    InNMonths=dummy
else
    InNMonths=ClosestFeasibleDate(dummy)
endif
end function InNMonths

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function Tomorrow(date)

!#**********************************************************************
!#* Purpose: Tomorrow's date (same hour)
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:23/02/2010
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.Tomorrow
!#**********************************************************************
type(DateType), intent(in)::date
type(DateType):: Tomorrow

Tomorrow%second=date%second
Tomorrow%minute=date%minute
Tomorrow%hour=date%hour
select case (date%month)
case(1,3,5,7,8,10)
    if (date%day==31) then
        Tomorrow%day=1;Tomorrow%month=date%month+1;Tomorrow%year=date%year
    else
        Tomorrow%day=date%day+1;Tomorrow%month=date%month;Tomorrow%year=date%year
    endif
case(4,6,9,11)
    if (date%day==30) then
        Tomorrow%day=1;Tomorrow%month=date%month+1;Tomorrow%year=date%year
    else
        Tomorrow%day=date%day+1;Tomorrow%month=date%month;Tomorrow%year=date%year
    endif
case(2)
    if(IsLeapYear(date%year)) then
        if (date%day==29) then
            Tomorrow%day=1;Tomorrow%month=date%month+1;Tomorrow%year=date%year
        else
            Tomorrow%day=date%day+1;Tomorrow%month=date%month;Tomorrow%year=date%year
        endif
    else
        if (date%day==28) then
            Tomorrow%day=1;Tomorrow%month=date%month+1;Tomorrow%year=date%year
        else
            Tomorrow%day=date%day+1;Tomorrow%month=date%month;Tomorrow%year=date%year
        endif
    endif
case(12)
    if (date%day==31) then
        Tomorrow%day=1;Tomorrow%month=1;Tomorrow%year=date%year+1
    else
        Tomorrow%day=date%day+1;Tomorrow%month=date%month;Tomorrow%year=date%year
    endif
case default ! should never happen
    Tomorrow%year=undefIN
end select

end function Tomorrow

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function Yesterday(date)

!#**********************************************************************
!#* Purpose: Yesterday's date (same hour)
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:23/02/2010
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.Yesterday
!#**********************************************************************
type(DateType), intent(in)::date
type(DateType):: Yesterday

Yesterday%second=date%second
Yesterday%minute=date%minute
Yesterday%hour=date%hour
if(date%day==1) then
    select case (date%month-1)
    case(1,3,5,7,8,10,12)
        Yesterday%day=31;Yesterday%month=date%month-1;Yesterday%year=date%year
    case(4,6,9,11)
        Yesterday%day=30;Yesterday%month=date%month-1;Yesterday%year=date%year
    case(2)
        if(IsLeapYear(date%year)) then
            Yesterday%day=29;Yesterday%month=date%month-1;Yesterday%year=date%year
        else
            Yesterday%day=28;Yesterday%month=date%month-1;Yesterday%year=date%year
        endif
    case(0)
        Yesterday%day=31;Yesterday%month=12;Yesterday%year=date%year-1
    case default ! should never happen
        Yesterday%year=undefIN
    end select
else
    Yesterday%day=date%day-1;Yesterday%month=date%month;Yesterday%year=date%year
endif

end function Yesterday

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function InNDays(date,N)

!#**********************************************************************
!#* Purpose: what will be the date in N days? (same hour)
!#* (or what was the date N days ago if N<0)
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:23/02/2010
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.N, >0 or <0
!#* OUT
!#*  1.InNDays
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in)::N
type(DateType):: InNDays
!locals
type(DateType):: dummy
integer(mik)::i

dummy=date
if(N>0) then
    do i=1,N
        dummy=Tomorrow(dummy)
    enddo
elseif(N<0) then
    do i=1,abs(N)
        dummy=Yesterday(dummy)
    enddo
endif
InNDays=dummy
end function InNDays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InXDays(date,X)

!#**********************************************************************
!#* Purpose: what will be the date in X days? (X is real!)
!#* mostly useful to convert back from continuous time to date
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified:11/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.X, >0 or <0
!#* OUT
!#*  1.InNDays
!#**********************************************************************
type(DateType), intent(in)::date
real(mrk), intent(in)::X
type(DateType):: InXDays
!locals
type(DateType):: dummy
integer(mik)::d,h,m
real(mrk)::rem

dummy=date
d=floor(X);rem=X-real(d,mrk)
dummy=InNDays(dummy,d)
h=floor(24._mrk*rem); rem=24._mrk*rem-real(h,mrk)
dummy=InNHours(dummy,h)
m=floor(60._mrk*rem); rem=60._mrk*rem-real(m,mrk)
dummy=InNMinutes(dummy,m)
dummy=InNSeconds(dummy,nint(60._mrk*rem))
InXDays=dummy
end function InXDays

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function InOneHour(date)

!#**********************************************************************
!#* Purpose: what will be the date in 1 hour?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:28/07/2011
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.date InOneHour
!#**********************************************************************
type(DateType), intent(in)::date
type(DateType):: InOneHour

InOneHour=date
if(date%hour<23) then
    InOneHour%hour=date%hour+1
else
    InOneHour=tomorrow(date)
    InOneHour%hour=0
endif

end function InOneHour

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function OneHourAgo(date)

!#**********************************************************************
!#* Purpose: what was the date 1 hour ago?
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:28/07/2011
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.date InOneHour
!#**********************************************************************
type(DateType), intent(in)::date
type(DateType):: OneHourAgo

OneHourAgo=date
if(date%hour>0) then
    OneHourAgo%hour=date%hour-1
else
    OneHourAgo=yesterday(date)
    OneHourAgo%hour=23
endif

end function OneHourAgo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InNHours(date,N)

!#**********************************************************************
!#* Purpose: what will be the date in N hours? 
!#* (or what was the date N hours ago if N<0)
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified: 28/07/2011
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.N, >0 or <0
!#* OUT
!#*  1.InNHours
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in)::N
type(DateType):: InNHours
!locals
type(DateType):: dummy
integer(mik)::i

dummy=date
if(N>0) then
    do i=1,N
        dummy=InOneHour(dummy)
    enddo
elseif(N<0) then
    do i=1,abs(N)
        dummy=OneHourAgo(dummy)
    enddo
endif
InNHours=dummy
end function InNHours

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function InOneMinute(date)

!#**********************************************************************
!#* Purpose: what will be the date in 1 minute?
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified: 10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.date InOneMinute
!#**********************************************************************
type(DateType), intent(in)::date
type(DateType):: InOneMinute
! local
type(DateType)::dummy

dummy=date
if(date%minute<59) then
    dummy%minute=date%minute+1
else
    dummy=InOneHour(date)
    dummy%minute=0
endif
InOneMinute=dummy
end function InOneMinute

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function OneMinuteAgo(date)

!#**********************************************************************
!#* Purpose: what was the date 1 minute ago?
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified: 10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.OneMinuteAgo
!#**********************************************************************
type(DateType), intent(in)::date
type(DateType):: OneMinuteAgo
! local
type(DateType)::dummy

dummy=date
if(date%minute>0) then
    dummy%minute=date%minute-1
else
    dummy=OneHourAgo(date)
    dummy%minute=59
endif
OneMinuteAgo=dummy
end function OneMinuteAgo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InNMinutes(date,N)

!#**********************************************************************
!#* Purpose: what will be the date in N minutes? 
!#* (or what was the date N minutes ago if N<0)
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified: 10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.N, >0 or <0
!#* OUT
!#*  1.InNMinutes
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in)::N
type(DateType):: InNMinutes
!locals
type(DateType):: dummy
integer(mik)::i

dummy=date
if(N>0) then
    do i=1,N
        dummy=InOneMinute(dummy)
    enddo
elseif(N<0) then
    do i=1,abs(N)
        dummy=OneMinuteAgo(dummy)
    enddo
endif
InNMinutes=dummy
end function InNMinutes

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function InOneSecond(date)

!#**********************************************************************
!#* Purpose: what will be the date in 1 Second?
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified: 10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.InOneSecond
!#**********************************************************************
type(DateType), intent(in)::date
type(DateType):: InOneSecond
! local
type(DateType)::dummy

dummy=date
if(date%second<59) then
    dummy%second=date%second+1
else
    dummy=InOneMinute(date)
    dummy%second=0
endif
InOneSecond=dummy
end function InOneSecond

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function OneSecondAgo(date)

!#**********************************************************************
!#* Purpose: what was the date 1 second ago?
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified: 10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#* OUT
!#*  1.OneSecondAgo
!#**********************************************************************
type(DateType), intent(in)::date
type(DateType):: OneSecondAgo
! local
type(DateType)::dummy

dummy=date
if(date%second>0) then
    dummy%second=date%second-1
else
    dummy=OneMinuteAgo(date)
    dummy%second=59
endif
OneSecondAgo=dummy
end function OneSecondAgo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InNSeconds(date,N)

!#**********************************************************************
!#* Purpose: what will be the date in N seconds? 
!#* (or what was the date N seconds ago if N<0)
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified: 10/05/2015
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*  1.date
!#*  2.N, >0 or <0
!#* OUT
!#*  1.InNSeconds
!#**********************************************************************
type(DateType), intent(in)::date
integer(mik), intent(in)::N
type(DateType):: InNseconds
!locals
type(DateType):: dummy
integer(mik)::i

dummy=date
if(N>0) then
    do i=1,N
        dummy=InOneSecond(dummy)
    enddo
elseif(N<0) then
    do i=1,abs(N)
        dummy=OneSecondAgo(dummy)
    enddo
endif
InNSeconds=dummy
end function InNSeconds

!--------------------------------------
! lists

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DateList(date1,date2,scale,scaleMult,list,err,message)
!#**********************************************************************
!#* Purpose: Generic subroutine for creating a list of dates, starting
!#*          from date1, and ending on or after date2.
!#*          Exemple: for a 15-minute time stepping, scale='m', scaleMult=15   
!#*          if date1>date2, the list starts from date2   
!#**********************************************************************
!#* Programmer: Ben Renard, Irstea Lyon
!#**********************************************************************
!#* Last modified: 10/05/2015
!#**********************************************************************
!#* Comments: 
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List: add step
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#*  3.scale: 'Y', 'M', 'D', 'h', 'm', or 's'
!#*  4.[scaleMult]
!#* OUT
!#*  1.list
!#**********************************************************************
type(DateType), intent(in)::date1, date2
character(*), intent(in)::scale
integer(mik), intent(in),optional::scaleMult
type(DateType), pointer::list(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
!locals
character(*),parameter::procnam='DateList'
integer(mik),parameter::option=1,scaleMult_def=1
type(DateType)::d1,d2
integer(mik)::i,n,mult
real(mrk)::duration

err=EXIT_SUCCESS; message=procnam//'/ok';
if(present(scaleMult)) then; mult=scaleMult; else; mult=scaleMult_def;endif
if(date1<=date2) then
    d1=date1;d2=date2
else
    d1=date2;d2=date1
endif
! workout number of steps
duration=DaysFromRef(d2,d1)
select case(scale)
case('Y')
    n=ceiling(duration/(real(mult,mrk)*365._mrk))+1
case('M')
    n=ceiling(duration/(real(mult,mrk)*30.4167_mrk))+1
case('D')
    n=ceiling(duration/real(mult,mrk))+1
case('h')
    n=ceiling( (24._mrk*duration)/real(mult,mrk))+1
case('m')
    n=ceiling( (1440._mrk*duration)/real(mult,mrk))+1
case('s')
    n=ceiling( (86400._mrk*duration)/real(mult,mrk))+1
case default
    err=EXIT_FAILURE;message='f-'//procnam//'/unknown[scale='//trim(scale)//']';return
end select
allocate(list(n))
! step
list(1)=d1
select case(scale)
case('Y')
    do i=2,n
        list(i)=InNYears(list(1),(i-1)*mult,option)
    enddo
case('M')
    do i=2,n
        list(i)=InNMonths(list(1),(i-1)*mult,option)
    enddo
case('D')
    do i=2,n
        list(i)=InNDays(list(i-1),mult)
    enddo
case('h')
    do i=2,n
        list(i)=InNHours(list(i-1),mult)
    enddo
case('m')
    do i=2,n
        list(i)=InNMinutes(list(i-1),mult)
    enddo
case('s')
    do i=2,n
        list(i)=InNSeconds(list(i-1),mult)
    enddo
case default
    err=EXIT_FAILURE;message='f-'//procnam//'/unknown[scale='//trim(scale)//']';return
end select

end subroutine DateList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function DailyList(date1,date2)
!#**********************************************************************
!#* Purpose: create a list of dates between date1 and date2
!#* if date1>date2, the list goes reverse
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:23/02/2010
!#**********************************************************************
!#* Comments: Hours, minutes and seconds ignored
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List: add step
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.DailyList
!#**********************************************************************
type(DateType), intent(in)::date1, date2
type(DateType)::DailyList(int(anint(DaysBetween(date1,date2)))+1)
!locals
integer(mik)::i,n

n=size(DailyList)
DailyList(1)=date1
if(date1<date2) then
    do i=1,n-1
        DailyList(i+1)=Tomorrow(DailyList(i))
    enddo
else
    do i=1,n-1
        DailyList(i+1)=Yesterday(DailyList(i))
    enddo
endif

end function DailyList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function HourlyList(date1,date2)
!#**********************************************************************
!#* Purpose: create a list of hourly dates between date1 and date2
!#* if date1>date2, the list goes reverse
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified: 27/07/2011
!#**********************************************************************
!#* Comments: minutes and seconds ignored
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List: add step
!#**********************************************************************
!#* IN
!#*  1.date1
!#*  2.date2
!#* OUT
!#*  1.HourlyList
!#**********************************************************************
type(DateType), intent(in)::date1, date2
type(DateType)::HourlyList(int(anint(24*DaysBetween(date1,date2)))+1)
!locals
integer(mik)::i,n

n=size(HourlyList)
HourlyList(1)=date1
if(date1<date2) then
    do i=1,n-1
        HourlyList(i+1)=InOneHour(HourlyList(i))
    enddo
else
    do i=1,n-1
        HourlyList(i+1)=OneHourAgo(HourlyList(i))
    enddo
endif

end function HourlyList

!--------------------------------------
! Formatting

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine FormatDate(date,fType,fDate,err,mess)

!^**********************************************************************
!^* Purpose: convert a date into a formated string; available types:
!^* 
!^* 1 -> 100 : Daily
!^* 1. 'YYYYMMDD'
!^* 2. 'DD/MM/YYYY'
!^* 3. 'day month year [French]'
!^* 4. 'day month year [English]'
!^* 5. 'DDMMYYYY'
!^* 
!^* 101 -> 200 : Hourly
!^* 101. 'YYYYMMDDhh'
!^* 
!^* 201 -> 300 : Minutely
!^* 
!^* 301 -> 400 : Secondly
!^* 301. 'DD/MM/YYYY hh:mm:ss'
!^* 302. 'DD/MM/YY hh:mm:ss'
!^* 
!^* 401 -> 500 : Monthly
!^* 401. 'YYYYMM' 
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:25/07/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: extend available types
!^**********************************************************************
!^* IN
!^*        1. date
!^*        2. format type - see list above
!^* OUT
!^*        1.fdate; string containing the formatted date (advise: trim it before use)
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
type(DateType), intent(in)::date
integer(mik),intent(in)::ftype
character(*), intent(out)::fdate
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(10), parameter,dimension(12)::FrenchName=(/'Janvier   ','Fevrier   ','Mars      ',&
                                                     'Avril     ','Mai       ','Juin      ',&
                                                     'Juillet   ','Aout      ','Septembre ',&
                                                     'Octobre   ','Novembre  ','Decembre  '/)
character(10), parameter,dimension(12)::EnglishName=(/'January   ','February  ','March     ',&
                                                      'April     ','May       ','June      ',&
                                                      'July      ','August    ','September ',&
                                                      'October   ','November  ','December  '/)
character(2)::day,month
character(4)::year

err=0;mess='';fdate=''
if( (.not. IsDateValid(date)) ) then
    err=1;mess='FormatDate:Fatal:Invalid [date]';return
endif

select case(ftype)
case(1) !'yyyymmdd'
    write(fdate,'(I8)') date%year*10000+date%month*100+date%day
    fdate=ReplaceSpaceByZeros(trim(fdate))
case(2) !'dd/mm/yyyy'
    1002 format(I2,'/',I2,'/',I4)
    write(fdate,1002) date%day,date%month,date%year
    fdate=ReplaceSpaceByZeros(trim(fdate))
case(3) !'day month year [French]'
    write(year,'(I4)') date%year
    write(day,'(I2)') date%day
    write(fdate,*) day//' '//trim(FrenchName(date%month))//' '//year
case(4) !'day month year [English]'
    write(year,'(I4)') date%year
    write(day,'(I2)') date%day
    write(fdate,*) day//' '//trim(EnglishName(date%month))//' '//year
case(5) !'ddmmyyyy'
    write(year,'(I4)') date%year
    write(month,'(I2)') date%month
    write(day,'(I2)') date%day
    write(fdate,'(A8)') day//month//year
    fdate=ReplaceSpaceByZeros(trim(fdate))
case(101) !'yyyymmddhh'
    write(fdate,'(I10)') date%year*1000000+date%month*10000+date%day*100+date%hour
    fdate=ReplaceSpaceByZeros(trim(fdate))
case(301) !'dd/mm/yyyy hh:mm:ss'
    1301 format(I2,'/',I2,'/',I4,' ',I2,':',I2,':',I2)
    write(fdate,1301) date%day,date%month,date%year,date%hour,date%minute,date%second
    fdate=ReplaceSpaceByZeros(trim(fdate))
    write(fdate(11:11),'(A1)') ' '
case(302) !'dd/mm/yy hh:mm:ss'
    1302 format(I2,'/',I2,'/',I2,' ',I2,':',I2,':',I2)
    write(fdate,1302) date%day,date%month,&
                      date%year-100*(date%year/100),&
                      date%hour,date%minute,date%second
    fdate=ReplaceSpaceByZeros(trim(fdate))
    write(fdate(9:9),'(A1)') ' '
case(401) !'yyyymm'
    write(fdate,'(I6)') date%year*100+date%month
    fdate=ReplaceSpaceByZeros(trim(fdate))
case default
    err=1;mess='FormatDate:Fatal:Unavailable [fType]';return
end select

end subroutine FormatDate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine UnformatDate(fDate,fType,date,err,mess)

!^**********************************************************************
!^* Purpose: convert a formated string into a date object
!^*          See sub FormatDate for avalailable fType
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:14/05/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: extend available types
!^**********************************************************************
!^* IN
!^*        1. fdate; string containing the formatted date
!^*        2. format type - see list above
!^* OUT
!^*        1.date; date object
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:string_number

character(*), intent(in)::fdate
integer(mik),intent(in)::ftype
type(DateType), intent(out)::date
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::today(8),incompleteYear,currentYear

err=0;mess='';

select case(ftype)
case(1) !'yyyymmdd'
    call string_number(fDate(1:4),date%year,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(5:6),date%month,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(7:8),date%day,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
case(2) !'dd/mm/yyyy'
    call string_number(fDate(7:10),date%year,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(4:5),date%month,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(1:2),date%day,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
case(3) !'day month year [French]'
    err=1;mess='UnformatDate:Fatal:unsupported [fType==3]';return
case(4) !'day month year [English]'
    err=1;mess='UnformatDate:Fatal:unsupported [fType==4]';return
case(5) !'ddmmyyyy'
    call string_number(fDate(5:8),date%year,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(3:4),date%month,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(1:2),date%day,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
case(101) !'yyyymmddhh'
    call string_number(fDate(1:4),date%year,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(5:6),date%month,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(7:8),date%day,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(9:10),date%hour,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
case(301) !'dd/mm/yyyy hh:mm:ss'
    call string_number(fDate(7:10),date%year,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(4:5),date%month,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(1:2),date%day,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(12:13),date%hour,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(15:16),date%minute,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(18:19),date%second,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
case(302) !'dd/mm/yy hh:mm:ss'
    call DATE_AND_TIME(values=today);currentYear=today(1)
    call string_number(fDate(7:8),incompleteYear,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    if(incompleteYear <= currentYear-100*(currentYear/100) )then
        date%year=100*(currentYear/100)+incompleteYear
    else
        date%year=(100*(currentYear/100)-100)+incompleteYear
    endif
    call string_number(fDate(4:5),date%month,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(1:2),date%day,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(10:11),date%hour,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(13:14),date%minute,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(16:17),date%second,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
case(401) !'yyyymm'
    call string_number(fDate(1:4),date%year,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
    call string_number(fDate(5:6),date%month,err); if (err/=0) then;mess='UnformatDate:Fatal:conversion failed';return;endif
case default
    err=1;mess='UnformatDate:Fatal:Unavailable [fType]';return
end select

end subroutine UnformatDate

!--------------------------------------
! Misc.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine FindExtremalDates(DateList,oldest,newest,err,mess)
!^**********************************************************************
!^* Purpose: find oldest and newest dates
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
!^*        1. DateList, vector of dates
!^* OUT
!^*        1.oldest, newest
!^*        2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        3.mess, error message
!^**********************************************************************
type(DateType), intent(in)::DateList(:)
type(DateType), intent(out)::oldest,newest
integer(mik), intent(out)::err
character(*), intent(out)::mess
! local
integer(mik)::Nx,i

err=0;mess=''
Nx=size(DateList)
oldest=DateList(1);newest=DateList(1)
do i=2,Nx
    if(DateList(i)>newest) newest=DateList(i)
    if(DateList(i)<oldest) oldest=DateList(i)
enddo
end subroutine FindExtremalDates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetCommonDates(d1,d2,d,i1,i2,err,mess)
!^**********************************************************************
!^* Purpose: Get Dates common to 2 date lists
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 22/11/2013
!^**********************************************************************
!^* Comments: Dates are assumes ordered to improve efficiency
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*        1. d1,d2, vectors of dates
!^* OUT
!^*        1.d, common dates
!^*        2.i1,i2, indexes of common dates for both series
!^*        3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*        4.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:ifirstTrueLoc
type(DateType), intent(in)::d1(:),d2(:)
type(DateType),pointer::d(:)
integer(mik),pointer::i1(:),i2(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! local
integer(mik)::n1,n2,i,j,m,k,di1(size(d1)),di2(size(d2))
logical::gotit

err=0;mess=''
n1=size(d1);n2=size(d2)
k=0;j=1;gotit=.false.
do i=1,n1
    do m=j,n2
        if(d1(i)==d2(m)) then
            j=m;gotit=.true.;exit
        endif
    enddo
    if(gotit) then
        k=k+1
        di1(k)=i;di2(k)=j
    endif
    gotit=.false.
enddo
if(associated(i1)) nullify(i1);allocate(i1(k))
if(associated(i2)) nullify(i2);allocate(i2(k))
if(associated(d)) nullify(d);allocate(d(k))
i1=di1(1:k);i2=di2(1:k)
d=d1(i1)
end subroutine GetCommonDates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function ClosestFeasibleDate(date)
! look backward for the closest feasible date
type(DateType), intent(in)::date
type(DateType)::ClosestFeasibleDate
! locals
type(DateType)::dummy

dummy=date
select case (date%month)
case(1,3,5,7,8,10,12)
    if (date%day>31) then
        dummy%day=31;dummy%hour=23;dummy%minute=59;dummy%second=59;
    endif
case(4,6,9,11)
    if (date%day>30) then
        dummy%day=30;dummy%hour=23;dummy%minute=59;dummy%second=59;
    endif
case(2)
    if(IsLeapYear(date%year)) then
        if (date%day>29) then
            dummy%day=29;dummy%hour=23;dummy%minute=59;dummy%second=59;
        endif
    else
        if (date%day>28) then
            dummy%day=28;dummy%hour=23;dummy%minute=59;dummy%second=59;
        endif
    endif
end select
ClosestFeasibleDate=dummy
end function ClosestFeasibleDate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!--------------------------------------
! internal private subs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure function ReplaceSpaceByZeros(string)
character(*),intent(in)::string
character(len(string))::ReplaceSpaceByZeros
integer(mik)::i

ReplaceSpaceByZeros=string
i=scan(ReplaceSpaceByZeros,' ')
do while (i/=0)
    write(ReplaceSpaceByZeros(i:i),'(I1)') 0
    i=scan(ReplaceSpaceByZeros,' ')
enddo

endfunction ReplaceSpaceByZeros
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Dates_tools

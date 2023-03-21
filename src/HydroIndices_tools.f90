module HydroIndices_tools

!~**********************************************************************
!~* Purpose: Various tools for extracting indices (e.g. annual mean,
!~* annual max, etc.) from data
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified: 15/09/2009
!~**********************************************************************
!~* Comments: mashup Antoine's nlowd with Ben's volumes as millions m3
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

implicit none
Private
public :: AnnualIndices_daily, SPAIndices_daily, BFI_daily,EmpiricalQuantile,&
          GetEvents_daily,GetSPAThreshold,GetLowFlowThreshold

Contains

Subroutine SPAIndices_daily(xin,SPAparameter,& !NON-SPA-filtred times series, SPA parameter (qzero)
                            threshold,& ! threshold above which to consider a SPA bump - default 0
                            Maxi,MaxiDate,& ! Bumb's maximum and its date
                            D,Dup,Ddown,& ! total duration + durations of the upward/downward phases
                            V,& ! Volume of the Bump
                            COMstart,COMcenter,COMend,& ! centers of mass
                            COMstartP,COMcenterP,COMendP,& ! p for COM definition above - defaults 10%,50%,90%
                            IED,& !Inter-Event Duration - Read WARNING below, handling of missing values problematic!
                            SPAfile, & !Name of the file to write the SPA series if requested
                            err,mess)
!^**********************************************************************
!^* Purpose: extract SPA indices from daily time series
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:23/03/2009
!^**********************************************************************
!^* Comments: MV problematic for SPA indices, because when a MV run occurs
!^* during a SPA bump, one can not say wether the bump corresponds to one
!^* or more events. Strategy used here is to discard all bumps containing
!^* a MV, but then the Inter-Event Duration may be wrong...
!^* So IED may be a biased indice for a series with a lot of MV
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin, input time series (NON-SPA FILTRED - filtring is done within this sub);
!^*    2.SPAparameter, qzero for applying the SPA filter
!^*    3.<threshold>, threshold above which to consider a SPA bump - default 0
!^*    4.<COMstartP>, p for COM definition - default 10%
!^*    5.<COMcenterP>, p for COM definition - default 50%
!^*    6.<COMendP>, p for COM definition - default 90%
!^*    7.<SPAFile>, File where the SPA-series can be written
!^* OUT
!^*    1. <Maxi>
!^*    2. <MaxiDate>
!^*    3. <D>
!^*    4. <Dup>
!^*    5. <Ddown>
!^*    6. <V>
!^*    7. <COMstart>
!^*    8. <COMcenter>
!^*    9. <COMend>
!^*    10. <IED>
!^*    11. err, error code; <0:Warning, ==0:OK, >0: Error
!^*    12. mess, error message
!^**********************************************************************
use TimeSeries_tools
use Dates_tools

type(OneDTimeSeriesType), intent(in)::xin
real(mrk), intent(in)::SPAparameter
real(mrk), intent(in),optional:: threshold, COMstartP, COMcenterP,COMendP
character(*), intent(in), optional::SPAFile
real(mrk), pointer, optional:: Maxi(:),D(:),Dup(:),Ddown(:),V(:),IED(:)
type(DateType),pointer,optional :: MaxiDate(:),COMstart(:),COMcenter(:),COMend(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
real(mrk):: thresh, COMsp, COMcp, COMep, cf, objvalue, cumul
real(mrk), allocatable::QJ(:)
type(DateType), allocatable::QJDate(:), MaxoDate(:)
integer(mik)::Nevents,i,j, t0, tn, k(1),FileUnit
type(OneDTimeSeriesType)::xSPA
integer(mik), pointer:: IndxEvents(:,:)


err=0;mess='';
!convertion factor to compute volumes in 10^6 cubic meters
cf=0.000001_mrk*60._mrk*60._mrk*24._mrk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART I: Handle default values for optional inputs !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! threshold
if(present(threshold)) then
    thresh=threshold;
else
    thresh=0._mrk
endif
! COMstartP
if(present(COMstartP)) then
    COMsp=COMstartP;
else
    COMsp=0.1_mrk
endif
! COMcenterP
if(present(COMcenterP)) then
    COMcp=COMcenterP;
else
    COMcp=0.5_mrk
endif
! COMendP
if(present(COMendP)) then
    COMep=COMendP;
else
    COMep=0.9_mrk
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART II: Apply SPA and get number of events       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply SPA
call TSFiltering_daily(filter=SPA,xin=xin,param=(/SPAparameter/), xout=xSPA,err=err,mess=mess)
if(err>0) then
    mess='SPAindices_Daily:FATAL:'//trim(mess);return
endif
if(present(SPAfile)) then
    FileUnit=RandomUnitNumber()
    open(unit=FileUnit,file=SPAfile,status='REPLACE')
    write(FileUnit,*) 'Year  Month Day   SPA'
    do i=1,xSPA%n
        write(FileUnit,'(I6,I6,I6,F15.3)') xSPA%ts(i)%date%year,xSPA%ts(i)%date%month,&
                                           xSPA%ts(i)%date%day,xSPA%ts(i)%q
    enddo
    close(FileUnit)
endif
! Compute number of events
call GetEvents_daily(x=xSPA,threshold=thresh,IsAbove=.true.,&
                     Nevents=Nevents,IndxEvents=IndxEvents,err=err,mess=mess)
if(err>0) then
    mess='SPAindices_Daily:FATAL:'//trim(mess);return
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART III: Allocation of optional outputs    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(present(Maxi)) then
    if(associated(Maxi)) nullify(Maxi)
    allocate(Maxi(Nevents))
    Maxi=xin%mv
endif
if(present(MaxiDate)) then
    if(associated(MaxiDate)) nullify(MaxiDate)
    allocate(MaxiDate(Nevents))
endif
if(present(D)) then
    if(associated(D)) nullify(D)
    allocate(D(Nevents))
    D=xin%mv
endif
if(present(Dup)) then
    if(associated(Dup)) nullify(Dup)
    allocate(Dup(Nevents))
    Dup=xin%mv
endif
if(present(Ddown)) then
    if(associated(Ddown)) nullify(Ddown)
    allocate(Ddown(Nevents))
    Ddown=xin%mv
endif
if(present(Ddown)) then
    if(associated(Ddown)) nullify(Ddown)
    allocate(Ddown(Nevents))
    Ddown=xin%mv
endif
if(present(V)) then
    if(associated(V)) nullify(V)
    allocate(V(Nevents))
    V=xin%mv
endif
if(present(COMstart)) then
    if(associated(COMstart)) nullify(COMstart)
    allocate(COMstart(Nevents))
endif
if(present(COMcenter)) then
    if(associated(COMcenter)) nullify(COMcenter)
    allocate(COMcenter(Nevents))
endif
if(present(COMend)) then
    if(associated(COMend)) nullify(COMend)
    allocate(COMend(Nevents))
endif
if(present(IED)) then
    if(associated(IED)) nullify(IED)
    allocate(IED(Nevents-1))
    IED=xin%mv
endif

!!!!!!!!!!!!!!!!!
! PART IV: GO!  !
!!!!!!!!!!!!!!!!!

if(allocated(MaxoDate)) deallocate(MaxoDate)
allocate(MaxoDate(Nevents))

do i=1,Nevents
    t0=IndxEvents(i,1);tn=IndxEvents(i,2)
    if(allocated(QJ)) deallocate(QJ)
    allocate(QJ(tn-t0+1))
    if(allocated(QJDate)) deallocate(QJDate)
    allocate(QJDate(tn-t0+1))
    QJ=xSPA%ts(t0:tn)%q
    QJDate=xSPA%ts(t0:tn)%date
    k=maxloc(QJ)
    MaxoDate(i)=QJDate(k(1))

    if(present(Maxi)) maxi(i)=maxval(QJ)
    if(present(MaxiDate)) MaxiDate(i)=MaxoDate(i)
    if(present(D)) D(i)=tn-t0+1
    if(present(Dup)) Dup(i)=k(1)
    if(present(Ddown)) Ddown(i)=(tn-t0+1)-k(1)
    if(present(V)) V(i)=sum(QJ)*cf
    if(present(COMstart)) then
        objvalue=sum(QJ)*COMsp
        cumul=0._mrk
        do j=1,tn-t0+1
            cumul=cumul+QJ(j)
            if(cumul>=objvalue) exit
        enddo
        COMstart(i)=QJDate(j)
    endif
    if(present(COMcenter)) then
        objvalue=sum(QJ)*COMcp
        cumul=0._mrk
        do j=1,tn-t0+1
            cumul=cumul+QJ(j)
            if(cumul>=objvalue) exit
        enddo
        COMcenter(i)=QJDate(j)
    endif
    if(present(COMend)) then
        objvalue=sum(QJ)*COMep
        cumul=0._mrk
        do j=1,tn-t0+1
            cumul=cumul+QJ(j)
            if(cumul>=objvalue) exit
        enddo
        COMend(i)=QJDate(j)
    endif
enddo

if(present(IED)) then
    do i=2,Nevents
        IED(i-1)=Date2Days(MaxoDate(i)-MaxoDate(i-1))
        ! Correct by removing number of days with MV in between two successive SPA bumps
        ! this is precisely where the handling of MV is problematic...
        IED(i-1)=IED(i-1) - count(xSPA%ts(IndxEvents(i-1,2)+1:IndxEvents(i,1)-1)%q == xSPA%mv )
    enddo
endif


!clean up the mess
if(associated(IndxEvents)) nullify(IndxEvents)
if(allocated(QJ)) deallocate(QJ)
if(allocated(QJDate)) deallocate(QJDate)
if(allocated(MaxoDate)) deallocate(MaxoDate)

end Subroutine SPAIndices_daily



subroutine AnnualIndices_daily(xin,&
                               YearStart,YearEnd, & !define season of interest (IN)
                               AMin,AMinDate,& !Annual minima + dates
                               AMax,AMaxDate,& !Annual minima + dates
                            !Threshold Indices
                               lowT,lowD,lowV,& ! low flow thresholds (IN)+ duration & volume below the threshold
                               lowP,lowComD,lowComV,& ! p-values (IN) + related Centers of mass in duration & volume. Computed for all lowT values
                               hiT,hiD,hiV,& ! high flow thresholds (IN)+ duration & volume above the threshold
                               hiP,hiComD,hiComV,& ! p-values (IN) + related Centers of mass in duration & volume. Computed for all hiT values
                            !End Threshold Indices
                               Total, & ! Annual Volume (WARNING: there is a conversion coefficient that only make sense for daily runoff in m3/s)
                               Modul,& ! Module (i.e. average computed for non-missing days of the season)
                               pval, Aquant, & ! p-values (IN) + related quantiles of the annual Flow Duration Curve
                               nlowd,nhid,& ! number of periods under / above given threshold, with associated max & min number of days
                               Pmv,& !percentage of missing values (% of season length)
                               YearList,& !convention: year from 10/1998 to 03/1999 is 1998
                               err,mess)

!^**********************************************************************
!^* Purpose: extract annual indices from daily time series
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:06/02/2009
!^**********************************************************************
!^* Comments: change yearlist from real to integer pointer
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. xin, input time series;
!^*    2.<YearStart>, date of the start of the year; default 01/01 00:00:00
!^*    3.<YearEnd>, date of the end of the year; default 31/12 23:59:59
!^*    4.<lowT>, low flow thresholds; default is HugeRe
!^*    5.<hiT>, high flow thresholds; default is -HugeRe
!^*    6.<lowP>, pvalues for low flows Centers of Mass; default is 0.5
!^*    7.<hiP>, pvalues for high flows Centers of Mass; default is 0.5
!^*    8.<pval>, p-values used for annual Flow duration curve quantiles. Default is 0.5
!^* OUT
!^*    1. <Amin>
!^*    2. <AminDate>
!^*    3. <Amax>
!^*    4. <AmaxDate>
!^*    5. <lowD>
!^*    6. <lowV>
!^*    7. <lowComD>
!^*    8. <lowComV>
!^*    9. <hiD>
!^*    10. <hiV>
!^*    11. <hiComD>
!^*    12. <hiComV>
!^*    13. <Total>
!^*    14. <Modul>
!^*    15. <Aquant>
!^*    16. <Pmv>
!^*    17. <YearList>
!^*    18. err, error code; <0:Warning, ==0:OK, >0: Error
!^*    19. mess, error message
!^**********************************************************************
use TimeSeries_tools, only: OneDTimeSeriesType
use Dates_tools

type(OneDTimeSeriesType), intent(in)::xin
type(DateType), intent(in),optional::YearStart, YearEnd
real(mrk), intent(in),optional::lowT(:),hiT(:), lowP(:), hiP(:), pval(:)
real(mrk), pointer, optional:: AMin(:),AMax(:),lowD(:,:),lowV(:,:),&
                               hiD(:,:),hiV(:,:),Total(:),Modul(:),&
                               Aquant(:, :), Pmv(:)
integer(mik), pointer, optional:: YearList(:), nlowd(:,:,:), nhid(:,:,:)
type(DateType), pointer,optional :: AMinDate(:),AMaxDate(:),&
                                    lowComD(:,:,:),lowComV(:,:,:),&
                                    hiComD(:,:,:),hiComV(:,:,:)
integer(mik), intent(out)::err
 character(*),intent(out)::mess
!locals
type(OneDTimeSeriesType) :: QJ
type(DateType)::Ystart,Yend,minDate, maxDate
real(mrk)::YLength, objvalue, currentvalue
real(mrk), allocatable::hT(:),lT(:),hP(:),lP(:),&
                        lowDur(:),lowVol(:),&
                        hiDur(:),hiVol(:),&
                        pv(:),quant(:)
type(DateType), allocatable::lowComDur(:,:),lowComVol(:,:),hiComDur(:,:),hiComVol(:,:)
integer(mik)::n,i,j,k,counter,nSeason,ndata,ndata2,dmin(1),dmax(1),nT, &
              nLowT, nHiT, nLowP, nHiP, nPval
integer(mik), pointer::lper(:), nld(:,:), hper(:), nhd(:,:)
real(mrk)::mini,maxi,totVol,mod,dmv,cf ! conversion factor
integer(mik), allocatable::Ylist(:)
logical,allocatable::mask(:)

err=0;mess='';

!convertion factor to compute volumes in 10^6 cubic meters
cf=0.000001_mrk*60._mrk*60._mrk*24._mrk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART I: Handle default values for optional inputs !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Start & End of the season
if(present(YearEnd)) then
    YEnd=YearEnd;
else
    YEnd=DateType(2,12,31,23,59,59)
endif

if(present(YearStart)) then
    Ystart=YearStart
else
    Ystart=DateType(2,1,1,0,0,0)
endif

! List of low-flows thresholds
if(present(lowT)) then
    nLowT=size(lowT)
    if(allocated(lT)) deallocate(lT)
    allocate(lT(nLowT))
    lT=lowT
else
    nLowT=1
    if(allocated(lT)) deallocate(lT)
    allocate(lT(nLowT))
    lT=HugeRe
endif
if (allocated(lowDur)) deallocate(lowDur)
if (allocated(lowVol)) deallocate(lowVol)
allocate(lowDur(nLowT),lowVol(nLowT))
if (associated(lper)) nullify(lper)
if (associated(nld)) nullify(nld)
allocate(lper(nlowT), nld(floor(date2days(Yend .Yminus. Ystart)),nlowT))

! List of low-flows pvals for centers of mass
if(present(lowP)) then
    nLowP=size(lowP)
    if(allocated(lP)) deallocate(lP)
    allocate(lP(nLowP))
    lP=lowP
else
    nLowP=1
    if(allocated(lP)) deallocate(lP)
    allocate(lP(nLowP))
    lP=0.5_mrk
endif
if (allocated(lowComDur)) deallocate(lowComDur)
if (allocated(lowComVol)) deallocate(lowComVol)
allocate(lowComDur(nLowT,nLowP),lowComVol(nLowT,nLowP))

! List of high-flows thresholds
if(present(hiT)) then
    nHiT=size(HiT)
    if(allocated(hT)) deallocate(hT)
    allocate(hT(nHiT))
    hT=hiT
else
    nHiT=1
    if(allocated(hT)) deallocate(hT)
    allocate(hT(nHiT))
    hT=-HugeRe
endif
if (allocated(hiDur)) deallocate(hiDur)
if (allocated(hiVol)) deallocate(hiVol)
allocate(hiDur(nhiT),hiVol(nhiT))
if (associated(hper)) nullify(hper)
if (associated(nhd)) nullify(nhd)
allocate(hper(nHiT), nhd(floor(date2days(Yend .Yminus. Ystart)),nHiT))

! List of high-flows pvals for centers of mass
if(present(hiP)) then
    nhiP=size(hiP)
    if(allocated(hP)) deallocate(hP)
    allocate(hP(nhiP))
    hP=hiP
else
    nhiP=1
    if(allocated(hP)) deallocate(hP)
    allocate(hP(nhiP))
    hP=0.5_mrk
endif
if (allocated(hiComDur)) deallocate(hiComDur)
if (allocated(hiComVol)) deallocate(hiComVol)
allocate(hiComDur(nhiT,nhiP),hiComVol(nhiT,nhiP))

! list of p-values for annual quantiles of the flow duration curve
if(present(pval)) then
    nPval=size(pval)
    if(allocated(pv)) deallocate(pv)
    allocate(pv(nPval))
    pv=pval
else
    nPval=1
    if(allocated(pv)) deallocate(pv)
    allocate(pv(nPval))
    pv=0.5_mrk
endif
if (allocated(quant)) deallocate(quant)
allocate(quant(nPval))






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART II: Preliminaries. Evaluate number of studied years !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n=xin%n
if(allocated(mask)) deallocate(mask)
allocate(mask(n))
Ylength=Date2Days(Yend.Yminus.Ystart) ! season length
! list of years encountered in xin, extended one year before/after to deal with incomplete seasons
nSeason=xin%ts(n)%date%year-xin%ts(1)%date%year+3
if (allocated(Ylist)) deallocate(Ylist)
allocate(Ylist(nSeason))
Ylist=(/(xin%ts(1)%date%year+i,i=-1,xin%ts(n)%date%year-xin%ts(1)%date%year+1)/)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART III: Preliminaries. Allocation of optional outputs  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(present(AMin)) then
    if(associated(Amin)) nullify(Amin)
    allocate(Amin(nSeason))
    Amin=xin%mv
endif
if(present(AminDate)) then
    if(associated(AminDate)) nullify(AminDate)
    allocate(AminDate(nSeason))
endif
if(present(Amax)) then
    if(associated(Amax)) nullify(Amax)
    allocate(Amax(nSeason))
    Amax=xin%mv
endif
if(present(AmaxDate)) then
    if(associated(AmaxDate)) nullify(AmaxDate)
    allocate(AmaxDate(nSeason))
endif
if(present(lowD)) then
    if(associated(lowD)) nullify(lowD)
    allocate(lowD(nSeason,nLowT))
    LowD=xin%mv
endif
if(present(lowV)) then
    if(associated(lowV)) nullify(lowV)
    allocate(lowV(nSeason, nLowT))
    LowV=xin%mv
endif
if(present(lowComD)) then
    if(associated(lowComD)) nullify(lowComD)
    allocate(lowComD(nSeason,nLowT,nLowP))
endif
if(present(lowComV)) then
    if(associated(lowComV)) nullify(lowComV)
    allocate(lowComV(nSeason,nLowT,nLowP))
endif
if(present(HiD)) then
    if(associated(HiD)) nullify(HiD)
    allocate(HiD(nSeason, nHiT))
    HiD=xin%mv
endif
if(present(hiV)) then
    if(associated(hiV)) nullify(hiV)
    allocate(hiV(nSeason,nHiT))
    HiV=xin%mv
endif
if(present(hiComD)) then
    if(associated(hiComD)) nullify(hiComD)
    allocate(hiComD(nSeason,nhiT,nhiP))
endif
if(present(hiComV)) then
    if(associated(hiComV)) nullify(hiComV)
    allocate(hiComV(nSeason,nhiT,nhiP))
endif
if(present(Total)) then
    if(associated(Total)) nullify(Total)
    allocate(Total(nSeason))
    Total=xin%mv
endif
if(present(Modul)) then
    if(associated(Modul)) nullify(Modul)
    allocate(Modul(nSeason))
    Modul=xin%mv
endif
if(present(Aquant)) then
    if(associated(Aquant)) nullify(Aquant)
    allocate(Aquant(nSeason, nPval))
    Aquant=xin%mv
endif
if(present(pmv)) then
    if(associated(pmv)) nullify(pmv)
    allocate(pmv(nSeason))
    pmv=xin%mv
endif
if(present(YearList)) then
    if(associated(YearList)) nullify(YearList)
    allocate(YearList(nSeason))
    YearList=xin%mv
    YearList=Ylist
endif
if(present(nlowd)) then
    if (associated(nlowd)) nullify(nlowd)
    allocate(nlowd(nSeason,nlowT,3))
    nlowd=xin%mv
endif
if(present(nhid)) then
    if (associated(nhid)) nullify(nhid)
    allocate(nhid(nSeason,nHiT,3))
    nhid=xin%mv
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART IV: scan series and extract indices                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,nSeason
    mask=IsInSeason(xin%ts(:)%date,Ystart,Yend,Ylist(i))
    ndata=count(mask)
    if(ndata==0) then ! no data for current season - missing value for all indices
        mini=xin%mv;maxi=xin%mv;
        lowDur=xin%mv;lowVol=xin%mv;
        lowComDur=DateType(undefIN,1,1,0,0,0);lowComVol=DateType(undefIN,1,1,0,0,0);
        hiDur=xin%mv;hiVol=xin%mv;
        hiComDur=DateType(undefIN,1,1,0,0,0);hiComVol=DateType(undefIN,1,1,0,0,0);
        totVol=xin%mv;mod=xin%mv;
        quant=xin%mv
        nld(:,:)=undefIN;lper(:)=undefIN
        nhd(:,:)=undefIN;hper(:)=undefIN
        dmv=Ylength;minDate=DateType(undefIN,1,1,0,0,0);maxDate=DateType(undefIN,1,1,0,0,0)
    else
        if(associated(QJ%ts)) nullify(QJ%ts)
        allocate(QJ%ts(ndata))
        QJ%mv=xin%mv
        QJ%ts=pack(xin%ts,mask) !QJ=data for the current season
        ndata2=count(QJ%ts%q/=QJ%mv)
        if(ndata2==0) then ! whole season missing - missing value for all indices
            mini=xin%mv;maxi=xin%mv;
            lowDur=xin%mv;lowVol=xin%mv;
            lowComDur=DateType(undefIN,1,1,0,0,0);lowComVol=DateType(undefIN,1,1,0,0,0);
            hiDur=xin%mv;hiVol=xin%mv;
            hiComDur=DateType(undefIN,1,1,0,0,0);hiComVol=DateType(undefIN,1,1,0,0,0);
            totVol=xin%mv;mod=xin%mv;
            quant=xin%mv
            nld(:,:)=undefIN;lper(:)=undefIN
            nhd(:,:)=undefIN;hper(:)=undefIN
            dmv=Ylength;minDate=DateType(undefIN,1,1,0,0,0);maxDate=DateType(undefIN,1,1,0,0,0)
        else
            !Annual max + date
            maxi=maxval(QJ%ts%q,mask=(QJ%ts%q/=QJ%mv));
            dmax=maxloc(QJ%ts%q,mask=(QJ%ts%q/=QJ%mv));maxDate=QJ%ts(dmax(1))%date
            !Annual min + date
            mini=minval(QJ%ts%q,mask=(QJ%ts%q/=QJ%mv));
            dmin=minloc(QJ%ts%q,mask=(QJ%ts%q/=QJ%mv));minDate=QJ%ts(dmin(1))%date
            !Hi threshold indices
            do j=1, nHiT
                nT=count(QJ%ts%q>hT(j) .and. QJ%ts%q/=QJ%mv)
                if(nT==0) then
                    hiDur(j)=0._mrk;hiVol(j)=0._mrk
                else
                    hiDur(j)=nT;hiVol(j)=cf*sum(QJ%ts%q,mask=(QJ%ts%q>hT(j) .and. QJ%ts%q/=QJ%mv))
                endif
            enddo
            ! Consecutice period of highflow
            nhd(:,:)=0;hper(:)=0
            do j=1,nHiT
                If(QJ%ts(1)%q/=QJ%mv .and. QJ%ts(1)%q>hT(j))then
                hper(j)=1
                nhd(hper(j),j)=1
                endif
                do k=2,ndata
                    if(QJ%ts(k)%q/=QJ%mv .and. QJ%ts(k)%q>hT(j)) then
                        if (QJ%ts(k-1)%q==QJ%mv .and. hper(j)==0) then
                            hper(j)=1
                            nhd(hper(j),j)=1
                        elseif (QJ%ts(k-1)%q==QJ%mv .and. hper(j)>0) then
                            nhd(hper(j),j)=nhd(hper(j),j)+1
                        elseif (QJ%ts(k-1)%q/=QJ%mv .and. QJ%ts(k-1)%q<=hT(j)) then
                            hper(j)=hper(j)+1
                            nhd(hper(j),j)=1
                        elseif(QJ%ts(k-1)%q/=QJ%mv .and. QJ%ts(k-1)%q>hT(j)) then
                            nhd(hper(j),j)=nhd(hper(j),j)+1
                        endif
                    endif
                enddo
            enddo
            ! Center of mass in Volume for high flows
            do j=1, nHiT
                do k=1,nHiP
                    objValue=hP(k)*hiVol(j) ! objective value to be exceeded (e.g. 50% of total volume/duration)
                    if(objValue==0._mrk) then! ie threshold was not crossed this year - undefined CoM
                        HiComVol(j,k)=DateType(undefIN,1,1,0,0,0)
                    else
                        counter=1;
                        if(QJ%ts(counter)%q/=QJ%mv .and. QJ%ts(counter)%q>hT(j)) then
                            currentValue=cf*QJ%ts(counter)%q
                        else
                            currentValue=0._mrk
                        endif
                        do while(currentValue<objValue .and. counter<ndata)
                            counter=counter+1;
                            if(QJ%ts(counter)%q/=QJ%mv .and. QJ%ts(counter)%q>hT(j)) & !threshold crossed
                                currentValue=currentValue+cf*QJ%ts(counter)%q
                        enddo
                        HiComVol(j,k)=QJ%ts(counter)%date
                    endif
                enddo
            enddo
            ! Center of mass in Duration for high flows
            do j=1, nHiT
                do k=1,nHiP
                    objValue=hP(k)*hiDur(j) ! objective value to be exceeded (e.g. 50% of total volume/duration)
                    if(objValue==0._mrk) then! ie threshold was not crossed this year - undefined CoM
                        HiComDur(j,k)=DateType(undefIN,1,1,0,0,0)
                    else
                        counter=1;
                        if(QJ%ts(counter)%q/=QJ%mv .and. QJ%ts(counter)%q>hT(j)) then
                            currentValue=1._mrk
                        else
                            currentValue=0._mrk
                        endif
                        do while(currentValue<objValue .and. counter<ndata)
                            counter=counter+1;
                            if(QJ%ts(counter)%q/=QJ%mv .and. QJ%ts(counter)%q>hT(j)) & !threshold crossed
                                currentValue=currentValue+1._mrk
                        enddo
                        HiComDur(j,k)=QJ%ts(counter)%date
                    endif
                enddo
            enddo
            !Low threshold indices
            do j=1,nLowT
                nT=count(QJ%ts%q<lT(j) .and. QJ%ts%q/=QJ%mv)
                if(nT==0) then
                    lowDur(j)=0._mrk;lowVol(j)=0._mrk
                else
                    lowDur(j)=nT;lowVol(j)=cf*(nT*lT(j)-sum(QJ%ts%q,mask=(QJ%ts%q<lT(j) .and. QJ%ts%q/=QJ%mv)))
                endif
            enddo
            ! Consecutice period of low flow
            nld(:,:)=0;lper(:)=0
            do j=1,nLowT
                If(QJ%ts(1)%q/=QJ%mv .and. QJ%ts(1)%q<lT(j))then
                lper(j)=1
                nld(lper(j),j)=1
                endif
                do k=2,ndata
                    if(QJ%ts(k)%q/=QJ%mv .and. QJ%ts(k)%q<lT(j)) then
                        if (QJ%ts(k-1)%q==QJ%mv .and. lper(j)==0) then
                            lper(j)=1
                            nld(lper(j),j)=1
                        elseif (QJ%ts(k-1)%q==QJ%mv .and. lper(j)>0) then
                            nld(lper(j),j)=nld(lper(j),j)+1
                        elseif (QJ%ts(k-1)%q/=QJ%mv .and. QJ%ts(k-1)%q>=lT(j)) then
                            lper(j)=lper(j)+1
                            nld(lper(j),j)=1
                        elseif(QJ%ts(k-1)%q/=QJ%mv .and. QJ%ts(k-1)%q<lT(j)) then
                            nld(lper(j),j)=nld(lper(j),j)+1
                        endif
                    endif
                enddo
            enddo
            ! Center of mass in Volume for Low flows
            do j=1, nLowT
                do k=1,nLowP
                    objValue=lP(k)*lowVol(j) ! objective value to be exceeded (e.g. 50% of total volume/duration)
                    if(objValue==0._mrk) then! ie threshold was not crossed this year - undefined CoM
                        LowComVol(j,k)=DateType(undefIN,1,1,0,0,0)
                    else
                        counter=1;
                        if(QJ%ts(counter)%q/=QJ%mv .and. QJ%ts(counter)%q<lT(j)) then
                            currentValue=cf*(lT(j)-QJ%ts(counter)%q)
                        else
                            currentValue=0._mrk
                        endif
                        do while(currentValue<objValue .and. counter<ndata)
                            counter=counter+1;
                            if(QJ%ts(counter)%q/=QJ%mv .and. QJ%ts(counter)%q<lT(j)) & !threshold crossed
                                currentValue=currentValue+cf*(lT(j)-QJ%ts(counter)%q)
                        enddo
                        LowComVol(j,k)=QJ%ts(counter)%date
                    endif
                enddo
            enddo
            ! Center of mass in Duration for Low flows
            do j=1, nLowT
                do k=1,nLowP
                    objValue=lP(k)*lowDur(j) ! objective value to be exceeded (e.g. 50% of total volume/duration)
                    if(objValue==0._mrk) then! ie threshold was not crossed this year - undefined CoM
                        LowComDur(j,k)=DateType(undefIN,1,1,0,0,0)
                    else
                        counter=1;
                        if(QJ%ts(counter)%q/=QJ%mv .and. QJ%ts(counter)%q<lT(j)) then !season first day
                            currentValue=1._mrk
                        else
                            currentValue=0._mrk
                        endif
                        do while(currentValue<objValue .and. counter<ndata)
                            counter=counter+1;

                            if(QJ%ts(counter)%q/=QJ%mv .and. QJ%ts(counter)%q<lT(j)) & !threshold crossed
                                currentValue=currentValue+1._mrk
                        enddo
                        LowComDur(j,k)=QJ%ts(counter)%date
                    endif
                enddo
            enddo
            !Annual Volume and module
            nT=count(QJ%ts%q/=QJ%mv)
            if(nT==0) then
                totVol=0._mrk;mod=0._mrk
            else
                totVol=sum(QJ%ts%q,mask=(QJ%ts%q/=QJ%mv))
                mod=totVol/nT
                totVol=cf*totVol
            endif
            !Annual quantiles of the FDC
            do j=1, npVal
                quant(j)=EmpiricalQuantile(x=QJ%ts%q,p=pv(j),mask=(QJ%ts%q/=QJ%mv))
            enddo
            ! missing values
            dmv=max(0._mrk,Ylength-ndata)+count(QJ%ts%q==QJ%mv)
        endif
    endif
    ! put results in optional outputs
    if(present(AMin)) Amin(i)=mini
    if(present(AMinDate)) AminDate(i)=mindate
    if(present(AMax)) Amax(i)=maxi
    if(present(AMaxDate)) AmaxDate(i)=maxdate
    if(present(lowD)) lowD(i,:)=lowDur
    if(present(LowComD)) LowComD(i,:,:)=LowComDur
    if(present(lowV)) lowV(i,:)=lowVol
    if(present(lowComV)) lowComV(i,:,:)=lowComVol
    if(present(hiD)) HiD(i,:)=hiDur
    if(present(hiComD)) HiComD(i,:,:)=hiComDur
    if(present(hiV)) hiV(i,:)=hiVol
    if(present(hiComV)) HiComV(i,:,:)=hiComVol
    if(present(Total)) Total(i)=totVol
    if(present(Modul)) Modul(i)=mod
    if(present(nlowd)) then
        do j=1,nlowT
            nlowd(i,j,1)=lper(j);nlowd(i,j,2)=maxval(nld(:,j));nlowd(i,j,3)=minval(nld(:,j),mask=(nld(:,j)/=0))
            if (lper(j)==0) nlowd(i,j,3)=0
        enddo
    endif
    if(present(nhid)) then
        do j=1,nhiT
            nhid(i,j,1)=hper(j);nhid(i,j,2)=maxval(nhd(:,j));nhid(i,j,3)=minval(nhd(:,j),mask=(nhd(:,j)/=0))
            if (hper(j)==0) nhid(i,j,3)=0
        enddo
    endif
    if(present(Aquant)) Aquant(i,:)=quant
    if(present(pmv)) pmv(i)=(100._mrk*dmv)/(1._mrk*Ylength)
enddo

! Deallocate everything
deallocate(hT,lT,hP,lP,lowDur,lowVol,hiDur,hiVol,pv,quant,lowComDur,lowComVol,hiComDur,hiComVol,Ylist,mask,lper,nld)
if (associated(QJ%ts)) nullify(QJ%ts)

end subroutine AnnualIndices_daily

pure subroutine BFI_daily(xin,&
                          YearStart,YearEnd, & !define season of interest (IN)
                          blocLength,turningPoint,& ! optional tuning for BFS algorithm
                          BFI, & ! BFI series
                          Pmv,& !percentage of missing values (% of season length)
                          YearList,& !convention: year from 10/1998 to 03/1999 is 1998
                          err,mess)

!^**********************************************************************
!^* Purpose: computes the series of BFI from a daily time series
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 20/05/2009
!^**********************************************************************
!^* Comments:
! 02/07/09 modified YearList(:) from real pointer to integer pointer
!^**********************************************************************
!^* References: Van Lannen and Tallaksen, 2004
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.xin
!^*    2.<YearStart>, date of the start of the year; default 01/01 00:00:00
!^*    3.<YearEnd>, date of the end of the year; default 31/12 23:59:59
!^*    4.<blocLength>, bloc length for BFS algorithm
!^*    5.<turningpoint>, definition of a turning point in BFS algorithm
!^* OUT
!^*    1.<BFI>
!^*    2. <Pmv>
!^*    3. <YearList>
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************

use TimeSeries_tools, only: OneDTimeSeriesType,TSFiltering_daily, BFS
use Dates_tools

type(OneDTimeSeriesType), intent(in)::xin
type(DateType), intent(in),optional::YearStart, YearEnd
real(mrk), intent(in),optional::blocLength, TurningPoint
real(mrk), pointer, optional:: BFI(:),Pmv(:)
integer(mik), pointer, optional::YearList(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
type(OneDTimeSeriesType) :: QJ, QJBFS, xBFS
type(DateType)::Ystart,Yend
logical,allocatable::mask(:)
integer(mik), allocatable::Ylist(:)
real(mrk)::YLength, bfindx, dmv, vol, volBFS
integer(mik)::n,i,nSeason,ndata,ndata2

err=0;mess='';

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART I: Handle default values for optional inputs !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Start & End of the season
if(present(YearEnd)) then
    YEnd=YearEnd;
else
    YEnd=DateType(2,12,31,23,59,59)
endif

if(present(YearStart)) then
    Ystart=YearStart
else
    Ystart=DateType(2,1,1,0,0,0)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART II: Preliminaries. Evaluate number of studied years !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n=xin%n
if(allocated(mask)) deallocate(mask)
allocate(mask(n))
Ylength=Date2Days(Yend.Yminus.Ystart) ! season length
! list of years encountered in xin, extended one year before/after to deal with incomplete seasons
nSeason=xin%ts(n)%date%year-xin%ts(1)%date%year+3
if (allocated(Ylist)) deallocate(Ylist)
allocate(Ylist(nSeason))
Ylist=(/(xin%ts(1)%date%year+i,i=-1,xin%ts(n)%date%year-xin%ts(1)%date%year+1)/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART III: Preliminaries. Allocation of optional outputs  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(present(BFI)) then
    if(associated(BFI)) nullify(BFI)
    allocate(BFI(nSeason))
    BFI=xin%mv
endif
if(present(pmv)) then
    if(associated(pmv)) nullify(pmv)
    allocate(pmv(nSeason))
    pmv=xin%mv
endif
if(present(YearList)) then
    if(associated(YearList)) nullify(YearList)
    allocate(YearList(nSeason))
    YearList=Ylist
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PART IV: scan series and extract indices                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Apply BFS algorithm
if(present(blocLength) .and. present(TurningPoint)) then
    call TSFiltering_daily(filter=BFS,param=(/blocLength,TurningPoint/),&
                           xin=xin,xout=xBFS,err=err,mess=mess)
else
    call TSFiltering_daily(filter=BFS,xin=xin,xout=xBFS,err=err,mess=mess)
endif
if(err>0) then
    mess='subroutine BFI_daily:'//trim(mess);return
endif

!scan series
do i=1,nSeason
    mask=IsInSeason(xBFS%ts(:)%date,Ystart,Yend,Ylist(i))
    ndata=count(mask)
    if(ndata==0) then ! no data for current season - missing value for all indices
        bfindx=xin%mv;dmv=Ylength;
    else
        if(associated(QJ%ts)) nullify(QJ%ts)
        if(associated(QJBFS%ts)) nullify(QJBFS%ts)
        allocate(QJ%ts(ndata),QJBFS%ts(ndata))
        QJ%ts=pack(xin%ts,mask) !QJ=raw data for the current season
        QJBFS%ts=pack(xBFS%ts,mask) !QJ=BFSized data for the current season
        ndata2=count(QJBFS%ts%q/=QJBFS%mv)
        if(ndata2==0) then ! whole season missing - missing value for all indices
            bfindx=xin%mv;dmv=Ylength;
        else
            !Compute BFI
            Vol=sum(QJ%ts%q,mask=(QJBFS%ts%q/=QJBFS%mv))
            VolBFS=sum(QJBFS%ts%q,mask=(QJBFS%ts%q/=QJBFS%mv))
            if(Vol==0._mrk) then
                bfindx=xin%mv
            else
                bfindx=VolBFS/Vol
            endif
            ! missing values
            dmv=max(0._mrk,Ylength-ndata)+count(QJBFS%ts%q==QJBFS%mv)
        endif
    endif
    ! put results in optional outputs
    if(present(BFI)) BFI(i)=bfindx
    if(present(pmv)) pmv(i)=(100._mrk*dmv)/(1._mrk*Ylength)
enddo
! deallocate everything
if(associated(QJ%ts)) nullify(QJ%ts)
if(associated(QJBFS%ts)) nullify(QJBFS%ts)
if(associated(xBFS%ts)) nullify(xBFS%ts)
deallocate(mask,Ylist)

end subroutine BFI_daily


!!!!!!!!!!!
! PRIVATE !
!!!!!!!!!!!

subroutine sort(xin, mask,xout)
!^**********************************************************************
!^* Purpose: sort xin into xout
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:20/05/2009
!^**********************************************************************
!^ Comments: Some old piece code, might not be optimal
!^ Consider replacement with DMSL routines.
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:DMSLize
!^**********************************************************************
!^* IN
!^*    1. xin, input time series;
!^*    2. mask, logical mask. Values with mask=FALSE are excluded
!^* OUT
!^*    1. xout, sorted time series
!^**********************************************************************
real(mrk), intent(in)::xin(:)
logical, intent(in)::mask(:)
real(mrk), intent(out)::xout(:)
real(mrk), allocatable:: z(:)
integer(mik):: N, ind, k
real(mrk):: temp

N=count(mask)
if(allocated(z)) deallocate(z)
allocate(z(N))
z=pack(xin,mask)
do k=1, N
    ind=minval(minloc(z(k:N)))
    xout(k)=z(ind+k-1)
    temp=z(ind+k-1)
    z(ind+k-1)=z(k)
    z(k)=temp
end do
if(allocated(z)) deallocate(z)
end subroutine sort

function EmpiricalQuantile(x,p,mask)

!#**********************************************************************
!#* Purpose: empirical p-quantile of data x
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:19/05/2009
!#**********************************************************************
!#* Comments: Some old piece code, might not be optimal
!#* Consider replacement with DMSL routines
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List: Use DMSL
!#**********************************************************************
!#* IN
!#*    1.x, data vector
!#*    2.p, p-value
!#*    3.mask, logical mask
!#* OUT
!#*    1.the p-quantile
!#**********************************************************************

Real(mrk), intent(in)::x(:),p
logical, intent(in)::mask(:)
real(mrk), allocatable::x2(:)
real(mrk):: EmpiricalQuantile
integer(mik)::k, n

n=count(mask)
if(allocated(x2)) deallocate(x2)
allocate(x2(n))
call sort(x,mask,x2)
k= floor(n*p)
if(k<1) then
    EmpiricalQuantile=x2(1)
else if(k>n) then
    EmpiricalQuantile=x2(n)
else
    EmpiricalQuantile=x2(k)
endif
deallocate(x2)
end function EmpiricalQuantile


subroutine GetEvents_daily(x,threshold,IsAbove,Nevents,IndxEvents,err,mess)

!^**********************************************************************
!^* Purpose: Extract number of events and their start/end indexes from
!^* a daily time series
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:22/03/2010
!^**********************************************************************
!^* Comments: An event is defined as a continous series of days verifying:
!^*           x(t)>threshold if IsAbove==.true.
!^*           x(t)<threshold if IsAbove==.false.
!^* WARNING: This is NOT POT-sampling! in particular, no independence
!^* constraint is used. Mainly useful for SPA-filtred series
!^* HANDLING OF MISSING VALUES: any continous run of days verifying
!^* x(t)> (or <) threshold but with at least one missing value is NOT
!^* considered as an event
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, data vector
!^*    2.threshold, the threshold defining the event
!^*    3.IsAbove, event defined above or below the threshold?
!^* OUT
!^*    1.Nevents, number of events
!^*    2.IndxEvents, Nevents x 2 matrix of start/end indexes of the events
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use TimeSeries_tools, only: OneDTimeSeriesType

type(OneDTimeSeriesType), intent(in)::x
real(mrk), intent(in)::threshold
logical, intent(in)::IsAbove
integer(mik), pointer:: IndxEvents(:,:)
integer(mik), intent(out)::Nevents,err
character(*),intent(out)::mess
!locals
integer(mik)::n,i, startEvent, endEvent
integer(mik), allocatable:: DummyTab(:,:)
logical :: AnyMV, EventInProgress, condition

err=0;mess='';
n=x%n
If(Allocated(DummyTab)) deallocate(DummyTab)
allocate(DummyTab(n,2))
!Init
AnyMV=.false.
EventInProgress=.false.
Nevents=0

do i=1,n
    if(IsAbove) then
        condition=(x%ts(i)%q>threshold)
    else
        condition=(x%ts(i)%q<threshold)
    endif
    if( (x%ts(i)%q==x%mv) .and. EventInProgress) then
        AnyMV=.true.
        cycle
    endif
    if(condition) then !Event
        if(EventInProgress) then
            cycle
        else ! First day of the event
            if(i>=2) then !Check that the day before wasn't MV - otherwise event is discarded
                if(x%ts(i-1)%q==x%mv) then
                    AnyMV=.true.
                endif
            endif
            EventInProgress=.true.
            startEvent=i
            cycle
        endif
    else ! No event
        if(EventInProgress) then ! Previous time step was last day of the event
            endEvent=i-1
            EventInProgress=.false.
            if(AnyMV) then ! At least one MV during event - discard event
                AnyMV=.false.;cycle
            else ! no MV during event - genuine event
                Nevents=Nevents+1
                DummyTab(Nevents,1)=startEvent;DummyTab(Nevents,2)=endEvent;
                AnyMV=.false.
            endif
        else
            EventInProgress=.false.
            cycle
        endif
    endif
enddo

! store results in pointer
If(Associated(IndxEvents)) nullify(IndxEvents)
allocate(IndxEvents(Nevents,2))
IndxEvents(:,:)=DummyTab(1:Nevents,:)

deallocate(DummyTab)

end subroutine GetEvents_daily

subroutine GetSPAThreshold(x,Nobj,thresholdIN,LogFile,SPAthreshold, SPApercentile,err,mess)

!^**********************************************************************
!^* Purpose: Compute the threshold yielding a number of events as close
!^* as possible to Nobj, amongst 8 candidates percentiles
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:23/03/2010
!^**********************************************************************
!^* Comments: An event is defined as a continous series of days verifying:
!^*           x(t)>thresholdIN
!^* HANDLING OF MISSING VALUES: inherited from sub GetEvents_daily
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, data vector
!^*    2.Nobj, Objective number of events
!^*    3.[thresholdIN], threshold for defining an event. Default 0._mrk
!^*    4.[LogFile] File 2 write summary of the search
!^* OUT
!^*    1.SPAThreshold, SPA threshold yielding about Nobj events
!^*    2.SPApercentile, same as above but expressed in percentile
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use TimeSeries_tools
type(OneDTimeSeriesType), intent(in)::x
integer(mik), intent(in)::Nobj
real(mrk), intent(in), optional::thresholdIN
character(*), intent(in), optional:: LogFile
real(mrk), intent(out):: SPAthreshold, SPApercentile
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik), parameter::Ncandidates=8
real(mrk), dimension(Ncandidates), parameter:: candidatesP=&
    (/100._mrk/365.25_mrk, 1._mrk, 2.5_mrk, 5._mrk, 10._mrk, 15._mrk, 20._mrk, 25._mrk/)
real(mrk):: thresh, candidatesQ(Ncandidates)
integer(mik)::fileUnit, i, Nlist(Ncandidates), k(1)
type(OneDTimeSeriesType)::xSPA
integer(mik), pointer:: IndxEvents(:,:)

err=0;mess='';
! Handle defaults for optional arguments
if(present(thresholdIN)) then
    thresh=thresholdIN
else
    thresh=0._mrk
endif

!misc.
fileUnit=RandomUnitNumber()
if(present(LogFile)) open(unit=fileUnit, file=trim(LogFile), status='REPLACE')

! compute candidate percentiles
do i=1,Ncandidates
    candidatesQ(i)=EmpiricalQuantile(x=x%ts(:)%q,p=0.01_mrk*candidatesP(i),mask=(x%ts(:)%q /= x%mv))
enddo

! Get Number of events for each of them
do i=1, Ncandidates
    ! Apply SPA
    call TSFiltering_daily(filter=SPA,xin=x,param=(/candidatesQ(i)/), xout=xSPA,err=err,mess=mess)
    if(err>0) then
        mess='GetSPAThreshold:FATAL:'//trim(mess);return
    endif
    ! Compute number of events
    call GetEvents_daily(x=xSPA,threshold=thresh,IsAbove=.true.,&
                         Nevents=Nlist(i),IndxEvents=IndxEvents,err=err,mess=mess)
    if(err>0) then
        mess='GetSPAThreshold:FATAL:'//trim(mess);return
    endif
enddo

! Select threshold yielding N closest to Nobs
k=minloc(abs(Nlist-Nobj))
SPAthreshold=candidatesQ(k(1))
SPApercentile=candidatesP(k(1))

! write 2 logfile if supplied
if(present(LogFile)) then
    write(FileUnit,*) 'N objective = ',Nobj
    do i=1,Ncandidates
        write(FileUnit,*) 'Percentile (%): ',candidatesP(i),&
                          'Value : ',candidatesQ(i), &
                          'N events: ',Nlist(i)
    enddo
    write(FileUnit,*) 'Selected percentile (%)',SPApercentile,&
                      'Value : ',SPAthreshold, &
                      'Absolute diff. with Nobj: ',abs(Nlist(k(1))-Nobj)
endif

if(present(LogFile)) close(fileUnit)

end subroutine GetSPAThreshold

subroutine GetLowFlowThreshold(x,Nobj,pmv,LogFile,LFthreshold, LFpercentile,err,mess)

!^**********************************************************************
!^* Purpose: Compute the lowflow threshold yielding a number of years
!^* with (drought duration>0) AND (missing data<pmv) as close as possible
!^* to Nobj, amongst 8 candidates percentiles
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:23/03/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, data vector
!^*    2.Nobj, Objective number of events
!^*    3.pmv, limit percentage of mv before discarding the year
!^*    4.[LogFile] File 2 write summary of the search
!^* OUT
!^*    1.LFThreshold, LF threshold yielding about Nobj events
!^*    2.LFpercentile, same as above but expressed in percentile
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use TimeSeries_tools
type(OneDTimeSeriesType), intent(in)::x
integer(mik), intent(in)::Nobj
real(mrk), intent(in)::pmv
character(*), intent(in), optional:: LogFile
real(mrk), intent(out):: LFthreshold, LFpercentile
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik), parameter::Ncandidates=8
real(mrk), dimension(Ncandidates), parameter:: candidatesP=&
    (/100._mrk/365.25_mrk, 1._mrk, 2.5_mrk, 5._mrk, 10._mrk, 15._mrk, 20._mrk, 25._mrk/)
real(mrk):: thresh, candidatesQ(Ncandidates)
integer(mik)::fileUnit, i, Nlist(Ncandidates), k(1)
real(mrk), pointer:: LowD(:,:), pmvList(:)

err=0;mess='';

!misc.
fileUnit=RandomUnitNumber()
if(present(LogFile)) open(unit=fileUnit, file=trim(LogFile), status='REPLACE')

! compute candidate percentiles
do i=1,Ncandidates
    candidatesQ(i)=EmpiricalQuantile(x=x%ts(:)%q,p=0.01_mrk*candidatesP(i),mask=(x%ts(:)%q /= x%mv))
enddo

! Get Number of events for each of them
! Extract yearly drought duration
Call AnnualIndices_daily(xin=x,lowD=lowD,lowT=candidatesQ,&
                         pmv=pmvList, err=err,mess=mess)
if(err>0) then
    mess='GetLowFlowThreshold:FATAL:'//trim(mess);return
endif
! Compute number of events
do i=1,Ncandidates
    Nlist(i) = count( (pmvList<=pmv) .and. (lowD(:,i)>0) )
enddo

! Select threshold yielding N closest to Nobs
k=minloc(abs(Nlist-Nobj))
LFthreshold=candidatesQ(k(1))
LFpercentile=candidatesP(k(1))

! write 2 logfile if supplied
if(present(LogFile)) then
    write(FileUnit,*) 'N objective = ',Nobj
    do i=1,Ncandidates
        write(FileUnit,*) 'Percentile (%): ',candidatesP(i),&
                          'Value : ',candidatesQ(i), &
                          'N events: ',Nlist(i)
    enddo
    write(FileUnit,*) 'Selected percentile (%)',LFpercentile,&
                      'Value : ',LFThreshold, &
                      'Absolute diff. with Nobj: ',abs(Nlist(k(1))-Nobj)
endif

if(present(LogFile)) close(fileUnit)

end subroutine GetLowFlowThreshold

function RandomUnitNumber()

!#**********************************************************************
!#* Purpose: small private function to generate a random unit number to
!#* avoid overwriting when writing from within a subroutine
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified: 23/03/2009
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* OUT
!#*    1.RandomUnitNumber, Integer
!#**********************************************************************
integer(mik)::RandomUnitNumber
real(mrk)::harvest

call random_number ( harvest )
RandomUnitNumber=1+floor(harvest*250)

end function RandomUnitNumber

end module HydroIndices_tools

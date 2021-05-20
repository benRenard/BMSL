module Verification_tools

!~**********************************************************************
!~* Purpose: Utilities for verifying probabilistic predictions
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea/Columbia
!~**********************************************************************
!~* Last modified: 20/09/2012
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
!~*		4. 
!~*		5. 
!~*		6. 
!~*		7. 
!~*		8. 
!~*		9. 
!~*		10.
!~*		11.
!~*     12.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use Dates_tools, only:DateType
implicit none
Private
public :: Value2Occurence,GetReliabilityDiagram,GetPITDiagram,GetPredSummary

! Validation Derived types
type, public:: OneTimeValidType
    type(DateType)::date=DateType(undefIN,undefIN,undefIN,undefIN,undefIN,undefIN)
    real(mrk)::obs=undefRN
    real(mrk), allocatable::pred(:)
end type OneTimeValidType

Type,public::ValidType
    type(OneTimeValidType), allocatable:: time(:)
end Type ValidType
! Plot objects
type, public:: ReliabilityDiagramType
    character(250)::nickname='AintGotNoName'
    integer(mik)::nBin=UndefIN
    integer(mik)::nT=UndefIN
    real(mrk),allocatable::bin(:)
    real(mrk),allocatable::BinCenter(:)
    integer(mik),allocatable::BinCount(:)
    real(mrk),allocatable::BinFreq(:)
    integer(mik),allocatable::ObsCount(:)
    real(mrk),allocatable::ObsFreq(:)
    real(mrk)::CIlevel=0.95_mrk
    real(mrk),allocatable::CI(:,:)
    real(mrk),allocatable::CI_x(:)
end type ReliabilityDiagramType
!interfaces
interface Value2Occurence
  module procedure Value2Occurence_1,Value2Occurence_v
endinterface Value2Occurence
interface Value2PIT
  module procedure Value2PIT_1,Value2PIT_v
endinterface Value2PIT


contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetPredSummary(Valid,file,Doplot,main,err,mess)
!^**********************************************************************
!^* Purpose: summarizes predictive distributions
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea/Columbia
!^**********************************************************************
!^* Last modified: 23/09/2012
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.Valid
!^*		2.[file]
!^*		3.DoPlot
!^*		4.[main], main title for the plot
!^* OUT
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:GetSpareUnit
use EmpiricalStats_tools,only:GetEmpiricalStats
type(ValidType), intent(in)::Valid(:)
character(*),intent(in),optional::file,main
logical, intent(in)::doplot
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
real(mrk)::all15(15)
integer(mik)::i,station,k,unt
real(mrk), allocatable::Mat(:,:),z(:)
logical, allocatable::mask(:)

err=0;mess=''
if(present(file)) then
    call getSpareUnit(unt,err,mess)
    if(err>0) then;mess='GetPredSummary:'//trim(mess);return;endif
    open(unit=unt,file=File,status='replace')
    write(unt,'(21A16)') 'num','Station','Year','Month','Day','Obs','n','mini','maxi',&
        'range','mean','median','Q10','Q25','Q75','Q90','std','var','CV','skew','kurtosis'
    k=0
    do station=1,size(Valid)
        do i=1,size(Valid(station)%time)
            if(allocated(mask)) deallocate(mask);allocate(mask(size(Valid(station)%time(i)%pred)))
            mask=Valid(station)%time(i)%pred>=0._mrk
            if(allocated(z)) deallocate(z);allocate(z(count(mask)))
            z=pack(Valid(station)%time(i)%pred,mask)
            call GetEmpiricalStats(x=z,all15=all15,err=err,mess=mess)
            if(err>0) then;mess='GetPredSummary:'//trim(mess);return;endif
            k=k+1
            write(unt,'(5I16,16F16.3)') k,station,Valid(station)%time(i)%date%year,&
            Valid(station)%time(i)%date%month,Valid(station)%time(i)%date%day,&
            Valid(station)%time(i)%obs,all15
        enddo
    enddo
    close(unt)
endif
if(DoPlot) then
    ! Put everything in Matrix Form
    if(allocated(Mat)) deallocate(Mat);allocate(Mat(k,7))
    k=0
    do station=1,size(Valid)
        do i=1,size(Valid(station)%time)
            k=k+1
            if(allocated(mask)) deallocate(mask);allocate(mask(size(Valid(station)%time(i)%pred)))
            mask=Valid(station)%time(i)%pred>=0._mrk
            if(allocated(z)) deallocate(z);allocate(z(count(mask)))
            z=pack(Valid(station)%time(i)%pred,mask)
            call GetEmpiricalStats(x=z,all15=all15,err=err,mess=mess)
            if(err>0) then;mess='GetPredSummary:'//trim(mess);return;endif
            Mat(k,:)=(/real(station,mrk),real(Valid(station)%time(i)%date%year,mrk),&
                      Valid(station)%time(i)%obs,all15((/7,6,10,13/))/)
        enddo
    enddo
    call PlotPredSummary(Mat,main)
endif

end subroutine GetPredSummary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetReliabilityDiagram(Valid,bins,Diagram,file,Doplot,main,err,mess)

!^**********************************************************************
!^* Purpose: Compute x/ys of reliability diagram
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea/Columbia
!^**********************************************************************
!^* Last modified: 20/09/2012
!^**********************************************************************
!^* Comments: Assumes ValidIN is already transformed in an Occurence obs/prob
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.ValidIN, original obs/forecast
!^*		2.bins
!^* OUT
!^*		1.ValidOut, transformed obs/forecast
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
use Distribution_tools, only:GetQuantile
type(ValidType), intent(in)::Valid
real(mrk), intent(in)::bins(:)
character(*),intent(in),optional::file,main
logical, intent(in)::doplot
type(ReliabilityDiagramType), intent(out)::Diagram
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,b
real(mrk)::forecast(size(Valid%time))
logical::mask(size(Valid%time)),feas

err=0;mess=''

call InitReliabilityDiagram(bins,Diagram)
Diagram%nT=size(Valid%time)
do i=1,Diagram%nT
    forecast(i)=Valid%time(i)%pred(1)
enddo
do b=1,Diagram%nBin
    Diagram%BinCenter(b)=0.5_mrk*(bins(b)+bins(b+1))
    mask= ( forecast>=bins(b) .and. forecast<bins(b+1) )
    Diagram%BinCount(b)=count(mask)
    Diagram%BinFreq(b)=real(Diagram%BinCount(b),mrk)/real(Diagram%nT,mrk)
    Diagram%CI_x( ((b-1)*3+1):(b*3) )=(/bins(b),Diagram%BinCenter(b),bins(b+1)/)
    if(Diagram%BinCount(b)==0) then
        Diagram%ObsCount(b)=undefIN
        Diagram%ObsFreq(b)=undefRN
        Diagram%CI(((b-1)*3+1):(b*3),:)=undefRN
    else
        Diagram%ObsCount(b)=sum(Valid%time(:)%obs,mask=mask)
        Diagram%ObsFreq(b)=real(Diagram%ObsCount(b),mrk)/real(count(mask),mrk)
        ! binomial prediction interval
        ! @BinCenter
        call GetQuantile(DistId='Binomial',p=0.5_mrk*(1._mrk-Diagram%CIlevel),&
                par=(/Diagram%BinCenter(b),real(Diagram%BinCount(b),mrk)/),&
                q=Diagram%CI((b-1)*3+2,1),feas=feas,err=err,mess=mess)
        if(err>0) then;mess='GetReliabilityDiagram:'//trim(mess);return;endif
        if(.not.feas) then;err=1;mess='GetReliabilityDiagram:Unfeasible binomial par';return;endif
        call GetQuantile(DistId='Binomial',p=1._mrk-0.5_mrk*(1._mrk-Diagram%CIlevel),&
                par=(/Diagram%BinCenter(b),real(Diagram%BinCount(b),mrk)/),&
                q=Diagram%CI((b-1)*3+2,2),feas=feas,err=err,mess=mess)
        if(err>0) then;mess='GetReliabilityDiagram:'//trim(mess);return;endif
        if(.not.feas) then;err=1;mess='GetReliabilityDiagram:Unfeasible binomial par';return;endif
        ! @Binleft
        call GetQuantile(DistId='Binomial',p=0.5_mrk*(1._mrk-Diagram%CIlevel),&
                par=(/bins(b),real(Diagram%BinCount(b),mrk)/),&
                q=Diagram%CI((b-1)*3+1,1),feas=feas,err=err,mess=mess)
        if(err>0) then;mess='GetReliabilityDiagram:'//trim(mess);return;endif
        if(.not.feas) then;err=1;mess='GetReliabilityDiagram:Unfeasible binomial par';return;endif
        call GetQuantile(DistId='Binomial',p=1._mrk-0.5_mrk*(1._mrk-Diagram%CIlevel),&
                par=(/bins(b),real(Diagram%BinCount(b),mrk)/),&
                q=Diagram%CI((b-1)*3+1,2),feas=feas,err=err,mess=mess)
        if(err>0) then;mess='GetReliabilityDiagram:'//trim(mess);return;endif
        if(.not.feas) then;err=1;mess='GetReliabilityDiagram:Unfeasible binomial par';return;endif
        ! @BinRight
        call GetQuantile(DistId='Binomial',p=0.5_mrk*(1._mrk-Diagram%CIlevel),&
                par=(/bins(b+1),real(Diagram%BinCount(b),mrk)/),&
                q=Diagram%CI((b-1)*3+3,1),feas=feas,err=err,mess=mess)
        if(err>0) then;mess='GetReliabilityDiagram:'//trim(mess);return;endif
        if(.not.feas) then;err=1;mess='GetReliabilityDiagram:Unfeasible binomial par';return;endif
        call GetQuantile(DistId='Binomial',p=1._mrk-0.5_mrk*(1._mrk-Diagram%CIlevel),&
                par=(/bins(b+1),real(Diagram%BinCount(b),mrk)/),&
                q=Diagram%CI((b-1)*3+3,2),feas=feas,err=err,mess=mess)
        if(err>0) then;mess='GetReliabilityDiagram:'//trim(mess);return;endif
        if(.not.feas) then;err=1;mess='GetReliabilityDiagram:Unfeasible binomial par';return;endif
        Diagram%CI( ((b-1)*3+1):(b*3),:)=Diagram%CI(((b-1)*3+1):(b*3),:)/real(count(mask),mrk)
    endif
enddo

if(Doplot) then
    if(present(main)) then
        call PlotReliabilityDiagram(Diagram,main)
    else
        call PlotReliabilityDiagram(Diagram,'')
    endif
endif
if(present(file)) then
    call WriteReliabilityDiagram(Diagram,File,err,mess)
    if(err>0) then;mess='GetReliabilityDiagram:'//trim(mess);return;endif
endif
end subroutine GetReliabilityDiagram
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetPITDiagram(Valid,file,Doplot,main,DoPerSite,DoPerYear,err,mess)

!^**********************************************************************
!^* Purpose: Compute and plot PIT diagram
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea/Columbia
!^**********************************************************************
!^* Last modified: 21/09/2012
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.Valid, PIT
!^*		2.[file]
!^*		3.DoPlot
!^*		4.[main], main title for the plot
!^* OUT
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:GetSpareUnit
type(ValidType), intent(in)::Valid(:)
character(*),intent(in),optional::file,main
logical, intent(in)::doplot
logical, intent(in), optional::DoPerSite,DoPerYear
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
type(ValidType), pointer::PIT(:)
integer(mik)::station,i,k,unt
real(mrk), allocatable::Mat(:,:)

err=0;mess=''
call Value2PIT(Valid,PIT,err, mess)
if(err>0) then;mess='GetPITDiagram:'//trim(mess);return;endif
if(present(file)) then
    call getSpareUnit(unt,err,mess)
    if(err>0) then;mess='GetPITDiagram:'//trim(mess);return;endif
    open(unit=unt,file=File,status='replace')
    write(unt,'(7A16)') 'num','Station','Year','Month','Day','Obs','PIT'
    k=0
    do station=1,size(PIT)
        do i=1,size(PIT(station)%time)
            k=k+1
            write(unt,'(5I16,2F16.3)') k,station,PIT(station)%time(i)%date%year,&
            PIT(station)%time(i)%date%month,PIT(station)%time(i)%date%day,&
            PIT(station)%time(i)%obs,PIT(station)%time(i)%pred(1)
        enddo
    enddo
    close(unt)
endif
if(DoPlot) then
    ! Put everything in Matrix Form
    if(allocated(Mat)) deallocate(Mat);allocate(Mat(k,4))
    k=0
    do station=1,size(PIT)
        do i=1,size(PIT(station)%time)
            k=k+1
            Mat(k,:)=(/real(station,mrk),real(PIT(station)%time(i)%date%year,mrk),PIT(station)%time(i)%obs,PIT(station)%time(i)%pred(1)/)
        enddo
    enddo
    call PlotPIT(Mat,main,DoPerSite,DoPerYear)
endif
end subroutine GetPITDiagram
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Value2PIT_1(Valid,PIT,err, mess)
!^**********************************************************************
!^* Purpose: transform values/predictive into PIT
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea/Columbia
!^**********************************************************************
!^* Last modified: 21/09/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.ValidIN, original obs/forecast
!^* OUT
!^*		1.PIT, Vector of PIT values
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
use EmpiricalStats_tools,only:GetEmpiricalPval
type(ValidType), intent(in)::Valid
real(mrk),intent(out)::PIT(size(Valid%time))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::i

do i=1,size(Valid%time)
    call GetEmpiricalPval(obs=Valid%time(i)%obs,x=Valid%time(i)%pred,pv=PIT(i),err=err,mess=mess)
    if(err>0) then;mess='Value2PIT_1:'//trim(mess);return;endif
enddo

end subroutine Value2PIT_1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Value2PIT_v(Valid,PIT,err, mess)
!^**********************************************************************
!^* Purpose: transform values/predictive into PIT
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea/Columbia
!^**********************************************************************
!^* Last modified: 21/09/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.Valid, original obs/forecast
!^* OUT
!^*		1.PIT, a ValidType object where %pred contains the PIT value
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
use EmpiricalStats_tools,only:GetEmpiricalPval
type(ValidType), intent(in)::Valid(:)
type(ValidType), pointer::PIT(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::i,station

if(associated(PIT)) deallocate(PIT);allocate(PIT(size(Valid)))
do station=1,size(Valid)
    if(allocated(PIT(station)%time)) deallocate(PIT(station)%time);allocate(PIT(station)%time(size(Valid(station)%time)))
    do i=1,size(Valid(station)%time)
        if(allocated(PIT(station)%time(i)%pred)) deallocate(PIT(station)%time(i)%pred);allocate(PIT(station)%time(i)%pred(1))
        call GetEmpiricalPval(obs=Valid(station)%time(i)%obs,x=Valid(station)%time(i)%pred,pv=PIT(station)%time(i)%pred(1),err=err,mess=mess)
        if(err>0) then;mess='Value2PIT_v:'//trim(mess);return;endif
        PIT(station)%time(i)%obs=Valid(station)%time(i)%obs
        PIT(station)%time(i)%date=Valid(station)%time(i)%date
    enddo
enddo
end subroutine Value2PIT_v
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine Value2Occurence_1(OccurenceType,par,ValidIN,ValidOut, err, mess)

!^**********************************************************************
!^* Purpose: transform values/predictive into binary occurences/forecast prob
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea/Columbia
!^**********************************************************************
!^* Last modified: 20/09/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.OccurenceType, e.g. '>=' for exceddance of some value
!^*		2.par, e.g. the value to be exceeded
!^*		3.ValidIN, original obs/forecast
!^* OUT
!^*		1.ValidOut, transformed obs/forecast
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
character(*), intent(in)::OccurenceType
real(mrk), intent(in)::par
type(ValidType), intent(in)::ValidIN
type(ValidType), intent(out)::ValidOut
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::nT,i

! Allocation
nT=size(ValidIn%time)
if(allocated(ValidOut%time)) deallocate(ValidOut%time);allocate(ValidOut%time(nT))
do i=1,nT
    if(allocated(ValidOut%time(i)%pred)) deallocate(ValidOut%time(i)%pred)
    allocate(ValidOut%time(i)%pred(1))
enddo
! Transform to occurence
do i=1,nT
    select case(trim(OccurenceType))
    case('>=')
        ValidOut%time(i)%obs=real(count((/ValidIN%time(i)%obs>=par/)),mrk)
        ValidOut%time(i)%pred=real(count(ValidIN%time(i)%pred>=par),mrk)/real(size(ValidIN%time(i)%pred),mrk)
    case default
        err=1;mess='Value2Occurence:Fatal:Unavailable OccurenceType'
    end select
enddo


end subroutine Value2Occurence_1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine Value2Occurence_v(OccurenceType,par,ValidIN,ValidOut, err, mess)

!^**********************************************************************
!^* Purpose: vectorial overload of Value2Occurence (eg multi-site)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea/Columbia
!^**********************************************************************
!^* Last modified: 20/09/2012
!^**********************************************************************
!^* Comments: merge results from all sites together
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.OccurenceType, e.g. '>=' for exceddance of some value
!^*		2.par, e.g. the value to be exceeded
!^*		3.ValidIN, original obs/forecast
!^* OUT
!^*		1.ValidOut, transformed obs/forecast
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************
character(*), intent(in)::OccurenceType
real(mrk), intent(in)::par(:)
type(ValidType), intent(in)::ValidIN(:)
type(ValidType), intent(out)::ValidOut
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::nS,i,nT,nTot,j,k
type(ValidType)::DummyVOut

! Allocation
nS=size(ValidIN)
nTot=0
do i=1,nS
    nT=size(ValidIn(i)%time)
    nTot=nTot+nT
enddo
if(allocated(ValidOut%time)) deallocate(ValidOut%time);allocate(ValidOut%time(nTot))
do i=1,nTot
    if(allocated(ValidOut%time(i)%pred)) deallocate(ValidOut%time(i)%pred)
    allocate(ValidOut%time(i)%pred(1))
enddo

k=0
do i=1,nS
    call Value2Occurence(OccurenceType,par(i),ValidIN(i),DummyVOut, err, mess)
    if(err>0) then; mess='Value2Occurence_v:'//trim(mess);return;endif
    do j=1,size(DummyVOut%time)
        k=k+1
        ValidOut%time(k)%obs=DummyVOut%time(j)%obs
        ValidOut%time(k)%pred=DummyVOut%time(j)%pred
    enddo
enddo

end subroutine Value2Occurence_v
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pure subroutine InitReliabilityDiagram(bins,Diagram)
real(mrk),intent(in)::bins(:)
type(ReliabilityDiagramType), intent(out)::Diagram

Diagram%nBin=size(bins)-1
if(allocated(Diagram%BinCenter)) deallocate(Diagram%BinCenter);allocate(Diagram%BinCenter(Diagram%nBin))
if(allocated(Diagram%BinCount)) deallocate(Diagram%BinCount);allocate(Diagram%BinCount(Diagram%nBin))
if(allocated(Diagram%BinFreq)) deallocate(Diagram%BinFreq);allocate(Diagram%BinFreq(Diagram%nBin))
if(allocated(Diagram%ObsCount)) deallocate(Diagram%ObsCount);allocate(Diagram%ObsCount(Diagram%nBin))
if(allocated(Diagram%ObsFreq)) deallocate(Diagram%ObsFreq);allocate(Diagram%ObsFreq(Diagram%nBin))
if(allocated(Diagram%CI)) deallocate(Diagram%CI);allocate(Diagram%CI(3*Diagram%nBin,2))
if(allocated(Diagram%CI_x)) deallocate(Diagram%CI_x);allocate(Diagram%CI_x(3*Diagram%nBin))
if(allocated(Diagram%bin)) deallocate(Diagram%bin);allocate(Diagram%bin(Diagram%nBin+1))
Diagram%bin=bins
end subroutine InitReliabilityDiagram
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PlotReliabilityDiagram(Diagram,main)
use RFortran
type(ReliabilityDiagramType), intent(in)::Diagram
character(*),intent(in)::main
!locals
logical::ok,mask(size(Diagram%BinCount))

mask=Diagram%BinCount>0
ok=Rput('mask',mask)
ok=Rput('bin',Diagram%Bin)
ok=Rput('x',Diagram%BinCenter)
ok=Rput('y',Diagram%ObsFreq)
ok=Rput('BinFreq',Diagram%BinFreq)
ok=Rput('CI',Diagram%CI)
ok=Rput('CI_x',Diagram%CI_x)
ok=Rput('main',trim(main))
ok=Reval('CImask=(CI[,1]>=0)')
ok=Reval('X11()')
ok=Reval('plot(x[mask],y[mask],xlim=c(0,1),ylim=c(0,1),col="white",xlab="Predicted probability",ylab="observed frequency",main=main)')
ok=Reval('polygon(c(CI_x[CImask],CI_x[CImask][sum(CImask):1]), c(CI[CImask,1],CI[CImask,2][sum(CImask):1]),border=NA,col="papayawhip")')
ok=Reval('lines(c(0,1),c(0,1))')
!ok=Reval('lines(x[mask],y[mask],type="b",lwd=2,col="red")')
ok=Reval('for(i in 1:length(bin)){lines(c(bin[i],bin[i]),c(0,1),type="l",lty=2,lwd=0.5,col="lightgrey")}')
ok=Reval('for(i in 1:(length(bin)-1)){lines(c(bin[i],bin[i+1]),c(y[i],y[i]),col="red",lwd=4)}')
ok=Reval('lines(x,BinFreq,type="b",lwd=0.5,col="steelblue")')
ok=Reval('legend("topleft",col=c("red","steelblue","papayawhip"),lty=1,lwd=c(5,1,20),legend=c("Event frequency","Bin frequency ","95% interval"))')
end subroutine PlotReliabilityDiagram
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WriteReliabilityDiagram(Diagram,File,err,mess)
use utilities_dmsl_kit,only:getSpareUnit
type(ReliabilityDiagramType), intent(in)::Diagram
character(*),intent(in)::File
integer(mik),intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::unt,i

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err>0) then;mess='WriteReliabilityDiagram:'//trim(mess);return;endif
open(unit=unt,file=File,status='replace')
write(unt,'(14A16)') 'num','BinCenter','BinLeft','BinRight','BinCount','BinFreq','EventCount',&
            'EventFreq','loCI@BinLeft','loCI@BinCenter','loCI@BinRight','hiCI@BinLeft',&
            'hiCI@BinCenter','hiCI@BinRight'
do i=1,Diagram%nBin
    write(unt,'(I16,3F16.3,I16,F16.3,I16,7F16.3)') i,Diagram%BinCenter(i),Diagram%bin(i),Diagram%bin(i+1),&
            Diagram%BinCount(i),Diagram%BinFreq(i),Diagram%ObsCount(i),Diagram%ObsFreq(i),&
            Diagram%CI( (i-1)*3+1 ,1),Diagram%CI( (i-1)*3+2 ,1),Diagram%CI( (i-1)*3+3 ,1),&
            Diagram%CI( (i-1)*3+1 ,2),Diagram%CI( (i-1)*3+2 ,2),Diagram%CI( (i-1)*3+3 ,2)
enddo
close(unt)

end subroutine WriteReliabilityDiagram
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PlotPIT(Mat,main,DoPerSite,DoPerYear)
use RFortran
real(mrk),intent(in)::Mat(:,:)
character(*),intent(in),optional::main
logical, intent(in), optional::DoPerSite,DoPerYear
! locals
logical::ok,DPS,DPY
integer(mik)::n,i

! handle optionals
if(present(DoPerSite)) then;DPS=DoPerSite;else;DPS=.true.;endif
if(present(DoPerYear)) then;DPY=DoPerYear;else;DPY=.true.;endif
! Send to R
n=size(Mat,dim=1)
ok=Rput('Mat',Mat)
ok=Rput('x',(/( (real(i,mrk)-0.5_mrk)/real(n,mrk),i=1,n )/))
ok=Rput('y',Mat(:,4))
ok=Reval('mask=Mat[,3]>0')
if(present(main)) then
    ok=Rput('main',trim(main))
else
    ok=Rput('main','')
endif
! Plot all PIT values
ok=Reval('X11()')
ok=Reval('plot(x[mask],sort(y[mask]),main=main,type="l",lwd=3,col="red",xlim=c(0,1),ylim=c(0,1),xlab="frequency",ylab="PIT")')
if(DPS) then
    ok=Reval('for(i in 1:max(Mat[,1])){mask=(Mat[,3]>0)&(Mat[,1]==i);k=sum(mask);if(k>1){lines((seq(1,k)-0.5)/k,sort(y[mask]),lwd=1,col="steelblue")}}')
endif
if(DPY) then
    ok=Reval('for(i in min(Mat[,2]):max(Mat[,2])){mask=(Mat[,3]>0)&(Mat[,2]==i);k=sum(mask);if(k>1){lines((seq(1,k)-0.5)/k,sort(y[mask]),lwd=1,col="lightgray")}}')
endif
ok=Reval('lines(c(0,1),c(0,1),lwd=1,col="black")')
ok=Reval('mask=Mat[,3]>0')
ok=Reval('lines(x[mask],sort(y[mask]),lwd=3,col="red")')
ok=Reval('legend("topleft",col=c("red","steelblue","lightgray"),lty=1,lwd=c(3,1,1),legend=c("PIT","PIT per station","PIT per year"))')
end subroutine PlotPIT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PlotPredSummary(Mat,main)
use RFortran
real(mrk),intent(in)::Mat(:,:)
character(*),intent(in),optional::main
! locals
logical::ok

ok=Rput('Mat',Mat)
if(present(main)) then
    ok=Rput('main',trim(main))
else
    ok=Rput('main','')
endif
ok=Reval('for(i in 1:max(Mat[,1])){if((i-1)%%9==0){X11();par(mfrow=c(3,3))};mask=((Mat[,1]==i)&(Mat[,3]>=0));plot(Mat[mask,2],Mat[mask,3],main=i,xlab="year",ylab="value");lines(Mat[mask,2],Mat[mask,4],col="green");lines(Mat[mask,2],Mat[mask,5],col="red");lines(Mat[mask,2],Mat[mask,6],col="blue");}')
!ok=Reval('lines(Mat[mask,2],Mat[mask,4],col="green");')
!ok=Reval('lines(Mat[mask,2],Mat[mask,5],col="red");')
!ok=Reval('lines(Mat[mask,2],Mat[mask,6],col="steelblue")}')

end subroutine PlotPredSummary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module Verification_tools

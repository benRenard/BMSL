program TestTSAggregation

use kinds_dmsl_kit
use Dates_tools
use TimeSeries_tools
use numerix_dmsl_kit,only:normaldev
use Aggregation_tools

implicit none
type(DateType), parameter::d1=DateType(1960,1,1,0,0,0),d2=DateType(2010,1,1,0,0,0)
character(250), parameter::funk='mean'
character(1), parameter::baseScale='h',targetScale='D'
integer(mik), parameter::scaleMult=1,option=666
real(mrk), parameter::mu=0._mrk,sigma=1._mrk
type(OneDTimeSeriesType)::TSin,TSout
type(DateType), pointer::list(:)
integer(mik)::err,i,n
character(250)::message,fDate
real(mrk)::dummy,tstart,tend
real(mrk), pointer::tx(:),tgrid(:)
Type(dateType)::FirstDate

! generate TS data
call DateList(d1,d2,baseScale,1,list,err,message)
n=size(list)
allocate(TSin%ts(n))
TSin%NameSeries='Input';TSin%n=n
TSin%ts(:)%date=list
do i=1,n
    call normaldev(mean=mu,sdev=sigma,gdev=TSin%ts(i)%q,err=err,message=message)
enddo
dummy=DaysFromRef(d2,d1)
! Perform aggregation
call CPU_TIME(tstart)

!! engine version
!call TSAggregate(operateur=basic_op,opID=(/LINEAR_interpol,AVE_funk,INFECT_qc,LIBERAL_cc/),&
!                TSin=TSin,scale=targetScale,scaleMult=scaleMult,DateListOption=option,&
!                TSout=TSout,tx=tx,tgrid=tgrid,err=err,message=message)

!! EZ version
!call TSAggregate(operateur=funk,&
!                TSin=TSin,scale=targetScale,scaleMult=scaleMult,DateListOption=option,&
!                TSout=TSout,tx=tx,tgrid=tgrid,err=err,message=message)

! nearly-minimal version
FirstDate=TSin%ts(1)%date
FirstDate%hour=6;FirstDate%minute=0;FirstDate%second=0
call TSAggregate(operateur=funk,&
                TSin=TSin,scale=targetScale,FirstDate=FirstDate ,&
                TSout=TSout,tx=tx,tgrid=tgrid,err=err,message=message)

!! minimal version
!call TSAggregate(operateur=funk,&
!                TSin=TSin,scale=targetScale,&
!                TSout=TSout,err=err,message=message)

call CPU_TIME(tend)
write(*,*) 'time:', tend-tstart
open(unit=666,file='Input.txt',status='replace')
do i=1,n
    call FormatDate(TSin%ts(i)%date,101,fDate,err,message)
    write(666,'(I14,A14,e24.14,e14.6)') i,trim(fDate),tx(i),TSin%ts(i)%q
enddo
close(666)

open(unit=666,file='Output.txt',status='replace')
do i=1,TSout%n
    call FormatDate(TSout%ts(i)%date,101,fDate,err,message)
    write(666,'(I14,A14,e24.14,e14.6)') i,trim(fDate),tgrid(i),TSout%ts(i)%q
enddo
close(666)

write(*,*)'Finito!'
read(*,*)

end program TestTSAggregation
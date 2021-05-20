program TestDates
use kinds_dmsl_kit
use dates_tools
implicit none

integer(mik), parameter::n=10000
type(DateType)::d0,d,dlist(n)
integer(mik)::i,err
real(mrk)::told(n),debut,fin
character(250)::mess,fDate
type(DateType), pointer::list(:)
real(mrk),allocatable::t(:)

d0=dateType(1950,1,31,12,30,10)
d=dateType(2010,3,5,12,0,0)
dlist=HourlyList(d0,InNHours(d0,n-1))
!do i=1,36; write(*,*) InNMonths(d0,-1*i,1); enddo
!write(*,*) DaysFromRef(d0,d)

call CPU_TIME(debut)
call DateList(d0,d,'m',15,list,err,mess)
call CPU_TIME(fin)
write(*,*) 'time:', fin-debut
write(*,*) '------------------------'
do i=1,min(size(list),20)
    write(*,*) list(i)
enddo
write(*,*) '...'
do i=(size(list)-3),size(list)
    write(*,*) list(i)
enddo
write(*,*) '------------------------'

call CPU_TIME(debut)
allocate(t(size(list)))
t=AbsoluteTime(list)
call CPU_TIME(fin)
write(*,*) 'time:', fin-debut
write(*,*) '------------------------'
do i=1,min(size(t),20)
    write(*,*) t(i)
enddo
write(*,*) '...'
do i=(size(t)-3),size(t)
    write(*,*) t(i)
enddo
write(*,*) '------------------------'

!call CPU_TIME(debut)
!do i=1,n
!    t(i)=DaysBetween(dlist(i),dlist(1))
!enddo
!call CPU_TIME(fin)
!write(*,*) 'time:', fin-debut
!do i=1,25;write(*,*) t(i);enddo
!
!call CPU_TIME(debut)
!do i=1,n
!    t(i)=DaysFromRef(dlist(i),dlist(1))
!enddo
!call CPU_TIME(fin)
!write(*,*) 'time:', fin-debut
!do i=1,25;write(*,*) t(i);enddo
!
!call CPU_TIME(debut)
!t=DaysFromRef(dlist,dlist(1))
!call CPU_TIME(fin)
!write(*,*) 'time:', fin-debut
!do i=1,25;write(*,*) t(i);enddo
!
!call CPU_TIME(debut)
!t(1)=0._mrk
!do i=2,n
!    t(i)=t(i-1)+DaysFromRef(dlist(i),dlist(i-1))
!enddo
!call CPU_TIME(fin)
!write(*,*) 'time:', fin-debut
!do i=1,25;write(*,*) t(i);enddo

write(*,*)'Finito!'
read(*,*)

end program TestDates
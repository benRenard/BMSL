program Test_RegionalTest

use Distribution_tools
use Test_tools
use Rfortran

implicit none

integer(mik),parameter::n=50,p=10,nsim=1000,mvOption=1
real(mrk), parameter::ro=0.8_mrk,mvprob=0._mrk,mv=-9999._mrk,level=0.95_mrk,trend=0.5_mrk
!locals
integer(mik)::i,j,k,err
real(mrk)::y(n),x(n,p),dev
logical::feas,ok
character(250)::mess
type(TestResultType)::out(nsim)

ok=Rinit()
do i=1,nsim
    ! Generate data
    call GenerateSample(DistId='Gaussian',par=(/0._mrk,1._mrk/),gen=y,feas=feas,err=err,mess=mess)
    do j=1,n
        call GenerateSample(DistId=AR1,par=(/0._mrk,1._mrk,ro/),gen=x(j,:),feas=feas,err=err,mess=mess)
        x(j,:)=x(j,:)+trend*y(j)
    enddo
    ! Add missing values
    do j=1,n
        do k=1,p
            call Generate(DistId=UNIF,par=(/0._mrk,1._mrk/),gen=dev,feas=feas,err=err,mess=mess)
            if(dev<=mvprob) x(j,k)=mv
        enddo
    enddo
    ! apply test
    call RegionalTest(y,x,level,mv,MVoption,out(i),err,mess)
    if(err>0) then;write(*,*) i,trim(mess);read(*,*);endif
enddo

write(*,*) 'rejact rate:', 1._mrk-real(count(out(:)%H0),mrk)/real(nsim,mrk)
! plot data
ok=Rput('y',y);ok=Rput('x',x);ok=Rput('mv',mv)
ok=Reval('matplot(y,x)')
ok=Rput('pval',out(:)%pval)
ok=Reval('X11()')
ok=Reval('hist(pval[pval>=0])')

read(*,*)

end program Test_RegionalTest
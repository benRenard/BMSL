program source1

use kinds_dmsl_kit
use Distribution_tools
implicit none

integer(mik)::err
character(100)::mess

character(50)::DistID,name
integer(mik)::npar,n_gen,i
logical::loga,feas,isnull
real(mrk),allocatable::par(:)
real(mrk)::val
real(mrk),allocatable::gen(:)

err=0;mess=''

!Bernoulli,GEOM,Poisson,Binomial,NegBinomial
            
!call GetParNumber(DistID, npar=npar, err=err, mess=mess)
!call GetParName(DistID, name, err, mess)

! Test Bernoulli Distribution
write (*,*) "Bernoulli Test"
write (*,*) "Pdf"
if(allocated(par)) deallocate(par)
allocate(par(1))
par=0.8_mrk
call GetPdf(DistId=Bernoulli,x=1._mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
if (val==par(1)) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

call GetPdf(DistId=Bernoulli,x=0.99_mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
if (val==par(1)) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

call GetPdf(DistId=Bernoulli,x=0._mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
if (val==1-par(1)) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

call GetPdf(DistId=Bernoulli,x=0.05_mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
if (val==1-par(1)) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

write (*,*) "Cdf"
call GetCdf(DistId=Bernoulli,x=-1._mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
if (val==0) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

call GetCdf(DistId=Bernoulli,x=0._mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
if (val==1-par(1)) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

call GetCdf(DistId=Bernoulli,x=0.9_mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
if (val==1-par(1)) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

call GetCdf(DistId=Bernoulli,x=1._mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
if (val==1) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

write (*,*) "Generate"
n_gen=2000
if(allocated(gen)) deallocate(gen)
allocate(gen(n_gen))

do i=1,n_gen
    call Generate(DistId=Bernoulli,par=par,gen=gen(i),feas=feas,err=err,mess=mess)
enddo
write(*,*) sum(gen)/n_gen,par(1)

! Test Geometric distribution
write (*,*) "GEOM Test"
write (*,*) "Pdf"
if(allocated(par)) deallocate(par)
allocate(par(1))
par=0.8_mrk
call GetPdf(DistId=GEOM,x=3._mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
if (val==(1-par(1))**2_mik*par(1)) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

call GetPdf(DistId=GEOM,x=4._mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
if (val==(1-par(1))**3_mik*par(1)) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif


write (*,*) "Cdf"
call GetCdf(DistId=GEOM,x=3._mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
if (val==1._mrk-(1-par(1))**3_mik) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif

call GetCdf(DistId=GEOM,x=4._mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
if (val==1._mrk-(1-par(1))**4_mik) then
    write (*,*) "OK"
else
    write (*,*) "error"
endif


write (*,*) "Generate"
n_gen=2000
if(allocated(gen)) deallocate(gen)
allocate(gen(n_gen))

do i=1,n_gen
    call Generate(DistId=GEOM,par=par,gen=gen(i),feas=feas,err=err,mess=mess)
enddo
write(*,*) sum(gen)/n_gen,1._mrk/par(1)

! Test Poisson Distribution
write (*,*) "Poisson Test"
write (*,*) "Pdf"
if(allocated(par)) deallocate(par)
allocate(par(1))
par=2._mrk
call GetPdf(DistId=Poisson,x=5._mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
write (*,*) val, 0.03608941

par=5._mrk
call GetPdf(DistId=Poisson,x=10._mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
write (*,*) val, 0.01813279


write (*,*) "Cdf"
par=2._mrk
call GetCdf(DistId=Poisson,x=5._mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
write (*,*) val, 0.9834364


par=5._mrk
call GetCdf(DistId=Poisson,x=10._mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
write (*,*) val, 0.9863047

write (*,*) "Generate"
n_gen=2000
if(allocated(gen)) deallocate(gen)
allocate(gen(n_gen))

do i=1,n_gen
    call Generate(DistId=Poisson,par=par,gen=gen(i),feas=feas,err=err,mess=mess)
enddo
write(*,*) sum(gen)/n_gen,par(1)

! Test Binomial Distribution
write (*,*) "Binomial Test"
write (*,*) "Pdf"
if(allocated(par)) deallocate(par)
allocate(par(2))
par(1)=0.8_mrk
par(2)=10._mrk
call GetPdf(DistId=Binomial,x=5._mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
write (*,*) val, 0.02642412


write (*,*) "Cdf"
call GetCdf(DistId=Binomial,x=5._mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
write (*,*) val, 0.0327935


write (*,*) "Generate"
n_gen=2000
if(allocated(gen)) deallocate(gen)
allocate(gen(n_gen))

do i=1,n_gen
    call Generate(DistId=Binomial,par=par,gen=gen(i),feas=feas,err=err,mess=mess)
enddo
write(*,*) sum(gen)/n_gen,par(1)*par(2)

! Test NegBinomial Distribution
write (*,*) "NegBinomial Test"
write (*,*) "Pdf"
if(allocated(par)) deallocate(par)
allocate(par(2))
par(1)=0.8_mrk
par(2)=10._mrk
call GetPdf(DistId=NegBinomial,x=5._mrk,par=par,loga=.false.,pdf=val,feas=feas,isnull=isnull,err=err,mess=mess)
write (*,*) val, 0.0687882


write (*,*) "Cdf"
call GetCdf(DistId=NegBinomial,x=5._mrk,par=par,cdf=val,feas=feas,err=err,mess=mess)
write (*,*) val, 0.9389486


write (*,*) "Generate"
n_gen=2000
if(allocated(gen)) deallocate(gen)
allocate(gen(n_gen))

do i=1,n_gen
    call Generate(DistId=NegBinomial,par=par,gen=gen(i),feas=feas,err=err,mess=mess)
enddo
write(*,*) sum(gen)/n_gen,par(1)*par(2)/(1._mrk-par(1))
read (*,*)

end program source1


program TestAggregation
use kinds_dmsl_kit
use aggregation_tools
use utilities_dmsl_kit,only:twopi
implicit none

integer(mik),parameter::n=8766,n_int=8766/5 !366
real(mrk)::x(n),y(n),xout,yout(1),xout_full(n_int),yout_full(n_int,1),intervals(n_int,2)
integer(mik)::qc(n),cc(n),qcout(1),ccout(1),qcout_full(n_int,1),ccout_full(n_int,1)
integer(mik)::i,err
character(250)::message
real(mrk)::tstart,tend

! define x,y,qc, etc.
x=(/(i*1._mrk,i=1,n)/)
y=cos(x*twopi/30._mrk) !y=2._mrk*x
qc=QC_PRESENT; qc(100:150)=QC_MISSVAL !qc(50)=QC_MISSVAL
cc=CC_CONT; cc(120:200)=CC_DISCONT !cc(3)=CC_DISCONT
intervals(:,1)=(/(x(1)+real(i-1,mrk)*(x(n)-x(1))/real(n_int,mrk),i=1,n_int)/)
intervals(:,2)=(/(x(1)+real(i,mrk)*(x(n)-x(1))/real(n_int,mrk),i=1,n_int)/)
! Apply operator on a single interval
call basic_op(opID=(/LINEAR_interpol,AVE_funk,INFECT_qc,CONSERV_cc/),&
              y=y,x=x,xmin=5.5_mrk,xmax=16.50_mrk,&
              qc=qc,&
              cc=cc,&
              !aux,& !IN
              yOut=yout,xOut=xout,&
              qcOut=qcout,&
              ccOut=ccOut,&
              err=err,message=message)
write(*,*)trim(message)

call CPU_TIME(tstart)
call aggregator(operateur=basic_op,&
                opID=(/LINEAR_interpol,SUM_funk,INFECT_qc,LIBERAL_cc/),&
                y=y,x=x,&
                grid=(/intervals(:,1),intervals(n_int,2)/),&
                !intervals=intervals,&
                qc=qc,cc=cc,&
                !aux,& !IN
                yOut=yout_full,xOut=xout_full,qcOut=qcout_full,ccOut=ccout_full,&
                err=err,message=message)
write(*,*)trim(message)
call CPU_TIME(tend)
write(*,*)'time:',tend-tstart

! write to file
if(.true.)then
  open(666,file='IN.txt',status='unknown')
  do i=1,n
    write(666,'(2E14.6,2I14)') x(i),y(i),qc(i),cc(i)
  enddo
  close(666)
  open(666,file='OUT.txt',status='unknown')
  do i=1,n_int
    write(666,'(4E14.6,2I14)') intervals(i,:),xout_full(i),yout_full(i,1),qcout_full(i,1),ccout_full(i,1)
  enddo
  close(666)
endif

write(*,*)'Finito!'; !read(*,*)

endprogram TestAggregation

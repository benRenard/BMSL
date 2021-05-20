program main
use kinds_dmsl_kit
use MCMCStrategy_tools
use distribution_tools,only:Generate, getPdf,GAUSS
use InOutTools
use omp_lib

implicit none

integer(mik)::N_x,i,j,N_par,N_chain
real(mrk),pointer::x(:),in_x(:,:),fx2(:),GR(:,:)
real(mrk)::par(2),fx,fAux(1,4),fAux2(1)
character(100)::distId

real(mrk)::MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor,scale
real(mrk),allocatable::start(:),startSTD(:),startM(:,:),startStdM(:,:),covarM(:,:,:),scaleM(:),covar(:,:)
integer(mik)::nAdapt,nCycles	

integer(mik)::showpar(3)
real(mrk)::ranges(3,2)

real(mrk),pointer::outMat1(:,:),outMat2(:,:),outMat3b(:,:),outMat3(:,:,:)

character(250)::mess
integer(mik)::err
logical::feas,isnull


err=0;mess='';feas=.true.;isnull=.false.
BurnFactor=0.5_mrk
DistId=GAUSS
N_par=2
N_x=100_mik
par=(/0._mrk,10._mrk/)

MinMoveRate=0.1_mrk;MaxMoveRate=0.5_mrk;DownMult=0.9_mrk;UpMult=1.1_mrk
nAdapt=200;nCycles=100

if(associated(x)) deallocate(x)
allocate(x(N_x))


do i=1,N_x
    call Generate(DistId=trim(DistID),par=par,gen=x(i),feas=feas,err=err,mess=mess)
    If(err>0) then
        mess='MeanProg error: '//trim(mess)
        write (*,*) mess
        read (*,*)
    endif
enddo

call mat_output(Filename="testdate.txt", vect=x, err=err, mess=mess)
If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif

call mat_input(Filename="testdate.txt", column=1, header=.false., tab=in_x, err=err, mess=mess)
If(err>0) then 
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif

x=in_x(:,1)

if(allocated(start)) deallocate(start)
allocate(start(N_par))
if(allocated(startStd)) deallocate(startStd)
allocate(startStd(N_par))


start(1)=20._mrk
start(2)=50._mrk

startStd(1)=2._mrk
startStd(2)=5._mrk

!call Adaptive_Metro_OAAT(f=XXX,x=start,&
!				fx=fx,std=startStd,&
!				nAdapt=nAdapt,nCycles=nCycles,&
!				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
!				DownMult=DownMult,UpMult=UpMult,&
!				OutFile="C:\Mydoc\XunProg\TestSubroutine\TestMCMC\testMCMC.txt",&
!				UseRFortran=.true.,&
!				err=err,mess=mess)
 !   If(err>0) then
!	    mess='Main:'//trim(mess)
!	    read(*,*)
 !   endif

N_chain=4

if(allocated(startM)) deallocate(startM)
allocate(startM(N_par,N_chain))
if(allocated(startStdM)) deallocate(startStdM)
allocate(startStdM(N_par,N_chain))
if(allocated(covarM)) deallocate(covarM)
allocate(covarM(N_par,N_par,N_chain))
if(allocated(scaleM)) deallocate(scaleM)
allocate(scaleM(n_chain))
if(allocated(covar)) deallocate(covar)
allocate(covar(N_par,N_par))
covar=0
covar(1,1)=1._mrk
covar(2,2)=1._mrk

startM(1,1)=20._mrk
startM(1,2)=40._mrk
startM(1,3)=-20._mrk
startM(1,4)=-40._mrk
startM(2,1)=50._mrk
startM(2,2)=60._mrk
startM(2,3)=70._mrk
startM(2,4)=80._mrk

startstdM(1,1)=2._mrk
startstdM(1,2)=4._mrk
startstdM(1,3)=2._mrk
startstdM(1,4)=4._mrk
startstdM(2,1)=5._mrk
startstdM(2,2)=6._mrk
startstdM(2,3)=7._mrk
startstdM(2,4)=8._mrk

startstdM=0.001*startstdM
startM(1,:)=20._mrk
startM(2,:)=40._mrk


covarM=0._mrk
do i=1,n_par
    covarM(i,i,:)=1._mrk
enddo
scaleM=2.4_mrk

showpar(1)=1
showpar(2)=2
showpar(3)=1

ranges(1,1)=-10._mrk
ranges(1,2)=50._mrk

ranges(2,1)=0._mrk
ranges(2,2)=100._mrk

ranges(3,1)=-10._mrk
ranges(3,2)=50._mrk

scaleM=0.001*scaleM
scale=2.4_mrk
scale=scale*0.001
startstd=startstd*0.001

!call MetroMS_Mix2(f=XXX,x=start,fx=fx2,fAux=fAux,std=startstd,&
!                nAdapt1=nAdapt,nCycles1=nCycles,BurnFactor=BurnFactor,&
!                nAdapt2=nAdapt,nCycles2=nCycles,&
!                nSim=50000,&
!!               scaleGJ,GRBlocSize,&
!				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
!				DownMult=DownMult,UpMult=UpMult,&
!				n_chain=n_chain,&
!!				UseRFortran=.false.,&
!				showPar=showpar,&
!!				ranges,&
!				GRIndex=GR, &
!!				GRBurnfactor,&
!!				OutDoc="",&
!!				Filename1="testMCMC_multi1.txt",&
!!				Filename2="testMCMC_multi2.txt",&
!!				Filename3="testMCMC_multi3.txt",&
!				outMat1=outMat1,&
!				outMat2=outMat2,&
!				outMat3=outMat3,&
!				err=err,mess=mess)
				


!call mat_output(Filename="GR.txt", mat=GR, err=err, mess=mess)


!call AdaptiveMS_Metro(f=XXX,x=startM,fx=fx2,covar=covarM,scale=scaleM,&
!                nAdapt=nAdapt,nCycles=nCycles,&
!				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
!				DownMult=DownMult,UpMult=UpMult,&
!				showPar=showPar,&
!				!ranges=ranges,&
!				OutDoc="",&
!				Filename="testMCMC_multi.txt",&
!				UseRFortran=.true.,&
!				err=err,mess=mess)

!call AdaptiveMS_Metro_OAAT(f=XXX,x=startM,&
!                fx=fx2,std=startStdM,&
!				nAdapt=nAdapt,nCycles=nCycles,&
!				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
!				DownMult=DownMult,UpMult=UpMult,&
!				showPar=showPar,&
!				!ranges=ranges,&
!				OutDoc="",&
!				Filename="testMCMC_multi.txt",&
!				UseRFortran=.true.,&
!				err=err,mess=mess)
!				
!call MetroMS_Mix(f=XXX,x=startM,&
!                fx=fx2,std=startStdM,&
!                nAdapt1=nAdapt,nCycles1=nCycles,BurnFactor=BurnFactor,&
!                nAdapt2=nAdapt,nCycles2=nCycles,&
!                nSim=40000,&
!!                fAux=fAux, &
!!                scaleGJ,&
!				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
!				DownMult=DownMult,UpMult=UpMult,&
!!				UseRFortran,&
!				showPar=showPar,&
!				!ranges=ranges,&
!				OutDoc="",&
!				Filename1="testMCMC_multi1.txt",&
!!				Filename2="testMCMC_multi2.txt",&
!				Filename3="testMCMC_multi3.txt",&
!				err=err,mess=mess)
!
!call Metro_Mix(f=XXX,x=start,fx=fx,fAux=fAux2,std=startstd,&
!                nAdapt1=nAdapt,nCycles1=nCycles,BurnFactor=BurnFactor,&
!                nAdapt2=nAdapt,nCycles2=nCycles,&
!                nSim=5000,&
!!                scaleGJ,&
!				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
!				DownMult=DownMult,UpMult=UpMult,&
!!				UseRFortran,&
!                showPar=showPar,&
!				!ranges=ranges,&
!				OutFile1="testMCMC1.txt",&
!				OutFile2="testMCMC2.txt",&
!				OutFile3="testMCMC3.txt",&
!				outMat1=outMat1,&
!				outMat2=outMat2,&
!				outMat3=outMat3b,&
!				err=err,mess=mess)
!
!call Adaptive_Metro(f=XXX,x=start,fx=fx,fAux=fAux2,&
!                covar=covar,scale=scale,&
!                nAdapt=nAdapt,nCycles=nCycles,&
!                MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
!				DownMult=DownMult,UpMult=UpMult,&
!!				UseRFortran,&
!				showPar=showpar,&
!!				ranges,&
!				OutFile="testAdaM.txt",&
!				outMat=outMat1,&
!				err=err,mess=mess)
!				
!call Adaptive_Metro_OAAT(f=XXX,x=start,fx=fx,&
!                std=startstd,&
!                nAdapt=nAdapt,nCycles=nCycles,&
!                MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
!				DownMult=DownMult,UpMult=UpMult,&
!!				UseRFortran,&
!				showPar=showpar,&
!!				ranges,&
!				OutFile="testMetro.txt",&
!				err=err,mess=mess)
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine XXX(theta,feas,isnull,fx,fAux,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:  /  /2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* OUT
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* INOUT
!^*     1.
!^*     2.
!^**********************************************************************
use kinds_dmsl_kit
real(mrk),intent(in)::theta(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
!local
real(mrk)::temp
integer(mik)::i

fx=0._mrk
do i=1,size(x)
call GetPdf(DistId=DistId,x=x(i),par=theta,loga=.true.,pdf=temp,feas=feas,isnull=isnull,err=err,mess=mess)
if (err>0) then
    mess='XXXerror'//trim(mess)
    return
endif
fx=temp+fx
enddo

if (present(fAux)) fAux(1)=fx+10._mrk

end subroutine XXX




end program main

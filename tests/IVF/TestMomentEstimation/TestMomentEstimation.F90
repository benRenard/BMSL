Program TestMomentEstimation
use Distribution_tools
use kinds_dmsl_kit
use MomentEstimation_tools
use MLEstimation_tools
implicit none

character(250), parameter::Dist=GEV,Mtype='LMOM'
integer(mik), parameter::n=200,Nboot=100,IterMax=10000
real(mrk), dimension(2), parameter::teta0_Gaussian=(/10.,1./),&
                                    teta0_Gumbel=(/10.,1./)
real(mrk), dimension(3), parameter::teta0_GEV=(/10.,1.,-0.2/),&
                                    teta0_PearsonIII=(/10.,1.,4./)    
integer(mik)::npar,err,i
character(250)::mess
real(mrk), allocatable::teta0(:),par(:),BootPar(:,:),Hess(:,:),StdError(:),MLstart(:)
real(mrk)::X(n)
logical::feas

call GetParNumber(Dist, npar, err, mess)
if(err>0) then; write(*,*) mess;STOP;endif
allocate(teta0(npar),par(npar),BootPar(Nboot,npar),Hess(npar,npar),StdError(npar),MLstart(npar))
select case(Dist)
case(Gaussian)
    teta0=teta0_Gaussian
case(Gumbel)
    teta0=teta0_Gumbel
case(GEV)
    teta0=teta0_GEV
case(PearsonIII)
    teta0=teta0_PearsonIII
case default
    write(*,*) 'Unavailable Dist';STOP
end select
! Generate data
call GenerateSample(DistID=Dist,par=teta0,gen=X,feas=feas,err=err,mess=mess)
if(err>0) then; write(*,*) mess;STOP;endif
if(.not.feas) then; write(*,*) 'Unfeasible teta0';STOP;endif
! Estimate
call GetMomentEstimate(Dist,X,Mtype,par,err,mess)
if(err>0) then; write(*,*) mess;STOP;endif
write(*,*) npar,teta0
write(*,*) par
!STOP
! Bootstrap
call GetUncertainty_Bootstrap(DistID=Dist,X=X,Mtype=Mtype,&
                    Nsim=Nboot,par=BootPar,err=err,mess=mess)
if(err>0) then; write(*,*) mess;STOP;endif
do i=1,Nboot
    write(*,'(<npar>F8.2)') BootPar(i,:)
enddo

! Try ML estimation
MLstart=par
Call ML_EstimPar(X=X,distID=Dist,mv=-99._mrk,&
                        teta0=MLstart,IterMax=IterMax,&
                        teta=par,Hessian=Hess,err=err,mess=mess)
write(*,*) trim(mess)
if(err>0) STOP
write(*,*) '---------------------------------------'
write(*,'(<npar>F8.2)') par
do i=1,npar
    StdError(i)=sqrt(Hess(i,i))
enddo
write(*,'(<npar>F8.2)') StdError
do i=1,Npar
    write(*,'(<npar>F8.2)') Hess(i,:)
enddo
call GetUncertainty_Fisher(MLestimate=par,Fisher=Hess,SkipUnfeas=.true.,&
                 distID=Dist,Nsim=Nboot,par=BootPar,err=err,mess=mess)
if(err>0) then; write(*,*) mess;STOP;endif
write(*,*) '---------------------------------------'
do i=1,Nboot
    write(*,'(<npar>F8.2)') BootPar(i,:)
enddo
end Program TestMomentEstimation
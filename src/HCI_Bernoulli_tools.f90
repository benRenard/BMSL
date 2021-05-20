module HCI_Bernoulli_tools

!~**********************************************************************
!~* Purpose: Hidden Climate Indices with Bernoulli 0/1 data
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified:23/01/2017
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~* 1. Consider replacing sep-files by tabulated files (easier to read  
!~*    for very wide data tables) 
!~* 2. Implement formula for predicting  
!~* 3. Improve and modularize type "SpaceTimeDataType"  
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. 
!~*		2. 
!~*		3. 
!~*		4. 
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: PerformHCI

! Site-specific properties needed to implement the speed-up strategy
type:: SpeedUpPropType
	integer(mik), allocatable::indx(:) ! indices defining the restriction
	real(mrk), allocatable::sigmainv(:,:) ! inverse of the restricted covariance matrix
	real(mrk)::sigmadet=undefRN ! determinant of the restricted covariance matrix
end type SpeedUpPropType

! Super-speed-up strategy: evaluate full gaussian pdf only on a subset of stations
type:: SuperSpeedUpType
	integer(mik)::Nslim=1 ! slimming factor (keep one site every Nslim)
	integer(mik)::N=undefIN ! Size after slimming
	integer(mik), allocatable::indx(:) ! indices defining the restriction
	real(mrk), allocatable::sigmainv(:,:) ! inverse of the restricted covariance matrix
	real(mrk)::sigmadet=undefRN ! determinant of the restricted covariance matrix
end type SuperSpeedUpType

! Speed Up object, used to pass properties of the speed-up strategy
type:: SpeedUpType
    integer(mik)::option=undefIN ! 0= no speed-up, 1=fixed-N strategy, 2=fixed-correl strategy
    integer(mik)::N=undefIN ! N in strategy 1 (fixed-N strategy)
    real(mrk)::corr=undefRN ! correlation in strategy 2 (fixed-correl strategy)
    type(SpeedUpPropType),allocatable::prop(:) ! Site-specific properties; size Nx.
    type(SuperSpeedUpType)::super
end type SpeedUpType

! MCMC object, used to store MCMC simulations
type:: MCMCType
    integer(mik)::Nsim=undefIN,Nx=undefIN,Nt=undefIN,Nbeta=undefIN,Ngamma=undefIN ! size of the case study
	real(mrk), allocatable::sam(:,:) ! samples, Nsim*(Nx+Nt+Nbeta+Ngamma)
	real(mrk), allocatable::post(:) ! Nsim
end type MCMCType

! MCMC tunings
type:: MCMCTunings
    integer(mik)::Nsim=undefIN,Nadapt=undefIN,StopAdapt=undefIN,Nburn=undefIN,Nslim=undefIN
    logical::DoBlockSampling=.true.
end type MCMCTunings

! Hypermodel object, used to store properties of lambda's hyperdistributions
type:: HyperModelType
	character(250)::Zmodel="undefined",Dmodel="undefined"
	integer(mik)::Nbeta=undefIN,Ngamma=undefIN ! number of beta/gamma parameters
	real(mrk), allocatable::Z(:,:) ! Nx*Nz
	character(250), allocatable::Znames(:) ! covariates names
	character(250)::Zfile=''
end type HyperModelType

! Y Data object, storing Y 0/1 values
type:: SpaceTimeDataType
	character(250)::file=''
	integer(mik)::Ns=undefIN,Nt=undefIN !dimensions
	integer(mik),allocatable::val(:,:) !data
	integer(mik)::mv=-99 !value denoting missing data
	character(250),allocatable::rownames(:) ! row names (stations)
	character(250),allocatable::colnames(:) ! column names (unused here)
	real(mrk),allocatable::sgrid(:,:) ! space grid (station/gridpoints coordinates)
	real(mrk),allocatable::tgrid(:) ! time grid
	real(mrk),allocatable::D(:,:) ! distance matrix
end type SpaceTimeDataType

! HCI object, storing all inference properties
type:: HCIType
	character(250)::workspace=''
	character(250)::Config_Model='',Config_Y='',Config_RunOptions='',Config_MCMC='',Config_W='' ! config files
    integer(mik)::K=undefIN ! number of components
    type(SpaceTimeDataType)::Y,W ! Y (predictand) and W (predictor) data
    type(HyperModelType),allocatable::hypermodel(:) ! hyper-model for lambdas
    type(MCMCtunings)::mcmc
    type(SpeedUpType)::speedUp
    logical::runOptions(3)
end type HCIType

! Formatting
character(1),parameter::sep=";"
character(250),parameter::lambda_file="Results_lambda.txt",tau_file="Results_tau.txt",&
                          mcmc_file_prefix="Results_MCMC_",mcmc_file_extension=".txt",&
                          forcings_file_prefix="Results_PHI_",forcings_file_extension=".txt",&
                          forcings_file_mean="Results_PHI_mean.txt",forcings_file_sdev="Results_PHI_sdev.txt",&
                          forcings_file_meanSE="Results_PHI_MSE.txt",forcings_file_design="Results_PHI_Design.txt",&
                          forcings_file_signif="Results_PHI_Signif.txt",&
                          SE_file_prefix="Results_PHI_SE_",SE_file_extension=".txt"

! MCMC hard-coded tunings
real(mrk), parameter:: minMR=0.1,maxMR=0.5,DownMult=0.9,UpMult=1.1,&
                       SDratio=0.5,SDzero=0.1
integer(mik),parameter::Tau0_option=1 ! 1 = smart, 2 = random

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PerformHCI(file,err,mess)
!^**********************************************************************
!^* Purpose: Perform HCI analysis
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:08/02/2017
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file, master config file
!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
character(*),intent(in)::file
integer(mik),intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='PerformHCI'
type(HCIType)::HCI
type(MCMCtype),allocatable::mcmc(:)

err=0;mess=''

! General information, always loaded
call LoadHCI(trim(file),HCI,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
! Step 1: Bayesian-MCMC estimation of the hidden climate indices
if(HCI%RunOptions(1)) then
    write(*,*) "--------------------------------------------------------------"
    write(*,*) "Step 1: Bayesian-MCMC estimation of the hidden climate indices"
    write(*,*) "--------------------------------------------------------------"
    if(allocated(mcmc)) deallocate(mcmc)
    allocate(mcmc(HCI%K+1))
    call StepwiseInference(HCI)
endif
! Step 2: Identification of the climate forcings
if(HCI%RunOptions(2)) then
    write(*,*) "----------------------------------------------"
    write(*,*) "Step 2: Identification of the climate forcings"
    write(*,*) "----------------------------------------------"
    call Config_W_Read(trim(HCI%workspace)//trim(HCI%Config_W),HCI,err,mess) 
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    call IdentifyForcings(HCI)
endif
! Step 3: prediction
if(HCI%RunOptions(3)) then
endif
 
end subroutine PerformHCI
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoadHCI(file,HCI,err,mess)
!^**********************************************************************
!^* Purpose: Load HCI object containing all inference info
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:08/02/2017
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file, master config file
!^* OUT
!^*		1.HCI, HCI object
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
character(*),intent(in)::file
type(HCIType),intent(out)::HCI
integer(mik),intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='LoadHCI'
integer(mik)::i

err=0;mess=''

call Config_Read(trim(file),HCI,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
call Config_RunOptions_Read(trim(HCI%workspace)//trim(HCI%Config_RunOptions),HCI,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
call Config_Model_Read(trim(HCI%workspace)//trim(HCI%Config_Model),HCI,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
call Config_Y_Read(trim(HCI%workspace)//trim(HCI%Config_Y),HCI,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif   
do i=1,HCI%K+1
    if( size(HCI%hypermodel(i)%Z,dim=1)/=HCI%Y%Ns ) then
        err=1;mess=trim(procname)//':Fatal:size mismatch [Z,Y]';return
    endif
enddo
call Config_MCMC_Read(trim(HCI%workspace)//trim(HCI%Config_MCMC),HCI,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

end subroutine LoadHCI
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine StepwiseInference(HCI)
!^**********************************************************************
!^* Purpose: Stepwise inference of K components
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:24/01/2017
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.HCI object
!^**********************************************************************
use numerix_dmsl_kit,only:normaldev !,uniran
use EmpiricalStats_tools,only:GetEmpiricalStats
type(HCIType),intent(in)::HCI
!locals
type(MCMCtype)::mcmc(HCI%K+1)
real(mrk),parameter::epsilon=0.000001_mrk
integer(mik)::i,j,n1,ntot,err,t,s
real(mrk)::lambda(HCI%Y%Ns,HCI%K+1),tau(HCI%Y%Nt,HCI%K),freq,u,expect,mean,std,ct
real(mrk),allocatable::beta(:),gamma(:)
type(MCMCtype)::kmcmc
logical::feas
character(250)::mess

! preliminaries
lambda=0._mrk;tau=0._mrk

!!!!!!!!!!!!!!!!!!!
! Component 0
! define starting point
write(*,*) 'Component:',0,'/',HCI%K
do j=1,HCI%Y%Ns
    n1=count(HCI%Y%val(j,:)==1)
    ntot=count(HCI%Y%val(j,:)>=0)
    if(n1==0) then
        freq=epsilon
    else if(n1==ntot) then
        freq=1._mrk-epsilon
    else
        freq=real(n1,mrk)/real(ntot,mrk)
    endif
    lambda(j,1)=log(freq/(1._mrk-freq))
enddo
if(allocated(beta)) deallocate(beta)
if(allocated(gamma)) deallocate(gamma)
allocate(beta(HCI%hypermodel(1)%Nbeta),gamma(HCI%hypermodel(1)%Ngamma))
beta=0._mrk;gamma=0.1_mrk
! do mcmc
call DoMCMC(HCI%Y%val,HCI%Y%D,HCI%hypermodel(1),lambda(:,1:1),&
            tau,beta,gamma,&
            HCI%mcmc%Nsim,HCI%mcmc%Nadapt,HCI%mcmc%StopAdapt,&
            HCI%mcmc%Nburn,HCI%mcmc%Nslim,&
            HCI%mcmc%DoBlockSampling,HCI%speedUp,kmcmc)
call CopyMCMC(kmcmc,mcmc(1))
call WriteSummary(trim(HCI%workspace),HCI%Y%val,mcmc(1:1))
!!!!!!!!!!!!!!!!!!!
! Component 1 to K
do i=1,HCI%K
    write(*,*) 'Component:',i,'/',HCI%K
    ! define starting point
    call GetMaxpost(mcmc(i)%sam(:,1:HCI%Y%Ns),mcmc(i)%post,lambda(:,i))
    lambda(:,i+1)=0._mrk
    if(i>1) call GetMaxpost(mcmc(i)%sam(:,(HCI%Y%Ns+1):(HCI%Y%Ns+HCI%Y%Nt)),&
                            mcmc(i)%post,tau(:,i-1))    
    ! Get starting point for tau
    if(Tau0_option==1) then ! "smart" Approach 
        do t=1,HCI%Y%Nt            
            ! Compute expected number of events
            expect=0._mrk
            do s=1,HCI%Y%Ns
                if(HCI%Y%val(s,t)<0) cycle
                ct=lambda(s,1)
                if(i>1) ct=ct+dot_product(lambda(s,2:i),tau(t,:))
                expect=expect+1._mrk/(1._mrk+exp(-1._mrk*ct))
            enddo
            ! actual number of events
            n1=count(HCI%Y%val(:,t)==1)
            ! difference
            tau(t,i)=n1-expect
        enddo
        ! center and scale tau
        call GetEmpiricalStats(x=tau(:,i),mean=mean,std=std,err=err,mess=mess)
        tau(:,i)=(tau(:,i)-mean)/std
        call applyConstraints(tau(:,i),feas)
        if(.not.feas) write(*,*) 'Infeasible starting point for tau, "random" approach is used'
    endif
    if(Tau0_option==2 .or. .not.feas) then ! "random" approach
        feas=.false.
        do while(.not.feas)
            do j=1,HCI%Y%Nt
                call normaldev(mean=0._mrk,sdev=1._mrk,gdev=tau(j,i),err=err,message=mess)
            enddo
            call applyConstraints(tau(:,i),feas)
            if(.not.feas) then
                write(*,*) 'Trying another starting point for tau...'
            else
                write(*,*) 'Starting point OK'
            endif            
        enddo
    endif
    ! beta and gamma
    if(allocated(beta)) deallocate(beta)
    if(allocated(gamma)) deallocate(gamma)
    allocate(beta(HCI%hypermodel(i+1)%Nbeta),gamma(HCI%hypermodel(i+1)%Ngamma))
    beta=0._mrk;gamma=0.1_mrk
    ! do mcmc
    call DoMCMC(HCI%Y%val,HCI%Y%D,HCI%hypermodel(i+1),&
                lambda(:,1:(i+1)),tau(:,1:i),beta,gamma,&
                HCI%mcmc%Nsim,HCI%mcmc%Nadapt,HCI%mcmc%StopAdapt,&
                HCI%mcmc%Nburn,HCI%mcmc%Nslim,&
                HCI%mcmc%DoBlockSampling,HCI%speedUp,kmcmc)
    call CopyMCMC(kmcmc,mcmc(i+1))
    call WriteSummary(trim(HCI%workspace),HCI%Y%val,mcmc(1:(i+1)))
enddo
 
end subroutine StepwiseInference
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine IdentifyForcings(HCI)
!^**********************************************************************
!^* Purpose: Identify climate forcings
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:25/02/2017
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.HCI object
!^**********************************************************************
use utilities_dmsl_kit, only:number_string
use DataRW_tools,only:ReadSeparatedFile
use linalg_dmsl_kit,only:choles_invrt
use EmpiricalSTats_tools,only:GetEmpiricalStats
use numerix_dmsl_kit,only:tdist_icdf
type(HCIType),intent(in)::HCI
!locals
real(mrk)::T(HCI%Y%Nt,HCI%K+1)
real(mrk),allocatable::bigT(:,:,:),Forcings(:,:,:),SE(:,:,:),M(:,:),saveTT(:,:),ttest(:,:,:)
integer(mik)::i,j,err,ncol,nsim
character(250)::mess,file
real(mrk),pointer::values(:,:)
real(mrk)::TT(HCI%K+1,HCI%K+1),sigmadet,mean(HCI%W%Ns,HCI%K+1),sdev(HCI%W%Ns,HCI%K+1),&
           res(HCI%W%Ns,HCI%Y%Nt),resStd(HCI%W%Ns),meanSE(HCI%W%Ns,HCI%K+1),meanTtest(HCI%W%Ns,HCI%K+1),tlim
logical::feas

! preliminaries
tlim=tdist_icdf(0.95_mrk,real(HCI%Y%Nt-2,mrk)) ! rejection limit for t-test


! Read MCMC files and populate bigT
do i=1,HCI%K+1
    ncol=HCI%Y%Ns+HCI%Y%Nt+HCI%hypermodel(i)%nbeta+HCI%hypermodel(i)%ngamma+1
    file=trim(HCI%workspace)//trim(mcmc_file_prefix)//trim(number_string(i-1))//trim(mcmc_file_extension)
    call ReadSeparatedFile(file=file,&
                       sep=sep,nhead=0,ncol=ncol,y=values,&
                       err=err,mess=mess)
    nsim=size(values,dim=1)
    if(i==1) then
        allocate(bigT(nsim,HCI%Y%Nt,HCI%K+1))
        allocate(Forcings(nsim,HCI%W%Ns,HCI%K+1))
        allocate(SE(nsim,HCI%W%Ns,HCI%K+1))
        allocate(saveTT(nsim*(HCI%K+1),HCI%K+1))
        allocate(ttest(nsim,HCI%W%Ns,HCI%K+1));ttest=0._mrk
    endif
    bigT(:,:,i)=values(:,(HCI%Y%Ns+1):(HCI%Y%Ns+HCI%Y%Nt))
enddo

! Apply identification formula for all mcmc samples
do i=1,nsim
    T=bigT(i,:,:)
    TT=matmul(transpose(T),T)
    call choles_invrt(ainv=TT,logDet=sigmadet,posDefinite=feas,err=err,message=mess)
    saveTT( ((i-1)*(HCI%K+1)+1) : (i*(HCI%K+1)) ,:)=TT ! for covariance matrix of estimated forcings
    if(.not.feas) then
        Forcings(i,:,:)=HCI%W%mv
        SE(i,:,:)=HCI%W%mv
        ttest(i,:,:)=HCI%W%mv
    else
        Forcings(i,:,:)=matmul(HCI%W%val,matmul(T,TT))
        ! computation of standard errors
        ! first get residuals
        res=matmul(Forcings(i,:,:),transpose(T))-HCI%W%val
        ! compute residual stdev
        do j=1,HCI%W%Ns
            call GetEmpiricalStats(x=res(j,:),std=resStd(j),err=err,mess=mess)
        enddo
        ! get standard errors
        do j=1,HCI%K+1
            SE(i,:,j)=sqrt(TT(j,j))*resStd
        enddo
        ! Perform t-test
        do j=1,HCI%K+1
            where(abs(Forcings(i,:,j)/SE(i,:,j))>tlim) ttest(i,:,j)=1._mrk
        enddo        
    endif
    write(*,*) 'Done:',100*real(i)/real(Nsim),'%'
enddo
 
! Write results to file
do i=1,HCI%K+1
    if(allocated(M)) deallocate(M)
    allocate(M(Nsim,HCI%W%Ns))
    ! estimates
    M=Forcings(:,:,i)
    call WriteMatrix(trim(HCI%workspace)//trim(forcings_file_prefix)//trim(number_string(i-1))//trim(forcings_file_extension),M,"\t")
    do j=1,HCI%W%Ns
        call GetEmpiricalStats(x=pack(M(:,j),M(:,j)/=HCI%W%mv),mean=mean(j,i),std=sdev(j,i),err=err,mess=mess)
    enddo
    ! standard errors
    M=SE(:,:,i)
    call WriteMatrix(trim(HCI%workspace)//trim(SE_file_prefix)//trim(number_string(i-1))//trim(SE_file_extension),M,"\t")
    do j=1,HCI%W%Ns
        call GetEmpiricalStats(x=pack(M(:,j),M(:,j)/=HCI%W%mv),mean=meanSE(j,i),err=err,mess=mess)
    enddo
    ! t-test
    M=ttest(:,:,i)
    do j=1,HCI%W%Ns
        call GetEmpiricalStats(x=pack(M(:,j),M(:,j)/=HCI%W%mv),mean=meanTtest(j,i),err=err,mess=mess)
    enddo
enddo
! Forcings: mean and stdev of point-estimate accross simulations
call WriteMatrix(trim(HCI%workspace)//trim(forcings_file_mean),mean,"\t")
call WriteMatrix(trim(HCI%workspace)//trim(forcings_file_sdev),sdev,"\t")
! SE: average of standard errors accross simulations
call WriteMatrix(trim(HCI%workspace)//trim(forcings_file_meanSE),meanSE,"\t")
! T-test: frequency of significant results accross simulations
call WriteMatrix(trim(HCI%workspace)//trim(forcings_file_signif),meanTtest,"\t")
! Save disign matrices accross simulations
call WriteMatrix(trim(HCI%workspace)//trim(forcings_file_design),saveTT,"\t")
end subroutine IdentifyForcings
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WriteSummary(folder,y,mcmc)
!^**********************************************************************
!^* Purpose: Write summary files of stepwise inference
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:24/01/2017
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.folder, where files are written
!^*		2.y, data matrix (Nt*Nx)
!^*		3.mcmc (K+1)
!^* OUT
!^**********************************************************************
use utilities_dmsl_kit, only:number_string
character(*),intent(in)::folder
type(MCMCtype), intent(in)::mcmc(:)
integer(mik),intent(in)::y(:,:)
!locals
integer(mik)::Nt,Nx,K,Nsim,i,j,p
real(mrk)::lambda(size(y,dim=1),size(mcmc)),&
           tau(size(y,dim=2),size(mcmc)-1)
real(mrk),allocatable::M(:,:)

! preliminaries
Nt=size(y,dim=2);Nx=size(y,dim=1)
K=size(mcmc)-1;Nsim=size(mcmc(1)%post)

! assemble maxpost estimate
call GetMaxpost(mcmc(1)%sam(:,1:Nx),mcmc(1)%post,lambda(:,1)) ! lambda0
do i=1,K
    call GetMaxpost(mcmc(i+1)%sam(:,1:Nx),mcmc(i+1)%post,lambda(:,i+1)) ! lambda
    call GetMaxpost(mcmc(i+1)%sam(:,(Nx+1):(Nx+Nt)),mcmc(i+1)%post,tau(:,i)) ! tau
enddo
! write result files
call WriteMatrix(trim(folder)//trim(lambda_file),lambda,sep)
if(allocated(M)) deallocate(M)
allocate(M(Nt,K+1))
M(:,1)=1._mrk
M(:,2:(K+1))=tau
call WriteMatrix(trim(folder)//trim(tau_file),M,sep)

! write mcmc files
do j=1,K+1
    p=size(mcmc(j)%sam,dim=2)
    if(allocated(M)) deallocate(M)
    allocate(M(Nsim,p+1))
    M(:,1:p)=mcmc(j)%sam
    if(j==1) then
        M(:,(Nx+1):(Nx+Nt))=1._mrk
    endif
    M(:,p+1)=mcmc(j)%post
    call WriteMatrix(trim(folder)//trim(mcmc_file_prefix)//trim(number_string(j-1))//trim(mcmc_file_extension),M,sep)
enddo

deallocate(M)

end subroutine WriteSummary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetLogLik(y,lambda,tau,mu,sigmainv,sigmadet,isDiag,t,x,LL,speedUp)
!^**********************************************************************
!^* Purpose: compute full LL
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:23/01/2017
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.y, data matrix (Nt*Nx)
!^*		2.lambda, lambda matrix (Nx*K+1)
!^*		3.tau, tau matrix (Nt*K)
!^*		4.mu, lambda-hypermean (Nx)
!^*		5.sigmainv, inverse of lambda-hypervariance (Nx*Nx)
!^*		6.sigmadet, determinant of lambda-hypervariance (1)
!^*		7.isDiag, is hypervariance diagonal?
!^*		8.t: LL restricted to time step t (t==0 => all time steps)
!^*		9.x: LL restricted to site x (x==0 => all sites)
!^*		10.[speedUp], optional "speed up" object
!^* OUT
!^*		1.LL, log-likelihood
!^**********************************************************************
use numerix_dmsl_kit,only:normal_logp
real(mrk), intent(in)::lambda(:,:),tau(:,:),mu(:),sigmainv(:,:),sigmadet
integer(mik),intent(in)::y(:,:),t,x
logical,intent(in)::isDiag
type(speedUpType),intent(in),optional::speedUp
real(mrk), intent(out)::LL
!locals
integer(mik)::Nx,Nt,K,x1,x2,t1,t2,i,j
real(mrk)::ct,theta0,theta

LL=0._mrk
! preliminaries
Nt=size(y,dim=2);Nx=size(y,dim=1);K=size(lambda,dim=2)-1
if(x<=0) then;x1=1;x2=Nx;else;x1=x;x2=x;endif
if(t<=0) then;t1=1;t2=Nt;else;t1=t;t2=t;endif

if(x>=0 .and. t>=0) then !otherwise only hyperdist is needed
    ! add up lkh contributions
    do j=t1,t2
        do i=x1,x2
            if(y(i,j)<0) cycle ! missing value
            if(K>0) then;ct=dot_product(lambda(i,2:),tau(j,:));else;ct=0._mrk;endif
            !theta0=lambda(j,1)*(1._mrk+ct)
            theta0=lambda(i,1)+ct
            theta=1._mrk/(1._mrk+exp(-1._mrk*theta0))
            if(y(i,j)==1) then
                LL=LL+log(theta)
            else
                LL=LL+log(1._mrk-theta)
            endif
        enddo
    enddo
endif

! Hyper-likelihood
! Note: only the contribution of the last component is computed
! LL is therefore unnormalized, OK only as long as stepwise inference is performed!
if(t==0 .and. x==0) then ! full logpost is required
    if(speedUp%super%Nslim>1) then ! restricted evaluation of multivariate Gaussian
        LL=LL+mnorm_logp(lambda(speedUp%super%indx,K+1),mu(speedUp%super%indx),&
                    speedUp%super%sigmainv,speedUp%super%sigmadet)
    else ! no speed-up: compute full multivariate Gaussian
        LL=LL+mnorm_logp(lambda(:,K+1),mu,sigmainv,sigmadet)
    endif
else if(t<=0) then ! otherwise a tau is being proposed, and hyperdist will cancel out in metropolis ratio
    if(isDiag) then ! no matrix computation, just a product
        do j=x1,x2
            LL=LL+normal_logp(x=lambda(j,K+1), mean=mu(j), var=1._mrk/sigmainv(j,j))
        enddo
    else if(present(speedUp)) then ! speed-up is possible 
        if(x>0) then ! computation for one lambda, use speedup
            if(speedUp%option>0) then ! restricted evaluation of multivariate Gaussian
                LL=LL+mnorm_logp(lambda(speedUp%prop(x)%indx,K+1),mu(speedUp%prop(x)%indx),&
                            speedUp%prop(x)%sigmainv,speedUp%prop(x)%sigmadet)
            else ! no speed-up: compute full multivariate Gaussian
                LL=LL+mnorm_logp(lambda(:,K+1),mu,sigmainv,sigmadet)
            endif
        else ! computation for beta or gamma, use super-speed-up
            if(speedUp%super%Nslim>1) then ! restricted evaluation of multivariate Gaussian
                LL=LL+mnorm_logp(lambda(speedUp%super%indx,K+1),mu(speedUp%super%indx),&
                            speedUp%super%sigmainv,speedUp%super%sigmadet)
            else ! no speed-up: compute full multivariate Gaussian
                LL=LL+mnorm_logp(lambda(:,K+1),mu,sigmainv,sigmadet)
            endif
        endif       
    else ! no speed-up: compute full multivariate Gaussian
        LL=LL+mnorm_logp(lambda(:,K+1),mu,sigmainv,sigmadet)
    endif
endif

end subroutine GetLogLik
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetLogPost(y,lambda,tau,beta,gamma,mu,sigmainv,sigmadet,isDiag,t,x,logp,speedUp)
!^**********************************************************************
!^* Purpose: compute unnormalized log-posterior
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:25/01/2017
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.y, data matrix (Nt*Nx)
!^*		2.lambda, lambda matrix (Nx*K+1)
!^*		3.tau, tau matrix (Nt*K)
!^*		4.beta, parameters for hypermean (Nbeta)
!^*		5.gamma, parameters for hypervariance (Ngamma)
!^*		6.mu, lambda-hypermean (Nx)
!^*		7.sigmainv, inverse of lambda-hypervariance (Nx*Nx)
!^*		8.sigmadet, determinant of lambda-hypervariance (1)
!^*		9.isDiag, is hypervariance diagonal?
!^*		10.t: LL restricted to time step t (t<=0 => all time steps)
!^*		11.x: LL restricted to site x (x<=0 => all sites)
!^*		12.[speedUp], optional "speed up" object
!^* OUT
!^*		1.logp, unnormalized log-posterior
!^**********************************************************************
use numerix_dmsl_kit,only:normal_logp,invChi2_logp
real(mrk), intent(in)::lambda(:,:),tau(:,:),beta(:),gamma(:),&
                       mu(:),sigmadet,sigmainv(:,:)
integer(mik),intent(in)::y(:,:),t,x
logical,intent(in)::isDiag
type(speedUpType),intent(in),optional::speedUp
real(mrk), intent(out)::logp
!locals
integer(mik)::Nx,Nt,K,x1,x2,t1,t2,i,j


! preliminaries
Nt=size(y,dim=2);Nx=size(y,dim=1);K=size(lambda,dim=2)-1
if(x<=0) then;x1=1;x2=Nx;else;x1=x;x2=x;endif
if(t<=0) then;t1=1;t2=Nt;else;t1=t;t2=t;endif
! Compute log-likelihood
call GetLogLik(y,lambda,tau,mu,sigmainv,sigmadet,isDiag,t,x,logp,speedUp)
! !!!!!!! BOTCH !!!!!!!
!add up priors
logp=logp+normal_logp(x=beta(1), mean=0._mrk, var=(2.5_mrk)**2)
do i=1,size(gamma)
    logp=logp+invChi2_logp(x=gamma(i),v=1._mrk)
enddo

!do i=1,K+1
!    do j=x1,x2
!        logp=logp+normal_logp(x=lambda(j,i), mean=0._mrk, var=(2.5_mrk)**2)
!    enddo
!enddo

end subroutine GetLogPost

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine DoMCMC(y,D,hypermodel,lambda,tau,beta,gamma,&
                  Nsim,Nadapt,StopAdapt,Nburn,Nslim,&
                  DoBlockSampling,speedUp,&
                  mcmc)
!^**********************************************************************
!^* Purpose: MCMC for the kth component, conditionally on all previous components
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:26/01/2017
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.y, data matrix (Nt*Nx)
!^*		2.D, distance matrix (Nx*Nx)
!^*		3.hypermodel
!^*		4.lambda, lambda matrix (Nx*K+1); last column interpreted as MCMC starting point
!^*		5.tau, tau matrix (Nt*K); last column interpreted as MCMC starting point
!^*		6.beta, hypermean parameters
!^*		7.gamma, hypervariance parameters
!^*		8.Nsim, number of MCMC simulations
!^*		9.Nadapt, adaptation period
!^*		10.StopAdapt, iteration of last adaptation
!^*		11.Nburn, number of iterations to burn
!^*		12.Nslim, slim factor (only save one iteration every Nslim)
!^*		13.DoBlockSampling, sample all lambdas in a single block?
!^*		14.SpeedUp, properties of speed-up strategy
!^* OUT
!^*		1.mcmc, MCMC samples
!^**********************************************************************
use numerix_dmsl_kit,only:normaldev,uniran
integer(mik),intent(in)::y(:,:),Nsim,Nadapt,StopAdapt,Nburn,Nslim
real(mrk), intent(in)::D(:,:),lambda(:,:),tau(:,:),beta(:),gamma(:)
type(HyperModelType),intent(in)::hypermodel
logical,intent(in)::DoBlockSampling
type(SPeedUpType),intent(in)::speedup
type(MCMCtype), intent(out)::mcmc
!locals
character(250)::mess
integer(mik)::K,Nt,Nx,iter,i,j,err,compt
logical::feas,isDiag
type(SPeedUpType)::beepbeep
real(mrk)::pcurrent,pcandid,dev,ratio,u,&
           lsd(size(lambda,dim=1)),tsd(size(tau,dim=1)),&
           bsd(size(beta)),gsd(size(gamma)),&
           lmr(size(lambda,dim=1)),tmr(size(tau,dim=1)),&
           bmr(size(beta)),gmr(size(gamma)),&
           lcurrent(size(lambda,dim=1),size(lambda,dim=2)),&
           lcandid(size(lambda,dim=1),size(lambda,dim=2)),&
           tcurrent(size(tau,dim=1),size(tau,dim=2)),&
           tcandid(size(tau,dim=1),size(tau,dim=2)),&
           bcurrent(size(beta)),bcandid(size(beta)),&           
           gcurrent(size(gamma)),gcandid(size(gamma)),&  
           mucurrent(size(lambda,dim=1)),mucandid(size(lambda,dim=1)),&     
           invcurrent(size(lambda,dim=1),size(lambda,dim=1)),&
           invcandid(size(lambda,dim=1),size(lambda,dim=1)),&
           detcurrent,detcandid,&
           choles(size(lambda,dim=1),size(lambda,dim=1)),&
           cov(size(lambda,dim=1),size(lambda,dim=1))
real(mrk)::saveDet
real(mrk),allocatable::saveSigmaInv(:,:),saveMcmc(:)

! preliminaries
Nt=size(y,dim=2);Nx=size(y,dim=1);K=size(lambda,dim=2)-1
lcurrent=lambda;tcurrent=tau;bcurrent=beta;gcurrent=gamma;
lsd=SDratio*abs(lambda(:,K+1));where(lsd==0._mrk) lsd=SDzero
bsd=SDratio*abs(beta);where(bsd==0._mrk) bsd=SDzero
gsd=SDratio*abs(gamma);where(gsd==0._mrk) gsd=SDzero
if(K>0) then
    tsd=SDratio*abs(tau(:,K));where(tsd==0._mrk) tsd=SDzero
endif
lmr=0._mrk;tmr=0._mrk;bmr=0._mrk;gmr=0._mrk
isDiag=(hypermodel%Dmodel=='inde' .or. hypermodel%Dmodel=='Inde')
! Compute mu,sigmainv and sigmadet
call GetHypermean(hypermodel,bcurrent,mucurrent)
if(isDiag .or. speedup%super%Nslim==1) then
    call GetHypervariance(D,hypermodel,gcurrent,invcurrent,detcurrent,feas)
endif
! Choleskize covariance matrix if block sampling is activated
if(DoBlockSampling) call GetCholes(D,hypermodel,gcurrent,choles,feas)
! prepare MCMC object
call InitMCMC(Nsim,Nburn,Nslim,Nx,Nt,hypermodel%Nbeta,hypermodel%Ngamma,mcmc)
allocate(saveMcmc(Nx+Nt+size(beta)+size(gamma)))
compt=0
! prepare speed-up object
beepbeep%option=speedup%option
beepbeep%N=speedup%N
beepbeep%corr=speedup%corr
beepbeep%super%Nslim=speedup%super%Nslim
call SetSpeedup(beepbeep,D,hypermodel,gcurrent)
allocate(saveSigmaInv(beepbeep%super%N,beepbeep%super%N))

do iter=1,Nsim
    ! update lambda's
    if(DoBlockSampling) then
        ! generate candidate
        lcandid=lcurrent
        cov=(lsd(1)**2)*choles
        call normaldev(mean=lcurrent(:,K+1),do_dcmp=.false.,covar_dcmp=cov,&
                       gdev=lcandid(:,K+1),err=err,message=mess)
        ! compute post
        call GetLogPost(y,lcurrent,tcurrent,bcurrent,gcurrent,&
                        mucurrent,invcurrent,detcurrent,isDiag,0,0,pcurrent)
        call GetLogPost(y,lcandid,tcurrent,bcurrent,gcurrent,&
                        mucurrent,invcurrent,detcurrent,isDiag,0,0,pcandid)
        !apply Metropolis acceptance rule
        ratio=exp(min(max(-100.0_mrk, pcandid-pcurrent), 0.0_mrk)) 
        call uniran(u)
        if(u<=ratio) then
            lcurrent=lcandid !accept candidate
            lmr(:)=lmr(:)+1._mrk/real(Nadapt,mrk) ! increment move rate
        endif
    else
        do j=1,Nx
            ! generate candidate
            lcandid=lcurrent
            call normaldev(mean=lcurrent(j,K+1),sdev=lsd(j),gdev=dev,err=err,message=mess)
            lcandid(j,K+1)=dev
            ! compute post
            call GetLogPost(y,lcurrent,tcurrent,bcurrent,gcurrent,&
                            mucurrent,invcurrent,detcurrent,isDiag,0,j,pcurrent,beepbeep)
            call GetLogPost(y,lcandid,tcurrent,bcurrent,gcurrent,&
                            mucurrent,invcurrent,detcurrent,isDiag,0,j,pcandid,beepbeep)
            !apply Metropolis acceptance rule
            ratio=exp(min(max(-100.0_mrk, pcandid-pcurrent), 0.0_mrk)) 
            call uniran(u)
            if(u<=ratio) then
                lcurrent=lcandid !accept candidate
                lmr(j)=lmr(j)+1._mrk/real(Nadapt,mrk) ! increment move rate
            endif
        enddo
    endif
    saveMCMC(1:Nx)=lcurrent(:,K+1)
    
    ! update tau's
    if(K>0) then
        call applyConstraints(tcurrent(:,K),feas) ! make sure i-constraints are respected
        do i=1,Nt-2
            ! generate candidate
            tcandid=tcurrent
            call normaldev(mean=tcurrent(i,K),sdev=tsd(i),gdev=dev,err=err,message=mess)
            tcandid(i,K)=dev
            call applyConstraints(tcandid(:,K),feas)
            if(.not.feas) cycle
            ! compute post
            call trickedLogPost(y,lcurrent,tcurrent,bcurrent,gcurrent,&
                                mucurrent,invcurrent,detcurrent,isDiag,i,pcurrent)
            call trickedLogPost(y,lcurrent,tcandid,bcurrent,gcurrent,&
                                mucurrent,invcurrent,detcurrent,isDiag,i,pcandid)
            ratio=exp(min(max(-100.0_mrk, pcandid-pcurrent), 0.0_mrk)) 
            call uniran(u)
            if(u<=ratio) then
                tcurrent=tcandid !accept candidate
                tmr(i)=tmr(i)+1._mrk/real(Nadapt,mrk) ! increment move rate
            endif
        enddo
        saveMCMC((Nx+1):(Nx+Nt))=tcurrent(:,K)
    endif
    
    ! update beta's
    do j=1,size(beta)
        ! generate candidate
        bcandid=bcurrent
        call normaldev(mean=bcurrent(j),sdev=bsd(j),gdev=dev,err=err,message=mess)
        bcandid(j)=dev
        call GetHypermean(hypermodel,bcandid,mucandid)
        ! compute post
        call GetLogPost(y,lcurrent,tcurrent,bcurrent,gcurrent,&
                        mucurrent,invcurrent,detcurrent,isDiag,-1,-1,pcurrent,beepbeep)
        call GetLogPost(y,lcurrent,tcurrent,bcandid,gcurrent,&
                        mucandid,invcurrent,detcurrent,isDiag,-1,-1,pcandid,beepbeep)
        !apply Metropolis acceptance rule
        ratio=exp(min(max(-100.0_mrk, pcandid-pcurrent), 0.0_mrk)) 
        call uniran(u)
        if(u<=ratio) then !accept candidate
            bcurrent=bcandid 
            mucurrent=mucandid
            bmr(j)=bmr(j)+1._mrk/real(Nadapt,mrk) ! increment move rate
        endif
    enddo
    saveMCMC((Nx+Nt+1):(Nx+Nt+size(beta)))=bcurrent

    ! update gamma's
    do j=1,size(gamma)
        ! generate candidate
        gcandid=gcurrent
        call normaldev(mean=gcurrent(j),sdev=gsd(j),gdev=dev,err=err,message=mess)
        gcandid(j)=dev
        ! compute current post
        call GetLogPost(y,lcurrent,tcurrent,bcurrent,gcurrent,&
                        mucurrent,invcurrent,detcurrent,isDiag,-1,-1,pcurrent,beepbeep)
        ! Get hypervariance with candidate parameters
        if(isDiag .or. speedup%super%Nslim==1) then
            call GetHypervariance(D,hypermodel,gcandid,invcandid,detcandid,feas)
            if(.not.feas) cycle
        else ! only reduced variance is needed
            saveDet=beepbeep%super%sigmadet
            saveSigmaInv=beepbeep%super%sigmainv
            call GetHypervariance(D(beepbeep%super%indx,beepbeep%super%indx),&
                                  hypermodel,gcandid,beepbeep%super%sigmainv,&
                                  beepbeep%super%sigmadet,feas)
            if(.not.feas) then
                beepbeep%super%sigmadet=saveDet
                beepbeep%super%sigmainv=saveSigmaInv
                cycle
            endif
        endif
        call GetLogPost(y,lcurrent,tcurrent,bcurrent,gcandid,&
                        mucurrent,invcandid,detcandid,isDiag,-1,-1,pcandid,beepbeep)
        !apply Metropolis acceptance rule
        ratio=exp(min(max(-100.0_mrk, pcandid-pcurrent), 0.0_mrk)) 
        call uniran(u)
        if(u<=ratio) then !accept candidate
            gcurrent=gcandid 
            invcurrent=invcandid
            detcurrent=detcandid
            gmr(j)=gmr(j)+1._mrk/real(Nadapt,mrk) ! increment move rate
        else ! get back to old sigma/sigmadet in speedup%super
            if(.not.isDiag .and. speedup%super%Nslim>1) then
                beepbeep%super%sigmadet=saveDet
                beepbeep%super%sigmainv=saveSigmaInv
            endif
        endif
    enddo
    saveMCMC((Nx+Nt+size(beta)+1):(Nx+Nt+size(beta)+size(gamma)))=gcurrent
    ! update inverse matrixes in Speed-Up
    call SetSpeedUp(beepbeep,D,hypermodel,gcurrent)
    
    if(iter>Nburn .and. mod(iter-(Nburn+1),Nslim)==0) then ! save this iteration
        compt=compt+1        
        mcmc%sam(compt,:)=saveMCMC        
        ! compute full post at the end of the iteration
        call GetLogPost(y,lcurrent,tcurrent,bcurrent,gcurrent,&
                        mucurrent,invcurrent,detcurrent,isDiag,0,0,mcmc%post(compt),beepbeep)
    endif
    
    ! update jump sdev
    if(mod(iter,Nadapt)==0) then
        write(*,*) 'Done:',100*real(iter)/real(Nsim),'%'
        if(iter<=StopAdapt) then 
            call updateJumpSD(lmr,lsd)
            if(DoBlockSampling)  call GetCholes(D,hypermodel,gcurrent,choles,feas) !update choles
            call updateJumpSD(tmr(1:(Nt-2)),tsd(1:(Nt-2)))
            call updateJumpSD(bmr,bsd)
            call updateJumpSD(gmr,gsd)
        endif
        ! reinitialize move rates
        lmr=0._mrk;tmr=0._mrk;bmr=0._mrk;gmr=0._mrk
    endif
enddo

! clean up
deallocate(saveSigmaInv,saveMCMC)

end subroutine DoMCMC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine applyConstraints(tau,feas)
real(mrk),intent(inout)::tau(:)
logical,intent(out)::feas
! locals
integer(mik)::n
real(mrk)::u,v,delta

feas=.true.
n=size(tau)
u=sum(tau(1:(n-2)))
v=dot_product(tau(1:(n-2)),tau(1:(n-2)))
delta=2._mrk*n-2._mrk*v-u**2
if(delta<0._mrk) then;feas=.false.;return;endif
tau(n-1)=-0.5_mrk*(u+sqrt(delta))
tau(n)=-0.5_mrk*(u-sqrt(delta))
end subroutine applyConstraints
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine trickedLogPost(y,lambda,tau,beta,gamma,mu,sigmainv,sigmadet,isDiag,t,logp)
real(mrk), intent(in)::lambda(:,:),tau(:,:),beta(:),gamma(:),&
                       mu(:),sigmainv(:,:),sigmadet
integer(mik),intent(in)::y(:,:),t
logical,intent(in)::isDiag
real(mrk), intent(out)::logp
! locals
integer(mik)::Nt
real(mrk)::p1,p2,p3

Nt=size(y,dim=2)
call GetLogPost(y,lambda,tau,beta,gamma,mu,sigmainv,sigmadet,isDiag,t,0,p1)
call GetLogPost(y,lambda,tau,beta,gamma,mu,sigmainv,sigmadet,isDiag,Nt-1,0,p2)
call GetLogPost(y,lambda,tau,beta,gamma,mu,sigmainv,sigmadet,isDiag,Nt,0,p3)
logp=p1+p2+p3
end subroutine trickedLogPost
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetMaxpost(mcmc,post,maxpost)
real(mrk),intent(in)::mcmc(:,:),post(:)
real(mrk), intent(out)::maxpost(size(mcmc,dim=2))
! locals
real(mrk)::ml(1)

ml=Maxloc(post)
maxpost=mcmc(ml(1),:)
end subroutine GetMaxpost
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WriteMatrix(file,M,sep)
use utilities_dmsl_kit, only:getspareunit
character(*),intent(in)::file,sep
real(mrk),intent(in)::M(:,:)
! locals
integer(mik)::nrow,ncol,unt,err,i,j
character(250)::mess

!preliminaries
nrow=size(M,dim=1);ncol=size(M,dim=2)
call getspareunit(unt,err,mess)
open(unit=unt,file=trim(file),status='replace')

! write
do i=1,nrow
    if(trim(sep)=="\t") then ! tabulated
        write(unt,'(<ncol>e14.6,A)') M(i,:)
    else ! user-defined separator
        do j=1,ncol
            write(unt,'(e14.6,A)',advance='NO') M(i,j),trim(sep)
        enddo
        write(unt,'(A)',advance='YES') ""
    endif
enddo
close(unt)
end subroutine WriteMatrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mnorm_logp(x,mu,sigmainv,sigmadet)
use utilities_dmsl_kit,only:lnTwoPi
real(mrk)::mnorm_logp
real(mrk), intent(in)::x(:),mu(:),sigmainv(:,:),sigmadet
!locals
integer(mik)::n
real(mrk)::centred(size(x)),z

n=size(x)
centred=x-mu
z=dot_product(matmul(centred,sigmainv),centred)
mnorm_logp=-0.5_mrk*(real(n,mrk)*lnTwoPi+sigmadet+z)
end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine InitMCMC(Nsim,Nburn,Nslim,Nx,Nt,Nbeta,Ngamma,mcmc)
integer(mik),intent(in)::Nsim,Nburn,Nslim,Nx,Nt,Nbeta,Ngamma
type(MCMCtype),intent(out)::mcmc
! local
integer(mik)::n,i

n=size((/(i,i= Nburn+1,Nsim,Nslim)/))
mcmc%Nsim=n
mcmc%Nt=Nt
mcmc%Nx=Nx
mcmc%Nbeta=Nbeta
mcmc%Ngamma=Ngamma
allocate(mcmc%sam(n,Nx+Nt+Nbeta+Ngamma))
allocate(mcmc%post(n))

end subroutine InitMCMC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CopyMCMC(original,copy)
type(MCMCtype),intent(in)::original
type(MCMCtype),intent(out)::copy

copy%Nsim=original%Nsim
copy%Nt=original%Nt
copy%Nx=original%Nx
copy%Nbeta=original%Nbeta
copy%Ngamma=original%Ngamma
allocate(copy%sam(copy%Nsim,copy%Nx+copy%Nt+copy%Nbeta+copy%Ngamma))
allocate(copy%post(copy%Nsim))
copy%sam=original%sam
copy%post=original%post

end subroutine CopyMCMC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetHypermean(hypermodel,beta,mu)
type(HyperModelType),intent(in)::hypermodel
real(mrk),intent(in)::beta(:)
real(mrk),intent(out)::mu(:)
mu=matmul(hypermodel%Z,beta)
end subroutine GetHypermean
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetHypervariance(D,hypermodel,gamma,sigmainv,sigmadet,feas)
use linalg_dmsl_kit,only:choles_invrt,LU_invrt
real(mrk),intent(in)::D(:,:)
type(HyperModelType),intent(in)::hypermodel
real(mrk),intent(in)::gamma(:)
real(mrk),intent(out)::sigmainv(:,:),sigmadet
logical,intent(out)::feas
! locals
real(mrk),parameter:: epsilon=0.000001_mrk
integer(mik)::err,i
character(250)::mess

feas=.true.
if(any(gamma<=0._mrk)) then;feas=.false.;return;endif
select case(trim(hypermodel%Dmodel))
case("inde")
    sigmainv=0._mrk;forall(i=1:size(D,dim=1)) sigmainv(i,i)=1._mrk/(gamma(1)**2)
    sigmadet=real(2*size(D,dim=1),mrk)*log(gamma(1))
case("exponential")
    ! get correlation matrix
    sigmainv=exp(-1._mrk*D/gamma(2))
    ! invert it
    call choles_invrt(ainv=sigmainv,logDet=sigmadet,posDefinite=feas,err=err,message=mess)
    ! COMMENTED OUT - Regularization appears to mess up the log-determinant
    ! If Choleski doesn't work, then try regularization
    !if(.not.feas) then
        !forall(i=1:size(D,dim=1)) sigmainv(i,i)=1._mrk+epsilon
        !call choles_invrt(ainv=sigmainv,logDet=sigmadet,posDefinite=feas,err=err,message=mess)
    !endif
    ! If nothing works, then try LU
    if(.not.feas) then
        call LU_invrt(ainv=sigmainv,logDet=sigmadet,sing=feas,err=err,message=mess)
        feas=.not.feas
    endif
    ! Get sigmainv
    sigmainv=sigmainv/(gamma(1)**2)
    sigmadet=sigmadet+real(2*size(D,dim=1),mrk)*log(gamma(1))
end select
end subroutine GetHypervariance
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetHyperCorr(D,hypermodel,gamma,corr,feas)
real(mrk),intent(in)::D(:,:)
type(HyperModelType),intent(in)::hypermodel
real(mrk),intent(in)::gamma(:)
real(mrk),intent(out)::corr(:,:)
logical,intent(out)::feas
! locals
integer(mik)::err,i
character(250)::mess

feas=.true.;Corr=undefRN
if(any(gamma<=0._mrk)) then;feas=.false.;return;endif
select case(trim(hypermodel%Dmodel))
case("inde")
    corr=0._mrk;forall(i=1:size(D,dim=1)) corr(i,i)=1._mrk
case("exponential")
    ! get correlation matrix
    corr=exp(-1._mrk*D/gamma(2))
end select
end subroutine GetHyperCorr
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetCholes(D,hypermodel,gamma,choles,feas)
use linalg_dmsl_kit,only:choles_dcmp
real(mrk),intent(in)::D(:,:)
type(HyperModelType),intent(in)::hypermodel
real(mrk),intent(in)::gamma(:)
real(mrk),intent(out)::choles(:,:)
logical,intent(out)::feas
! locals
integer(mik)::err,i
character(250)::mess

feas=.true.
if(any(gamma<=0._mrk)) then;feas=.false.;return;endif
! get covariance matrix
select case(trim(hypermodel%Dmodel))
case("inde")
    choles=0._mrk;forall(i=1:size(D,dim=1)) choles(i,i)=1._mrk/(gamma(1)**2)
case("exponential")
    choles=exp(-1._mrk*D/gamma(2))
end select
call choles_dcmp(A=choles, posDefinite=feas, err=err, message=mess)
end subroutine GetCholes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine updateJumpSD(mr,sd)
real(mrk),intent(in)::mr(:)
real(mrk),intent(inout)::sd(:)
! locals
integer(mik)::j

do j=1,size(mr)
    if(mr(j)<=minMR) sd(j)=sd(j)*DownMult
    if(mr(j)>=maxMR) sd(j)=sd(j)*UpMult
enddo
end subroutine updateJumpSD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SetSpeedUp(speedup,D,hypermodel,gamma)
use numerix_dmsl_kit,only:indexx_qsort
use utilities_dmsl_kit,only:trueLocIndx_f
Type(SpeedUpType),intent(inout)::speedup
real(mrk),intent(in)::D(:,:)
type(HyperModelType),intent(in)::hypermodel
real(mrk),intent(in)::gamma(:)
! locals
integer(mik)::Nx,i,indx(size(D,dim=1)),err,Nn
character(250)::mess
logical::feas,mask(size(D,dim=1))
real(mrk)::corr(size(D,dim=1),size(D,dim=2))

Nx=size(D,dim=1)
if(.not.allocated(speedup%prop)) then ! initialize
    ! super-speedup strategy
    speedup%super%N=size((/(i,i=1,Nx,speedup%super%Nslim)/))
    allocate(speedup%super%indx(speedup%super%N))
    speedup%super%indx=(/(i,i=1,Nx,speedup%super%Nslim)/)
    if(speedup%super%Nslim>1) then
        allocate(speedup%super%sigmainv(speedup%super%N,speedup%super%N))
        call GetHypervariance(D(speedup%super%indx,speedup%super%indx),&
                            hypermodel,gamma,&
                            speedup%super%sigmainv,speedup%super%sigmadet,feas)
    endif
    ! speed-up strategy
    allocate(speedup%prop(Nx))
    if(speedup%option==0) return ! no speed up, nothing more to do
    if(speedup%option==1) then ! fixed-N strategy: allocate properties once and for all
        do i=1,Nx
            ! allocate arrays
            allocate(speedup%prop(i)%indx(speedup%N))
            allocate(speedup%prop(i)%sigmainv(speedup%N,speedup%N))
            ! find the N closest neighbors of site i
            call indexx_qsort(arr=D(i,:),indx=indx,err=err,message=mess)
            speedup%prop(i)%indx=indx(1:speedup%N)
        enddo  
    endif
endif

! computation of indexes, inverse matrixes and determinant
if(speedup%option==2) call GetHyperCorr(D,hypermodel,gamma,corr,feas) ! needed for fixed-corr option only
do i=1,Nx
    if(speedup%option==2) then ! fixed-corr option: select sites with correlation larger than threshold
        mask=corr(i,:)>=speedup%corr
        Nn=count(mask)
        if(allocated(speedup%prop(i)%indx)) deallocate(speedup%prop(i)%indx)
        allocate(speedup%prop(i)%indx(Nn))
        if(allocated(speedup%prop(i)%sigmainv)) deallocate(speedup%prop(i)%sigmainv)
        allocate(speedup%prop(i)%sigmainv(Nn,Nn))
        speedup%prop(i)%indx=trueLocIndx_f(mask)
    endif
    ! get inverse matrix and determinant
    call GetHypervariance(D(speedup%prop(i)%indx,speedup%prop(i)%indx),&
                            hypermodel,gamma,&
                            speedup%prop(i)%sigmainv,speedup%prop(i)%sigmadet,feas)
enddo

end subroutine SetSpeedUp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_Read(file,HCI,err,mess)                                
!^**********************************************************************
!^* Purpose: Read main config file
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 07/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file
!^* INOUT
!^*		1.HCI object
!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
type(HCIType),intent(inout)::HCI
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Read'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
! workspace
read(unt,*,iostat=err) HCI%workspace
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! run options
read(unt,*,iostat=err) HCI%Config_runoptions
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! config file for model
read(unt,*,iostat=err) HCI%Config_Model
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! config file for Y data
read(unt,*,iostat=err) HCI%Config_Y
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! config file for MCMC
read(unt,*,iostat=err) HCI%Config_MCMC
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! config file for W data
read(unt,*,iostat=err) HCI%Config_W
if(err/=0) then;mess=trim(procname)//':read error';return;endif
close(unt)

end subroutine Config_Read  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_RunOptions_Read(file,HCI,err,mess)                                
!^**********************************************************************
!^* Purpose: Read config file for RunOptions
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 25/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file
!^* INOUT
!^*		1. HCI object
!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
type(HCItype),intent(inout)::HCI
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_RunOptions_Read'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
! read
read(unt,*,iostat=err) HCI%RunOptions(1)
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) HCI%RunOptions(2)
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) HCI%RunOptions(3)
if(err/=0) then;mess=trim(procname)//':read error';return;endif
close(unt)

end subroutine Config_RunOptions_Read  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_Model_Read(file,HCI,err,mess)                                
!^**********************************************************************
!^* Purpose: Read config file for model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 07/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file
!^* INOUT
!^*		1. HCI object
!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
type(HCItype),intent(inout)::HCI
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Model_Read'
integer(mik)::unt,K,i
character(250):: Zfile

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
! Number of components
read(unt,*,iostat=err) K
if(err/=0) then;mess=trim(procname)//':read error';return;endif
HCI%K=K
! Hypermodel for each component
if(allocated(HCI%hypermodel)) deallocate(HCI%hypermodel)
allocate(HCI%hypermodel(K+1))
do i=1,K+1
    ! not used
    read(unt,*,iostat=err) 
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    ! Z model
    read(unt,*,iostat=err) HCI%hypermodel(i)%Zmodel
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    ! Z config file
    read(unt,*,iostat=err) Zfile
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    ! read Z data
    call Config_Z_Read(file=trim(HCI%workspace)//trim(Zfile),&
                       workspace=trim(HCI%workspace),&
                       hypermodel=HCI%hypermodel(i),&
                       err=err,mess=mess) 
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    ! number of parameters for Z-model
    read(unt,*,iostat=err) HCI%hypermodel(i)%nbeta
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    ! D-model
    read(unt,*,iostat=err) HCI%hypermodel(i)%Dmodel
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
    ! number of parameters for D-model
    read(unt,*,iostat=err) HCI%hypermodel(i)%ngamma
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
enddo
close(unt)

end subroutine Config_Model_Read  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_Z_Read(file,workspace,hypermodel,err,mess)                                
!^**********************************************************************
!^* Purpose: Read config file for Z covariates
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 08/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file
!^*		2.workspace
!^* INOUT
!^*		1. hypermodel object
!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use DataRW_tools, only:ReadSeparatedFile
character(*), intent(in)::file,workspace
type(HyperModeltype),intent(inout)::hypermodel
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Z_Read'
integer(mik)::unt,Nrow,Ncol
real(mrk),pointer::values(:,:)

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
! Z data file
read(unt,*,iostat=err) hypermodel%Zfile
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Z dimensions
read(unt,*,iostat=err) Nrow,Ncol
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Read Z data
call ReadSeparatedFile(file=trim(workspace)//trim(hypermodel%Zfile),&
                       sep=sep,nhead=0,ncol=Ncol,y=values,&
                       err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(allocated(hypermodel%Z)) deallocate(hypermodel%Z)
allocate(hypermodel%Z(size(values,dim=1),Ncol))
hypermodel%Z=values
! Z names
if(allocated(hypermodel%Znames)) deallocate(hypermodel%Znames)
allocate(hypermodel%Znames(Ncol))
read(unt,*,iostat=err) hypermodel%Znames
if(err/=0) then;mess=trim(procname)//':read error';return;endif
close(unt)

end subroutine Config_Z_Read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_Y_Read(file,HCI,err,mess)                                
!^**********************************************************************
!^* Purpose: Read config file for Y data
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 08/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file
!^* INOUT
!^*		1. HCI object
!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use DataRW_tools, only:ReadSeparatedFile
use Geodesy_tools,only:GetDistance
character(*), intent(in)::file
type(HCItype),intent(inout)::HCI
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Y_Read'
integer(mik)::unt,unt2,i
character(250)::gridfile,Dmethod
real(mrk),pointer::values(:,:)

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
! Y data file
read(unt,*,iostat=err) HCI%Y%file
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Y dimensions
read(unt,*,iostat=err) HCI%Y%Ns,HCI%Y%Nt
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Read data
if(allocated(HCI%Y%val)) deallocate(HCI%Y%val)
allocate(HCI%Y%val(HCI%Y%Ns,HCI%Y%Nt))
call ReadSeparatedFile(file=trim(HCI%workspace)//trim(HCI%Y%file),&
                       sep=sep,nhead=0,ncol=HCI%Y%Nt,y=values,&
                       err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(size(values,dim=2)/=HCI%Y%Nt) then
    err=1;mess=trim(procname)//':incorrect size [Y]';return
endif
HCI%Y%val=nint(values)
! Read station names
if(allocated(HCI%Y%rownames)) deallocate(HCI%Y%rownames)
allocate(HCI%Y%rownames(HCI%Y%Ns))
read(unt,*,iostat=err) gridfile
if(err/=0) then;mess=trim(procname)//':read error';return;endif
call getSpareUnit(unt2,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt2,file=trim(HCI%workspace)//trim(gridfile), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
do i=1,HCI%Y%Ns
    read(unt2,*,iostat=err) HCI%Y%rownames(i)
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
enddo
close(unt2)
! Read time grid
read(unt,*,iostat=err) gridfile
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(allocated(HCI%Y%tgrid)) deallocate(HCI%Y%tgrid)
allocate(HCI%Y%tgrid(HCI%Y%Nt))
call ReadSeparatedFile(file=trim(HCI%workspace)//trim(gridfile),&
                       sep=sep,nhead=0,ncol=1,y=values,&
                       err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(size(values,dim=1)/=HCI%Y%Nt) then
    err=1;mess=trim(procname)//':incorrect size [tgrid]';return
endif
HCI%Y%tgrid=values(:,1)
! Read space grid
read(unt,*,iostat=err) gridfile
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(allocated(HCI%Y%sgrid)) deallocate(HCI%Y%sgrid)
allocate(HCI%Y%sgrid(HCI%Y%Ns,2))
call ReadSeparatedFile(file=trim(HCI%workspace)//trim(gridfile),&
                       sep=sep,nhead=0,ncol=2,y=values,&
                       err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(size(values,dim=1)/=HCI%Y%Ns) then
    err=1;mess=trim(procname)//':incorrect size [sgrid]';return
endif
HCI%Y%sgrid=values
! Distance computation
read(unt,*,iostat=err) Dmethod
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(allocated(HCI%Y%D)) deallocate(HCI%Y%D)
allocate(HCI%Y%D(HCI%Y%Ns,HCI%Y%Ns))
call GetDistance(formula=Dmethod,pts=HCI%Y%sgrid,D=HCI%Y%D,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
! MV flag
read(unt,*,iostat=err) HCI%Y%mv
if(err/=0) then;mess=trim(procname)//':read error';return;endif
close(unt)

end subroutine Config_Y_Read  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_MCMC_Read(file,HCI,err,mess)                                
!^**********************************************************************
!^* Purpose: Read config file for MCMC
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 08/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file
!^* INOUT
!^*		1. HCI object
!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
type(HCItype),intent(inout)::HCI
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_MCMC_Read'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
read(unt,*,iostat=err) HCI%mcmc%Nsim
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) HCI%mcmc%Nadapt
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) HCI%mcmc%StopAdapt
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) HCI%mcmc%Nburn
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) HCI%mcmc%Nslim
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) HCI%mcmc%DoBlockSampling
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) HCI%speedup%option
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(HCI%speedup%option==1) then
    read(unt,*,iostat=err) HCI%speedup%N
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
else if(HCI%speedup%option==2) then
    read(unt,*,iostat=err) HCI%speedup%corr
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
else
    read(unt,*,iostat=err) 
    if(err/=0) then;mess=trim(procname)//':read error';return;endif
endif
read(unt,*,iostat=err) HCI%speedup%super%Nslim
if(err/=0) then;mess=trim(procname)//':read error';return;endif

close(unt)

end subroutine Config_MCMC_Read  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_W_Read(file,HCI,err,mess)                                
!^**********************************************************************
!^* Purpose: Read config file for W data
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 25/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file
!^* INOUT
!^*		1. HCI object
!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use DataRW_tools, only:ReadSeparatedFile
character(*), intent(in)::file
type(HCItype),intent(inout)::HCI
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_W_Read'
integer(mik)::unt
character(250)::gridfile
real(mrk),pointer::values(:,:)

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
! W data file
read(unt,*,iostat=err) HCI%W%file
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! W dimensions
read(unt,*,iostat=err) HCI%W%Ns,HCI%W%Nt
if(err/=0) then;mess=trim(procname)//':read error';return;endif
! Read data
if(allocated(HCI%W%val)) deallocate(HCI%W%val)
allocate(HCI%W%val(HCI%W%Ns,HCI%W%Nt))
call ReadSeparatedFile(file=trim(HCI%workspace)//trim(HCI%W%file),&
                       sep=sep,nhead=0,ncol=HCI%W%Nt,y=values,&
                       err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(size(values,dim=2)/=HCI%W%Nt) then
    err=1;mess=trim(procname)//':incorrect size [W]';return
endif
HCI%W%val=nint(values)
! Read space grid
read(unt,*,iostat=err) gridfile
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(allocated(HCI%W%sgrid)) deallocate(HCI%W%sgrid)
allocate(HCI%W%sgrid(HCI%W%Ns,2))
call ReadSeparatedFile(file=trim(HCI%workspace)//trim(gridfile),&
                       sep=sep,nhead=0,ncol=2,y=values,&
                       err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(size(values,dim=1)/=HCI%W%Ns) then
    err=1;mess=trim(procname)//':incorrect size [sgrid]';return
endif
HCI%W%sgrid=values
close(unt)

end subroutine Config_W_Read  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end module HCI_Bernoulli_tools
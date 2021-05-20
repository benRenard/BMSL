program TestRegCop
use kinds_dmsl_kit
use RegionalRegMCMC
use RegionalSpatialType
use RegionalModel_tools2
use Distribution_tools, only:UNIF,GAUSS, GEV, Generate
use BayesianEstimation_tools, only: PriorListType
use InOutTools
use Copula_tools
use utilities_dmsl_kit, only:number_string
implicit none


character(500)::mess
integer(mik)::err,i,j,temp,loop
logical::feas,isnull

! Model
type(ModelType)::Model
character(10)::distId
character(250)::DistName

! Data
real(mrk),pointer::X(:,:),seaDist(:,:),distMat(:,:)
integer(mik)::N_time,N_site,filenumber

! Spatial
type::SpatType
type(SpatialParType),pointer::SpatialPar(:)=>NULL()
end type SpatType
type(SpatType),allocatable::SpatPar(:)
real(mrk),allocatable::s_cov(:,:)
character(100)::spatName,SpatReg
real(mrk),allocatable::mcmcPar(:)

! Dpar
type(DparType),allocatable::Dpar(:)
character(250),allocatable::Name(:),Ivslink(:),regression(:)
real(mrk),allocatable::bounds(:,:)
type::Index_regionalParType
integer(mik),pointer::Indice(:)=>NULL()
end type Index_regionalParType
type(Index_regionalParType),allocatable::Index_regionalPar(:)
!integer(mik),allocatable::Index_regionalPar1(:),Index_regionalPar2(:),Index_regionalPar3(:)

! Covariates
type(covariateType),pointer::cov(:)
real(mrk),allocatable::t_cov(:,:)

type(tCovType),allocatable::alltCov(:)
type(sCovType),allocatable::allsCov(:)
type(stCovType),allocatable::allstCov(:)
type(allCovariateType),pointer::allcov(:)

! Copula
type(CopulaType)::copula
character(250)::CopName,dependID
real(mrk),pointer:: CoSiteInfo(:,:),CoSiteInfoTemp(:,:)

! RegMCMC
integer(mik)::Ntotpar
type(PriorListType),allocatable:: PL(:)
real(mrk),allocatable::start(:),startSTD(:)
character(1)::outdoc

!MCMC parameters
character(100)::OutFile,inFile
real(mrk)::mv,MinMoveRate,MaxMoveRate,DownMult,UpMult
integer(mik)::nAdapt,nCycles


! Main
err=0;mess='';feas=.true.

N_time=50_mik;N_site=20_mik
filenumber=45_mik
inFile="C:\Mydoc\Forge\hydrostat\trunk\Fortran\test\Tests_IVF\TestRegionalCopula\Data\"

do loop=1,filenumber

call mat_input(Filename=trim(inFile)//"SimuData_"//trim(number_string(loop))//".txt", row=N_time, column=N_site, header=.true., tab=x, err=err, mess=mess)
If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif


!X(1:10,1:5)=-9999._mrk
!x(11:20,6:8)=-9999._mrk
!x(21:30,9:10)=-9999._mrk

!##################################################################################
!############### 
! Initialize settings
distId=GAUSS
DistName='GAUSS'
OutFile='C:\Mydoc\Forge\hydrostat\trunk\Fortran\test\Tests_IVF\TestRegionalCopula\MCMC\'
mv=-9999._mrk;MinMoveRate=0.1_mrk;MaxMoveRate=0.5_mrk;DownMult=0.9_mrk;UpMult=1.1_mrk
nAdapt=100;nCycles=200
!###############
!################################################################################## 


! Initialization
! Model
call CompleteMT(DistName=DistName,DistID=distID,MT=Model,err=err,mess=mess) !OK
If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif

!##################################################################################
!############### 
! Spatial
if(allocated(SpatPar)) deallocate(SpatPar)
allocate(SpatPar(Model%N))

if(associated(SpatPar(1)%SpatialPar)) deallocate(SpatPar(1)%SpatialPar)
allocate(SpatPar(1)%SpatialPar(2))
if(associated(SpatPar(2)%SpatialPar)) deallocate(SpatPar(2)%SpatialPar)
allocate(SpatPar(2)%SpatialPar(1))

SpatName=undefCH
SpatReg='Exp'

if(allocated(s_cov)) deallocate(s_cov)
allocate(s_cov(1,N_site))
s_cov(1,:)=1._mrk


call InitializeSpatialpar(Name=SpatName, SpatReg='Linear_L1', SpatialPar=SpatPar(1)%SpatialPar(1),s_cov=s_cov(1:1,:),err=err, mess=mess)
If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif
call InitializeSpatialpar(Name=SpatName, SpatReg='Linear_L1', SpatialPar=SpatPar(1)%SpatialPar(2),s_cov=s_cov(1:1,:),err=err, mess=mess)
If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif
call InitializeSpatialpar(Name=SpatName, SpatReg='Linear_L1', SpatialPar=SpatPar(2)%SpatialPar(1),s_cov=s_cov(1:1,:),err=err, mess=mess)
If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif


!###############
!################################################################################## 



! Dpar
if(allocated(Dpar)) deallocate(Dpar)
allocate(Dpar(Model%N))

if(allocated(Name)) deallocate(Name)
allocate(Name(Model%N))
if(allocated(Ivslink)) deallocate(Ivslink)
allocate(Ivslink(Model%N))
if(allocated(regression)) deallocate(regression)
allocate(regression(Model%N))
if(allocated(bounds)) deallocate(bounds)
allocate(bounds(Model%N,2))
!##################################################################################
!############### 
!Name=(/'location','scale','shape'/)
!!Bounds(2,:)=(/0._mrk,HugeRe/)
!Ivslink=(/'Identity','Identity','Identity'/)
!regression=(/'Regional_1','Identity','Identity'/)
Name=(/'mean','stddev'/)
!Bounds(2,:)=(/0._mrk,HugeRe/)
Ivslink=(/'Identity','Identity'/)
regression=(/'Regional_1','Identity'/)


!###############
!################################################################################## 


if(allocated(Index_regionalPar)) deallocate(Index_regionalPar)
allocate(Index_regionalPar(Model%N))

!##################################################################################
!############### 
if(associated(Index_regionalPar(1)%Indice)) deallocate(Index_regionalPar(1)%Indice)
allocate(Index_regionalPar(1)%Indice(2))
Index_regionalPar(1)%Indice=(/1_mik,2_mik/)

if(associated(Index_regionalPar(2)%Indice)) deallocate(Index_regionalPar(2)%Indice)
allocate(Index_regionalPar(2)%Indice(1))
Index_regionalPar(2)%Indice=(/1_mik/)
!###############
!################################################################################## 


do i=1, Model%N
    if (associated(Index_regionalPar(i)%Indice)) then
        call InitializeDpar(Name=Name(i), Ivslink=Ivslink(i), regression=regression(i), Index_regionalPar=Index_regionalPar(i)%Indice,&
            SpatialPar=SpatPar(i)%SpatialPar, N_site=N_site, Dpar=Dpar(i), feas=feas, err=err, mess=mess)
        If(err>0) then
            mess='MeanProg error: '//trim(mess)
            write (*,*) mess
            read (*,*)
        endif

    else
        call InitializeDpar(Name=Name(i), Ivslink=Ivslink(i), regression=regression(i), &
            SpatialPar=SpatPar(i)%SpatialPar, N_site=N_site, Dpar=Dpar(i), feas=feas, err=err, mess=mess)
        If(err>0) then
            mess='MeanProg error: '//trim(mess)
            write (*,*) mess
            read (*,*)
        endif
    endif
enddo

!Copula
!##################################################################################
!############### 
CopName="Gaussian"
dependID='Power'

if (associated(distMat)) deallocate(distMat)
allocate(distMat(N_site,N_site))

forall(i=1:N_site,j=1:N_site) distMat(i,j)=abs(i-j)

!###############
!##################################################################################  
call InitializeCopula(CopName, N_time, dependID, distMat, Copula, err, mess)
If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif


! covariates
!##################################################################################
!############### 
if(allocated(t_cov)) deallocate(t_cov)
allocate(t_cov(Dpar(1)%N_cov(1),N_time))
t_cov(1,:)=(/(i,i=1,N_time)/)


if(allocated(alltCov)) deallocate(alltCov)
allocate(alltCov(Model%N))
if(allocated(allsCov)) deallocate(allsCov)
allocate(allsCov(Model%N))
if(allocated(allstCov)) deallocate(allstCov)
allocate(allstCov(Model%N))

if(associated(alltCov(1)%t_covS)) deallocate(alltCov(1)%t_covS)
allocate(alltCov(1)%t_covS(1,N_time))
alltCov(1)%t_covS=t_cov

if(associated(allsCov(1)%s_covS)) deallocate(allsCov(1)%s_covS)
allocate(allsCov(1)%s_covS(1,N_site))
allsCov(1)%s_covS(1,:)=1

!if(associated(allsCov(3)%s_covS)) deallocate(allsCov(3)%s_covS)
!allocate(allsCov(3)%s_covS(1,N_site))
!allsCov(3)%s_covS(1,:)=1
!!temp try
!if(associated(alltCov(3)%t_covS)) deallocate(alltCov(3)%t_covS)
!allocate(alltCov(3)%t_covS(1,N_time))
!alltCov(3)%t_covS=t_cov

!###############
!################################################################################## 


call allCovInput(tCovs=alltCov, sCovs=allsCov, stCovs=allstCov, Nmodel=model%N,&
    N_time=N_time, Nsite=N_site, covariate=allcov, err=err, mess=mess)
If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif 


!##################################################################################
!############### 
!prior
Ntotpar=0
do i=1, Model%N
    Ntotpar=Ntotpar+Dpar(i)%N_totpar
enddo
Ntotpar=Ntotpar+Copula%N_par
if(allocated(PL)) deallocate(PL)
allocate(PL(Ntotpar))
do i=1, Ntotpar
    PL(i)%dist=Unif
    if(allocated(PL(i)%par)) deallocate(PL(i)%par)
    allocate(PL(i)%par(2))
    PL(i)%par=(/0._mrk,500._mrk/)
enddo


PL(1)%par=(/20._mrk,20.001_mrk/)
PL(3)%par=(/5._mrk,5.001_mrk/)
PL(Ntotpar)%par=(/0._mrk,1._mrk/)

if(allocated(start)) deallocate(start)
allocate(start(Ntotpar))
if(allocated(startStd)) deallocate(startStd)
allocate(startStd(Ntotpar))

start(1)=20._mrk
start(2)=1._mrk
start(3)=5._mrk
start(4)=0.5_mrk

startStd(1)=0.2_mrk
startStd(2)=0.1_mrk
startStd(3)=0.01_mrk
startStd(4)=0.05_mrk
!###############
!################################################################################## 
!
call RegMCMC(X=x, Model=Model, Dpar=Dpar, Copula=Copula, Covariates=allcov, PriorList=PL,&
            start=start, startStd=startStd, &
            nAdapt=nAdapt,nCycles=nCycles,&
            MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
            DownMult=DownMult,UpMult=UpMult,&
            OutFile=trim(OutFile)//"SimulationCop"//trim(number_string(loop))//".txt", &
            ! error handling
            err=err,mess=mess)

If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif



Copula%N_par=0

!prior
Ntotpar=0
do i=1, Model%N
    Ntotpar=Ntotpar+Dpar(i)%N_totpar
enddo
Ntotpar=Ntotpar+Copula%N_par
if(allocated(PL)) deallocate(PL)
allocate(PL(Ntotpar))
do i=1, Ntotpar
    PL(i)%dist=Unif
    if(allocated(PL(i)%par)) deallocate(PL(i)%par)
    allocate(PL(i)%par(2))
    PL(i)%par=(/0._mrk,500._mrk/)
enddo


PL(1)%par=(/20._mrk,20.001_mrk/)
PL(3)%par=(/5._mrk,5.001_mrk/)

!PL(Ntotpar)%par=(/0._mrk,1._mrk/)

if(allocated(start)) deallocate(start)
allocate(start(Ntotpar))
if(allocated(startStd)) deallocate(startStd)
allocate(startStd(Ntotpar))

start(1)=20._mrk
start(2)=1._mrk
start(3)=5._mrk

startStd(1)=0.2_mrk
startStd(2)=0.1_mrk
startStd(3)=0.01_mrk
!###############
!################################################################################## 
!Copula=Copula,
call RegMCMC(X=x, Model=Model, Dpar=Dpar,  Covariates=allcov, PriorList=PL,&
            start=start, startStd=startStd, &
            nAdapt=nAdapt,nCycles=nCycles,&
            MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
            DownMult=DownMult,UpMult=UpMult,&
            OutFile=trim(OutFile)//"Simulation"//trim(number_string(loop))//".txt", &
            ! error handling
            err=err,mess=mess)

If(err>0) then
    mess='MeanProg error: '//trim(mess)
    write (*,*) mess
    read (*,*)
endif






enddo
end program TestRegCop
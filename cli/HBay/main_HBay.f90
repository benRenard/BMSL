program Hbay

use kinds_dmsl_kit
use BayesianEstimation_tools
use DataRW_tools, only:DatRead,DatWrite
Use Distribution_tools, only:GetRoughEstimate,GetQuantile
use EmpiricalStats_tools, only:GetEmpiricalDistSummary,p2T
use HistoricBay_tools, only:HistoricBayType, ReadHistoricData,HistoricBay_EstimPar,HistoricBay_Predictive,FatalExit
use Estimator_tools
implicit none

!character(250), parameter::exepath="..\..\Core\HBay\" ! CORE version
!character(250), parameter::exepath="..\..\..\Core\HBay\"! DEBUG/RELEASE
character(250), parameter::exepath="" ! exe-only version
character(250), parameter::ConfigFile="Config.txt"
character(250), parameter::ConfigFolder="" ! "Config\"
character(250), parameter::Config_Data_File="Config_Data.txt"
character(250), parameter::Config_Inference_File="Config_Inference.txt"
character(250), parameter::Config_MCMC_File="Config_MCMC.txt"
character(250), parameter::Config_ResultFiles_File="Config_ResultFiles.txt"
character(250), parameter::Config_ResultOptions_File="Config_ResultOptions.txt"
character(250), parameter::Config_SystematicErrors_File="Config_SystematicErrors.txt"
character(250), parameter::MCMC_tempFile="MCMC.txt"
integer(mik),parameter::ncol_DATfile=3,obscol=2

character(250):: workspace,Cfolder,mess,DataFile,DataType,&
                 mvOption,O2_cond,dist,headers(ncol_DATfile),line
integer(mik)::err,npar,i,errcode,nrow,nSystErr,Ninfer
real(mrk)::O1_value,O2_value,xmini,xmaxi
type(PriorListType), pointer:: priors(:)
character(250), pointer:: parnames(:),SystErrNames(:)
character(250), allocatable:: MCMChead(:)
real(mrk), pointer:: y(:,:),x(:,:),mcmc(:,:),LogPost(:),pred(:)
real(mrk), allocatable::teta0(:),teta_std0(:),mode(:),xlist(:),plist(:),&
    pdf_modal(:),pdf_pred(:),cdf_modal(:),cdf_pred(:),Qp_modal(:),Qp_pred(:),&
    IC_modal(:,:),matrice(:,:),SlimMCMC(:,:)
logical::feas
type(HistoricBayType)::Hb
!!!!!!!!!!!!!!!!!!!
! MCMC settings
integer(mik):: Nadapt,Ncycles,Nslim,Nrep
real(mrk)::BurnFactor,MinMoveRate,MaxMoveRate,DownMult,UpMult,MultFactor,Eps
! END MCMC settings
!!!!!!!!!!!!!!!!!!!
! ResultFiles
character(250)::mc2_File,par_File,dat_File,pdf_File,&
                           cdf_File,prd_File,qtl_File,emp_File
logical::mc2_Save,par_Save,dat_Save,pdf_Save,cdf_Save,&
                          prd_Save,qtl_Save,emp_Save
! END ResultFiles
!!!!!!!!!!!!!!!!!!!
! ResultOptions
real(mrk)::xmin,xmax,xn,pmin,pmax,pn,level,PtoT,step
character(250)::pp,kernel
! END ResultOptions
!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                      STEP 1: READ CONFIG FILES                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write(*,*) 'Read configuration and data files...'
! Read workspace into config file found in exe folder
open(unit=1, file=trim(exepath)//trim(ConfigFile),status="OLD",iostat=err)
if(err/=0) call FatalExit('HBay:FATAL:Problem opening config file')
read(1,'(A250)', iostat=err) workspace
if(err/=0) call FatalExit('HBay:FATAL:Problem reading config file')
close(1)

Cfolder=trim(workspace)//trim(ConfigFolder)

! Read Config files
call CRead_data_OLD(file=trim(Cfolder)//trim(Config_Data_File),&
                DataFile=DataFile,DataType=DataType,&
                mvOption=mvOption,O1_value=O1_value,&
                O2_cond=O2_cond,O2_value=O2_value,&
                err=err)
if(err/=0) call FatalExit('HBay:FATAL:Problem reading Config_Data file')

call CRead_Inference(file=trim(Cfolder)//trim(Config_Inference_File),&
                dist=dist,parnames=parnames,priors=Hb%DistPrior,err=err)
if(err/=0) call FatalExit('HBay:FATAL:Problem reading Config_Inference file')
npar=size(parnames)


call CRead_MCMC(trim(Cfolder)//trim(Config_MCMC_File),&
                Nadapt,Ncycles,BurnFactor,Nslim,Nrep,MinMoveRate,&
                MaxMoveRate,DownMult,UpMult,MultFactor,Eps,err)
if(err/=0) call FatalExit('HBay:FATAL:Problem reading Config_MCMC file')

call CRead_ResultFiles_HBay(file=trim(Cfolder)//trim(Config_ResultFiles_File),&
                mc2_File=mc2_File,mc2_Save=mc2_Save,&
                par_File=par_File,par_Save=par_Save,&
                pdf_File=pdf_File,pdf_Save=pdf_Save,&
                cdf_File=cdf_File,cdf_Save=cdf_Save,&
                prd_File=prd_File,prd_Save=prd_Save,&
                qtl_File=qtl_File,qtl_Save=qtl_Save,&
                err=err)
if(err/=0) call FatalExit('HBay:FATAL:Problem reading Config_ResultFiles file')

call CRead_ResultOptions_OLD(file=trim(Cfolder)//trim(Config_ResultOptions_File),&
                       xmin=xmin,xmax=xmax,xn=xn,& ! options for pdf/cdf
                       pp=pp,kernel=kernel,& ! options for empirical stats
                       pmin=pmin,pmax=pmax,pn=pn,level=level,PtoT=PtoT,& ! options for quantiles
                       err=err)
if(err/=0) call FatalExit('HBay:FATAL:Problem reading Config_ResultOptions file')

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!            STEP 2: Read Data and handle missing values               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! first read to get nrow
open(unit=1,file=trim(DataFile),status='old')
read(1,*) line !headers
nrow=0;errcode=0
do while(errcode==0) !just count number of days
    nrow=nrow+1
    read(1,*,iostat=errcode) line
    if(trim(line)=='') errcode=1
enddo
nrow=nrow-1
close(1)
! Read histo data
call ReadHistoricData(file=trim(DataFile),& ! path to file
                  nskip=1,& ! number of headers row to be skipped
                  nrow=nrow,ncol=9,& ! number of row/columns (excluding headers)
                  IndxCol=1,& ! column containing the indexing variable (eg years)
                  DataTypeCol=2,& ! column containing the data type (0 to 3)
                  IsEqualToCol=3,& ! column containing value for type-0 data
                  IsSmallerThanCol=4,& ! column containing value for type-1 data
                  IsLargerThanCol=5,& ! column containing value for type-2 data
                  IsBetweenCol=(/6,7/),& ! 2 columns containing value for type-3 data
                  QDcol=8,& ! column containing QD values (only works if corresponding Qpeak is of type 0)
                  SystErrCol=9,& ! column containing the systematic error number (0 for no error)
                  mv=O1_value,& ! value denoting missing data
                  nickname='Histo',& ! name for this dataset
                  X=Hb%X,& ! data stored in an object HistoSampleType
                  err=err,mess=mess)
if(err/=0) call FatalExit('HBay:FATAL:Problem reading Histo-Data')
Hb%nickname='Histo'
Hb%distId=dist
Hb%npar=npar
do i=1,2
    Hb%DToPeakPrior(i)%dist='Uniform'
    allocate(Hb%DToPeakPrior(i)%par(2))
    Hb%DToPeakPrior(i)%par=(/0._mrk,1._mrk/)
enddo
if(associated(y)) nullify(y);allocate(y(nrow,2));
y(:,1)=Hb%X%Hsample(:)%Indx
y(:,2)=Hb%X%Hsample(:)%IsEqualTo
call RemoveMV(yin=y,mvOption='Option1',O1_value=O1_value,&
              O2_cond=O2_cond,O2_value=O2_value,&
              yout=x,err=err)
if(err/=0) call FatalExit('HBay:FATAL:Problem removing MV from Data')
call CRead_SystematicErrors(file=trim(Cfolder)//trim(Config_SystematicErrors_File),&
                npar=Hb%X%nSystErr,parnames=SystErrNames,priors=Hb%SystErrPrior,err=err)
if(err/=0) call FatalExit('HBay:FATAL:Problem reading Config_SystematicErrors file')
write(*,*) 'Reading OK!'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                        STEP 3: Run MCMC sampling                     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compute starting points
Ninfer=npar+2+Hb%X%nSystErr
allocate(teta0(Ninfer))
allocate(teta_std0(Ninfer))
! initial values for D-pars
call GetRoughEstimate(DistID=trim(dist), X=x(:,obscol), par=teta0(1:npar), err=err, mess=mess)
! Initial values for D2peak-par
teta0( (npar+1):(npar+2) )=(/0.6_mrk, 0.2_mrk/)
! Initial values for SystErr-par
if(Hb%X%nSystErr>0) teta0( (npar+3):)=1._mrk
where(teta0/=0._mrk)
  teta_std0=MultFactor*abs(teta0)
elsewhere
  teta_std0=eps
end where
write(*,*) '--------------------------------------------'
write(*,*) 'Start MCMC sampling (might take a while...)'
call HistoricBay_EstimPar(Hbay=Hb,& ! Object "HistoricBayType", containing data, model, priors, etc...
                ! Tuning of the MCMC sampler
                start=teta0,startStd=teta_std0,&
                nAdapt=nAdapt,nCycles=nCycles,&
                MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
                DownMult=DownMult,UpMult=UpMult,&
                OutFile=trim(exepath)//trim(MCMC_tempFile), &
                ! error handling
                err=err,mess=mess)
if(err/=0) call FatalExit('HBay:FATAL:Problem during MCMC sampling')
write(*,*) 'MCMC sampling OK!'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                    STEP 4: Post-process MCMC samples                 !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write(*,*) '--------------------------------------------'
write(*,*) 'Post-process MCMC samples...'


! BAY_MC2: Load MCMC samples
call LoadMCMCsample(File=trim(exepath)//trim(MCMC_tempFile),&
                    BurnFactor=BurnFactor,Nslim=Nslim,distID=dist,npar=Ninfer,&
                                    mcmc=mcmc,LogPost=LogPost,err=err,mess=mess)
if(err/=0) call FatalExit('HBay:FATAL:Problem during loading MCMC samples')

! BAY_MC2: write burned & slimmed MCMC samples
if(allocated(MCMChead)) deallocate(MCMChead)
allocate(MCMChead(Ninfer+1))
MCMChead(1:npar)=parnames
MCMChead((npar+1):(npar+2))=(/'Peak2Day_mean','Peak2Day_sdev'/)
if(Hb%X%nSystErr>0) MCMChead((npar+3):Ninfer)=SystErrNames
MCMChead(Ninfer+1)='LogPosterior'

if(MC2_save) then
  if(allocated(SlimMCMC)) deallocate(SlimMCMC)
  allocate(SlimMCMC(size(mcmc,dim=1),size(mcmc,dim=2)+1))
  SlimMCMC(:,1:size(mcmc,dim=2))=mcmc
  SlimMCMC(:,1+size(mcmc,dim=2))=LogPost
  call DatWrite(y=SlimMCMC,file=trim(workspace)//trim(mc2_File),&
                headers=MCMChead,err=err,mess=mess)
  if(err/=0) call FatalExit('HBay:FATAL:Problem writing MCMC samples')
endif

! BAY_PAR: get MCMC summary
if(par_Save) then
    call GetMCMCSummary(mcmc=mcmc,LogPost=LogPost,&
                        OutFile=trim(workspace)//trim(par_File),&
                        parnames=MCMChead(1:Ninfer),&
                        err=err,mess=mess)
    if(err/=0) call FatalExit('HBay:FATAL:Problem summarizing MCMC samples')
endif

! BAY_PRD
if(prd_Save .or. pdf_Save .or. cdf_Save .or. qtl_Save) then
    call GeneratePred(mcmc=mcmc(:,1:npar),distID=dist,nrep=nrep,pred=pred,&
                      err=err,mess=mess)
    if(err/=0) call FatalExit('HBay:FATAL:Problem Generating values from predictive')
endif

! BAY_PDF, BAY_CDF and BAY_QTL
if(pdf_Save .or. cdf_Save .or. qtl_Save) then
    if(allocated(mode)) deallocate(mode)
    allocate(mode(npar))

    if(allocated(xlist)) deallocate(xlist)
    allocate(xlist(nint(xn)))
    call GetMode(mcmc=mcmc(:,1:npar),LogPost=LogPost,mode=mode,err=err,mess=mess)
    if(err/=0) call FatalExit('HBay:FATAL:Problem Getting modal estimate')
    call GetQuantile(DistId=dist,p=xmin,par=mode,q=xmini,feas=feas,err=err,mess=mess)
    if(err/=0) call FatalExit('HBay:FATAL:Problem computing xlist')
    call GetQuantile(DistId=dist,p=xmax,par=mode,q=xmaxi,feas=feas,err=err,mess=mess)
    if(err/=0) call FatalExit('HBay:FATAL:Problem computing xlist')
    step=(xmaxi-xmini)/real(xn-1,mrk)
    forall(i=1:nint(xn)) xlist(i)=xmini+(i-1)*step

    if(allocated(plist)) deallocate(plist)
    allocate(plist(nint(pn)))
    step=(10._mrk**pmax-10._mrk**pmin)/real(pn-1,mrk)
    forall(i=1:nint(pn)) plist(i)=log10(10._mrk**pmin+(i-1)*step)

    if(allocated(pdf_modal)) deallocate(pdf_modal)
    allocate(pdf_modal(nint(xn)))
    if(allocated(pdf_pred)) deallocate(pdf_pred)
    allocate(pdf_pred(nint(xn)))
    if(allocated(cdf_modal)) deallocate(cdf_modal)
    allocate(cdf_modal(nint(xn)))
    if(allocated(cdf_pred)) deallocate(cdf_pred)
    allocate(cdf_pred(nint(xn)))
    if(allocated(Qp_modal)) deallocate(Qp_modal)
    allocate(Qp_modal(nint(pn)))
    if(allocated(Qp_pred)) deallocate(Qp_pred)
    allocate(Qp_pred(nint(pn)))
    if(allocated(IC_modal)) deallocate(IC_modal)
    allocate(IC_modal(nint(pn),2))

    call GetModalEstimates(mcmc=mcmc(:,1:npar),LogPost=LogPost,distID=dist,&
                          mode=mode,&
                          x=xlist,pdf=pdf_modal,cdf=cdf_modal,&
                          p=plist,Qp=Qp_modal,IC=IC_modal,&
                          level=level,err=err,mess=mess)
    if(err/=0) call FatalExit('HBay:FATAL:Problem computing modal estimates')

    call GetPredictiveEstimates(pred=pred,&
                          x=xlist,pdf=pdf_pred,cdf=cdf_pred,&
                          p=plist,Qp=Qp_pred,&
                          ppFormula=pp, &! optional, formula for plotting position [default Hazen]
                          kernelID=kernel,& !optional, parameters of kde [defaults Gaussian + built_in h value]
                          err=err,mess=mess)
    if(err/=0) call FatalExit('HBay:FATAL:Problem computing modal estimates')
    ! BAY_CDF
    if(cdf_Save) then
        if(allocated(matrice)) deallocate(matrice)
        allocate(matrice(nint(xn),3))
        matrice(:,1)=xlist
        matrice(:,2)=cdf_modal
        matrice(:,3)=cdf_pred
        call DatWrite(y=matrice,file=trim(workspace)//trim(cdf_File),&
                headers=(/'x         ','modal_F(x)','pred_F(x) '/),err=err,mess=mess)
        if(err/=0) call FatalExit('HBay:FATAL:Problem writing bay_cdf file')
    endif
    ! BAY_PDF
    if(pdf_Save) then
        if(allocated(matrice)) deallocate(matrice)
        allocate(matrice(nint(xn),3))
        matrice(:,1)=xlist
        matrice(:,2)=pdf_modal
        matrice(:,3)=pdf_pred
        call DatWrite(y=matrice,file=trim(workspace)//trim(pdf_File),&
                headers=(/'x         ','modal_f(x)','pred_f(x) '/),err=err,mess=mess)
        if(err/=0) call FatalExit('HBay:FATAL:Problem writing bay_pdf file')
    endif
    ! BAY_QTL
    if(qtl_Save) then
        if(allocated(matrice)) deallocate(matrice)
        allocate(matrice(nint(pn),6))
        matrice(:,1)=plist
        matrice(:,2)=p2T(plist,pToT)
        matrice(:,3)=Qp_modal
        matrice(:,4)=IC_modal(:,1)
        matrice(:,5)=IC_modal(:,2)
        matrice(:,6)=Qp_pred
        call DatWrite(y=matrice,file=trim(workspace)//trim(qtl_File),&
                headers=(/'p         ','T         ','modal_Q(p)','IC_lower  ','IC_upper  ','pred_Q(p) '/),&
                err=err,mess=mess)
        if(err/=0) call FatalExit('HBay:FATAL:Problem writing bay_pdf file')
    endif
    ! BAY_PRD
    if(prd_Save) then
        if(allocated(matrice)) deallocate(matrice)
        allocate(matrice(size(pred),1))
        matrice(:,1)=pred
        call DatWrite(y=matrice,file=trim(workspace)//trim(prd_File),&
                headers=(/'Predictive replicates'/),err=err,mess=mess)
        if(err/=0) call FatalExit('HBay:FATAL:Problem writing bay_prd file')
    endif
endif
write(*,*) 'Post-processing OK!'

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write(*,*) ' '
write(*,*) ' '
write(*,*) 'All done!!!';read(*,*)
end program Hbay

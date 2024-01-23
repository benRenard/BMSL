program Estimator

use kinds_dmsl_kit
use BayesianEstimation_tools
use DataRW_tools, only:DatRead,DatWrite,ReadSeparatedFile,WriteSeparatedFile
Use Distribution_tools, only:GetRoughEstimate,GetQuantile,GetParNumber,GetCdf
use EmpiricalStats_tools, only:GetEmpiricalDistSummary,p2T,T2p,GetEmpiricalQuantile
use Estimator_tools
implicit none

!character(250), parameter::exepath="..\..\Core\Estimator\" ! CORE version
!character(250), parameter::exepath="..\..\..\Core\Estimator\"! DEBUG/RELEASE
character(250), parameter::exepath="" ! CORE version
character(250), parameter::ConfigFile="Estimator_Config.txt"
character(250), parameter::ConfigFolder="" !"Config\"
character(250), parameter::Config_Data_File="Config_Data.txt"
character(250), parameter::Config_Inference_File="Config_Inference.txt"
character(250), parameter::Config_MCMC_File="Config_MCMC.txt"
character(250), parameter::Config_ResultFiles_File="Config_ResultFiles.txt"
character(250), parameter::Config_ResultOptions_File="Config_ResultOptions.txt"
character(250), parameter::Config_GetQT_File="Config_GetQT.txt"
character(250), parameter::MCMC_tempFile="MCMC.txt"
character(1),parameter::sep=";"
integer(mik),parameter::ncol_DATfile=3,obscol=2
integer(mik),parameter::err_ok=0,err_openConfig=1,err_readConfig=2,err_readConfigData=3,&
                        err_readConfigInference=4,err_readConfigMCMC=5,err_readConfigResultFiles=6,&
                        err_readConfigResultOptions=7,err_readData=8,err_MCMC=9,&
                        err_loadMCMC=10,err_writeMCMC=11,err_summarizeMCMC=12,&
                        err_emp=13,err_pred=14,err_modal=15,err_xlist=16,err_writeCdf=17,&
                        err_writePdf=18,err_writeQtl=19,err_unknownDataFormat=20,err_writeEmp=21,&
                        err_writePrd=22,err_writePar=23,err_readConfigGetQT=24,err_GetQT=25,&
                        err_misc=666

character(250):: workspace,Cfolder,mess,DataFile,DataType,&
                 mvOption,O2_cond,dist,headers(ncol_DATfile)
character(1)::QorT
integer(mik)::err,npar,i
real(mrk)::O1_value,O2_value,xmini,xmaxi,val,prob,maxpost,lo,hi
type(PriorListType), pointer:: priors(:)
character(250), pointer:: parnames(:)
character(250), allocatable:: MCMChead(:)
real(mrk), pointer:: y(:,:),x(:,:),mcmc(:,:),LogPost(:),pred(:)
real(mrk), allocatable::teta0(:),teta_std0(:),mode(:),xlist(:),plist(:),&
    pdf_modal(:),pdf_pred(:),cdf_modal(:),cdf_pred(:),Qp_modal(:),Qp_pred(:),&
    IC_modal(:,:),matrice(:,:),SlimMCMC(:,:),par_modal(:),par_mean(:),par_median(:),&
    par_std(:),par_IC(:,:),par_cor(:,:)
logical::feas,invertT,ok
!!!!!!!!!!!!!!!!!!!
! MCMC settings
integer(mik):: Nadapt,Ncycles,Nslim,Nrep,N
real(mrk)::BurnFactor,MinMoveRate,MaxMoveRate,DownMult,UpMult,MultFactor,Eps
! END MCMC settings
!!!!!!!!!!!!!!!!!!!
! ResultFiles
character(250)::mc2_File,par_File,dat_File,pdf_File,&
                cdf_File,prd_File,qtl_File,emp_File,QT_File
logical::mc2_Save,par_Save,dat_Save,pdf_Save,cdf_Save,&
         prd_Save,qtl_Save,emp_Save,QTonly
! END ResultFiles
!!!!!!!!!!!!!!!!!!!
! ResultOptions
real(mrk)::xmin,xmax,xn,pmin,pmax,pn,level,PtoT,step 
character(250)::pp,kernel
! END ResultOptions
!!!!!!!!!!!!!!!!!!!
real(mrk), allocatable::sortedX(:),Epdf(:),Ecdf(:),Et(:)

mc2_Save=.true.
par_Save=.true.
dat_Save=.true.
pdf_Save=.true.
cdf_Save=.true.
prd_Save=.true.
qtl_Save=.true.
emp_Save=.true.
                          
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                      STEP 1: READ CONFIG FILES                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Read workspace into config file found in exe folder
open(unit=1, file=trim(exepath)//trim(ConfigFile),status="OLD",iostat=err)
if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem opening config file';call exit(err_openConfig);endif
read(1,'(A250)', iostat=err) workspace
if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem reading config file';call exit(err_readConfig);endif
close(1)

Cfolder=trim(workspace)//trim(ConfigFolder)

inquire(file=trim(Cfolder)//trim(Config_GetQT_File),exist=QTonly)

if(.not.QTonly) then ! Full inference
    ! Read Config files
    call CRead_data(file=trim(Cfolder)//trim(Config_Data_File),&
                    DataFile=DataFile,DataType=DataType,&
        !                mvOption=mvOption,O1_value=O1_value,&
        !                O2_cond=O2_cond,O2_value=O2_value,&
                    err=err)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem reading Config_Data file';call exit(err_readConfigData);endif

    call CRead_Inference(file=trim(Cfolder)//trim(Config_Inference_File),&
                    dist=dist,parnames=parnames,priors=priors,err=err)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem reading Config_Inference file';call exit(err_readConfigInference);endif
    npar=size(parnames)

    call CRead_MCMC(trim(Cfolder)//trim(Config_MCMC_File),&
                    Nadapt,Ncycles,BurnFactor,Nslim,Nrep,MinMoveRate,&
                    MaxMoveRate,DownMult,UpMult,MultFactor,Eps,err)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem reading Config_MCMC file';call exit(err_readConfigMCMC);endif

    call CRead_ResultFiles(file=trim(Cfolder)//trim(Config_ResultFiles_File),&
                    mc2_File=mc2_File,& !mc2_Save=mc2_Save,&
                    par_File=par_File,& !par_Save=par_Save,&
                    dat_File=dat_File,& !dat_Save=dat_Save,&
                    pdf_File=pdf_File,& !pdf_Save=pdf_Save,&
                    cdf_File=cdf_File,& !cdf_Save=cdf_Save,&
                    prd_File=prd_File,& !prd_Save=prd_Save,&
                    qtl_File=qtl_File,& !qtl_Save=qtl_Save,&
                    emp_File=emp_File,& !emp_Save=emp_Save,&
                    err=err)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem reading Config_ResultFiles file';call exit(err_readConfigResultFiles);endif

    call CRead_ResultOptions(file=trim(Cfolder)//trim(Config_ResultOptions_File),&
                           xmin=xmin,xmax=xmax,xn=xn,& ! options for pdf/cdf
                           pp=pp,kernel=kernel,& ! options for empirical stats
                           pmin=pmin,pmax=pmax,pn=pn,level=level,PtoT=PtoT,invertT=invertT,& ! options for quantiles
                           err=err)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem reading Config_ResultOptions file';call exit(err_readConfigResultOptions);endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !            STEP 2: Read Data and handle missing values               !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    select case(DataType)
    case("bay_var")
        call ReadSeparatedFile(file=trim(DataFile),sep=sep,nhead=1,ncol=3,y=y,err=err,mess=mess)
        !call DatRead(file=DataFile,ncol=ncol_DATfile,y=y,&
        !          headers=headers,err=err,mess=mess)
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem reading Data file';call exit(err_readData);endif
    case default
        write(*,*) 'Estimator:FATAL:Unknown data format';call exit(err_unknownDataFormat);
    end select
    if(associated(x)) nullify(x);allocate(x(size(y,1),size(y,2)))
    x=y
    !call RemoveMV(yin=y,mvOption=mvOption,O1_value=O1_value,&
    !              O2_cond=O2_cond,O2_value=O2_value,&
    !              yout=x,err=err)
    !if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem removing MV from Data';call exit(err_removeMV);endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !                        STEP 3: Run MCMC sampling                     !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Compute starting points
    allocate(teta0(npar))
    allocate(teta_std0(npar))
    call GetRoughEstimate(DistID=trim(dist), X=x(:,obscol), par=teta0, err=err, mess=mess)
    where(teta0/=0._mrk)
      teta_std0=MultFactor*abs(teta0)
    elsewhere
      teta_std0=eps
    end where
    if(allocated(MCMChead)) deallocate(MCMChead)
    allocate(MCMChead(npar+1))
    MCMChead(1:npar)=parnames
    MCMChead(npar+1)='LogPosterior'
    call Bayes_EstimPar(X=x(:,obscol),distID=trim(dist),mv=O1_value,&
                        PriorList=priors,& ! Data & Model
				        ! Tuning of the MCMC sampler
				        start=teta0,startStd=teta_std0,&
				        nAdapt=nAdapt,nCycles=nCycles,&
				        MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				        DownMult=DownMult,UpMult=UpMult,&
				        OutFile=trim(exepath)//trim(MCMC_tempFile), &
				        ! error handling
				        headers=MCMChead,err=err,mess=mess)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem during MCMC sampling';call exit(err_MCMC);endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !                    STEP 4: Post-process MCMC samples                 !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !! BAY_DAT: Re-write data without MV
    !if(dat_Save) then
    !    call DatWrite(y=x,file=trim(workspace)//trim(dat_File),&
    !                    headers=headers,err=err,mess=mess)
    !    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem Re-writing data file';call exit(err_writeData);endif
    !endif

    ! BAY_MC2: Load MCMC samples
    call LoadMCMCsample(File=trim(exepath)//trim(MCMC_tempFile),&
                        BurnFactor=BurnFactor,Nslim=Nslim,distID=dist,&
								        mcmc=mcmc,LogPost=LogPost,err=err,mess=mess)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem during loading MCMC samples';call exit(err_loadMCMC);endif

    ! BAY_MC2: write burned & slimmed MCMC samples
    if(MC2_save) then
      if(allocated(SlimMCMC)) deallocate(SlimMCMC)
      allocate(SlimMCMC(size(mcmc,dim=1),size(mcmc,dim=2)+1))
      SlimMCMC(:,1:size(mcmc,dim=2))=mcmc
      SlimMCMC(:,1+size(mcmc,dim=2))=LogPost
      call WriteSeparatedFile(file=trim(workspace)//trim(mc2_File),sep=sep,y=SlimMCMC,headers=MCMChead,err=err,mess=mess)
      !call DatWrite(y=SlimMCMC,file=trim(workspace)//trim(mc2_File),&
      !              headers=MCMChead,err=err,mess=mess)
      if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem writing MCMC samples';call exit(err_writeMCMC);endif
    endif

    ! BAY_PAR: get MCMC summary
    if(par_Save) then
        !call GetMCMCSummary(mcmc=mcmc,LogPost=LogPost,&
        !                    OutFile=trim(workspace)//trim(par_File),&
        !                    parnames=parnames,&
        !                    err=err,mess=mess)
        allocate(par_modal(npar),par_median(npar),par_mean(npar),par_std(npar))
        allocate(par_IC(npar,2),par_cor(npar,npar))
        call GetMCMCSummary(mcmc=mcmc,LogPost=LogPost,&
                            parnames=parnames,&
                            modal=par_modal,median=par_median,mean=par_mean,&
                            sdev=par_std,IC90=par_IC,correl=par_cor,&
                            err=err,mess=mess)
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem summarizing MCMC samples';call exit(err_summarizeMCMC);endif
        if(allocated(matrice)) deallocate(matrice);allocate(matrice(npar,6))
        matrice(:,1)=par_modal
        matrice(:,2)=par_median
        matrice(:,3)=par_mean
        matrice(:,4)=par_std
        matrice(:,5:6)=par_IC
        call WriteSeparatedFile(file=trim(workspace)//trim(par_File),sep=sep,y=matrice,&
                            headers=(/'mode','median','mean','stdev','q05','q95'/),err=err,mess=mess)
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem writing bay_par file';call exit(err_writePar);endif
        call WriteSeparatedFile(file=trim(workspace)//trim(par_File),sep=sep,y=par_cor,&
                            headers=(/'correlation'/),append=.true.,err=err,mess=mess)
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem writing bay_par file';call exit(err_writePar);endif
    endif

    ! BAY_EMP: empirical cdf/pdf from data
    if(emp_Save) then
        allocate(sortedX(size(x(:,obscol))))
        allocate(Epdf(size(x(:,obscol))))
        allocate(Ecdf(size(x(:,obscol))))
        allocate(ET(size(x(:,obscol))))
        call GetEmpiricalDistSummary(X=x(:,obscol),& ! data
                                ppFormula=pp, &! optional,, formula for plotting position [default Hazen]
                                kernelID=kernel,& !optional, parameters of kde [defaults Gaussian + built_in h value]
                                p2Tfactor=PtoT,&! optional, factor for converting probabilities to return periods, default 1
                                invertT=invertT,&
                                !Xname=headers(2),& ! optional, name given to X, will be written in header of output file
                                !Outfile=trim(workspace)//trim(emp_File),&! optional, file to write results
                                SortedX=SortedX,Epdf=Epdf,Ecdf=Ecdf,ET=ET,&
                                err=err,mess=mess)
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem Computing Empirical estimates';call exit(err_Emp);endif
        if(allocated(matrice)) deallocate(matrice);allocate(matrice(size(ET),4))
            matrice(:,1)=SortedX
            matrice(:,2)=Epdf
            matrice(:,3)=Ecdf
            matrice(:,4)=ET
            call WriteSeparatedFile(file=trim(workspace)//trim(emp_File),sep=sep,y=matrice,&
                                headers=(/'x','Epdf(x)','Ecdf(x)','ET(x)'/),err=err,mess=mess)
            if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem writing bay_emp file';call exit(err_writeEmp);endif
    endif

    ! BAY_PRD
    if(prd_Save .or. pdf_Save .or. cdf_Save .or. qtl_Save) then
        call GeneratePred(mcmc=mcmc,distID=dist,nrep=nrep,pred=pred,&
                          err=err,mess=mess)   
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem Generating values from predictive';call exit(err_pred);endif
        if(allocated(matrice)) deallocate(matrice);allocate(matrice(size(pred),1))
            matrice(:,1)=pred
            call WriteSeparatedFile(file=trim(workspace)//trim(prd_File),sep=sep,y=matrice,&
                                headers=(/'replicates'/),err=err,mess=mess)
            if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem writing bay_prd file';call exit(err_writePrd);endif
    endif

    ! BAY_PDF, BAY_CDF and BAY_QTL
    if(pdf_Save .or. cdf_Save .or. qtl_Save) then
        if(allocated(mode)) deallocate(mode)
        allocate(mode(npar))
        
        if(allocated(xlist)) deallocate(xlist)
        allocate(xlist(nint(xn)))
        call GetMode(mcmc=mcmc,LogPost=LogPost,mode=mode,err=err,mess=mess)
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem Getting modal estimate';call exit(err_modal);endif
        call GetQuantile(DistId=dist,p=xmin,par=mode,q=xmini,feas=feas,err=err,mess=mess)
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem computing xlist';call exit(err_xlist);endif
        call GetQuantile(DistId=dist,p=xmax,par=mode,q=xmaxi,feas=feas,err=err,mess=mess)
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem computing xlist';call exit(err_xlist);endif
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
       
        call GetModalEstimates(mcmc=mcmc,LogPost=LogPost,distID=dist,&
                              mode=mode,&
                              x=xlist,pdf=pdf_modal,cdf=cdf_modal,&
                              p=plist,Qp=Qp_modal,IC=IC_modal,&
                              level=level,err=err,mess=mess)    
        if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem computing modal estimates';call exit(err_modal);endif
        
        call GetPredictiveEstimates(pred=pred,&
                              x=xlist,pdf=pdf_pred,cdf=cdf_pred,&
                              p=plist,Qp=Qp_pred,&
                              ppFormula=pp, &! optional, formula for plotting position [default Hazen]
                              kernelID=trim(kernel),& !optional, parameters of kde [defaults Gaussian + built_in h value]
                              err=err,mess=mess)
        if(err>0) then;write(*,*) 'Estimator:FATAL:Problem computing predictive estimates';call exit(err_pred);endif
        ! BAY_CDF
        if(cdf_Save) then
            if(allocated(matrice)) deallocate(matrice)
            allocate(matrice(nint(xn),3))
            matrice(:,1)=xlist
            matrice(:,2)=cdf_modal
            matrice(:,3)=cdf_pred
            call WriteSeparatedFile(file=trim(workspace)//trim(cdf_File),sep=sep,y=matrice,&
                                headers=(/'x','modal_F(x)','pred_F(x)'/),err=err,mess=mess)
            !call DatWrite(y=matrice,file=trim(workspace)//trim(cdf_File),&
            !        headers=(/'x','modal_F(x)','pred_F(x)'/),err=err,mess=mess)
            if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem writing bay_cdf file';call exit(err_writeCdf);endif
        endif
        ! BAY_PDF
        if(pdf_Save) then
            if(allocated(matrice)) deallocate(matrice)
            allocate(matrice(nint(xn),3))
            matrice(:,1)=xlist
            matrice(:,2)=pdf_modal
            matrice(:,3)=pdf_pred
            call WriteSeparatedFile(file=trim(workspace)//trim(pdf_File),sep=sep,y=matrice,&
                                headers=(/'x','modal_f(x)','pred_f(x)'/),err=err,mess=mess)
            !call DatWrite(y=matrice,file=trim(workspace)//trim(pdf_File),&
            !        headers=(/'x','modal_f(x)','pred_f(x)'/),err=err,mess=mess)
            if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem writing bay_pdf file';call exit(err_writePdf);endif
        endif
        ! BAY_QTL
        if(qtl_Save) then
            if(allocated(matrice)) deallocate(matrice)
            allocate(matrice(nint(pn),6))
            matrice(:,1)=plist
            matrice(:,2)=p2T(plist,pToT,invertT)
            matrice(:,3)=Qp_modal
            matrice(:,4)=IC_modal(:,1)
            matrice(:,5)=IC_modal(:,2)
            matrice(:,6)=Qp_pred
            call WriteSeparatedFile(file=trim(workspace)//trim(qtl_File),sep=sep,y=matrice,&
                                headers=(/'p','T','modal_Q(p)','IC_lower','IC_upper','pred_Q(p)'/),err=err,mess=mess)
            !call DatWrite(y=matrice,file=trim(workspace)//trim(qtl_File),&
            !        headers=(/'p','T','modal_Q(p)','IC_lower','IC_upper','pred_Q(p)'/),&
            !        err=err,mess=mess)
            if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem writing bay_qtl file';call exit(err_writeQtl);endif
        endif
    endif
else ! Only compute Q or T
    ! Read Config file
    call CRead_GetQT(file=trim(Cfolder)//trim(Config_GetQT_File),&
                    dist=dist,QorT=QorT,val=val,&
                    PtoT=PtoT,invertT=invertT,level=level,&
                    mc2_File=mc2_File,QT_File=QT_File,err=err)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem reading Config_GetQT file';call exit(err_readConfigGetQT);endif
    ! get npar
    call GetParNumber(DistID=trim(dist),npar=npar,err=err,mess=mess)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:'//trim(mess);call exit(err_misc);endif
    ! Read MCMC file and get mode
    call LoadMCMCsample(File=trim(mc2_File),BurnFactor=0._mrk,Nslim=1,&
                        distID=trim(dist),npar=npar,sep=sep,&
                        mcmc=mcmc,LogPost=LogPost,err=err,mess=mess)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:Problem reading MCMC file';call exit(err_loadMCMC);endif
    N=size(mcmc,dim=1)
    if(allocated(xlist)) deallocate(xlist);allocate(xlist(N))
    if(allocated(mode)) deallocate(mode);allocate(mode(npar))
    call GetMode(mcmc,LogPost,mode,err,mess)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:'//trim(mess);call exit(err_misc);endif            
    ! Compute Q or T for mode and all MCMC parameters
    ok=.true.
    if(QorT=='Q') then
        prob=T2p(val,PtoT,invertT)
        if(abs(prob)/=undefRN) then
            call GetQuantile(DistId=dist,p=prob,par=mode,q=maxpost,feas=feas,err=err,mess=mess)
            if(err/=0 .or. (.not. feas)) then;write(*,*) 'Estimator:FATAL:'//trim(mess);call exit(err_misc);endif            
            do i=1,N
                call GetQuantile(DistId=dist,p=prob,par=mcmc(i,:),q=xlist(i),feas=feas,err=err,mess=mess)
                if(err/=0 .or. (.not. feas)) then;write(*,*) 'Estimator:FATAL:'//trim(mess);call exit(err_misc);endif            
            enddo  
        else
            ok=.false.
        endif      
    else if(QorT=='T') then
        call GetCdf(DistId=dist,x=val,par=mode,cdf=prob,feas=feas,err=err,mess=mess)
        if(err/=0 .or. (.not. feas)) then;write(*,*) 'Estimator:FATAL:'//trim(mess);call exit(err_misc);endif            
        if(prob>0._mrk .and. prob<1.0_mrk) then
            maxpost=p2T(prob,PtoT,invertT)
            do i=1,N
                call GetCdf(DistId=dist,x=val,par=mcmc(i,:),cdf=prob,feas=feas,err=err,mess=mess)
                if(err/=0 .or. (.not. feas)) then;write(*,*) 'Estimator:FATAL:'//trim(mess);call exit(err_misc);endif            
                xlist(i)=p2T(prob,PtoT,invertT)
            enddo
        else
            ok=.false.
        endif      
    else
        write(*,*) 'Estimator:FATAL: in Config_GetQT, invalid QorT value';call exit(err_misc)
    endif
    if(.not.ok) then;write(*,*) 'Estimator:FATAL:Problem computing '//QorT;call exit(err_GetQT);endif
    ! Compute uncertainty interval
    call GetEmpiricalQuantile(p=0.5_mrk*(1._mrk-level),x=xlist,q=lo,err=err,mess=mess)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:'//trim(mess);call exit(err_misc);endif            
    call GetEmpiricalQuantile(p=1._mrk-0.5_mrk*(1._mrk-level),x=xlist,q=hi,err=err,mess=mess)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:'//trim(mess);call exit(err_misc);endif  
   ! write result
   call WriteSeparatedFile(file=trim(workspace)//trim(QT_File),sep=sep,&
                           y=reshape((/maxpost,lo,hi/),(/1,3/)),&
                           headers=(/'maxpost','UC_low','UC_high'/),err=err,mess=mess)
    if(err/=0) then;write(*,*) 'Estimator:FATAL:'//trim(mess);call exit(err_misc);endif  
endif

write(*,*) 'Estimator ran successfully'
call exit(err_ok)

end program Estimator

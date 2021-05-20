module LocalLikelihood_tools

!~**********************************************************************
!~* Purpose: Local Likelihood estimation of Hydro-model 
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 12/02/2014
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1.
!~*		2.
!~*		3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: LoL_EstimPar !LoL_LoadHmodel,LoL_runHmodel,LoL_GetLogLikelihood_Baseline

type,public::HmodelType ! contains properties of the Hydrological model
    integer(mik)::modelID=undefIN
    integer(mik)::ninput=undefIN,nstate=undefIN,npar=undefIN
    character(250)::modelName='AintGotNoName',indxName='AintGotNoName'
    character(250),allocatable::inputName(:),stateName(:),parName(:)
    real(mrk),allocatable::stateLo(:),stateHi(:),parLo(:),parHi(:)
    real(mrk),allocatable::inScal(:),stateScal(:),parScal(:)
    real(mrk),allocatable::stateDef(:),parDef(:),parSD(:)
    integer(mik),allocatable::parTran(:)
    logical, allocatable::parFit(:)
end type HmodelType

type,public::LoLType ! contains properties of the local-likelihood setup
    type(HmodelType)::Hmodel
    real(mrk),allocatable::input(:,:),output(:) ! input/output time series. Single output is assumed for the moment
    real(mrk),allocatable::baseline(:) ! baseline parameter vector
    real(mrk),allocatable::BaselineSim(:) ! simulation corresponding to the baseline
    real(mrk),allocatable::BaselineResPar(:) ! parameters of the residual model corresponding to the baseline
    real(mrk),allocatable::covariate(:,:) ! covariates time series
    real(mrk)::smooth=undefRN ! smoothness parameter
    integer(mik)::targetpar=undefIN ! index of the parameter that will be bumped
    real(mrk)::mvcode=-99.99_mrk ! code used to denote mv in input/output time series
    logical::EstimBaseline=.true. ! if false, will instead estimate the bump
    real(mrk),allocatable::dist(:) ! Nt-sized vector of distances (in covariate space) between t and target time step
end type LoLType
type(LoLType)::LL

Contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoL_LoadHmodel(modelID,myDMSL_path,setupCmd,infoCmd,Hmodel,err,mess)

!^**********************************************************************
!^* Purpose: Load properties of the selected Hydro-model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:12/02/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Hmodel
!^*		2. myDMSL_path - Hmodel info is stored in a subfolder
!^*		3. [infoCmd]generally unused, see dynamicModelLibrary
!^*		4. [setupCmd]generally unused, see dynamicModelLibrary
!^* OUT
!^*		1.Hmodel, object storing model information
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use dynamicModelLibrary,only:DMDL_setModel,DMDL_getModelInfo
use intercom_dmsl_kit,only:setFilePaths_dmsl
integer(mik), intent(in)::modelID
character(*), intent(in)::myDMSL_path
character(*), intent(in), optional::setupCmd,infoCmd
type(HmodelType), intent(out)::Hmodel
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250), parameter::setupCmd_default='',infoCmd_default=''
character(250)::Scmd,Icmd

err=0;mess='';Hmodel%modelID=modelID
! handle optional arguments
if(present(setupCmd)) then;Scmd=setupCmd;else;Scmd=setupCmd_default;endif
if(present(infoCmd)) then;Icmd=infoCmd;else;Icmd=infoCmd_default;endif
! Set DMSL folder - needed to know where to read the PCO files
call setFilePaths_dmsl(myDMSL_path=myDMSL_path,err=err,message=mess)
if(err>0) then;mess="LoL_LoadHmodel: "//trim(mess);return;endif
! Choose the H-model
call DMDL_setModel(modelID=(/modelID/),setupCmd=Scmd,err=err,message=mess) 
if(err>0) then;mess="LoL_LoadHmodel: "//trim(mess);return;endif
! Get dimensions of input/states/par
call DMDL_getModelInfo(modelID=(/modelID/),infoCmd=Icmd,modelName=Hmodel%modelName,indxName=Hmodel%indxName,&
    ninput=Hmodel%ninput,nstate=Hmodel%nstate,npar=Hmodel%npar,err=err,message=mess)
if(err>0) then;mess="LoL_LoadHmodel: "//trim(mess);return;endif
if(allocated(Hmodel%inputName)) deallocate(Hmodel%inputName);allocate(Hmodel%inputName(Hmodel%ninput))
if(allocated(Hmodel%stateName)) deallocate(Hmodel%stateName);allocate(Hmodel%stateName(Hmodel%nstate))
if(allocated(Hmodel%parName)) deallocate(Hmodel%parName);allocate(Hmodel%parName(Hmodel%npar))
if(allocated(Hmodel%stateLo)) deallocate(Hmodel%stateLo);allocate(Hmodel%stateLo(Hmodel%nstate))
if(allocated(Hmodel%stateHi)) deallocate(Hmodel%stateHi);allocate(Hmodel%stateHi(Hmodel%nstate))
if(allocated(Hmodel%parLo)) deallocate(Hmodel%parLo);allocate(Hmodel%parLo(Hmodel%npar))
if(allocated(Hmodel%parHi)) deallocate(Hmodel%parHi);allocate(Hmodel%parHi(Hmodel%npar))
if(allocated(Hmodel%inScal)) deallocate(Hmodel%inScal);allocate(Hmodel%inScal(Hmodel%ninput))
if(allocated(Hmodel%stateScal)) deallocate(Hmodel%stateScal);allocate(Hmodel%stateScal(Hmodel%nstate))
if(allocated(Hmodel%parScal)) deallocate(Hmodel%parScal);allocate(Hmodel%parScal(Hmodel%npar))
if(allocated(Hmodel%stateDef)) deallocate(Hmodel%stateDef);allocate(Hmodel%stateDef(Hmodel%nstate))
if(allocated(Hmodel%parDef)) deallocate(Hmodel%parDef);allocate(Hmodel%parDef(Hmodel%npar))
if(allocated(Hmodel%parSD)) deallocate(Hmodel%parSD);allocate(Hmodel%parSD(Hmodel%npar))
if(allocated(Hmodel%parTran)) deallocate(Hmodel%parTran);allocate(Hmodel%parTran(Hmodel%npar))
if(allocated(Hmodel%parFit)) deallocate(Hmodel%parFit);allocate(Hmodel%parFit(Hmodel%npar))
! Get infos for input/states/par
call DMDL_getModelInfo(modelID=(/modelID/),infoCmd=Icmd,&
    inputName=Hmodel%inputName,stateName=Hmodel%stateName,parName=Hmodel%parName,&
    stateLo=Hmodel%stateLo,stateHi=Hmodel%stateHi,parLo=Hmodel%parLo,parHi=Hmodel%parHi,&
    inScal=Hmodel%inScal,stateScal=Hmodel%stateScal,parScal=Hmodel%parScal,&
    stateDef=Hmodel%stateDef,parDef=Hmodel%parDef,parSD=Hmodel%parSD,&
    parTranDef=Hmodel%parTran,parFitDef=Hmodel%parFit,&
    err=err,message=mess)
if(err>0) then;mess="LoL_LoadHmodel: "//trim(mess);return;endif

end subroutine LoL_LoadHmodel
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoL_runHmodel(Hmodel,input,par,output,feas,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.Hmodel
!^*		2.input, Nt*Ninput matrix
!^*		3.par,Nt*Npar matrix
!^* OUT
!^*		1.output, streamflow series
!^*		2.feas, is run feasible?
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************
use dynamicModelLibrary,only:DMDL_controlModel,DMDL_runModel 
type(HmodelType), intent(in)::Hmodel
real(mrk), intent(in)::input(:,:),par(:,:)
real(mrk), intent(out)::output(:)
logical,intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::t,Nt
real(mrk)::state(Hmodel%nstate)

err=0;mess='';feas=.true.
! size checks
Nt=size(input,dim=1)
if(size(par,dim=1)/=Nt .or. size(output)/=Nt) then
    err=1;mess='LoL_runHmodel:Fatal:Nt size mismatch [input,par,output]';return
endif
if(size(par,dim=2)/=Hmodel%npar .or. size(input,dim=2)/=Hmodel%ninput) then
    err=1;mess='LoL_runHmodel:Fatal:Hmodel size mismatch [input,par]';return
endif
! Initialisation
call DMDL_controlModel(modelID=(/Hmodel%modelID/),parIn=par(1,:),setS0in=.true.,feas=feas,err=err,message=mess)
if(err>0) then;mess="LoL_runHmodel: "//trim(mess);return;endif
if(.not.feas) return
do t=1,Nt
    call DMDL_controlModel(modelID=(/Hmodel%modelID/),parIn=par(t,:),feas=feas,err=err,message=mess)
    if(err>0) then;mess="LoL_runHmodel: "//trim(mess);return;endif
    if(.not.feas) return
    call DMDL_runModel(modelID=(/Hmodel%modelID/),&
            runitCmd='',iT=undefIN,dataProps=(/1._mrk/),& ! BOTCH: works for GR4J but may need to be changed, not sure what those variables are doing
            input=input(t,:),state=state,feas=feas,err=err,message=mess)
    if(err>0) then;mess="LoL_runHmodel: "//trim(mess);return;endif
    if(.not.feas) return
    output(t)=state(1) ! BOTCH: not sure what's the general rule to decide which state is streamflow - I assume it's the first one
enddo

end subroutine LoL_runHmodel
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoL_GetLogLikelihood_Baseline(teta,LoL,logp,feas,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments: the residual error model is hard-coded for the moment
!^* Using sd.res=a+bQsim
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.teta, Hmodel_npar+2 vector of parameters to be inferred
!^*		2.LoL: object storing all info for local likelihood implementation
!^* OUT
!^*		1.logp: log-lkh
!^*		2.feas, feasability
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************
use numerix_dmsl_kit, only:normal_logp
real(mrk), intent(in)::teta(:)
type(LoLtype), intent(in)::LoL
real(mrk), intent(out)::logp
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik),parameter::npar_ResModel=2
integer(mik)::npar_Hmodel,t,Nt,i,k
real(mrk)::par(size(LoL%input,dim=1),LoL%Hmodel%npar),dummypar(LoL%Hmodel%npar)
real(mrk)::a,b

err=0;mess=''
logp=undefRN
feas=.true.
Nt=size(LoL%input,dim=1)
! Check size
npar_Hmodel=count(LoL%Hmodel%ParFit)
if(size(teta)/=npar_Hmodel+npar_ResModel) then
    mess='LoL_GetLogLikelihood_Baseline:Fatal:wronf size [teta]';err=1;return
endif
! Run model
dummypar=LoL%Hmodel%ParDef;k=0
do i=1,LoL%Hmodel%npar
    if(LoL%Hmodel%ParFit(i)) then;k=k+1;dummypar(i)=teta(k);endif
enddo
forall(t=1:Nt) par(t,:)=dummypar
! compute lkh
a=teta(npar_Hmodel+1);b=teta(npar_Hmodel+2)
call lkh_engine(par,a,b,LoL,logp,feas,err,mess)
if(err>0) then;mess="LoL_GetLogLikelihood_Baseline: "//trim(mess);return;endif

end subroutine LoL_GetLogLikelihood_Baseline
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoL_GetLogLikelihood(teta,LoL,dist,logp,feas,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments: 1/ Hard-coded bump function for the moment
!^* bump(d)=(-delta/w)*d+delta if d<=w, 0 otherwise
!^* 2/ Bump applies to a single parameter for the moment
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.teta, Hmodel_npar+2 vector of parameters to be inferred
!^*		2.LoL: object storing all info for local likelihood implementation
!^*		3.dist, Nt-sized vector of distances (in covariate space) between t and target time step
!^* OUT
!^*		1.logp: log-lkh
!^*		2.feas, feasability
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************
use numerix_dmsl_kit, only:normal_logp
real(mrk), intent(in)::teta(:),dist(:)
type(LoLtype), intent(in)::LoL
real(mrk), intent(out)::logp
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik),parameter::npar_ResModel=2
integer(mik)::t,Nt
real(mrk)::par(size(LoL%input,dim=1),LoL%Hmodel%npar)
real(mrk)::bump,a,b

err=0;mess='';feas=.true.;logp=undefRN
Nt=size(LoL%input,dim=1)
! Check size
if(size(dist)/=Nt) then
    mess='LoL_GetLogLikelihood:Fatal:wrong size [dist]';err=1;return
endif
if(size(teta)/=1+npar_ResModel) then
    mess='LoL_GetLogLikelihood:Fatal:wrong size [teta]';err=1;return
endif
! Run model
do t=1,Nt
    par(t,:)=LoL%baseline
    if(dist(t)<=LoL%smooth) then
        bump=teta(1)-(teta(1)/LoL%smooth)*dist(t)
    else
        bump=0._mrk
    endif
    par(t,LoL%targetpar)=LoL%baseline(LoL%targetpar)+bump
enddo
! compute lkh
a=teta(1+1);b=teta(1+2)
call lkh_engine(par,a,b,LoL,logp,feas,err,mess)
if(err>0) then;mess="LoL_GetLogLikelihood: "//trim(mess);return;endif

end subroutine LoL_GetLogLikelihood
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine lkh_engine(par,a,b,LoL,logp,feas,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. par: Nt*Npar matrix
!^*		2.a,b: residual error model parameters
!^*		3.LoL: object storing all info for local likelihood implementation
!^* OUT
!^*		1.logp: log-lkh
!^*		2.feas, feasability
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************
use numerix_dmsl_kit,only:normal_logp
real(mrk),intent(in)::par(:,:),a,b
type(LoLtype), intent(in)::LoL
real(mrk),intent(out)::logp
logical,intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk)::Qsim(size(LoL%input,dim=1)),lkh(size(LoL%input,dim=1))
integer(mik)::t,Nt
real(mrk)::sd

err=0;mess='';feas=.true.;logp=undefRN
if(a<=0._mrk .or. b<0._mrk) then;feas=.false.;return;endif
Nt=size(LoL%input,dim=1)
call LoL_runHmodel(LoL%Hmodel,LoL%input,par,Qsim,feas,err,mess)
if(.not.feas) return
if(err>0) then;mess="lkh_engine: "//trim(mess);return;endif
! Get lkh
do t=1,Nt
    if(LoL%output(t)==LoL%mvcode) then; lkh(t)=0._mrk;cycle;endif
    sd=a+b*Qsim(t)
    if(sd<=0._mrk) then;feas=.false.;return;endif
    lkh(t)=normal_logp(LoL%output(t),Qsim(t),sd**2)
enddo
logp=sum(lkh)

end subroutine lkh_engine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoL_LogLikelihoodWrapper (dataIN, dataOUT, x, feas, fx, gradFx, hessFx, err, message)

!^**********************************************************************
!^* Purpose: Wrapper to Log-Likelihood sub, compatible with D's optimizer
!^**********************************************************************
!^* Programmer: Ben Renard
!^**********************************************************************
!^* Last modified: 12/02/2014
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. dataIN, unused
!^*		2. dataOUT, unused
!^*		3. x, -> teta
!^* OUT
!^*		1. feas, feasible?
!^*		2. fx, Log-lok(teta)
!^*		3. GradFx & hessFx: unused
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type

TYPE(data_ricz_type),INTENT(in),OPTIONAL::dataIN
TYPE(data_ricz_type),INTENT(inout),OPTIONAL::dataOUT
REAL(mrk), INTENT(in) :: x(:)
LOGICAL, INTENT(out) :: feas
REAL(mrk), INTENT(out), OPTIONAL :: fx, gradFx(:), hessFx(:,:)
INTEGER(mik),INTENT(out)::err
CHARACTER(*),INTENT(out)::message
! locals
real(mrk)::logp

if(LL%EstimBaseline) then
    call LoL_GetLogLikelihood_Baseline(teta=x,LoL=LL,logp=logp,feas=feas,err=err,mess=message)
else
    call LoL_GetLogLikelihood(teta=x,LoL=LL,dist=LL%dist,logp=logp,feas=feas,err=err,mess=message)
endif
if(present(fx)) fx=-1*logp
end SUBROUTINE LoL_LogLikelihoodWrapper
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoL_EstimPar(modelID,myDMSL_path,setupCmd,infoCmd,& ! properties of the chosen hydro-model
                             ParFit,ParVal,& ! optional properties - if missing will be set to the default as defined in PCO files
                             input,output,& ! input/output time series
                             TargetPar,covariate,smooth,slim,& ! properties of the LoL case study
                             EstimatedBumps,EstimatedPar,BaselinePar,&!results-parameters
                             EstimatedSim,BaselineSim,& !results-simulated reponses
                             err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use Optimisation_tools
integer(mik), intent(in)::modelID
character(*), intent(in)::myDMSL_path
character(*), intent(in), optional::setupCmd,infoCmd
real(mrk),intent(in),optional::ParVal(:)
logical,intent(in),optional::ParFit(:)
real(mrk),intent(in)::input(:,:),output(:),covariate(:,:),smooth
integer(mik),intent(in)::TargetPar
integer(mik), intent(in),optional::slim
real(mrk),intent(out)::EstimatedBumps(size(output)),EstimatedPar(size(output)),&
    BaselinePar(:),EstimatedSim(size(output)),BaselineSim(size(output))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
! BOTCH: make the following optional input arguments of the procedure
integer(mik), parameter::smeth=qndmsl_smeth,mometh=qndmsl_smeth
integer(mik),parameter::nopt=10,itermax=10000
character(250),parameter::msFile="Optim_MultiStart_Results.txt",hisFile="Optim_History.txt",resFile="Optim_Results.txt"
real(mrk),parameter::NAval=UndefRN
! END BOTCH
integer(mik),parameter::npar_ResModel=2,slim_default=1
integer(mik)::t,Nt,Ndim,npar,Ninfer,i,step
integer(mik)::fcalls,gcalls,hcalls
real(mrk),pointer::xLo(:),xHi(:),msLo(:),msHi(:),xScale(:),x0(:),parmatrix(:,:)
real(mrk),allocatable::tetaHat(:),Hessian(:,:),bumps(:)
integer(mik),allocatable::activeSet(:)
real(mrk)::fmin,EstimParOut(size(output),size(BaselinePar))
logical::feas

err=0;mess=''
! handle optional arguments
if(present(slim)) then;step=slim;else;step=slim_default;endif
! Load hydrologic model
call LoL_LoadHmodel(modelID,myDMSL_path,setupCmd,infoCmd,LL%Hmodel,err,mess)
if(err>0) then;mess="LoL_EstimPar: "//trim(mess);return;endif
! Start populating LL object
Nt=size(input,dim=1)
if(allocated(LL%input)) deallocate(LL%input);allocate(LL%input(Nt,size(input,dim=2)))
if(allocated(LL%output)) deallocate(LL%output);allocate(LL%output(size(output)))
if(allocated(LL%baseline)) deallocate(LL%baseline);allocate(LL%baseline(LL%Hmodel%npar))
if(allocated(LL%dist)) deallocate(LL%dist);allocate(LL%dist(Nt))
if(allocated(LL%BaselineSim)) deallocate(LL%BaselineSim);allocate(LL%BaselineSim(Nt))
if(allocated(LL%covariate)) deallocate(LL%covariate);allocate(LL%covariate(Nt,size(covariate,dim=2)))
if(allocated(LL%BaselineResPar)) deallocate(LL%BaselineResPar);allocate(LL%BaselineResPar(npar_ResModel))
LL%input=input;LL%output=output;LL%smooth=smooth;LL%targetPar=TargetPar;LL%covariate=covariate
if(present(ParFit)) then
    if(size(ParFit)/=LL%Hmodel%npar) then;err=1;mess="LoL_EstimPar:Fatal:wrong size [ParFit]";return;endif
    LL%HModel%Parfit=parFit
endif
if(present(ParVal)) then
    if(size(ParVal)/=LL%Hmodel%npar) then;err=1;mess="LoL_EstimPar:Fatal:wrong size [ParVal]";return;endif
    LL%HModel%ParDef=parVal
endif
! Start by estimating the baseline
LL%EstimBaseline=.true.
Ndim=LL%Hmodel%npar+npar_ResModel
Ninfer=count(LL%Hmodel%ParFit)+npar_ResModel
if(allocated(activeSet)) deallocate(activeSet);allocate(activeSet(Ninfer))
if(allocated(tetaHat)) deallocate(tetaHat);allocate(tetaHat(Ninfer))
if(allocated(Hessian)) deallocate(Hessian);allocate(Hessian(Ninfer,Ninfer))
call StitchIt(LL%Hmodel%parLo,LL%Hmodel%ParFit,(/(0._mrk,i=1,npar_ResModel)/),xLo)
call StitchIt(LL%Hmodel%parHi,LL%Hmodel%ParFit,(/(HugeRe,i=1,npar_ResModel)/),xHi)
call StitchIt(LL%Hmodel%parLo,LL%Hmodel%ParFit,(/(0._mrk,i=1,npar_ResModel)/),msLo)
call StitchIt(LL%Hmodel%parHi,LL%Hmodel%ParFit,(/(HugeRe,i=1,npar_ResModel)/),msHi)
call StitchIt(LL%Hmodel%parScal,LL%Hmodel%ParFit,(/(10._mrk,i=1,npar_ResModel)/),xScale)
call StitchIt(LL%Hmodel%parDef,LL%Hmodel%ParFit,(/(1._mrk,i=1,npar_ResModel)/),x0)
activeSet=0
call optimiserWrap(smeth=smeth,nopt=nopt,evalFunc=LoL_LogLikelihoodWrapper,&
                          nDim=Ninfer,&
                          xLo=xLo,xHi=xHi,&
                          msLo=msLo,msHi=msHi,&
                          xscale=xscale,&
                          activeSet=activeSet,&
                          msFile=msFile, &
                          imeth_qn=5,gmeth=1,hmeth_qn=6,himeth_qn=1,&  ! QN options
                          mMem_lbfgs=10,&                      ! LBFGS method
                          mometh=mometh,&
                          x0=x0, &
                          xMin=tetaHat, fMin=fMin, &
                          chkGrad=.false.,chkHess=.false.,&
                          inverseHessian=.true., hess=Hessian,&
                          iterMax=iterMax, hisFile=hisFile, &
                          iternfo=0,&
                          resFile=resFile,&
                          fcalls=fcalls,gcalls=gcalls,hcalls=hcalls,&
                          err=err, message=mess)
if(err>0) then;mess="LoL_EstimPar: "//trim(mess);return;endif
call BlendIt(LL%Hmodel%ParDef,LL%Hmodel%ParFit,tetaHat(1:count(LL%Hmodel%ParFit)),LL%Baseline)
LL%BaselineResPar=tetaHat((count(LL%Hmodel%ParFit)+1):)
call UnrollIt(LL%Baseline,Nt,parmatrix)
call LoL_runHmodel(LL%Hmodel,LL%input,parmatrix,LL%BaselineSim,feas,err,mess)
if(err>0) then;mess="LoL_EstimPar: "//trim(mess);return;endif
BaselinePar=LL%Baseline;BaselineSim=LL%BaselineSim
! Now proceed to estimating bumps
LL%EstimBaseline=.false.
if(err>0) then;mess="LoL_EstimPar: "//trim(mess);return;endif
Ndim=LL%Hmodel%npar+npar_ResModel
Ninfer=1+npar_ResModel
call StitchIt((/-HugeRe/),(/.true./),(/(0._mrk,i=1,npar_ResModel)/),xLo)
call StitchIt((/HugeRe/),(/.true./),(/(HugeRe,i=1,npar_ResModel)/),xHi)
call StitchIt((/-HugeRe/),(/.true./),(/(0._mrk,i=1,npar_ResModel)/),msLo)
call StitchIt((/HugeRe/),(/.true./),(/(HugeRe,i=1,npar_ResModel)/),msHi)
call StitchIt((/10._mrk/),(/.true./),(/(10._mrk,i=1,npar_ResModel)/),xScale)
call StitchIt((/0._mrk/),(/.true./),LL%BaselineResPar,x0)
activeSet=0
if(allocated(activeSet)) deallocate(activeSet);allocate(activeSet(Ninfer))
if(allocated(tetaHat)) deallocate(tetaHat);allocate(tetaHat(Ninfer))
if(allocated(Hessian)) deallocate(Hessian);allocate(Hessian(Ninfer,Ninfer))
if(allocated(bumps)) deallocate(bumps);allocate(bumps(Nt))
bumps=NAval
do t=1,Nt,step
    write(*,*) 'Processing time step:',t,'/',Nt
    call LoL_GetDist(covariate=LL%covariate,TargetTime=t,DoStdz=.true.,dist=LL%dist,err=err,mess=mess)
    call optimiserWrap(smeth=smeth,nopt=nopt,evalFunc=LoL_LogLikelihoodWrapper,&
                          nDim=Ninfer,&
                          xLo=xLo,xHi=xHi,&
                          msLo=msLo,msHi=msHi,&
                          xscale=xscale,&
                          activeSet=activeSet,&
                          msFile=msFile, &
                          imeth_qn=5,gmeth=1,hmeth_qn=6,himeth_qn=1,&  ! QN options
                          mMem_lbfgs=10,&                      ! LBFGS method
                          mometh=mometh,&
                          x0=x0, &
                          xMin=tetaHat, fMin=fMin, &
                          chkGrad=.false.,chkHess=.false.,&
                          inverseHessian=.true., hess=Hessian,&
                          iterMax=iterMax, hisFile=hisFile, &
                          iternfo=0,&
                          resFile=resFile,&
                          fcalls=fcalls,gcalls=gcalls,hcalls=hcalls,&
                          err=err, message=mess)
    if(err>0) then
        mess="LoL_EstimPar: "//trim(mess);
        write(*,*) trim(mess)
        cycle
        bumps(t)=NAval;x0(1)=0._mrk
    else
        bumps(t)=tetaHat(1);x0(1)=tetaHat(1)
    endif
enddo
EstimatedBumps=bumps
! Try a last run with bumped parameters
call Interpolate(LL%covariate,EstimatedBumps,LL%Baseline,LL%TargetPar,NAval,EstimParOut)
EstimatedPar=EstimParOut(:,LL%TargetPar)
call LoL_runHmodel(LL%Hmodel,LL%input,EstimParOut,EstimatedSim,feas,err,mess)
if(err>0 .or. (.not.feas)) then
        mess="LoL_EstimPar: Warning:unfeasible last run with EstimPar";
        write(*,*) trim(mess)
        EstimatedSim=NAval
        err=-1
endif
end subroutine LoL_EstimPar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine StitchIt(u1,mask,u2,out)
real(mrk),intent(in)::u1(:),u2(:)
logical,intent(in)::mask(:)
real(mrk),pointer::out(:)
integer(mik)::n
n=count(mask)+size(u2)
if(associated(out)) deallocate(out);allocate(out(n))
out=(/pack(u1,mask),u2/)
end subroutine StitchIt
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine BlendIt(u1,mask,u2,out)
real(mrk),intent(in)::u1(:),u2(:)
logical,intent(in)::mask(:)
real(mrk),intent(out)::out(:)
integer(mik)::k,i
out=u1;k=0
do i=1,size(mask)
    if(mask(i)) then;k=k+1;out(i)=u2(k);endif
enddo
end subroutine BlendIt
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine UnrollIt(u1,n,out)
real(mrk),intent(in)::u1(:)
integer(mik),intent(in)::n
real(mrk),pointer::out(:,:)
integer(mik)::i
if(associated(out)) deallocate(out);allocate(out(n,size(u1)))
do i=1,n;out(i,:)=u1;enddo
end subroutine UnrollIt
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Interpolate(covariates,EstimParIn,BaselinePar,TargetPar,NAval,EstimParOut)
! Fills in EstimPar with the nearest (in covariate space) neighbours' EstimPar value
real(mrk),intent(in)::covariates(:,:),EstimParIn(:),BaselinePar(:),NAval
integer(mik),intent(in)::TargetPar
real(mrk),intent(out)::EstimParOut(size(covariates,dim=1),size(BaselinePar))
! locals
integer(mik)::t,Nt,N,mini(1),err
real(mrk)::dist(size(covariates,dim=1))
logical::mask(size(covariates,dim=1))
character(250)::mess

Nt=size(covariates,dim=1)
mask=(EstimParIn/=NAval)
N=count(mask)
do t=1,Nt
    EstimParOut(t,:)=BaselinePar
    if(mask(t)) then
        EstimParOut(t,TargetPar)=BaselinePar(TargetPar)+EstimParIn(t)
    else
        call LoL_GetDist(covariate=covariates,TargetTime=t,DoStdz=.true.,dist=dist,err=err,mess=mess)
        mini=minloc(dist,mask) !Nearest neighbour
        EstimParOut(t,TargetPar)=BaselinePar(TargetPar)+EstimParIn(mini(1))
    endif
enddo
end subroutine Interpolate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoL_GetDist(covariate,TargetTime,DoStdz,dist,err,mess)

!^**********************************************************************
!^* Purpose: Defines a distance in covariate space
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use numerix_dmsl_kit,only:getmeanvar
real(mrk),intent(in)::covariate(:,:)
integer(mik),intent(in)::TargetTime
logical,intent(in)::DoStdz
real(mrk),intent(out)::dist(size(covariate,dim=1))
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::t,Nt
real(mrk)::m(size(covariate,dim=2)),v(size(covariate,dim=2)),sd(size(covariate,dim=2))

err=0;mess='';dist=undefRN
Nt=size(covariate,dim=1)
if(TargetTime<=0 .or. TargetTime>Nt) then;err=1;mess='LoL_GetDist:Fatal:impossible value [TargetTime]';return;endif
if(DoStdz) then
    call getmeanvar(x=covariate,package=2,mean=m,var=v,method='c',err=err,message=mess)
    if(err>0) then;mess="LoL_GetDist: "//trim(mess);return;endif
    sd=sqrt(v)
else
    sd=1._mrk
endif
do t=1,Nt
    v=(1._mrk/sd)*(covariate(t,:)-covariate(TargetTime,:))
    dist(t)=dot_product(v,v)
enddo

end subroutine LoL_GetDist
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module LocalLikelihood_tools
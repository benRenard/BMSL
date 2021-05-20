module ExtraFlo_RegionalModels

!~**********************************************************************
!~* Purpose: 
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 25/03/2012
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

public :: Regional_IF,Regional_IR

type,public:: RegionalEstimateType
    real(mrk),pointer::Kpred(:)=>NULL()
    real(mrk),pointer::Kpredictive(:,:)=>NULL()
end type RegionalEstimateType

Contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Index flood regional model
subroutine Regional_IF(IF_C,IF_V,& ! Index flood values 
                    xy_C,xy_V,& ! Stations coordinates
                    covC,covV,& ! covariates
                    NormalizedData,& ! Normalized data for regional distribution
                    dist,smoothing,& ! Model
                    FilePrefix,path2Rscript,path2Results,& ! Paths
                    nAdapt,nCycles,Nslim,Nrep,& ! MCMC
                    MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor,& !MCMC
                    LogSpace,modal,pred,err,mess)
! USE
use utilities_dmsl_kit, only:number_string,GetSpareUnit
use BayesianEstimation_tools
use ExtraFlo_utilities
use Distribution_tools
use EmpiricalStats_tools, only:GetEmpiricalStats
use RegressionFunk, only:identity
! IN
real(mrk), intent(in)::IF_C(:),IF_V(:),covC(:,:),covV(:,:),xy_C(:,:),xy_V(:,:)
type(VoVtype),intent(in)::NormalizedData(:)
character(*),intent(in)::FilePrefix,dist,smoothing,path2Rscript,path2Results
integer(mik), intent(in)::nAdapt,nCycles,Nslim,Nrep
real(mrk), intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor
Logical, intent(in)::LogSpace
! OUT
real(mrk),pointer::modal(:,:),pred(:,:)
character(*), intent(out)::mess
integer(mik),intent(out)::err
! locals
integer(mik)::nc,nv,nregion,foo(1),regType,nDistPar,site,reg,p,i,j,nsim_pred,unt1,unt2              
real(mrk), allocatable::WLS(:), Y(:),Y_V(:),mode_RegionalDist(:,:),&
                        modal_dummy(:),par_dummy(:),par_dummy2(:,:),pmean(:),pstd(:)
logical::feas
character(250)::LinkFunk,PredictiveFile,ModalPredFile,ResidualFile,MCMCfile,&
                PriorStdFile,PriorMeanFile
type(RegionalEstimateType), allocatable::Regine(:)

! Preliminaries
if(allocated(Regine)) deallocate(Regine)
allocate(Regine(1))
nc=size(IF_C);nv=size(IF_V)
if(allocated(WLS)) deallocate(WLS);allocate(WLS(nc))
if(allocated(Y)) deallocate(Y);allocate(Y(nc))
if(allocated(Y_V)) deallocate(Y_V);allocate(Y_V(nv))
foo=maxval(nint(covC(:,1),mik));nregion=foo(1)
call GetParNumber(DistID=dist, npar=nDistPar, err=err, mess=mess)
If(err>0) then
  mess='Regional_IF:'//trim(mess);return
endif
if(allocated(mode_RegionalDist)) deallocate(mode_RegionalDist);allocate(mode_RegionalDist(nregion,nDistPar))    
if(allocated(modal_dummy)) deallocate(modal_dummy);allocate(modal_dummy(size(Regine)))
if(allocated(par_dummy)) deallocate(par_dummy);allocate(par_dummy(nDistPar))

! Index Flood regression preliminaries
Y=IF_C;Y_V=IF_V
WLS=0._mrk ! no weighting for IF-regression
regType=1 ! Additive (1) or Power (2) regression?
if(LogSpace) then
    LinkFunk= identity
else
    LinkFunk= 'Log' 
endif
PredictiveFile=trim(FilePrefix)//"_IF_Predictive.txt"
ModalPredFile=trim(FilePrefix)//"_IF_ModalPred.txt"
ResidualFile=trim(FilePrefix)//"_IF_Residual.txt"
MCMCfile=trim(FilePrefix)//"_IF_MCMC.txt"
! Index Flood regression 
call BigRegression(Y,Y_V,WLS,regType,LinkFunk,smoothing,path2Rscript,path2Results,&
                     nregion,covC,covV,xy_C,xy_V,&
                     nAdapt,nCycles,Nslim,Nrep,&
                     MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor,&
                     MCMCfile,ResidualFile,ModalPredFile,PredictiveFile,&
                     Regine(1)%Kpred,Regine(1)%Kpredictive,err,mess)
If(err>0) then
  mess='Regional_IF:'//trim(mess);return
endif
! Estimate regional distribution in each region
do reg=1,nregion
    ! Estimate regional distribution
    call EstimateRegionalDist(X=NormalizedData(reg)%val,dist=dist,&
                          path2MCMC=trim(path2Results)//trim(FilePrefix)//'_RD_MCMC_'//trim(number_string(reg))//'.txt',&
                          nAdapt=nAdapt,nCycles=nCycles,Nslim=Nslim,&
                          MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
                          DownMult=DownMult,UpMult=UpMult,BurnFactor=BurnFactor,&
                          modal=mode_RegionalDist(reg,:),err=err,mess=mess)
    If(err>0) then
        mess='Regional_IF:'//trim(mess);return
    endif
enddo 
! Estimate at validation sites
PriorMeanFile=trim(FilePrefix)//"_PriorMean.txt"
PriorStdFile=trim(FilePrefix)//"_PriorStd.txt"
call GetSpareUnit(unt1,err,mess)
open(unit=unt1, file=trim(path2Results)//trim(PriorMeanFile), status='REPLACE',iostat=err)
call GetSpareUnit(unt2,err,mess)
open(unit=unt2, file=trim(path2Results)//trim(PriorStdFile), status='REPLACE',iostat=err)
nsim_pred=size(Regine(1)%Kpredictive(site,:))
if(associated(modal)) deallocate(modal);allocate(modal(nV,nDistPar))
if(associated(pred)) deallocate(pred);allocate(pred(nV,nsim_pred))
if(allocated(par_dummy2)) deallocate(par_dummy2);allocate(par_dummy2(nsim_pred,nDistPar))
if(allocated(pmean)) deallocate(pmean);allocate(pmean(nDistPar))
if(allocated(pstd)) deallocate(pstd);allocate(pstd(nDistPar))
do site=1,nV
    do j=1,size(Regine)
        modal_dummy(j)=Regine(j)%Kpred(site)
    enddo
    reg=nint(covV(site,1),mik)
    modal(site,1)=modal_dummy(1)*mode_RegionalDist(reg,1) ! location
    modal(site,2)=mode_RegionalDist(reg,2) ! cv
    if(nDistPar>=3) modal(site,3)=mode_RegionalDist(reg,3) !shape
    do i=1,nsim_pred
        par_dummy(1)=Regine(1)%Kpredictive(site,i)*mode_RegionalDist(reg,1) ! location
        par_dummy(2)=mode_RegionalDist(reg,2) ! cv
        if(nDistPar>=3) par_dummy(3)=mode_RegionalDist(reg,3) ! shape
        par_dummy2(i,:)=par_dummy
        call Generate(DistId=dist,par=par_dummy,&
		        gen=pred(site,i),feas=feas,err=err,mess=mess)
        If(err>0) then
            mess='Regional_IF:'//trim(mess);return
        endif  
    enddo
    do i=1,nDistPar
        call GetEmpiricalStats(x=par_dummy2(:,i), & 
                mean=pmean(i),std=pstd(i),err=err,mess=mess)
        If(err>0) then
            mess='Regional_IF:'//trim(mess);return
        endif  
    enddo
    write(unt1,'(<nDistPar>e14.6)') pmean
    write(unt2,'(<nDistPar>e14.6)') pstd
enddo
close(unt1);close(unt2)                      
end subroutine Regional_IF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Individual regressions regional model
subroutine Regional_IR(par_C,par_V,wls_weights,& ! Parameter values 
                    xy_C,xy_V,& ! Stations coordinates
                    covC,covV,& ! covariates
                    dist,smoothing,regression,DoShapeRegression,& ! Model
                    FilePrefix,path2Rscript,path2Results,& ! Paths
                    nAdapt,nCycles,Nslim,Nrep,& ! MCMC
                    MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor,& !MCMC
                    LogSpace,modal,pred,err,mess)
! USE
use utilities_dmsl_kit, only:GetSpareUnit
use BayesianEstimation_tools
use RegressionFunk, only:identity
use Distribution_tools
use EmpiricalStats_tools, only:GetEmpiricalStats
! IN
real(mrk), intent(in)::par_C(:,:),par_V(:,:),wls_weights(:,:),&
                       covC(:,:),covV(:,:),xy_C(:,:),xy_V(:,:)
character(*),intent(in)::FilePrefix,dist,smoothing,regression,path2Rscript,path2Results
integer(mik), intent(in)::nAdapt,nCycles,Nslim,Nrep
logical, intent(in)::DoShapeRegression,LogSpace
real(mrk), intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor
! OUT
real(mrk),pointer::modal(:,:),pred(:,:)
character(*), intent(out)::mess
integer(mik),intent(out)::err
! locals
integer(mik)::nc,nv,regType,nDistPar,site,nsim_pred,parx,nregion,foo(1),i,j,unt1,unt2              
real(mrk), allocatable::WLS(:), Y(:),Y_V(:),modal_dummy(:),par_dummy(:),par_dummy2(:,:),pmean(:),pstd(:)
logical::feas
character(250)::LinkFunk,PredictiveFile,ModalPredFile,ResidualFile,MCMCfile,&
                PriorStdFile,PriorMeanFile
type(RegionalEstimateType), allocatable::Regine(:)

! Preliminaries
call GetParNumber(DistID=dist, npar=nDistPar, err=err, mess=mess)
If(err>0) then
  mess='Regional_IR:'//trim(mess);return
endif
if(allocated(Regine)) deallocate(Regine)
allocate(Regine(nDistPar))
nc=size(par_C,dim=1);nv=size(par_V,dim=1)
foo=maxval(nint(covC(:,1),mik));nregion=foo(1)
if(allocated(WLS)) deallocate(WLS);allocate(WLS(nc))
if(allocated(Y)) deallocate(Y);allocate(Y(nc))
if(allocated(Y_V)) deallocate(Y_V);allocate(Y_V(nv))
if(allocated(modal_dummy)) deallocate(modal_dummy);allocate(modal_dummy(size(Regine)))
if(allocated(par_dummy)) deallocate(par_dummy);allocate(par_dummy(nDistPar))
! Location regression preliminaries
parx=1
Y=par_C(:,parx);Y_V=par_V(:,parx)
regType=1 ! Additive (1) or Power (2) regression?
if(LogSpace) then
    LinkFunk=identity
    WLS=SetWLS(regression,wls_weights(:,parx)) 
else
    LinkFunk= 'Log'
    WLS=SetWLS(regression,sqrt(log(1._mrk+(wls_weights(:,parx)/Y)**2))) 
endif
PredictiveFile=trim(FilePrefix)//"_LOC_Predictive.txt"
ModalPredFile=trim(FilePrefix)//"_LOC_ModalPred.txt"
ResidualFile=trim(FilePrefix)//"_LOC_Residual.txt"
MCMCfile=trim(FilePrefix)//"_LOC_MCMC.txt"
! Location regression 
call BigRegression(Y,Y_V,WLS,regType,LinkFunk,smoothing,path2Rscript,path2Results,&
                     nregion,covC,covV,xy_C,xy_V,&
                     nAdapt,nCycles,Nslim,Nrep,&
                     MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor,&
                     MCMCfile,ResidualFile,ModalPredFile,PredictiveFile,&
                     Regine(parx)%Kpred,Regine(parx)%Kpredictive,err,mess)
If(err>0) then
  mess='Regional_IR:'//trim(mess);return
endif
! CV regression preliminaries
parx=2
Y=par_C(:,parx);Y_V=par_V(:,parx)
regType=1 ! Additive (1) or Power (2) regression?
if(LogSpace) then
    LinkFunk=identity
    WLS=SetWLS(regression,wls_weights(:,parx)) 
else
    LinkFunk= 'Log'
    WLS=SetWLS(regression,sqrt(log(1._mrk+(wls_weights(:,parx)/Y)**2))) 
endif
PredictiveFile=trim(FilePrefix)//"_CV_Predictive.txt"
ModalPredFile=trim(FilePrefix)//"_CV_ModalPred.txt"
ResidualFile=trim(FilePrefix)//"_CV_Residual.txt"
MCMCfile=trim(FilePrefix)//"_CV_MCMC.txt"
! CV regression 
call BigRegression(Y,Y_V,WLS,regType,LinkFunk,smoothing,path2Rscript,path2Results,&
                     nregion,covC,covV,xy_C,xy_V,&
                     nAdapt,nCycles,Nslim,Nrep,&
                     MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor,&
                     MCMCfile,ResidualFile,ModalPredFile,PredictiveFile,&
                     Regine(parx)%Kpred,Regine(parx)%Kpredictive,err,mess)
If(err>0) then
  mess='Regional_IR:'//trim(mess);return
endif
! SHAPE regression preliminaries
if(Ndistpar>=3) then
    parx=3
    Y=par_C(:,parx);Y_V=par_V(:,parx)
    if(DoShapeRegression) then
        regType=1 ! Additive (1) or Power (2) regression? (0) for botchy constant model
    else
        regType=0 ! Additive (1) or Power (2) regression? (0) for botchy constant model
    endif
    LinkFunk= identity
    WLS=SetWLS(regression,wls_weights(:,parx)) 
    PredictiveFile=trim(FilePrefix)//"_SHAPE_Predictive.txt"
    ModalPredFile=trim(FilePrefix)//"_SHAPE_ModalPred.txt"
    ResidualFile=trim(FilePrefix)//"_SHAPE_Residual.txt"
    MCMCfile=trim(FilePrefix)//"_SHAPE_MCMC.txt"
! SHAPE regression 
    call BigRegression(Y,Y_V,WLS,regType,LinkFunk,smoothing,path2Rscript,path2Results,&
                         nregion,covC,covV,xy_C,xy_V,&
                         nAdapt,nCycles,Nslim,Nrep,&
                         MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor,&
                         MCMCfile,ResidualFile,ModalPredFile,PredictiveFile,&
                         Regine(parx)%Kpred,Regine(parx)%Kpredictive,err,mess)
    If(err>0) then
      mess='Regional_IR:'//trim(mess);return
    endif
endif
! Estimate at validation sites
PriorMeanFile=trim(FilePrefix)//"_PriorMean.txt"
PriorStdFile=trim(FilePrefix)//"_PriorStd.txt"
call GetSpareUnit(unt1,err,mess)
open(unit=unt1, file=trim(path2Results)//trim(PriorMeanFile), status='REPLACE',iostat=err)
call GetSpareUnit(unt2,err,mess)
open(unit=unt2, file=trim(path2Results)//trim(PriorStdFile), status='REPLACE',iostat=err)
nsim_pred=size(Regine(1)%Kpredictive(site,:))
if(associated(modal)) deallocate(modal);allocate(modal(nV,nDistPar))
if(associated(pred)) deallocate(pred);allocate(pred(nV,nsim_pred))
if(allocated(par_dummy2)) deallocate(par_dummy2);allocate(par_dummy2(nsim_pred,nDistPar))
if(allocated(pmean)) deallocate(pmean);allocate(pmean(nDistPar))
if(allocated(pstd)) deallocate(pstd);allocate(pstd(nDistPar))
do site=1,nV
    do j=1,size(Regine)
        modal_dummy(j)=Regine(j)%Kpred(site)
    enddo
    modal(site,:)=modal_dummy
    do i=1,nsim_pred
        do j=1,size(Regine)
            par_dummy(j)=Regine(j)%Kpredictive(site,i)
        enddo
        call Generate(DistId=dist,par=par_dummy,&
		        gen=pred(site,i),feas=feas,err=err,mess=mess)
        If(err>0) then
            mess='Regional_IR:'//trim(mess);return
        endif  
        par_dummy2(i,:)=par_dummy
    enddo
    do i=1,nDistPar
        call GetEmpiricalStats(x=par_dummy2(:,i), & 
                mean=pmean(i),std=pstd(i),err=err,mess=mess)
        If(err>0) then
            mess='Regional_IR:'//trim(mess);return
        endif  
    enddo
    write(unt1,'(<nDistPar>e14.6)') pmean
    write(unt2,'(<nDistPar>e14.6)') pstd
enddo          
close(unt1);close(unt2)            
end subroutine Regional_IR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PRIVATE SUBROUTINES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine BigRegression(Y,Y_V,WLS,regType,LinkFunk,smoothing,path2Rscript,path2Results,&
                         nregion,covC,covV,xy_C,xy_V,&
                         nAdapt,nCycles,Nslim,Nrep,&
                         MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor,&
                         MCMCfile,ResidualFile,ModalPredFile,PredictiveFile,&
                         Kpred,Kpredictive,err,mess)
use Distribution_tools
use RegressionFunk, only:identity, XtraFlo_Reg
use RegBay_tools
use EmpiricalStats_tools, only:GetEmpiricalStats,GetEmpiricalQuantile
use BayesianEstimation_tools, only:PriorListType
use RFortran
real(mrk), intent(in)::Y(:),Y_V(:),WLS(:),covC(:,:),covV(:,:),xy_C(:,:),xy_V(:,:)
integer(mik), intent(in)::nregion,nAdapt,nCycles,Nslim,Nrep,regType
real(mrk), intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor
character(*), intent(in)::MCMCfile,ResidualFile,ModalPredFile,PredictiveFile,&
                          path2Rscript,path2Results,smoothing,LinkFunk
integer(mik), intent(out)::err
character(*), intent(out)::mess
real(mrk), pointer::Kpred(:),Kpredictive(:,:)
!locals
Type(PriorListType), allocatable::PriorList_teta(:), PriorList_RemnantSigma(:)
real(mrk),allocatable::teta0(:),teta_std0(:)
real(mrk)::RemnantSigma0,RemnantSigma_std0
real(mrk), allocatable::modalPar(:), modalPred_C(:),modalPred_V(:),Res(:),KRes(:),&
            KRes_var(:),stdRes(:),gen(:)
real(mrk)::MarginalLkh
real(mrk), pointer::predictive(:,:)
integer(mik)::npar,nc,nv,i,j,npredi
logical::ok,feas

! Allocate
npar=nregion*size(covC,dim=2)
nc=size(covC,dim=1)
nv=size(covV,dim=1)
if(allocated(teta0)) deallocate(teta0)
allocate(teta0(npar))
if(allocated(teta_std0)) deallocate(teta_std0)
allocate(teta_std0(npar))
if(allocated(PriorList_teta)) deallocate(PriorList_teta)
allocate(PriorList_teta(npar))
if(allocated(PriorList_RemnantSigma)) deallocate(PriorList_RemnantSigma)
allocate(PriorList_RemnantSigma(1))
if(allocated(modalPar)) deallocate(modalPar)
allocate(modalPar(npar+1))
if(allocated(modalPred_C)) deallocate(modalPred_C)
allocate(modalPred_C(nc))
if(allocated(modalPred_V)) deallocate(modalPred_V)
allocate(modalPred_V(nv))
if(allocated(Res)) deallocate(Res)
allocate(Res(nc))
if(allocated(KRes)) deallocate(KRes)
allocate(KRes(nv))
if(allocated(KRes_var)) deallocate(KRes_var)
allocate(KRes_var(nv))
if(associated(Kpred)) deallocate(Kpred)
allocate(Kpred(nv))
if(allocated(gen)) deallocate(gen)
allocate(gen(nv))
if(allocated(stdRes)) deallocate(stdRes)
allocate(stdRes(nc))
! Settings
CALL SetDefaultPriors(Y=Y,X=covC,nregion=nregion,regType=regType,LinkFunk=LinkFunk,&
                   PriorList_teta=PriorList_teta,PriorList_RemnantSigma=PriorList_RemnantSigma,& ! Priors
                   teta0=teta0,RemnantSigma0=RemnantSigma0,teta_std0=teta_std0,RemnantSigma_std0=RemnantSigma_std0,&
                   err=err,mess=mess)
if(err>0) then
    mess="BigRegression: "//trim(mess);return
endif

!Launch IF regression
call RegBay_Fit(Y=Y,X=covC,& ! observations for predictand (Y) and predictors (X)
                   RegFunk=XtraFlo_reg,Link=LinkFunk,& ! Regression function, link function
                   Xtra=1._mrk*(/nregion,regType/),& ! optional arguments to pass to the regression funk
                   PriorList_teta=PriorList_teta,PriorList_RemnantSigma=PriorList_RemnantSigma,& ! Priors
                   WLS_sdev=WLS,& ! optional vector (1:n) containing standard deviations of individuals yi's (for WLS-like inferrence). Default is zero (ie SLS) 
                   !!!!!!! Tuning of the MCMC sampler !!!!!!!!
                   teta0=teta0,RemnantSigma0=RemnantSigma0, &!initial values for teta and remnant std
			       teta_std0=teta_std0,RemnantSigma_std0=RemnantSigma_std0,& ! initial values for the std of jump distribution
			       nAdapt=nAdapt,nCycles=nCycles,& 
			       MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
			       DownMult=DownMult,UpMult=UpMult,&
                   !!!!!!! END Tuning of the MCMC sampler !!!!!!!!
			       OutFile=trim(path2Results)//trim(MCMCfile), & ! Output file (for MCMC samples)
			       err=err,mess=mess)! error handling
if(err>0) then
    mess="BigRegression: "//trim(mess);return
endif

! Predict on Calib data
call RegBay_Predict(File=trim(path2Results)//trim(MCMCfile),&
                   Nburn=NINT(BurnFactor*nAdapt*nCycles),&
                   Nread=NINT((1._mrk-BurnFactor)*nAdapt*nCycles),&
                   Nslim=Nslim,Nrep=Nrep,& ! IN: Read properties
                   X=covC,& ! Matrix of covariates where prediction is made
                   RegFunk=XtraFlo_reg,Link=LinkFunk,& ! Regression function, link function
                   Xtra=1._mrk*(/nregion,regType/),& ! optional arguments to pass to the regression funk
                   modalPar=modalPar,modalPred=modalPred_C,& ! OUT: modal parameter estimates & correponding predictions
                   predictive=predictive,& !OUT: replicates from the predictive distribution
                   MarginalLkh=MarginalLkh,& !OUT: marginal likelihood (crude) approximation
				   err=err, mess=mess)
if(err>0) then
    mess="BigRegression: "//trim(mess);return
endif
write(*,*)"---------------------------------------"
write(*,*) "Modal parameters:", modalPar
write(*,*)"---------------------------------------"
write(*,*) "Marginal lkh:", MarginalLkh
write(*,*)"---------------------------------------"


!residuals analysis
call RegBay_Residuals(Y=Y,X=covC,& ! observations for predictand (Y) and predictors (X)
                   RegFunk=XtraFlo_reg,Link=LinkFunk,& ! Regression function, link function
                   Xtra=1._mrk*(/nregion,regType/),& ! optional arguments to pass to the regression funk
                   RegPar=modalPar,& ! Par used to apply regression (e.g. modal)
                   WLS_sdev=WLS,& ! optional vector (1:n) containing standard deviations of individuals yi's (for WLS-like inferrence). Default is zero (ie SLS) 
                   Res=Res,& ! OUT: residuals=link(y)-regression(x|RegPar)
                   standardizedRes=stdRes,& ! OUT: residuals=link(y)-regression(x|RegPar)
                   ResFile=trim(path2Results)//trim(ResidualFile),&! optional, adress to write result file
				   err=err, mess=mess)
if(err>0) then
    mess="BigRegression: "//trim(mess);return
endif

! Residual kriging using R
if(smoothing/='NOK') then
    ok = Rinit()
    ! Send XY_C, XY_V,Standardized residuals
    ok=Rput('XY',xy_C/1000)
    ok=Rput('XY_V',xy_V/1000)
    ok=Rput('res',stdRes)
    ok=Reval(trim("source('")//trim(path2Rscript)//trim("')"))
    ok=Rget('k$predict',Kres)
    ok=Rget('k$krige.var',Kres_var)
    ok=Rclose()
else
    Kres=0._mrk
    Kres_var=0._mrk
endif
! Predict at validation sites
call RegBay_Predict(File=trim(path2Results)//trim(MCMCfile),&
                   Nburn=NINT(BurnFactor*nAdapt*nCycles),&
                   Nread=NINT((1._mrk-BurnFactor)*nAdapt*nCycles),&
                   Nslim=Nslim,Nrep=Nrep,& ! IN: Read properties
                   X=covV,& ! Matrix of covariates where prediction is made
                   RegFunk=XtraFlo_reg,Link=LinkFunk,& ! Regression function, link function
                   Xtra=1._mrk*(/nregion,regType/),& ! optional arguments to pass to the regression funk
                   PropagateSigmaU=.false.,& ! Set to false since this going to be handled in the regression kriging step
                   modalPar=modalPar,modalPred=modalPred_V,& ! OUT: modal parameter estimates & correponding predictions
                   predictive=predictive,& !OUT: replicates from the predictive distribution
                   MarginalLkh=MarginalLkh,& !OUT: marginal likelihood (crude) approximation
				   err=err, mess=mess)
if(err>0) then
    mess="BigRegression: "//trim(mess);return
endif

! Get Modal Prediction
call CorrectPred(Pred=modalPred_V,Kres=modalPar(npar+1)*Kres,&
                 RegFunk=XtraFlo_reg,Link=LinkFunk,& ! Regression function, link function
                 Xtra=1._mrk*(/nregion,regType/),&
                 Kpred=Kpred,err=err,mess=mess)
if(err>0) then
    mess="BigRegression: "//trim(mess);return
endif
! Write results
open(unit=1,file=trim(path2Results)//trim(ModalPredFile),status='replace')
write(1,'(3A15)') 'yobs','ypred','yKpred'
do i=1,nv
    write(1,'(3F15.3)') Y_V(i),modalPred_V(i),Kpred(i)
enddo
close(1)
! Get predictive
if(associated(Kpredictive)) deallocate(Kpredictive)
allocate(Kpredictive(size(predictive,dim=1),size(predictive,dim=2)))

do i=1,size(predictive,dim=2)
    ! Generate kriged residual
    do j=1,size(Kres)
        if(Kres_var(j)>0._mrk) then
            call Generate(DistId=GAUSS,par=(/Kres(j),sqrt(Kres_var(j))/),gen=gen(j),feas=feas,err=err,mess=mess)
            if(err>0) then
                mess="BigRegression: "//trim(mess);return
            endif
        else
            gen(j)=Kres(j)
        endif
    enddo
    call CorrectPred(Pred=predictive(:,i),Kres=modalPar(npar+1)*gen,&
                 RegFunk=XtraFlo_reg,Link=LinkFunk,& ! Regression function, link function
                 Xtra=1._mrk*(/nregion,regType/),&
                 Kpred=Kpredictive(:,i),err=err,mess=mess)
    if(err>0) then
        mess="BigRegression: "//trim(mess);return
    endif
enddo
! Write results
open(unit=1,file=trim(path2Results)//trim(PredictiveFile),status='replace')
write(1,'(A90)') "Replicates (columns) from predictive at each validation site (rows)"
npredi=size(predictive,dim=2)
do i=1,nv
    write(1,'(<npredi>F8.2)') Kpredictive(i,:)
enddo
close(1)
end subroutine BigRegression
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine EstimateRegionalDist(X,dist,path2MCMC,&
                              nAdapt,nCycles,Nslim,&
                              MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor,&
                              modal,err,mess)
use BayesianEstimation_tools
use Distribution_tools
character(*), intent(in)::dist,path2MCMC
real(mrk), intent(in)::X(:)
integer(mik), intent(in)::nAdapt,nCycles,Nslim
real(mrk), intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor
real(mrk), intent(out)::modal(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
Type(PriorListType), allocatable::PriorList_RegionalDist(:)
real(mrk),allocatable::start_RegionalDist(:),startStd_RegionalDist(:)
real(mrk),pointer::mcmc(:,:),LogPost(:)
integer(mik)::nDistPar,i
real(mrk), parameter::MultFactor=0.1_mrk, BigRe=1000000._mrk

call GetParNumber(DistID=dist, npar=nDistPar, err=err, mess=mess)
If(err>0) then
  mess='EstimateRegionalDist:'//trim(mess);return
endif
! Set priors and starting point/stdev
if(allocated(PriorList_RegionalDist)) deallocate(PriorList_RegionalDist)
if(allocated(start_RegionalDist)) deallocate(start_RegionalDist)
if(allocated(startStd_RegionalDist)) deallocate(startStd_RegionalDist)
allocate(PriorList_RegionalDist(nDistPar),start_RegionalDist(nDistPar),startStd_RegionalDist(nDistPar))
do i=1,nDistPar ! large unif for everybody
    if(allocated(PriorList_RegionalDist(i)%par)) deallocate(PriorList_RegionalDist(i)%par)
    allocate(PriorList_RegionalDist(i)%par(2))
    PriorList_RegionalDist(i)%dist=UNIF
    PriorList_RegionalDist(i)%par=(/-BigRe,BigRe/)
enddo
call GetRoughEstimate(DistID=dist, X=X, par=start_RegionalDist, err=err, mess=mess)
If(err>0) then
  mess='EstimateRegionalDist:'//trim(mess);return
endif
!start_RegionalDist(1)=1._mrk
!start_RegionalDist(2)=0.5_mrk
!if(nDistPar>2) start_RegionalDist(3)=0._mrk
startStd_RegionalDist=abs(start_RegionalDist)*MultFactor
where(start_RegionalDist==0._mrk) startStd_RegionalDist=MultFactor

call Bayes_EstimPar(X=X,distID=dist,mv=-99._mrk,&
        PriorList=PriorList_RegionalDist,&
		! Tuning of the MCMC sampler
		start=start_RegionalDist,startStd=startStd_RegionalDist,&
		nAdapt=nAdapt,nCycles=nCycles,&
		MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
		DownMult=DownMult,UpMult=UpMult,&
		OutFile=trim(path2MCMC), &
		! error handling
		err=err,mess=mess)
If(err>0) then
  mess='EstimateRegionalDist:'//trim(mess);return
endif
! Load slimmed MCMC and get mode
call LoadMCMCsample(File=trim(path2MCMC),BurnFactor=BurnFactor,&
                        Nslim=Nslim,distID=dist,& 
						mcmc=mcmc,LogPost=LogPost,err=err,mess=mess)
If(err>0) then
  mess='EstimateRegionalDist:'//trim(mess);return
endif
call GetMCMCSummary(mcmc=mcmc,LogPost=LogPost,modal=modal,&
                          err=err,mess=mess)
If(err>0) then
  mess='EstimateRegionalDist:'//trim(mess);return
endif

end subroutine EstimateRegionalDist
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pure subroutine SetDefaultPriors(Y,X,nregion,regType,LinkFunk,&
                       PriorList_teta,PriorList_RemnantSigma,& ! Priors
                       teta0,RemnantSigma0,teta_std0,RemnantSigma_std0,&
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
use BayesianEstimation_tools
use Distribution_tools
real(mrk), intent(in)::Y(:),X(:,:)
integer(mik), intent(in)::nregion,regType
character(*), intent(in)::LinkFunk
real(mrk), intent(out)::teta0(:),RemnantSigma0,teta_std0(:),RemnantSigma_std0
integer(mik), intent(out)::err
character(*),intent(out)::mess
Type(PriorListType), intent(out)::PriorList_teta(:), PriorList_RemnantSigma(:)
! locals
real(mrk), parameter::MultFactor=0.1_mrk, TetaFactor=0.1_mrk, RSigma0=1._mrk,BigRe=10000._mrk
integer(mik)::i,npar

err=0;mess=''

do i=1,nregion
    ! intercept = region average
    if(trim(LinkFunk)=="Log") then
        teta0(i)=sum(log(Y),NINT(X(:,1))==i)/count(NINT(X(:,1))==i)
    else
        teta0(i)=sum(Y,NINT(X(:,1))==i)/count(NINT(X(:,1))==i)
    endif
    teta_std0(i)=MultFactor*abs(teta0(1))
enddo
! other teta = starts at zero
teta0(nregion+1:)=0._mrk
teta_std0(nregion+1:)=TetaFactor
RemnantSigma0=RSigma0
RemnantSigma_std0=MultFactor*RemnantSigma0
!priors
npar=nregion*size(X,dim=2)
do i=1,npar ! large unif for everybody
    if(allocated(PriorList_teta(i)%par)) deallocate(PriorList_teta(i)%par)
    allocate(PriorList_teta(i)%par(2))
    PriorList_teta(i)%dist=UNIF
    PriorList_teta(i)%par=(/-BigRe,BigRe/)
enddo
if(allocated(PriorList_RemnantSigma(1)%par)) deallocate(PriorList_RemnantSigma(1)%par)
allocate(PriorList_RemnantSigma(1)%par(2))
PriorList_RemnantSigma%dist=UNIF
PriorList_RemnantSigma(1)%par=(/-BigRe,BigRe/)

end subroutine SetDefaultPriors
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine CorrectPred(Pred,Kres,&
                     RegFunk,Link,& ! Regression function, link function
                     Xtra,&
                     Kpred,err,mess)

!^**********************************************************************
!^* Purpose: Correct raw prediction with kriged residuals
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
use RegressionFunk
real(mrk), intent(in)::Pred(:),Kres(:),Xtra(:)
character(*), intent(in)::RegFunk,Link
integer(mik), intent(out)::err
character(*),intent(out)::mess
real(mrk), intent(out)::Kpred(size(Pred))
!locals
integer(mik)::i,n
real(mrk)::z
logical::feas

err=0;mess=''
n=size(Pred)
do i=1,n
    ! apply Link
    call LinkFunk(LFunk=Link,x=pred(i),fx=z,feas=feas,err=err,mess=mess)
    if(err>0) then
        mess="CorrectPred: "//trim(mess);return
    endif
    if(.not.feas) then
        Kpred(i)=undefRN
        cycle
    endif
    ! add kriged residual
    z=z+Kres(i)
    ! re-apply inverse link
    call InverseLinkFunk(IvsLinkFunk=Link,x=z,fx=KPred(i),feas=feas,err=err,mess=mess)
    if(err>0) then
        mess="CorrectPred: "//trim(mess);return
    endif
    if(.not.feas) then
        Kpred(i)=undefRN
        cycle
    endif
enddo

end subroutine CorrectPred
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function SetWLS(regression,V)
  character(*), intent(in)::regression
  real(mrk), intent(in)::V(:)
  real(mrk)::WLS(size(V))
  real(mrk)::SetWLS(size(V))

  select case(regression)
    case('SLS')
      SetWLS=0._mrk
    case('WLS')
      SetWLS=V
    case default
      Write(*,*) "ERROR in SetWLS!!!"
      pause
  end select
end function SetWLS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end module ExtraFlo_RegionalModels
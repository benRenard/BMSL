module RegBay_tools

!~**********************************************************************
!~* Purpose: Multiple regression à la Bayesian - SLS and WLS for now
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified: 17/05/2011
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
use BayesianEstimation_tools, only:PriorListType

implicit none
Private
public :: RegBay_Fit,RegBay_Predict,RegBay_Residuals

! variables globally available to this module
type, public:: RegModType ! the "regression model" object
	character(100)::RegFunk='AintGotNoName' ! Regression formula - see RegressionFunk.f90 for the catalogue
	character(100)::Link='AintGotNoName' !  Link function - see RegressionFunk.f90 for the catalogue
	integer(mik)::nteta=undefIN ! number of regression parameters
	integer(mik)::nobs=undefIN ! number of observations used to fit the regression
	real(mrk), allocatable::Y(:) ! predictand
	real(mrk), allocatable::X(:,:) ! predictors
	real(mrk), allocatable::WLS_sdev(:) ! sstandard deviations of individuals yi's (for WLS-like inferrence)
    type(PriorListType), allocatable:: PriorList_RemnantSigma(:) ! priors for remnant std - see BayesianEstimation_tools module for PriorListType definition
    type(PriorListType), allocatable:: PriorList_teta(:) ! priors for regression parameters
    real(mrk), allocatable::Xtra(:) ! optional arguments to pass to the regression funk
end type RegModType
type(RegModType):: MoMo


Contains

subroutine RegBay_Fit(Y,X,& ! observations for predictand (Y) and predictors (X)
                       RegFunk,Link,& ! Regression function, link function
                       Xtra,& ! optional arguments to pass to the regression funk
                       PriorList_teta,PriorList_RemnantSigma,& ! Priors
                       WLS_sdev,& ! optional vector (1:n) containing standard deviations of individuals yi's (for WLS-like inferrence). Default is zero (ie SLS) 
                       !!!!!!! Tuning of the MCMC sampler !!!!!!!!
                       teta0,RemnantSigma0, &!initial values for teta and remnant std
				       teta_std0,RemnantSigma_std0,& ! initial values for the std of jump distribution
				       nAdapt,nCycles,& 
				       MinMoveRate,MaxMoveRate,&
				       DownMult,UpMult,&
                       !!!!!!! END Tuning of the MCMC sampler !!!!!!!!
				       OutFile, & ! Output file (for MCMC samples)
				       err,mess)! error handling
!^**********************************************************************
!^* Purpose: compute post(teta,RemnantSigma)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 30/01/2012
!^**********************************************************************
!^* Comments: New Xtra variable, used to pass optional arguments
!^* (eg number of regions, regression type, atc)
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Y, predictand
!^*		2. X, predictor
!^*		3. RegFunk, regression function (see RegressionFunk.f90)
!^*		4. Link, link function (see RegressionFunk.f90)
!^*		5. [Xtra], optional arguments to pass to regression funk
!^*		6. PriorList_teta, priors for teta
!^*		7. PriorList_RemnantSigma, prior for remnant error
!^*		9. [WLS_sdev], standard deviations of individuals yi's. Default is zero (ie SLS) 
!^*		10 to 18. tuning of MCMC sampler. See MCMC_strategy_tools.f90 for details. 
!^*		19. Output File. 
!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
use RegressionFunk, only:GetRegParNumber
use MCMCStrategy_tools, only:Adaptive_Metro_OAAT
use utilities_dmsl_kit, only:number_string

! INPUTS
real(mrk), intent(in)::Y(:),X(:,:)
real(mrk), intent(in),optional::Xtra(:)
character(*), intent(in)::RegFunk,Link
Type(PriorListType), intent(in)::PriorList_teta(:), PriorList_RemnantSigma(:)
real(mrk), intent(in),optional:: WLS_sdev(:)
real(mrk), intent(in)::teta0(:),RemnantSigma0
real(mrk), intent(in)::teta_std0(:),RemnantSigma_std0
integer(mik), intent(in)::nAdapt,nCycles
real(mrk), intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in)::OutFile
! OUTPUTS
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,n,p,nteta
real(mrk), allocatable:: WLS(:)
real(mrk), allocatable:: start(:), startStd(:)
character(100), allocatable:: headers(:)
real(mrk)::fx

! Init
err=0;mess=''
n=size(Y);p=size(X,dim=2)
CALL GetRegParNumber(RegressionFunk=RegFunk,covariate=X(1,:),Xtra=Xtra,npar=nteta,err=err,mess=mess)
if(err>0) then
    mess="RegBay_Fit: Fatal: "//trim(mess);return
endif

! Optional WLS_sdev varable
if(allocated(WLS)) deallocate(WLS)
if(present(WLS_sdev)) then
    allocate(WLS(size(WLS_sdev)))
    WLS=WLS_sdev
else
    allocate(WLS(n))
    WLS=0._mrk
endif

! Size checks
if( size(X,dim=1)/=n .or. size(PriorList_teta)/=nteta .or. size(PriorList_RemnantSigma)/=1 &
    .or. size(WLS)/=n) then
    err=1;mess='RegBay_Fit: Fatal: size mismatch in input data';return
endif

! Populate Regression Model object
MoMo%RegFunk=RegFunk ! Regression formula - see RegressionFunk.f90 for the catalogue
MoMo%Link=Link ! Link function - see RegressionFunk.f90 for the catalogue
MoMo%nobs=n ! number of observations used to fit the regression
MoMo%nteta=nteta ! number of regression parameters
if(allocated(MoMo%Y)) deallocate(MoMo%Y)
allocate(MoMo%Y(n))
MoMo%Y=Y ! predictand
if(allocated(MoMo%X)) deallocate(MoMo%X)
allocate(MoMo%X(n,p))
MoMo%X=X ! predictors
if(allocated(MoMo%WLS_sdev)) deallocate(MoMo%WLS_sdev)
allocate(MoMo%WLS_sdev(n))
MoMo%WLS_sdev=WLS ! standard deviations of individuals yi's (for WLS-like inferrence)
if(allocated(MoMo%PriorList_RemnantSigma)) deallocate(MoMo%PriorList_RemnantSigma)
allocate(MoMo%PriorList_RemnantSigma(1))
MoMo%PriorList_RemnantSigma=PriorList_RemnantSigma ! priors for remnant std - see BayesianEstimation_tools module for PriorListType definition
if(allocated(MoMo%PriorList_teta)) deallocate(MoMo%PriorList_teta)
allocate(MoMo%PriorList_teta(nteta))
MoMo%PriorList_teta=PriorList_teta ! priors for regression parameters
if(present(Xtra)) then
    if(allocated(MoMo%Xtra)) deallocate(MoMo%Xtra)
    allocate(MoMo%Xtra(size(Xtra)))
    MoMo%Xtra=Xtra
endif

! MCMC sampling
if(allocated(start)) deallocate(start)
allocate(start(MoMo%nteta+1))
if(allocated(startStd)) deallocate(startStd)
allocate(startStd(MoMo%nteta+1))
if(allocated(headers)) deallocate(headers)
allocate(headers(MoMo%nteta+1+1))
start(1:MoMo%nteta+1)=(/teta0,RemnantSigma0/)
startStd(1:MoMo%nteta+1)=(/teta_std0,RemnantSigma_std0/)
do i=1,MoMo%nteta
    headers(i)="Teta"//trim(number_string(i))
enddo
headers(MoMo%nteta+1)="RemnantSTD"
headers(MoMo%nteta+1+1)="LogPost"

call Adaptive_Metro_OAAT(f=Posterior_wrapper,x=start,&
				fx=fx,std=startStd,&
				nAdapt=nAdapt,nCycles=nCycles,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				OutFile=OutFile,headers=headers, err=err,mess=mess)		
if(err>0) then
    mess="RegBay_Fit: "//trim(mess);return
endif

end subroutine RegBay_Fit


subroutine RegBay_Residuals(Y,X,& ! observations for predictand (Y) and predictors (X)
                       RegFunk,Link,& ! Regression function, link function
                       Xtra,& ! optional arguments to pass to the regression funk
                       RegPar,& ! Par used to apply regression (e.g. modal)
                       WLS_sdev,& ! optional vector (1:n) containing standard deviations of individuals yi's (for WLS-like inferrence). Default is zero (ie SLS) 
                       Res,& ! OUT: residuals=link(y)-regression(x|RegPar)
                       StandardizedRes,& ! optional, OUT: standardized residuals = (y-InverseLink(regression(x|RegPar)))/sqrt(RemnantSigma**2+WLS**2)
                       InverseLinkRes,& ! optional, OUT: transformed residuals = y-InverseLink(regression(x|RegPar))
                       ResFile,&! optional, adress to write result file
					   err, mess)

!^**********************************************************************
!^* Purpose: Residual Analysis
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 31/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.Y
!^*		2.X
!^*		3.RegFunk
!^*		4.Link
!^*		5.[Xtra]
!^*		6.RegPar
!^*		7.WLS_sdev
!^*		8.ResFile
!^* OUT
!^*		1.Res
!^*		2.[StandardizedRes]
!^*		3.[InverseLinkRes]
!^*		4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use RegressionFunk

real(mrk), intent(in)::Y(:),X(:,:)
character(*),intent(in)::RegFunk, Link
real(mrk), intent(in), optional::Xtra(:),WLS_sdev(:)
real(mrk), intent(in)::regPar(:)
character(*),intent(in),optional::ResFile
real(mrk), intent(out)::Res(size(Y))
real(mrk), intent(out), optional::StandardizedRes(size(Y)),InverseLinkRes(size(Y))
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::npred,npar,nteta,unt,ncov,i
real(mrk)::yhat(size(Y)),LY(size(Y)),ILyhat(size(Y))
real(mrk), allocatable::WLS(:)
logical::feas_L,feas_IL
real(mrk)::InvLinkRes(size(Y)),StdRes(size(Y))

!init
err=0;mess='';Res=undefRN
if(present(StandardizedRes)) StandardizedRes=undefRN
if(present(InverseLinkRes)) InverseLinkRes=undefRN

! Optional WLS_sdev varable
if(allocated(WLS)) deallocate(WLS)
if(present(WLS_sdev)) then
    allocate(WLS(size(WLS_sdev)))
    WLS=WLS_sdev
else
    allocate(WLS(size(Y)))
    WLS=0._mrk
endif

! get number of parameters
CALL GetRegParNumber(RegressionFunk=RegFunk,covariate=X(1,:),Xtra=Xtra,npar=nteta,err=err,mess=mess)
npar=nteta+1
if(err>0) then
    mess="RegBay_Residuals: "//trim(mess);return
endif
npred=size(Y)
ncov=size(X,dim=2)

! get residuals
do i=1,npred
    ! apply regression
    call ApplyRegressionFunk(RegressionFunk=RegFunk,covariate=X(i,:),&
                    regressionPar=regPar(1:nteta),Xtra=Xtra,out=yhat(i),err=err,mess=mess)
    if(err>0) then
        mess="RegBay_Residuals: "//trim(mess);return
    endif
    ! Apply link to Y
    call LinkFunk(LFunk=Link,x=Y(i),fx=LY(i),feas=feas_L,err=err,mess=mess)
    if(err>0) then
        mess="RegBay_Residuals: "//trim(mess);return
    endif
    ! apply Inverse Link to regression
    call InverseLinkFunk(IvsLinkFunk=Link,x=yhat(i),fx=ILyhat(i),feas=feas_IL,err=err,mess=mess)
    if(err>0) then
        mess="RegBay_Residuals: "//trim(mess);return
    endif
    ! Compute Residual
    if(.not.feas_L) then
        Res(i)=undefRN
        StdRes(i)=undefRN
    else
        Res(i)=LY(i)-yhat(i)
        StdRes(i)=(LY(i)-yhat(i))/sqrt(regPar(nteta+1)**2+WLS(i)**2)
    endif
    if(.not.feas_IL) then
        InvLinkRes(i)=undefRN
    else
        InvLinkRes(i)=Y(i)-ILyhat(i)
    endif
    if(present(StandardizedRes)) StandardizedRes(i)=StdRes(i)
    if(present(InverseLinkRes)) InverseLinkRes(i)=InvLinkRes(i)
enddo

! Write 2 file if requested
if(present(resfile)) then
    call getSpareUnit(unt,err,mess)
    if(err/=0) then
        mess="RegBay_Residuals: "//trim(mess);return
    endif
    open(unit=unt,file=ResFile,status="REPLACE",iostat=err)
    if(err/=0) then
        mess="RegBay_Residuals: problem opening result file";return
    endif
    write(unt,'(<7+1>A15)') "Y","LinkY","Yhat","InvLinkYhat","Res","StdRes","InvLinkRes","covariates-->"
    do i=1,npred
        write(unt,'(<7+ncov>e15.6)') Y(i),LY(i),yhat(i),ILyhat(i),Res(i),StdRes(i),InvLinkRes(i),X(i,:)
    enddo
    close(unt)
endif

end subroutine RegBay_Residuals


subroutine RegBay_Predict(File,Nburn,Nread,Nslim,Nrep,& ! IN: Read properties
                       X,& ! Matrix of covariates where prediction is made
                       RegFunk,Link,& ! Regression function, link function
                       Xtra,& ! optional arguments to pass to the regression funk
                       PropagateTetaU,PropagateSigmaU,& !IN, options for computation of predictive, default true 
                       modalPar,modalPred,& ! OUT: modal parameter estimates & correponding predictions
                       predictive,& !OUT: replicates from the predictive distribution
                       MarginalLkh,& !OUT: marginal likelihood (crude) approximation
					   err, mess)
!^**********************************************************************
!^* Purpose: Use the regression to predict
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 01/02/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Think about propagating WLS uncertainty in predictive???
!^**********************************************************************
!^* IN
!^*		1. File, containing MCMC samples
!^*		2. Nburn, number of lines to discard
!^*		3. Nread, number of lines to read
!^*		4. Nslim, only one line every Nslim will be used
!^*		5. Nrep, for each line, Nrep deviates will be generated
!^*(advise: keep Nrep ~ Nslim - Nrep>>Nslim will yield a huge number of replicates!)
!^*		6. X, Matrix (npred x ncov) of covariates where prediction is made (1 prediction per line)
!^*		7. RegFunk, regression function (see RegressionFunk.f90)
!^*		8. Link, link function (see RegressionFunk.f90)
!^*		9. [Xtra],optional arguments to pass to the regression funk)
!^*		10. [PropagateTetaU], propagate teta uncertainty when computing predictive? default TRUE
!^*		11. [PropagateSigmaU], propagate remnant sigma uncertainty when computing predictive? default TRUE
!^* OUT
!^*		1.modalPar, modal parameter estimate
!^*		2.modalPred, prediction corresponding to modal parameters (vector 1 x npred)
!^*		3.predictive, deviates from the predictive (matrix npred x ndeviates)
!^*		4.MarginalLkh, marginal likelihood (crude) approximation (in log space)
!^*		5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		6.mess, error message
!^**********************************************************************
use RegressionFunk, only:GetRegParNumber,ApplyRegressionFunk,InverseLinkFunk
use Distribution_tools, only: GenerateSample, GAUSS 
use numerix_dmsl_kit, only:getmeanvar,mnormal_logp_sub

integer(mik), intent(in)::Nburn, Nread, Nslim, Nrep
character(*),intent(in)::File, RegFunk, Link
real(mrk), intent(in)::X(:,:)
real(mrk), intent(in), optional::Xtra(:)
logical, intent(in), optional::PropagateTetaU,PropagateSigmaU
real(mrk), intent(out)::modalPar(:), modalPred(size(X,dim=1)),MarginalLkh
real(mrk), pointer::predictive(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::nteta,npar,i,j,k,n,errcode,npred,npredi
real(mrk), allocatable::Temp(:,:), mcmc(:,:),mu(:),sigma(:,:),var_invrs(:,:)
real(mrk):: ml(1),yhat,gen(nrep),gen2(nrep)
logical::feas,PropTetaU,PropSigmaU
real(mrk)::numerator,denominator,logDet

! Init
err=0;mess='';MarginalLkh=undefRN;modalPred=undefRN;modalPar=undefRN;
if(present(PropagateTetaU)) then
    PropTetaU=PropagateTetaU
else
    PropTetaU=.true.
endif
if(present(PropagateSigmaU)) then
    PropSigmaU=PropagateSigmaU
else
    PropSigmaU=.true.
endif

CALL GetRegParNumber(RegressionFunk=RegFunk,covariate=X(1,:),Xtra=Xtra,npar=nteta,err=err,mess=mess)
npar=nteta+1
if(err>0) then
    mess="RegBay_Predict: "//trim(mess);return
endif
npred=size(X,dim=1)

!***************************
!********* STEP 1: READ MCMC
!***************************
! Read file
open(unit=1,file=File,status='old')
read(1,*) !headers
do i=1,Nburn ! burn
	read(1,*,iostat=errcode)
	if(errcode/=0) then ! problem reading file
		err=1
		mess='RegBay_Predict:FATAL:problem reading file during burnin - cannot proceed!'
		return
	endif
enddo
if(allocated(Temp)) deallocate(Temp)
allocate(Temp(Nread,npar+1))
do i=1,Nread ! read used lines
	k=i
	read(1,*,iostat=errcode) Temp(i,:)
	if(errcode/=0) then ! problem reading file
		err=-1;k=k-1
		mess='RegBay_Predict:WARNING:problem reading file before [Nread] lines!'
		EXIT
	endif
enddo
close(1)
! Slim
n=k/Nslim
if(allocated(mcmc)) deallocate(mcmc)
allocate(mcmc(n,Npar+1))
mcmc=Temp(1:n:nSlim,:)

!****************************************************
!********* STEP 2: Marginal likelihood approximation
!****************************************************
! get modal estimate
ml=Maxloc(mcmc(:,npar+1))
if(size(modalPar)/=npar) then
	err=1; mess='RegBay_Predict:SEVERE:size mismath [modalPar]'
	modalPar=undefRN
else
	modalPar=mcmc(ml(1),1:npar)
endif

numerator=mcmc(ml(1),npar+1) ! unnormalized log-posterior
! Gaussian approximation of the posterior
! Allocate
if(allocated(mu)) deallocate(mu)
allocate(mu(npar))
if(allocated(sigma)) deallocate(sigma)
allocate(sigma(npar,npar))
if(allocated(var_invrs)) deallocate(var_invrs)
allocate(var_invrs(npar,npar))

! Compute posterior mean/covar from MCMC samples
call getmeanvar(x=mcmc(:,1:npar),package=2,mean=mu,covar=sigma,method="c",err=err,message=mess)
if(err>0) then
    mess="RegBay_Predict: "//trim(mess);return
endif
! Compute gaussian log-pdf at mode
call mnormal_logp_sub(x=modalPar,mean=mu,var=sigma,do_invrs=.true.,&
            var_invrs=var_invrs,logDet=logDet,posdef=feas,logp=denominator,&
            err=err,message=mess)
if(err>0) then
    write(*,*) "RegBay_Predict: Problem while computing Marginal Likelihood"
    write(*,*) "Error is non-fatal, execution will continue"
else ! get MarginalLkh
    if(feas) then
        MarginalLkh=numerator-denominator    
    else
        MarginalLkh=undefRN    
    endif
endif 
!***************************
!********* STEP 3: MODAL prediction
!***************************
! get modal prediction
do i=1,npred
    ! apply regression
    call ApplyRegressionFunk(RegressionFunk=RegFunk,covariate=X(i,:),&
                    regressionPar=modalPar(1:nteta),Xtra=Xtra,out=yhat,err=err,mess=mess)
    if(err>0) then
        mess="RegBay_Predict: "//trim(mess);return
    endif
    ! apply Inverse Link
    call InverseLinkFunk(IvsLinkFunk=Link,x=yhat,fx=modalPred(i),feas=feas,err=err,mess=mess)
    if(err>0) then
        mess="RegBay_Predict: "//trim(mess);return
    endif
    if(.not.feas) modalPred(i)=undefRN
enddo

!***************************
!********* STEP 4: predictive
!***************************
npredi=n*nrep
if(associated(predictive)) deallocate(predictive)
allocate(predictive(npred,npredi))
do i=1,npred
    do j=1,n
        ! apply regression
        if(PropTetaU) then ! teta from MCMC samples
            call ApplyRegressionFunk(RegressionFunk=RegFunk,covariate=X(i,:),&
                        regressionPar=mcmc(j,1:nteta),Xtra=Xtra,out=yhat,err=err,mess=mess)
        else ! teta fixed to modal estimate
            call ApplyRegressionFunk(RegressionFunk=RegFunk,covariate=X(i,:),&
                        regressionPar=modalPar(1:nteta),Xtra=Xtra,out=yhat,err=err,mess=mess)
        endif
        if(err>0) then
            mess="RegBay_Predict: "//trim(mess);return
        endif
        ! Generate deviates
        if(PropSigmaU) then ! noisify regression
           call GenerateSample(DistId=GAUSS,par=(/yhat,mcmc(j,nteta+1)/),&
                            gen=gen,feas=feas,err=err,mess=mess)
        else ! do not add noise to regression
            gen=yhat
        endif
        if(err>0) then
            mess="RegBay_Predict: "//trim(mess);return
        endif
        ! Inverse-link transformation
        do k=1,nrep
            call InverseLinkFunk(IvsLinkFunk=Link,x=gen(k),fx=gen2(k),feas=feas,err=err,mess=mess)
            if(err>0) then
                mess="RegBay_Predict: "//trim(mess);return
            endif
        enddo
        predictive(i,((j-1)*nrep+1):(j*nrep))=gen2
    enddo
enddo

end subroutine RegBay_Predict

!!!!!!!!!!!!!!!!
! Private Subs !
!!!!!!!!!!!!!!!!

subroutine GetLogPost(teta,RemnantSigma,& ! inferred quantities: regression parameters + remnant sdev
                      Y,X,& ! Y(1:n)=predictand, X(1:n,1:p)=predictor
                      RegFunk,Link,& ! Regression function, link function- Link(Y)=RegFunk(X|teta)+eps
                      Xtra,& ! optional arguments to pass to the regression funk
                      WLS,& ! vector (1:n) containing standard deviations of individuals yi's (for WLS-like inferrence). Default is zero (ie SLS) 
                      PriorList_teta,PriorList_RemnantSigma,& ! Priors
                      lp,& ! log-posterior 
                      feas, isnull,err,mess) ! flags & messages

!^**********************************************************************
!^* Purpose: compute post(teta,RemnantSigma)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 17/05/2011
!^**********************************************************************
!^* Comments: Missing data forbidden for now
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. teta, regression parameters
!^*		2. RemnantSigma, std of remnant errors
!^*		3. Y, predictand
!^*		4. X, predictor
!^*		5. RegFunk, regression function (see RegressionFunk.f90)
!^*		6. Link, link function (see RegressionFunk.f90)
!^*		7. Xtra, optional arguments to pass to the regression funk
!^*		8. WLS, standard deviations of individuals yi's. Default is zero (ie SLS) 
!^*		9. PriorList_teta, priors for teta
!^*		10. PriorList_RemnantSigma, prior for remnant error
!^* OUT
!^*		1.lp, log-posterior
!^*		2.feas, feasible?
!^*		3.is null, is (natural) posterior = zero?
!^*		4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5.mess, error message
!^**********************************************************************
use BayesianEstimation_tools, only:GetLogPrior, PriorListType
use Distribution_tools, only: GetPdf, GAUSS
use RegressionFunk, only:ApplyRegressionFunk,LinkFunk

real(mrk), intent(in)::teta(:),RemnantSigma,Y(:),X(:,:)
real(mrk), intent(in):: WLS(:)
real(mrk), intent(in), optional::Xtra(:)
Type(PriorListType), intent(in)::PriorList_teta(:), PriorList_RemnantSigma(:)
character(*), intent(in)::RegFunk,Link
real(mrk), intent(out)::lp
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::n,p !(number of obs, number of covariates)
integer(mik)::nteta ! number of regression parameters
real(mrk)::lkh,lik,prior,pr,yhat0,yhat
integer(mik)::i

!Init
err=0;mess='';feas=.true.;isnull=.false.;lp=undefRN
n=size(Y);p=size(X,dim=2);nteta=size(teta)

! Size checks
if( size(X,dim=1)/=n .or. size(PriorList_teta)/=nteta .or. size(PriorList_RemnantSigma)/=1 &
    .or. size(WLS)/=n) then
    err=1;mess='GetLogPost: Fatal: size mismatch in input data';feas=.false.;return
endif
! check Remnant Sigma
if(RemnantSigma<=0._mrk) then 
    feas=.false.;return
endif

! compute log-likelihood                      
lik=0._mrk
do i=1,n
    ! apply regression function
    call ApplyRegressionFunk(RegressionFunk=RegFunk,covariate=X(i,:),&
                            regressionPar=teta,Xtra=Xtra,out=yhat0,err=err,mess=mess)
    if(err>0) then
        mess="GetLogPost: "//trim(mess);feas=.false.;return
    endif
    ! apply link
    call LinkFunk(LFunk=Link,x=Y(i),fx=yhat,feas=feas,err=err,mess=mess)
    if(err>0) then
        mess="GetLogPost: "//trim(mess);feas=.false.;return
    endif
    if(.not. feas) return
    ! compute contribution to the likelihood
    call Getpdf(DistId=GAUSS,x=yhat,par=(/yhat0,sqrt(WLS(i)**2 + RemnantSigma**2)/),loga=.true.,&
                pdf=lkh,feas=feas,isnull=isnull,err=err,mess=mess)
    if(err>0) then
        mess="GetLogPost: "//trim(mess);feas=.false.;return
    endif
    if(.not. feas) return
    if(isnull) return
    lik=lik+lkh
enddo

! compute log-prior 
prior=0._mrk
! prior for Remnant sdev
call GetLogPrior(teta=(/RemnantSigma/),PriorList=PriorList_RemnantSigma,&
                lp=pr, feas=feas, isnull=isnull,err=err,mess=mess)
If(err>0) then
	mess='GetLogPost: '//trim(mess)
	feas=.false.;return
endif
if ( (.not. feas) .or. (isnull) ) return
prior=prior+pr
! prior for teta
call GetLogPrior(teta=teta,PriorList=PriorList_teta,&
                lp=pr, feas=feas, isnull=isnull,err=err,mess=mess)
If(err>0) then
	mess='GetLogPost: '//trim(mess)
	feas=.false.;return
endif
if ( (.not. feas) .or. (isnull) ) return
prior=prior+pr

! Log-post    
lp = prior + lik
                 
end subroutine GetLogPost


subroutine Posterior_wrapper(x,feas,isnull,fx,fAux,err,mess)

!^**********************************************************************
!^* Purpose: wrapper to GetLogPost to comply with MCMC interface
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:17/05/2011
!^**********************************************************************
!^* Comments: mv handled out of here
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.x, teta
!^* OUT
!^*		1.feas, is x feasible?
!^*		2.isnull, is fx=0?
!^*		3.fx (WARNING: LOG-posterior)
!^*		4.[fAux], by-product from the posterior computation that may be usefull
!^*		5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		6.mess, error message
!^**********************************************************************
real(mrk),intent(in)::x(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess

if(allocated(MoMo%Xtra)) then
    call GetLogPost(teta=x(1:MoMo%nteta),RemnantSigma=x(MoMo%nteta+1),& ! inferred quantities: regression parameters + remnant sdev
                      Y=MoMo%Y,X=MoMo%X,& ! Y(1:n)=predictand, X(1:n,1:p)=predictor
                      RegFunk=MoMo%RegFunk,Link=MoMo%Link,& ! Regression function, link function- Link(Y)=RegFunk(X|teta)+eps
                      Xtra=MoMo%Xtra,&
                      WLS=MoMo%WLS_sdev,& ! optional vector (1:n) containing standard deviations of individuals yi's (for WLS-like inferrence). Default is zero (ie SLS) 
                      PriorList_teta=MoMo%PriorList_teta,PriorList_RemnantSigma=MoMo%PriorList_RemnantSigma,& ! Priors
                      lp=fx,& ! log-posterior 
                      feas=feas, isnull=isnull,err=err,mess=mess) ! flags & messages
else
    call GetLogPost(teta=x(1:MoMo%nteta),RemnantSigma=x(MoMo%nteta+1),& ! inferred quantities: regression parameters + remnant sdev
                      Y=MoMo%Y,X=MoMo%X,& ! Y(1:n)=predictand, X(1:n,1:p)=predictor
                      RegFunk=MoMo%RegFunk,Link=MoMo%Link,& ! Regression function, link function- Link(Y)=RegFunk(X|teta)+eps
                      WLS=MoMo%WLS_sdev,& ! optional vector (1:n) containing standard deviations of individuals yi's (for WLS-like inferrence). Default is zero (ie SLS) 
                      PriorList_teta=MoMo%PriorList_teta,PriorList_RemnantSigma=MoMo%PriorList_RemnantSigma,& ! Priors
                      lp=fx,& ! log-posterior 
                      feas=feas, isnull=isnull,err=err,mess=mess) ! flags & messages
endif
                      
if(present(fAux)) fAux=UndefRN ! not handled for now

end subroutine Posterior_wrapper


end module RegBay_tools
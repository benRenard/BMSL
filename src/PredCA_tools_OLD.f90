module PredCA_tools

!~**********************************************************************
!~* Purpose: Implement Predictive Component Analysis
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified:13/15/2013
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
!~*		4. 
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: DoMCMC

!-------------------------------------------------------------------------------------------
! Definition of parameters used by PredCA
!-------------------------------------------------------------------------------------------
! Various inference types 
Character(3), parameter, PUBLIC:: FIX='FIX',& ! parameter is fixed to its default value
								  REG='REG',& ! parameter is inferred and assumed regional
                                  LOC='LOC',& ! parameter is inferred and assumed local
								  LAT='LAT'   ! parameter is inferred and considered as a latent variable in a hierachical model

!-------------------------------------------------------------------------------------------
! Definition of types
!-------------------------------------------------------------------------------------------
! Data Type
type, public::PredCA_DataType
    character(250)::nickname='AintGotNoName' ! name of the dataset
    integer(mik)::Nt=undefIN, Ns=undefIN     ! dimension (time*space) of data matrix
    real(mrk), allocatable::val(:,:)         ! values, size Nt*Ns
	real(mrk)::mv=-99._mrk                   ! value denoting missing data
	real(mrk)::typscale=undefRN              ! typical scale (i.e. order of magnitude) of the data
    real(mrk), allocatable::distance(:,:)    ! Matrix of spatial distances, size Ns*Ns
	character(250)::file='undefined'         ! file where data has been read
	character(250)::grdfile_t='undefined'    ! grd file for temporal grid
	character(250)::grdfile_s='undefined'    ! grd file for spatial grid
	character(250)::grdfile_f='undefined'    ! grd file for field grid
    real(mrk), allocatable::tgrid(:)         ! Temporal grid, size Nt
    real(mrk), allocatable::sgrid(:,:)       ! Spatial grid, size Ns*2
    integer(mik), allocatable::fgrid(:)      ! field grid, size Ns
end type PredCA_DataType

! Model for each inferred parameter
type, public::PredCA_InfType
    character(250)::nickname='AintGotNoName'    ! name of the inferred parameter
    real(mrk)::val=undefRN                      ! value
    character(100)::priorDist='undefined'       ! prior distribution
    real(mrk), allocatable::priorpar(:)         ! prior parameters
    integer(mik)::indx=undefIN                  ! index in the vector of all inferred parameters
    real(mrk)::JumpStd=undefRN                  ! initial value for the Gaussian Jump distribution
end type PredCA_InfType
! vectorial version
type, public::PredCA_InfVType
    type(PredCA_InfType),allocatable::par(:)    ! individual inferred parameters
    character(3)::typ                           ! type (FIX,REG,LOC,LAT)
end type PredCA_InfVType

! Model for each D-parameter
type, public::PredCA_DparType
    character(250)::nickname='AintGotNoName'            ! name of the D-par
	integer(mik)::Nk=undefIN                            ! number of components
    real(mrk)::defval=undefRN                           ! default value of the D-par
	character(250)::link='undefined'                    ! link function
    !character(3),allocatable::Lambda_Type(:)            ! Types for lambdas (size Nk+1). 'FIX' only allowed for lambda0
    !character(3),allocatable::Psi_Type(:)               ! Types for Psi (size Nk). 'LOC', 'REG' or 'LAT'
    !character(3),allocatable::Tau_Type(:)               ! Types for Tau (size Nk). 'LOC' or 'LAT'
    type(PredCA_InfVType), allocatable::lambda(:)       ! lambdas (size Nk+1)
    type(PredCA_InfVType), allocatable::psi(:)          ! psis (size Nk)
    type(PredCA_InfVType), allocatable::tau(:)          ! taus (size Nk)
    type(PredCA_InfType), allocatable::sigma_tau(:)     ! hyper-sdev of taus (if 'LAT') (size Nk)
    ! NOTE: need to add spatial variances when lanmbdas and psis are 'LAT' - not implemented yet
end type PredCA_DparType

! Information on the PredCA model (i.e. how Hyd & Clim are linked)
type, public::PredCA_Model
    character(250)::nickname='AintGotNoName'         ! name of the model
    character(100)::dist='undefined'                 ! Y: distribution
    integer(mik)::Nd=UndefIN                         ! Y: number of D-par
	type(PredCA_DparType), allocatable::Dpar(:)      ! Size Nd. properties of each Dpar 
	character(250), allocatable::DparFile(:)         ! Size Nd. Config file for each Dpar 
	logical::DoCenter=.true.,DoScale=.false.         ! Phi: Center and scale data?
    !character(3)::Mu_Phi_Type=LOC                    ! Phi: Use local or regional mean for centering?
    !character(3)::Sigma_Phi_Type=LOC                 ! Phi: Use local or regional sdev for scaling?
    type(PredCA_InfVType)::Mu_Phi                    ! Mean of Phi 
    type(PredCA_InfVType)::Sigma_Phi                 ! Sdev of Phi 
    type(PredCA_InfType)::sigma_eps                  ! sdev of Phi-residuals
end type PredCA_Model

! All required info for a case study
type, public::PredCA_CSType
    character(250)::nickname='AintGotNoName' ! name of the case study
    integer(mik)::Ninf=undefIN               ! number of parameters to be inferred
    integer(mik)::Nlat=undefIN               ! number of latent variables
    integer(mik)::NnonLat=undefIN            ! number of non-latent variables, ie pars+hyperpars
    integer(mik)::Nt=undefIN                 ! number of time steps
    integer(mik)::Nd=undefIN                 ! number of Dpar
    integer(mik)::Nx=undefIN                 ! number of columns in Y
    integer(mik)::Ns=undefIN                 ! number of columns in Phi
	type(PredCA_DataType)::Y                 ! Y data
	type(PredCA_DataType)::Phi               ! Phi data
	type(PredCA_Model)::mdl                  ! PredCA model
	logical::AnyLAT=.false.                  ! Any latent variable in the inference?
	character(250)::OUTfolder='undefined'    ! folder where results will be saved
	character(250)::cas_file='undefined'     ! file where case study spec is read
	character(250)::mdl_file='undefined'     ! file where model spec is read
	character(250)::dat_fileY='undefined'    ! file where Y spec is read
	character(250)::dat_filePhi='undefined'  ! file where Phi spec is read

end type PredCA_CSType            
! Object made globally available to this module
type(PredCA_CSType),target::Fred

contains

subroutine DoMCMC(file,priorOption,initOption,err,mess)
use MCMCStrategy_tools,only:MetroMS_Mix2
character(*),intent(in)::file,priorOption,initOption
integer(mik),intent(out)::err
character(*),intent(out)::mess
!~~~~~~~~~~~~~~~~~~~~
integer(mik),parameter::nAdapt=20,nCycles=100,nSim=10000,GRBlocSize=100,n_chain=4
real(mrk), parameter::BurnFactor=0.5_mrk,MinMoveRate=0.1_mrk,MaxMoveRate=0.3_mrk,DownMult=0.9_mrk,UpMult=1.1_mrk
integer(mik),dimension(4)::ShowPar=(/1,11,12,13/)
!~~~~~~~~~~~~~~~~~~~~
real(mrk),allocatable::start(:),startstd(:)
real(mrk),pointer::fx(:)

call LoadFred(file,priorOption,initOption,err,mess)
if(err/=0) then;mess='DoMCMC:'//trim(mess);return;endif
allocate(start(Fred%Ninf),startstd(Fred%Ninf))
call UnfoldStructure2Vector(val=start,JumpStd=startstd,err=err,mess=mess)
if(err/=0) then;mess='DoMCMC:'//trim(mess);return;endif
call MetroMS_Mix2(f=GetLogPost_full,x=start,fx=fx,std=startstd,&
                nAdapt1=nAdapt,nCycles1=nCycles,BurnFactor=BurnFactor,&
                nAdapt2=nAdapt,nCycles2=nCycles,&
                nSim=nSim,GRBlocSize=GRBlocSize,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				n_chain=n_chain,&
				UseRFortran=.true.,&
				showPar=showPar,&
				GRBurnfactor=BurnFactor,&
				!OutDoc,Filename1,Filename2,Filename3,&
				!OutMat1, OutMat2, OutMat3,&
				!headers,&
				err=err,mess=mess)
if(err/=0) then;mess='DoMCMC:'//trim(mess);return;endif
end subroutine DoMCMC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine UnfoldStructure2Vector(val,JumpStd,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:25/05/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		3. 
!^*		4. 
!^*		5. 
!^*		6. 
!^*		7. 
!^*		8. 
!^*		9. 
!^* OUT
!^*		1. [val], values
!^*		2. [JumpStd], stds of Jump distributions
!^*		3.err
!^*		4.mess
!^**********************************************************************
real(mrk), intent(out),optional::val(:),JumpStd(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::d,k,n,c
type(PredCA_DparType), pointer::dpar

err=0;mess='';n=0
if(present(val)) val=undefRN
if(present(JumpStd)) JumpStd=undefRN

do d=1,Fred%Nd ! loop on all D-par
    dpar=>Fred%mdl%Dpar(d) ! define shortcut
    !lambdas
    do k=1,dpar%Nk+1
        if(dpar%lambda(k)%typ/=FIX) then
            c=size(dpar%lambda(k)%par)
            if(present(val)) val((n+1):(n+c))=dpar%lambda(k)%par(:)%val 
            if(present(JumpStd)) JumpStd((n+1):(n+c))=dpar%lambda(k)%par(:)%JumpStd
            n=n+c
        endif 
        if(dpar%lambda(k)%typ==LAT) then
            err=1;mess='UnfoldStructure2Vector:Fatal:spatial latents not implemented'
        endif      
    enddo
    !tau
    do k=1,dpar%Nk
        if(dpar%tau(k)%typ/=FIX) then
            c=size(dpar%tau(k)%par)
            if(present(val)) val((n+1):(n+c))=dpar%tau(k)%par(:)%val
            if(present(JumpStd)) JumpStd((n+1):(n+c))=dpar%tau(k)%par(:)%JumpStd
            n=n+c
        endif 
        if(dpar%tau(k)%typ==LAT) then
            c=1
            if(present(val)) val((n+1):(n+c))=dpar%sigma_tau(k)%val
            if(present(JumpStd)) JumpStd((n+1):(n+c))=dpar%sigma_tau(k)%JumpStd
            n=n+c
        endif      
    enddo
    !psi
    do k=1,dpar%Nk
        if(dpar%psi(k)%typ/=FIX) then
            c=size(dpar%psi(k)%par)
            if(present(val)) val((n+1):(n+c))=dpar%psi(k)%par(:)%val
            if(present(JumpStd)) JumpStd((n+1):(n+c))=dpar%psi(k)%par(:)%JumpStd
            n=n+c
        endif 
        if(dpar%psi(k)%typ==LAT) then
            err=1;mess='UnfoldStructure2Vector:Fatal:spatial latents not implemented'
        endif      
    enddo
enddo

! Phi mean
if(Fred%mdl%DoCenter) then
    if(Fred%mdl%Mu_Phi%typ/=FIX) then
        c=size(Fred%mdl%Mu_Phi%par)
        if(present(val)) val((n+1):(n+c))=Fred%mdl%Mu_Phi%par(:)%val
        if(present(JumpStd)) JumpStd((n+1):(n+c))=Fred%mdl%Mu_Phi%par(:)%JumpStd
        n=n+c
    endif 
    if(Fred%mdl%Mu_Phi%typ==LAT) then
        err=1;mess='UnfoldStructure2Vector:Fatal:spatial latents not implemented'
    endif      
endif

! Phi Sdev
if(Fred%mdl%DoScale) then
    if(Fred%mdl%Sigma_Phi%typ/=FIX) then
        c=size(Fred%mdl%Sigma_Phi%par)
        if(present(val)) val((n+1):(n+c))=Fred%mdl%Sigma_Phi%par(:)%val
        if(present(JumpStd)) JumpStd((n+1):(n+c))=Fred%mdl%Sigma_Phi%par(:)%JumpStd
        n=n+c
    endif 
    if(Fred%mdl%Sigma_Phi%typ==LAT) then
        err=1;mess='UnfoldStructure2Vector:Fatal:spatial latents not implemented'
    endif      
endif

!! Phi residuals
!c=1
!if(present(val)) val((n+1):(n+c))=Fred%mdl%Sigma_eps%val
!if(present(JumpStd)) JumpStd((n+1):(n+c))=Fred%mdl%Sigma_eps%JumpStd

end subroutine UnfoldStructure2Vector      
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine FoldVector2Structure(parvector,priorlist,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:15/05/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. [parvector], vector to fold into structure
!^*		2. [priorList], list of prior distributions
!^*		3. 
!^*		4. 
!^*		5. 
!^*		6. 
!^*		7. 
!^*		8. 
!^*		9. 
!^* OUT
!^*		1.
!^*		2.err
!^*		3.mess
!^**********************************************************************
use BayesianEstimation_tools, only:PriorListType
real(mrk), intent(in),optional::parvector(:)
type(PriorListType), intent(in),optional:: priorlist(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::d,k,n
type(PredCA_DparType), pointer::dpar

err=0;mess='';n=0
do d=1,Fred%Nd ! loop on all D-par
    dpar=>Fred%mdl%Dpar(d) ! define shortcut
    ! lambda 0
    call FoldComponent(comp=dpar%lambda(1),priorlist=priorlist,n=n,parvector=parvector)
    if(dpar%lambda(1)%typ==LAT) then
        err=1;mess='FoldVector2Structure:Fatal:spatial latents not implemented'
    endif      
    ! loop on all components for the current Dpar
    do k=1,dpar%Nk 
        ! Lambda
        call FoldComponent(comp=dpar%lambda(k+1),priorlist=priorlist,n=n,parvector=parvector)
        if(dpar%lambda(k+1)%typ==LAT) then
            err=1;mess='FoldVector2Structure:Fatal:spatial latents not implemented'
        endif      
        ! tau
        call FoldComponent(comp=dpar%tau(k),priorlist=priorlist,n=n,parvector=parvector)
        if(dpar%tau(k)%typ==LAT) then
            n=n+1;if(present(parvector)) dpar%sigma_tau(k)%val=parvector(n)
            if(present(priorlist)) then
                dpar%sigma_tau(k)%priorDist=priorlist(n)%dist
                if(allocated(dpar%sigma_tau(k)%priorPar)) deallocate (dpar%sigma_tau(k)%priorPar)
                allocate(dpar%sigma_tau(k)%priorPar(size(priorlist(n)%par)))
                dpar%sigma_tau(k)%priorPar=priorlist(n)%par
            endif
        endif      
        ! psi
        call FoldComponent(comp=dpar%psi(k),priorlist=priorlist,n=n,parvector=parvector)
        if(dpar%psi(k)%typ==LAT) then
            err=1;mess='FoldVector2Structure:Fatal:spatial latents not implemented'
        endif
    enddo
enddo
! Phi mean
if(Fred%mdl%DoCenter) then
    call FoldComponent(comp=Fred%mdl%Mu_Phi,priorlist=priorlist,n=n,parvector=parvector)
endif
! Phi Sdev
if(Fred%mdl%DoScale) then
    call FoldComponent(comp=Fred%mdl%Sigma_Phi,priorlist=priorlist,n=n,parvector=parvector)
endif
!! Phi residuals
!n=n+1
!if(present(parvector)) Fred%mdl%sigma_eps%val=parvector(n)
!if(present(priorlist)) then
!    Fred%mdl%sigma_eps%priorDist=priorlist(n)%dist
!    if(allocated(Fred%mdl%sigma_eps%priorPar)) deallocate (Fred%mdl%sigma_eps%priorPar)
!    allocate(Fred%mdl%sigma_eps%priorPar(size(priorlist(n)%par)))
!    Fred%mdl%sigma_eps%priorPar=priorlist(n)%par
!endif

end subroutine FoldVector2Structure      
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine FoldComponent(comp,n,parvector,priorlist)
!^**********************************************************************
!^* Purpose: generic sub to fold a component (lambda, psi or tau) into structure
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:15/05/2014
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. comp, component to treat
!^*		2. n, counter
!^*		3. [parvector], vector to fold into structure
!^* OUT
!^*		1.
!^*		2.err
!^*		3.mess
!^**********************************************************************
use BayesianEstimation_tools, only:PriorListType
type(PredCA_InfVType), intent(inout)::comp
integer(mik), intent(inout)::n
real(mrk), intent(in),optional::parvector(:)
type(PriorListType), intent(in),optional:: priorlist(:)
!locals
integer(mik)::x

if(comp%typ/=FIX) then
    do x=1,size(comp%par)
        n=n+1
        if(present(parvector)) comp%par(x)%val=parvector(n)
        if(present(priorlist)) then
            if(comp%typ/=LAT) then
                comp%par(x)%priorDist=priorlist(n)%dist
                if(allocated(comp%par(x)%priorPar)) deallocate (comp%par(x)%priorPar)
                allocate(comp%par(x)%priorPar(size(priorlist(n)%par)))
                comp%par(x)%priorPar=priorlist(n)%par
            endif
        endif
    enddo        
endif      
end subroutine FoldComponent

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetLogPost_full(x,feas, isnull,fx,fAux,err,mess)
!^**********************************************************************
!^* Purpose: compute full log posterior
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:16/05/2014
!^**********************************************************************
!^* Comments: LL associated with the parameter values currently stored
!^*           in structure - folding should already have been done 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.x, value at which log-posterior is evaluated
!^* OUT
!^*		1.feas, log-posterior
!^*		2.isnull
!^*		3.fx,posterior value
!^*		3.[fAux],unused
!^*		4.err
!^*		5.mess
!^**********************************************************************
real(mrk),intent(in)::x(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk)::logprior,loglik, loghyper

err=0;mess='';fx=0._mrk
call FoldVector2Structure(parvector=x,err=err,mess=mess)
if(err/=0) then;mess='GetLogPost_full:Fatal:'//trim(mess);return;endif
call GetLogPrior_full(logprior,feas,isnull,err,mess)
if( (.not. feas) .or. isnull) return
if(err>0) then;mess='GetLogPost_full:Fatal:'//trim(mess);return;endif
call GetLogLik_full(loglik,feas,isnull,err,mess)
if( (.not. feas) .or. isnull) return
if(err>0) then;mess='GetLogPost_full:Fatal:'//trim(mess);return;endif
call GetLogHyper_full(loghyper,feas,isnull,err,mess)
if( (.not. feas) .or. isnull) return
if(err>0) then;mess='GetLogPost_full:Fatal:'//trim(mess);return;endif
fx=logprior+loglik+loghyper 
end subroutine GetLogPost_full

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetLogLik_full(LL,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: compute full LL
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:16/05/2014
!^**********************************************************************
!^* Comments: LL associated with the parameter values currently stored
!^*           in structure - folding should already have been done 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^* OUT
!^*		1.LL, log-likelihood
!^*		2.feas
!^*		3.isnull
!^*		4.err
!^*		5.mess
!^**********************************************************************
use RegressionFunk, only:InverseLinkFunk
use Distribution_tools, only:GetPdf,GAUSS
real(mrk), intent(out)::LL
logical,intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::t,x,d,k,s
type(PredCA_DparType), pointer::dpar
real(mrk)::theta(Fred%Nd)
real(mrk)::dummy,lambda,tau,psi,mu,sig,eps,pdf

err=0;mess='';LL=0._mrk
! Y likelihood
do t=1,Fred%Nt
    do x=1,Fred%Nx
        do d=1,Fred%Nd ! retrieve each D-par
            dpar=>Fred%mdl%dpar(d) ! define shortcut
            dummy=0._mrk ! reinitialise
            ! lambda0
            dummy=dummy+GetCompVal(dpar%lambda(1),x)
            ! other lambdas*tau
            do k=1,dpar%Nk
                lambda=GetCompVal(dpar%lambda(k+1),x)
                tau=GetCompVal(dpar%tau(k),t)
                dummy=dummy+lambda*tau
            enddo
            ! Apply inverse link function and store in theta
            call InverseLinkFunk(dpar%link,dummy,theta(d),feas,err,mess)
            if(.not. feas) return
            if(err>0) then;mess='GetLogLik_full:Fatal:'//trim(mess);return;endif
         enddo
         !Get LL contribution for Y(t,x)
         if(Fred%Y%val(t,x)/=Fred%Y%mv) then
            call GetPdf(DistId=Fred%mdl%dist,x=Fred%Y%val(t,x),par=theta,loga=.true.,&
                        pdf=pdf,feas=feas,isnull=isnull,err=err,mess=mess)    
            if( (.not. feas) .or. isnull) return
            if(err>0) then;mess='GetLogLik_full:Fatal:'//trim(mess);return;endif
            LL=LL+pdf
         endif
    enddo
enddo
! Phi likelihood
do s=1,Fred%Ns
    mu=GetCompVal(Fred%mdl%Mu_Phi,s)
    sig=GetCompVal(Fred%mdl%Sigma_Phi,s)
    do t=1,Fred%Nt
        dummy=0._mrk ! reinitialise
        do d=1,Fred%Nd ! D-loop 
            dpar=>Fred%mdl%dpar(d) ! define shortcut
            ! psi*tau
            do k=1,dpar%Nk
                psi=GetCompVal(dpar%psi(k),s)
                tau=GetCompVal(dpar%tau(k),t)
                dummy=dummy+psi*tau
            enddo
         enddo
         !Get residual for Phi(t,s) and compute its contribution
         if(Fred%Phi%val(t,s)/=Fred%Phi%mv) then
            eps= ((Fred%Phi%val(t,s)-mu)/sig)-dummy
            !call GetPdf(DistId=GAUSS,x=eps,par=(/0._mrk,Fred%mdl%sigma_eps%val/),loga=.true.,&
            !            pdf=pdf,feas=feas,isnull=isnull,err=err,mess=mess)    
            call GetPdf(DistId=GAUSS,x=eps,par=(/0._mrk,1._mrk/),loga=.true.,&
                        pdf=pdf,feas=feas,isnull=isnull,err=err,mess=mess)    
            if( (.not. feas) .or. isnull) return
            if(err>0) then;mess='GetLogLik_full:Fatal:'//trim(mess);return;endif
            LL=LL+pdf
         endif
    enddo
enddo
end subroutine GetLogLik_full

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetLogPrior_full(logprior,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: compute full logprior
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:19/05/2014
!^**********************************************************************
!^* Comments: LogPrior associated with the parameter values currently stored
!^*           in structure - folding should already have been done 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^* OUT
!^*		1.LogPrior, log-prior
!^*		2.feas
!^*		3.isnull
!^*		4.err
!^*		5.mess
!^**********************************************************************
use Distribution_tools, only:GetPdf
real(mrk), intent(out)::logprior
logical,intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::t,x,d,k,s
type(PredCA_DparType), pointer::dpar
type(PredCA_InfType), pointer::ipar
real(mrk)::pdf

err=0;mess='';logprior=0._mrk
do d=1,Fred%Nd
    dpar=>Fred%mdl%dpar(d) ! define shortcut
    ! lambdas
    do k=1,dpar%Nk+1
        if (dpar%lambda(k)%typ==LAT) then
            ! get prior for hyper-parameters
            err=1;mess='GetLogPrior_full:Fatal:spatial latents not implemented'
        else
            do x=1,size(dpar%lambda(k)%par)
                ipar=>dpar%lambda(k)%par(x) ! define shortcut
                call GetCompPrior(ipar,logprior,feas,isnull,err,mess)
                if( (.not. feas) .or. isnull) return
                if(err>0) then;mess='GetLogPrior_full:Fatal:'//trim(mess);return;endif
            enddo
        endif
     enddo
    ! taus
    do k=1,dpar%Nk
        if (dpar%tau(k)%typ==LAT) then
            ! get prior for hyper-sdev
             ipar=>dpar%sigma_tau(k) ! define shortcut
             call GetCompPrior(ipar,logprior,feas,isnull,err,mess)
             if( (.not. feas) .or. isnull) return
             if(err>0) then;mess='GetLogPrior_full:Fatal:'//trim(mess);return;endif
        else
            do t=1,size(dpar%tau(k)%par)
                ipar=>dpar%tau(k)%par(t) ! define shortcut
                call GetCompPrior(ipar,logprior,feas,isnull,err,mess)
                if( (.not. feas) .or. isnull) return
                if(err>0) then;mess='GetLogPrior_full:Fatal:'//trim(mess);return;endif
            enddo
        endif
     enddo
    ! psis
    do k=1,dpar%Nk
        if (dpar%psi(k)%typ==LAT) then
            ! get prior for hyper-parameters
            err=1;mess='GetLogPrior_full:Fatal:spatial latents not implemented'
        else
            do s=1,size(dpar%psi(k)%par)
                ipar=>dpar%psi(k)%par(s) ! define shortcut
                call GetCompPrior(ipar,logprior,feas,isnull,err,mess)
                if( (.not. feas) .or. isnull) return
                if(err>0) then;mess='GetLogPrior_full:Fatal:'//trim(mess);return;endif
            enddo
        endif
     enddo
enddo
! Mu Phi
if(Fred%mdl%DoCenter) then
    if (Fred%mdl%Mu_Phi%typ==LAT) then
        ! get prior for hyper-parameters
        err=1;mess='GetLogPrior_full:Fatal:spatial latents not implemented'
    else
        do s=1,size(Fred%mdl%Mu_Phi%par)
            ipar=>Fred%mdl%Mu_Phi%par(s) ! define shortcut
            call GetCompPrior(ipar,logprior,feas,isnull,err,mess)
            if( (.not. feas) .or. isnull) return
            if(err>0) then;mess='GetLogPrior_full:Fatal:'//trim(mess);return;endif
        enddo
    endif
endif

!Sigma Phi
if(Fred%mdl%DoScale) then
    if (Fred%mdl%Sigma_Phi%typ==LAT) then
        ! get prior for hyper-parameters
        err=1;mess='GetLogPrior_full:Fatal:spatial latents not implemented'
    else
        do s=1,size(Fred%mdl%Sigma_Phi%par)
            ipar=>Fred%mdl%Sigma_Phi%par(s) ! define shortcut
            call GetCompPrior(ipar,logprior,feas,isnull,err,mess)
            if( (.not. feas) .or. isnull) return
            if(err>0) then;mess='GetLogPrior_full:Fatal:'//trim(mess);return;endif
        enddo
    endif
endif
!!Sigma epsilon
!ipar=>Fred%mdl%sigma_eps ! define shortcut
!call GetCompPrior(ipar,logprior,feas,isnull,err,mess)
!if( (.not. feas) .or. isnull) return
!if(err>0) then;mess='GetLogPrior_full:Fatal:'//trim(mess);return;endif

end subroutine GetLogPrior_full

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetLogHyper_full(loghyper,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: compute full hyper-lik
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:19/05/2014
!^**********************************************************************
!^* Comments: LogHyper associated with the parameter values currently stored
!^*           in structure - folding should already have been done 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^* OUT
!^*		1.LogHyper, log-prior
!^*		2.feas
!^*		3.isnull
!^*		4.err
!^*		5.mess
!^**********************************************************************
use Distribution_tools, only:GetSamplePdf,GAUSS
real(mrk), intent(out)::loghyper
logical,intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::d,k
type(PredCA_DparType), pointer::dpar
real(mrk)::pdf

err=0;mess='';loghyper=0._mrk
if(.not.Fred%AnyLAT) return

do d=1,Fred%Nd
    dpar=>Fred%mdl%dpar(d) ! define shortcut
    ! taus
    do k=1,dpar%Nk
        if (dpar%tau(k)%typ==LAT) then
            call GetSamplePdf(DistId=GAUSS,x=dpar%tau(k)%par(:)%val,par=(/0._mrk,dpar%sigma_tau(k)%val/),&
                    loga=.true.,pdf=pdf,feas=feas,isnull=isnull,err=err,mess=mess)
            if( (.not. feas) .or. isnull) return
            if(err>0) then;mess='GetLogHyper_full:Fatal:'//trim(mess);return;endif
            loghyper=loghyper+pdf
        endif
     enddo
enddo
! NOTE: code for spatial latents missing... 

end subroutine GetLogHyper_full

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetCompPrior(comp,logprior,feas,isnull,err,mess)
use Distribution_tools, only:GetPdf
type(PredCA_InfType),intent(in)::comp
real(mrk),intent(inout)::logprior
logical,intent(out)::feas,isnull
integer(mik),intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::pdf

call GetPdf(DistId=comp%PriorDist,x=comp%val,par=comp%PriorPar,&
        loga=.true.,pdf=pdf,feas=feas,isnull=isnull,err=err,mess=mess)
if( (.not. feas) .or. isnull) return
if(err>0) then;mess='GetCompPrior:'//trim(mess);return;endif
logprior=logprior+pdf

end subroutine GetCompPrior

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GetCompVal(comp,indx)
type(PredCA_InfVType), intent(in)::comp
integer(mik)::indx
real(mrk)::GetCompVal
if(comp%typ==REG .or. comp%typ==FIX) then
    GetCompVal=comp%par(1)%val
else
    GetCompVal=comp%par(indx)%val
endif
end function GetCompVal

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine LoadFred(file,priorOption,initOption,err,mess)
!^**********************************************************************
!^* Purpose: Load Fred
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:19/05/2014
!^**********************************************************************
!^* Comments:  
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file (full path)
!^*		2.priorOption, only available option is 'FLAT' for the time being
!^*		3.initOption, only available option is 'DEFAULT' for the time being
!^* OUT
!^*		1.err
!^*		2.mess
!^**********************************************************************
use utilities_dmsl_kit, only:getspareunit
use BayesianEstimation_tools, only:CreateFlatPriorList,PriorListType

character(*), intent(in)::file,prioroption,initoption
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::unt,d,k,c,m,i
type(PredCA_DparType), pointer::dpar
type(PriorListType),allocatable::priorlist(:)

err=0;mess='';m=0;
call getspareunit(unt,err,mess)
if(err/=0) then;mess='LoadFred:'//trim(mess);return;endif
open(unit=unt,file=trim(file),status='old',iostat=err)
if (err/=0) then;mess='LoadFred:Error opening file ['//trim(file)//']';return;endif
c=c+1;read(unt,*,iostat=err) Fred%cas_file
if (err/=0) then;mess=ReadConfig_errmsg(file,err,c);return;endif
close(unt)

call ReadConfig_cas(Fred%cas_file,err,mess)
if (err/=0) then;mess='LoadFred:'//trim(mess);return;endif
call ReadConfig_dat(Fred%dat_fileY,'Y',err,mess)
if (err/=0) then;mess='LoadFred:'//trim(mess);return;endif
call ReadConfig_dat(Fred%dat_filePhi,'PHI',err,mess)
if (err/=0) then;mess='LoadFred:'//trim(mess);return;endif
call ReadConfig_mdl(Fred%mdl_file,err,mess)
if (err/=0) then;mess='LoadFred:'//trim(mess);return;endif

! Complete filling of Fred
Fred%Nt=Fred%Y%Nt
Fred%Nx=Fred%Y%Ns
Fred%Ns=Fred%Phi%Ns
Fred%Nd=Fred%mdl%Nd

! Compute Ninf & anyLAT, and do extra allocations
Fred%AnyLAT=.false.;Fred%Ninf=0;Fred%Nlat=0
do d=1,Fred%Nd
    dpar=>Fred%mdl%dpar(d) ! define shortcut
    ! lambdas
    do k=1,dpar%Nk+1
        select case (dpar%lambda(k)%typ)
        case(LOC,LAT);c=Fred%Nx
        case(REG);c=1
        end select
        Fred%Ninf=Fred%Ninf+c
        if(dpar%lambda(k)%typ==LAT) Fred%Nlat=Fred%Nlat+c
        if(allocated(dpar%lambda(k)%par)) deallocate(dpar%lambda(k)%par)
        allocate(dpar%lambda(k)%par(c))
        dpar%lambda(k)%par(:)%indx=(/(i,i=m+1,m+c)/);m=m+c
        if(dpar%lambda(k)%typ==LAT) then
            Fred%AnyLAT=.true.
            ! BOTCH: Need to add hyperparameters!!!!
        endif
    enddo    
    ! taus
    do k=1,dpar%Nk
        select case (dpar%tau(k)%typ)
        case(LOC,LAT);c=Fred%Nt
        case(REG);err=1;mess='LoadFred:Fatal:option REG not allowed for tau';return
        end select
        Fred%Ninf=Fred%Ninf+c
        if(dpar%tau(k)%typ==LAT) Fred%Nlat=Fred%Nlat+c
        if(allocated(dpar%tau(k)%par)) deallocate(dpar%tau(k)%par)
        allocate(dpar%tau(k)%par(c))
        dpar%tau(k)%par(:)%indx=(/(i,i=m+1,m+c)/);m=m+c
        if(dpar%tau(k)%typ==LAT) then
            Fred%AnyLAT=.true.
            Fred%Ninf=Fred%Ninf+1
            dpar%sigma_tau(k)%indx=m+1;m=m+1
        endif
    enddo
    ! psis
    do k=1,dpar%Nk
        select case (dpar%psi(k)%typ)
        case(LOC,LAT);c=Fred%Ns
        case(REG);c=1
        end select
        Fred%Ninf=Fred%Ninf+c
        if(dpar%psi(k)%typ==LAT) Fred%Nlat=Fred%Nlat+c
        if(allocated(dpar%psi(k)%par)) deallocate(dpar%psi(k)%par)
        allocate(dpar%psi(k)%par(c))
        dpar%psi(k)%par(:)%indx=(/(i,i=m+1,m+c)/);m=m+c
        if(dpar%psi(k)%typ==LAT) then
            Fred%AnyLAT=.true.
            ! BOTCH: Need to add hyperparameters!!!!
        endif
    enddo    
enddo

! Mu Phi
if(Fred%mdl%DoCenter) then
    select case (Fred%mdl%Mu_Phi%typ)
    case(LOC,LAT);c=Fred%Ns
    case(REG);c=1
    end select
    Fred%Ninf=Fred%Ninf+c
    if(Fred%mdl%Mu_Phi%typ==LAT) Fred%Nlat=Fred%Nlat+c
    if(allocated(Fred%mdl%Mu_Phi%par)) deallocate(Fred%mdl%Mu_Phi%par)
    allocate(Fred%mdl%Mu_Phi%par(c))
    Fred%mdl%Mu_Phi%par(:)%indx=(/(i,i=m+1,m+c)/);m=m+c
    if(Fred%mdl%Mu_Phi%typ==LAT) then
        Fred%AnyLAT=.true.
        ! BOTCH: Need to add hyperparameters!!!!
    endif
else
    if(allocated(Fred%mdl%Mu_Phi%par)) deallocate(Fred%mdl%Mu_Phi%par)
    allocate(Fred%mdl%Mu_Phi%par(Fred%Ns))
endif

!Sigma Phi
if(Fred%mdl%DoScale) then
    select case (Fred%mdl%Sigma_Phi%typ)
    case(LOC,LAT);c=Fred%Ns
    case(REG);c=1
    end select
    Fred%Ninf=Fred%Ninf+c
    if(Fred%mdl%Sigma_Phi%typ==LAT) Fred%Nlat=Fred%Nlat+c
    if(allocated(Fred%mdl%Sigma_Phi%par)) deallocate(Fred%mdl%Sigma_Phi%par)
    allocate(Fred%mdl%Sigma_Phi%par(c))
    Fred%mdl%Sigma_Phi%par(:)%indx=(/(i,i=m+1,m+c)/);m=m+c
    if(Fred%mdl%Sigma_Phi%typ==LAT) then
        Fred%AnyLAT=.true.
        ! BOTCH: Need to add hyperparameters!!!!
    endif
else
    if(allocated(Fred%mdl%Sigma_Phi%par)) deallocate(Fred%mdl%Sigma_Phi%par)
    allocate(Fred%mdl%Sigma_Phi%par(Fred%Ns))
endif

!!Sigma epsilon
!Fred%Ninf=Fred%Ninf+1
!Fred%mdl%Sigma_eps%indx=m+1;m=m+1

Fred%Nnonlat=Fred%Ninf-Fred%Nlat

! Handle priors
select case(priorOption)
case('FLAT','flat','Flat')
    if(allocated(priorlist)) deallocate(priorlist);allocate(priorlist(Fred%nNonLat))
    call CreateFlatPriorList(npar=Fred%nNonLat,PriorList=priorlist,err=err,mess=mess)
    if (err/=0) then;mess='LoadFred:'//trim(mess);return;endif
    call FoldVector2Structure(priorlist=priorlist,err=err,mess=mess)
case default
    err=1;mess='LoadFred:Fatal:unknown[priorOption]';return
end select

! Handle initial value
call SetInitValue(initOption,err,mess)
if (err/=0) then;mess='LoadFred:'//trim(mess);return;endif

end subroutine LoadFred

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine SetInitValue(initOption,err,mess)
!^**********************************************************************
!^* Purpose: Set initial value of inferred parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:25/05/2014
!^**********************************************************************
!^* Comments:  
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.initOption, only available option is 'DEFAULT' for the time being
!^* OUT
!^*		1.err
!^*		2.mess
!^**********************************************************************
use Distribution_tools, only:GetRoughEstimate
use EmpiricalStats_tools, only:GetEmpiricalStats
character(*), intent(in)::initoption
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
real(mrk),parameter::JumpRatio=0.1_mrk,comp_init=0._mrk,tauscale=1._mrk,eps=1._mrk
integer(mik)::d,k,x,t,s,n
type(PredCA_DparType), pointer::dpar
real(mrk)::par(Fred%Nd,Fred%Nx),foo,muphi(Fred%Ns),sigmaphi(Fred%Ns)
real(mrk), allocatable::Xpacked(:)
logical::mask(Fred%Nt)


! Get rough estimate of D-par
do x=1,Fred%Nx
    mask=Fred%Y%val(:,x)/=Fred%Y%mv;n=count(mask)
    if(allocated(Xpacked)) deallocate(Xpacked);allocate(Xpacked(n))
    Xpacked=pack(Fred%Y%val(:,x),mask)
    call GetRoughEstimate(DistID=Fred%mdl%dist, X=Xpacked, par=par(:,x), err=err, mess=mess)
    if(err/=0) then;mess='SetInitValue:'//trim(mess);return;endif
enddo

! Estimate muphi and sigmaphi
do s=1,Fred%Ns
    mask=Fred%Phi%val(:,s)/=Fred%Phi%mv;n=count(mask)
    if(allocated(Xpacked)) deallocate(Xpacked);allocate(Xpacked(n))
    Xpacked=pack(Fred%Phi%val(:,s),mask)
    call GetEmpiricalStats(x=Xpacked,mean=muphi(s),std=sigmaphi(s),err=err,mess=mess)
    if (err/=0) then;mess='SetInitValue:'//trim(mess);return;endif
enddo

do d=1,Fred%Nd ! loop on all D-par
    dpar=>Fred%mdl%Dpar(d) ! define shortcut
    ! lambda 0
    select case(initOption)
    case('DEFAULT','Default','default')
        select case(dpar%lambda(1)%typ)
        case(LAT,LOC)
            do x=1,Fred%Nx
                call Init_engine(dpar%lambda(1)%par(x),par(d,x),eps,JumpRatio)
            enddo
        case(REG)
            call Init_engine(dpar%lambda(1)%par(1),sum(par(d,:))/real(Fred%Nx,mrk),eps,JumpRatio)
        case(FIX)
            call Init_engine(dpar%lambda(1)%par(1),dpar%defval,eps,JumpRatio)
        end select
    case default
        err=1;mess='SetInitValue:Fatal:unknown[initOption]';return
    end select
    ! loop on all components for the current Dpar
    do k=1,dpar%Nk 
        ! Lambda
        select case(initOption)
        case('DEFAULT','Default','default')
            select case(dpar%lambda(k+1)%typ)
            case(LAT,LOC)
                do x=1,Fred%Nx
                    call Init_engine(dpar%lambda(k+1)%par(x),comp_init,max(abs(dpar%lambda(1)%par(x)%val),eps),JumpRatio)
                enddo
            case(REG)
                call Init_engine(dpar%lambda(k+1)%par(1),comp_init,max(abs(dpar%lambda(1)%par(1)%val),eps),JumpRatio)
            end select
        case default
            err=1;mess='SetInitValue:Fatal:unknown[initOption]';return
        end select
        ! tau
        select case(initOption)
        case('DEFAULT','Default','default')
            select case(dpar%tau(k)%typ)
            case(LAT,LOC)
                do t=1,Fred%Nt
                    call Init_engine(dpar%tau(k)%par(t),comp_init,tauscale,JumpRatio)
                enddo
            case(REG)
                err=1;mess='SetInitValue:Fatal:option REG not allowed for tau';return
            end select
        case default
            err=1;mess='SetInitValue:Fatal:unknown[initOption]';return
        end select
        if(dpar%tau(k)%typ==LAT) then
            select case(initOption)
            case('DEFAULT','Default','default')
                call Init_engine(dpar%sigma_tau(k),tauscale,tauscale,JumpRatio)
            case default
                err=1;mess='SetInitValue:Fatal:unknown[initOption]';return
            end select
        endif      
        ! psi
        if(Fred%mdl%DoScale) then;foo=tauscale;else;foo=Fred%Phi%typscale;endif
        select case(initOption)
        case('DEFAULT','Default','default')
            select case(dpar%psi(k)%typ)
            case(LAT,LOC)
                do s=1,Fred%Ns
                    call Init_engine(dpar%psi(k)%par(s),comp_init,foo,JumpRatio)
                enddo
            case(REG)
                call Init_engine(dpar%psi(k)%par(1),comp_init,foo,JumpRatio)
            end select
        case default
            err=1;mess='SetInitValue:Fatal:unknown[initOption]';return
        end select        
    enddo
enddo

! Phi mean
if(Fred%mdl%DoCenter) then
    select case(initOption)
    case('DEFAULT','Default','default')
        select case(Fred%mdl%Mu_Phi%typ)
        case(LAT,LOC)
            do s=1,Fred%Ns
                call Init_engine(Fred%mdl%Mu_Phi%par(s),muphi(s),Fred%Phi%typscale,JumpRatio)
            enddo
        case(REG)
            call Init_engine(Fred%mdl%Mu_Phi%par(1),sum(muphi)/real(Fred%Ns,mrk),Fred%Phi%typscale,JumpRatio)
        end select
    case default
        err=1;mess='SetInitValue:Fatal:unknown[initOption]';return
    end select        
else
    Fred%mdl%Mu_Phi%par(:)%val=0._mrk
endif

! Phi Sdev
if(Fred%mdl%DoScale) then
    select case(initOption)
    case('DEFAULT','Default','default')
        select case(Fred%mdl%Sigma_Phi%typ)
        case(LAT,LOC)
            do s=1,Fred%Ns
                call Init_engine(Fred%mdl%Sigma_Phi%par(s),sigmaphi(s),Fred%Phi%typscale,JumpRatio)
            enddo
        case(REG)
            call Init_engine(Fred%mdl%Sigma_Phi%par(1),sum(sigmaphi)/real(Fred%Ns,mrk),Fred%Phi%typscale,JumpRatio)
        end select
    case default
        err=1;mess='SetInitValue:Fatal:unknown[initOption]';return
    end select        
else
    Fred%mdl%Sigma_Phi%par(:)%val=1._mrk
endif

!! Phi residuals
!select case(initOption)
!case('DEFAULT','Default','default')
!    if(Fred%mdl%DoScale) then;foo=tauscale;else;foo=Fred%Phi%typscale;endif
!    call Init_engine(Fred%mdl%Sigma_eps,foo,foo,JumpRatio)
!case default
!    err=1;mess='SetInitValue:Fatal:unknown[initOption]';return
!end select        

end subroutine SetInitValue

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Init_engine(comp,val,typscale,JumpRatio)
type(PredCA_InfType), intent(inout)::comp
real(mrk), intent(in)::val,JumpRatio,typscale

comp%val=val
if(val/=0._mrk) then
    comp%JumpStd=JumpRatio*abs(val)
else
    comp%JumpStd=JumpRatio*abs(typscale)
endif

end subroutine Init_engine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadConfig_cas(file,err,mess)
!^**********************************************************************
!^* Purpose: read .cas file and start filling Fred
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:19/05/2014
!^**********************************************************************
!^* Comments:  
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file (full path)
!^* OUT
!^*		1.err
!^*		2.mess
!^**********************************************************************
use utilities_dmsl_kit, only:getspareunit
character(*), intent(in)::file
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::unt,k

err=0;mess='';k=0
call getspareunit(unt,err,mess)
if(err/=0) then;mess='ReadConfig_cas:'//trim(mess);return;endif
open(unit=unt,file=trim(file),status='old',iostat=err)
if (err/=0) then;mess='ReadConfig_cas:Error opening file ['//trim(file)//']';return;endif
k=k+1;read(unt,*,iostat=err) Fred%nickname
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Fred%dat_fileY
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Fred%dat_filePhi
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Fred%mdl_file
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Fred%OUTfolder
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
close(unt)

end subroutine ReadConfig_cas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadConfig_dat(file,YorPHI,err,mess)
!^**********************************************************************
!^* Purpose: read .dat file and start filling Fred
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:19/05/2014
!^**********************************************************************
!^* Comments:  
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file (full path)
!^*		2.YorPHI, 'Y' or 'PHI'
!^* OUT
!^*		1.err
!^*		2.mess
!^**********************************************************************
use utilities_dmsl_kit, only:getspareunit
character(*), intent(in)::file,YorPHI
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::unt,k,t,s
type(PredCA_DataType),pointer::Z

err=0;mess='';k=0
if(trim(YorPHI)=='Y'.or.trim(YorPHI)=='y') then
    Z=>Fred%Y
elseif(trim(YorPHI)=='PHI'.or.trim(YorPHI)=='phi'.or.trim(YorPHI)=='Phi') then
    Z=>Fred%Phi
else
    err=1;mess="ReadConfig_dat:Fatal:unknown [YorPHI] ['Y','PHI']";return;
endif
call getspareunit(unt,err,mess)
if(err/=0) then;mess='ReadConfig_dat:'//trim(mess);return;endif
open(unit=unt,file=trim(file),status='old',iostat=err)
if (err/=0) then;mess='ReadConfig_dat:Error opening file ['//trim(file)//']';return;endif
k=k+1;read(unt,*,iostat=err) Z%nickname
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Z%Nt,Z%Ns
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Z%file
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Z%grdfile_t
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Z%grdfile_s
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Z%grdfile_f
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Z%mv
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Z%typscale
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
close(unt)
! Read data
if(allocated(Z%val)) deallocate(Z%val);allocate(Z%val(Z%Nt,Z%Ns))
call getspareunit(unt,err,mess)
if(err/=0) then;mess='ReadConfig_dat:'//trim(mess);return;endif
open(unit=unt,file=trim(Z%file),status='old',iostat=err)
if (err/=0) then;mess='ReadConfig_dat:Error opening file ['//trim(Z%file)//']';return;endif
do t=1,Z%Nt
    read(unt,*,iostat=err) Z%val(t,:)
    if (err/=0) then;mess=ReadConfig_errmsg(Z%file,err,t);return;endif
enddo
close(unt)
! Read t-grid
if(allocated(Z%tgrid)) deallocate(Z%tgrid);allocate(Z%tgrid(Z%Nt))
call getspareunit(unt,err,mess)
if(err/=0) then;mess='ReadConfig_dat:'//trim(mess);return;endif
open(unit=unt,file=trim(Z%grdfile_t),status='old',iostat=err)
if (err/=0) then;mess='ReadConfig_dat:Error opening file ['//trim(Z%grdfile_t)//']';return;endif
do t=1,Z%Nt
    read(unt,*,iostat=err) Z%tgrid(t)
    if (err/=0) then;mess=ReadConfig_errmsg(Z%grdfile_t,err,t);return;endif
enddo
close(unt)
! Read s-grid
if(allocated(Z%sgrid)) deallocate(Z%sgrid);allocate(Z%sgrid(Z%Ns,2))
call getspareunit(unt,err,mess)
if(err/=0) then;mess='ReadConfig_dat:'//trim(mess);return;endif
open(unit=unt,file=trim(Z%grdfile_s),status='old',iostat=err)
if (err/=0) then;mess='ReadConfig_dat:Error opening file ['//trim(Z%grdfile_s)//']';return;endif
do s=1,Z%Ns
    read(unt,*,iostat=err) Z%sgrid(s,:)
    if (err/=0) then;mess=ReadConfig_errmsg(Z%grdfile_s,err,t);return;endif
enddo
close(unt)
! Read f-grid
if(allocated(Z%fgrid)) deallocate(Z%fgrid);allocate(Z%fgrid(Z%Ns))
call getspareunit(unt,err,mess)
if(err/=0) then;mess='ReadConfig_dat:'//trim(mess);return;endif
open(unit=unt,file=trim(Z%grdfile_f),status='old',iostat=err)
if (err/=0) then;mess='ReadConfig_dat:Error opening file ['//trim(Z%grdfile_f)//']';return;endif
do s=1,Z%Ns
    read(unt,*,iostat=err) Z%fgrid(s)
    if (err/=0) then;mess=ReadConfig_errmsg(Z%grdfile_f,err,t);return;endif
enddo
close(unt)

end subroutine ReadConfig_dat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadConfig_mdl(file,err,mess)
!^**********************************************************************
!^* Purpose: read .cas file and start filling Fred
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:19/05/2014
!^**********************************************************************
!^* Comments:  
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file (full path)
!^* OUT
!^*		1.err
!^*		2.mess
!^**********************************************************************
use utilities_dmsl_kit, only:getspareunit
use distribution_tools, only:GetParNumber,GetParName
character(*), intent(in)::file
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::unt,k,d
character(250), allocatable::parname(:)

err=0;mess='';k=0
call getspareunit(unt,err,mess)
if(err/=0) then;mess='ReadConfig_mdl:'//trim(mess);return;endif
open(unit=unt,file=trim(file),status='old',iostat=err)
if (err/=0) then;mess='ReadConfig_mdl:Error opening file ['//trim(file)//']';return;endif
k=k+1;read(unt,*,iostat=err) Fred%mdl%nickname
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Fred%mdl%dist
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
call GetParNumber(DistID=Fred%mdl%dist, npar=Fred%mdl%Nd, err=err, mess=mess)
if (err/=0) then;mess='ReadConfig_mdl:'//trim(mess);return;endif
k=k+1;read(unt,*,iostat=err) Fred%mdl%DoCenter,Fred%mdl%DoScale
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
k=k+1;read(unt,*,iostat=err) Fred%mdl%Mu_Phi%typ,Fred%mdl%Sigma_Phi%typ
if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
if(allocated(Fred%mdl%DparFile)) deallocate(Fred%mdl%DparFile);allocate(Fred%mdl%DparFile(Fred%mdl%Nd))
do d=1,Fred%mdl%Nd
    k=k+1;read(unt,*,iostat=err) Fred%mdl%DparFile(d)
    if (err/=0) then;mess=ReadConfig_errmsg(file,err,k);return;endif
enddo
close(unt)
if(allocated(Fred%mdl%Dpar)) deallocate(Fred%mdl%Dpar);allocate(Fred%mdl%Dpar(Fred%mdl%Nd))
do d=1,Fred%mdl%Nd
    call ReadConfig_dpr(Fred%mdl%DparFile(d),Fred%mdl%Dpar(d),err,mess)
    if (err/=0) then;mess='ReadConfig_mdl:'//trim(mess);return;endif
enddo
if(allocated(parname)) deallocate(parname);allocate(parname(Fred%mdl%Nd))
call GetParName(DistID=Fred%mdl%dist, name=parname, err=err, mess=mess)
if (err/=0) then;mess='ReadConfig_mdl:'//trim(mess);return;endif
do d=1,Fred%mdl%Nd
   Fred%mdl%Dpar(d)%nickname=parname(d)
enddo

end subroutine ReadConfig_mdl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadConfig_dpr(file,comp,err,mess)
!^**********************************************************************
!^* Purpose: read .cas file and start filling Fred
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:19/05/2014
!^**********************************************************************
!^* Comments:  
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.file (full path)
!^* INOUT
!^*		1.comp
!^* OUT
!^*		1.err
!^*		2.mess
!^**********************************************************************
use utilities_dmsl_kit, only:getspareunit
character(*), intent(in)::file
type(PredCA_DparType), intent(inout)::comp
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::unt,c,k

err=0;mess='';c=0
call getspareunit(unt,err,mess)
if(err/=0) then;mess='ReadConfig_dpr:'//trim(mess);return;endif
open(unit=unt,file=trim(file),status='old',iostat=err)
if (err/=0) then;mess='ReadConfig_dpr:Error opening file ['//trim(file)//']';return;endif
c=c+1;read(unt,*,iostat=err) comp%Nk
if (err/=0) then;mess=ReadConfig_errmsg(file,err,c);return;endif
c=c+1;read(unt,*,iostat=err) comp%defval
if (err/=0) then;mess=ReadConfig_errmsg(file,err,c);return;endif
c=c+1;read(unt,*,iostat=err) comp%link
if (err/=0) then;mess=ReadConfig_errmsg(file,err,c);return;endif
if(allocated(comp%lambda)) deallocate(comp%lambda);allocate(comp%lambda(comp%Nk+1))
c=c+1;read(unt,*,iostat=err) comp%lambda(1)%typ
if (err/=0) then;mess=ReadConfig_errmsg(file,err,c);return;endif
if(comp%Nk>0) then
    c=c+1;read(unt,*,iostat=err) comp%lambda(2:(comp%Nk+1))%typ
    if (err/=0) then;mess=ReadConfig_errmsg(file,err,c);return;endif
endif
if(allocated(comp%psi)) deallocate(comp%psi);allocate(comp%psi(comp%Nk))
if(comp%Nk>0) then
    c=c+1;read(unt,*,iostat=err) comp%psi(:)%typ
    if (err/=0) then;mess=ReadConfig_errmsg(file,err,c);return;endif
endif
if(allocated(comp%tau)) deallocate(comp%tau);allocate(comp%tau(comp%Nk))
if(allocated(comp%sigma_tau)) deallocate(comp%sigma_tau);allocate(comp%sigma_tau(comp%Nk))
if(comp%Nk>0) then
    c=c+1;read(unt,*,iostat=err) comp%tau(:)%typ
    if (err/=0) then;mess=ReadConfig_errmsg(file,err,c);return;endif
endif

end subroutine ReadConfig_dpr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ReadConfig_errmsg(file,errcode,line)
use utilities_dmsl_kit, only:number_string
character(*),intent(in)::file
integer(mik),intent(in)::errcode,line
character(250)::ReadConfig_errmsg

ReadConfig_errmsg='Error reading file['//trim(file)//'],line['//trim(number_string(line))//&
                   &'],error code:['//trim(number_string(errcode))//']'
end function ReadConfig_errmsg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module PredCA_tools

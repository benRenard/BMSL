module Optimisation_tools
! Purpose:
! Programmer: Dmitri Kavetski
! Last modified:
! Comments:
use kinds_dmsl_kit
implicit none
! type definitions
! public::qndmsl_smeth,lbfgs_smeth,multiStart_smeth
integer(mik),parameter::qndmsl_smeth=1,lbfgs_smeth=2,multiStart_smeth=3 ! algorithm options
! variable definitions
!----------------------------------------------------
contains
!----------------------------------------------------
SUBROUTINE optimiserWrap (smeth,nopt,evalFunc,dataIN,nDim,xLo,xHi,msLo,msHi,xscale,activeSet,msFile, &
                          imeth_qn,gmeth,hmeth_qn,himeth_qn,&  ! QN options
                          mMem_lbfgs,&                      ! LBFGS method
                          mometh,&
                          x0, xMin, fMin, chkGrad,chkHess,&
                          inverseHessian, hess, iterMax, hisFile, &
                          iternfo, resFile,&
                          fcalls,gcalls,hcalls,&
                          err, message)
! Estimates optimum of objective function and computes estimate of covariance matrix.
USE linalg_dmsl_kit,ONLY:getEigValsVecsRealSym,invertEigs,eig_makemat,covar_to_correl
USE utilities_dmsl_kit,ONLY:norm2,flip_UtoL,write_matrix,arthsi,ns=>number_string,trimv=>TRIM, &
  getSpareUnit,getHessFromFunc,getHessFromGrad,&
  fatalError=>fatalError_conDef,one,zero,putDiag,fmatmul_mm
USE optimiser_dmsl_kit,ONLY:qnewton,LBFGS,multiStart,qnewtonUnwise_type
use numerix_dmsl_kit,only:gradvec_ridders,uniran
use types_dmsl_kit, only:data_ricz_type
implicit none
! Dummy args
type(data_ricz_type),intent(in), optional::dataIN
integer(mik),intent(in)::smeth,mometh
real(mrk),intent(in)::xlo(:),xhi(:),mslo(:),mshi(:)
integer(mik),intent(out)::activeSet(:)
integer(mik), intent(in)::nDim,nopt,iternfo
! real(mrk), intent(in)::x0(:)
real(mrk), intent(inout)::x0(:)
real(mrk), intent(in)::xscale(:)
real(mrk), intent(inout)::xMin(:)
real(mrk), intent(out)::fMin
character(*),intent(in)::msFile
logical(mlk), intent(in)::inverseHessian, chkGrad,chkHess
real(mrk), intent(out)::hess(:,:)
integer(mik), intent(in)::iterMax,imeth_qn,gmeth,hmeth_qn,himeth_qn,mMem_lbfgs
character(*), intent(in) :: hisFile, resFile
integer(mik), intent(out)::fcalls,gcalls,hcalls
integer(mik), intent(out)::err
character(*), intent(out)::message
! interface to objective function
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! locals

! Tolerances
integer(mik),parameter::maxFev=100000000
REAL(mrk), PARAMETER :: gtol=1.e-15, ftol=1.0E-30_mrk, stol=1.e-15_mrk
REAL(mrk), PARAMETER ::ftol_sce=1.e-10_mrk
!
! computation cost
INTEGER(mik)::iter,fcalls0
integer(mik),parameter::itfile=666,itfile0=777,umult=666,untCheck=134 !,nopt=10
integer(mik)::useed
!
! QN variables
REAL(mrk):: grad(nDim),eigvals(nDim), fscale, stepmax,epsFmax,dumM(ndim,ndim)
real(mrk)::hessOpt(size(x0),size(x0))
REAL(mrk) :: trustRad,condMax,condEff
real(mrk)::memHess2, memadd
type(qnewtonUnwise_type)::qnewtonUnwise
! miscellaneous
character(100)::msg
INTEGER(mik) :: aSize(4), err2,fdigits,uOut, unt,iBad,nmodE
REAL(mrk) :: t0, t1, cn1,cn2
logical(mlk)::okE(size(x0)),checkEigs,checkFuncGrad

! choices within multistart
integer(mik),parameter::none_mometh=    1,&   ! random search
                        QN_DMSL_mometh= 2,&   ! native DMSL quasi-Newton (best for middle-D problems)
                        QN_IMSL_mometh= 3,&   ! IMSL-based quasi-Newton
                        LBFGS_mometh=   4,&   ! LBFGS scheme (best for huge-D problems)
                        SCE_mometh=     5     ! SCE search (multistart version probably redundant)
! integer(mik),parameter::mometh=QN_DMSL_mometh ! LBFGS_mometh ! 
! Start procedure here
! Initialise
iter=0; fcalls=0; gcalls=0; hcalls=0
aSize(1) = SIZE(x0)
aSize(2) = SIZE(xMin)
aSize(3) = SIZE(hess,1)
aSize(4) = SIZE(hess,2)
IF (ANY(aSize /= nDim)) THEN
  err = 1000; message = 'f-quasiNewTonSolver/Inconsistent dimensions'; RETURN
END IF
!
! History
CALL getSpareUnit(unt=uout,err=err,message=message)
IF (err /= 0) then
    message='f-optimerWrap/Unable to get unit for history file '//TRIM(message);return
endif
OPEN(uout, file=TRIM(hisFile), status="unknown", action="WRITE", iostat=err)
IF(err /= 0) then
    !CALL fatalError ('f-optimerWrap/Failed to open history file: '//TRIM(hisFile))
    message='Fatal error! Failed to open history file: '//TRIM(hisFile);return
!    return
endif
!
! Output
CALL getSpareUnit(unt=unt,err=err,message=message)
IF (err /= 0) then
    message='f-optimerWrap/Unable to get unit for results file '//TRIM(message);return
endif
OPEN(unt, file=TRIM(resFile), status="unknown", action="WRITE", iostat=err)
IF(err /= 0) then
    !CALL fatalError ('f-optimerWrap/Failed to open results file: '//TRIM(resFile))
    message='Fatal error! Failed to open results file: '//TRIM(resFile);return
!    return
endif

WRITE(unt,*)
CALL write_matrix(unt=unt,header=" Starting vector:",v=X0,nFig=12,&
  transp=.false.,vLabel="  var"//trimv(ns(arthsi(ndim)),3)//": ",err=err,message=msg)
WRITE(unt,*)

fscale=1._mrk; fdigits=-1; epsFmax=1.e-15_mrk
trustrad = 1.0e-2_mrk  ! <0 --> default initial trust radius
stepmax = 1.e5_mrk*MAX(norm2(x0/xscale),norm2(1.0_mrk/xscale)) ! Recommended
!       stepmax=1._mrk ! more recommended

qnewtonUnwise%iternfo=iternfo ! iteration progress
IF (chkGrad) THEN
! Gradient checking
  qnewtonUnwise%chkGrd= 4  ! full gradient check
  qnewtonUnwise%chkGrd_gmeth=2 ! 1=forward; 2=central differences for checking
ELSE
!   qnewtonUnwise%chkGrd = 0  ! check when fail (usual default)
  qnewtonUnwise%chkGrd = -1 ! never check
ENDIF

if(chkHess)then
! Hessian checking options
  qnewtonUnwise%chkHess=2 ! full Hessian check
  qnewtonUnwise%chkHess_hmeth=2 ! 1=forward; 2=central differences of gradient for checking
else
  qnewtonUnwise%chkHess=0 ! no Hessian check
endif


checkFuncGrad=.false.
! checkFuncGrad=.true.
! if(ndim>40)checkFuncGrad=.true.
! checkFuncGrad=.true.
if(checkFuncGrad)then
! test function
!   open(untCheck,file="bateacheck.dat",status="old",action="read",iostat=err)
!   read(untCheck,*)condEff
!   do iter=1,ndim
!     read(untCheck,*)xmin(iter),eigvals(iter)
!   enddo
!   call evalFunc(x=xmin,feas=checkEigs,fx=cn1,gradFx=grad,err=err,message=message)
!   call evalFunc(x=xmin,feas=checkEigs,fx=cn2,err=err,message=message)
!   call evalFunc(x=xmin,feas=checkEigs,gradFx=hess(:,1),err=err,message=message)
!   hess(:,2)=hess(:,1)-grad; hess(:,3)=grad-eigvals
!   x0=xmin
!   hess(:,4)=1.e-4_mrk
!   call gradvec_ridders(evalFunc,x=xmin,dx=hess(:,4),tol=1.e-8_mrk,fscale=fscale,&
!     epsF=epsRe,xscale=xscale,ncoarse=50,redCoarse=1.e-2_mrk,&
!     ntab=10,con=1.4_mrk,safe=1.e1_mrk,gradvec=eigvals,&
!     E_gradvec=hess(:,2),hx=hess(:,3),err=err,message=message)
!   goto 100
!DK: test correct reinitialisation (should produce identical results on each call)
  xmin=x0
  write(*,*)"checking function time"
  call cpu_time(t0)
  do iter=1,1000
!     if(modulo(iter,2)==0)then
    call uniran(harvest=xmin,low=msLo,high=msHi)

!       xmin=xmin+1._mrk
!     else
!       xmin=xmin-1._mrk
!     endif
!     do ibad=1,3
!       call evalFunc(x=xmin,feas=checkEigs,fx=cn1,err=err,message=message)
      call evalFunc(x=xmin,feas=checkEigs,fx=cn2,gradFx=grad,err=err,message=message)
!     enddo
    if(err/=0)then
      write(*,'(a,i0,a)')"error @",iter,"; "//trim(message)
    endif
  enddo
  call cpu_time(t1)
  write(*,*)"total time is",t1-t0
  stop
! 100 continue
endif
!
! Solver
CALL CPU_TIME(t0)
selectcase(smeth)
case(qndmsl_smeth)
  CALL qnewton(                                                  &
     evalfunc=evalFunc, dataIN=dataIN,                                         &
     x0=x0,xscale=xscale,fscale=fscale,fdigits=fdigits,          &
     xLo=xLo,xHi=xHi,activeSet=activeSet,                        &
     gtol=gtol,stol=stol,ftol=ftol,                              &
     imeth=imeth_qn,gmeth=gmeth,hmeth=hmeth_qn,                  & 
     himeth=himeth_qn,trustRad=trustRad,                         &
     stpmax=stepmax,                                             &
     maxIter=itermax,maxFev=maxFev,                              &
     xopt=xmin,fopt=fmin,                                        &
     iter=iter,fcalls=fcalls,gcalls=gcalls,hcalls=hcalls,        &
     gradopt=grad,hessopt=hessOpt,                               &
     qnewtonUnwise=qnewtonUnwise,                                &
     uout=uout,err=err,message=message, memhess2=memhess2)
case(lbfgs_smeth)
  call LBFGS(evalfunc=evalFunc,x0=x0,xLo=xLo,xHi=xHi,bType=spread(2,1,ndim), dataIN=dataIN,&
    gtol=gtol,ftol=ftol,xscale=xscale,fscale=fscale,fdigits=fdigits,        &
    stpmax=stepmax,gmeth=gmeth,mMem=mMem_lbfgs,FDCDswitch=.true.,           &
    activeSet=activeSet,&
    maxIter=itermax,maxFev=maxFev,itprint=2,itfile=itfile,itfile0=itfile0,  &
    xmin=xmin,fmin=fmin,gmin=grad,iter=iter,fcalls=fcalls,gcalls=gcalls,err=err,message=message)
case(multiStart_smeth)
  useed=-230981
  open(umult,file=msFile,action="write",status="unknown")
  call multiStart(evalfunc=evalFunc, dataIN=dataIN,                                            &
    xlo=xlo,xhi=xhi,                                                            &
    xscale=xscale,fscale=fscale,fdigits=fdigits,                                &
    x0=spread(x0,2,1),x0lo=msLo,x0hi=msHi,                                      &
    nopt=nopt,uout=umult,uouti=-abs(itfile),useed=useed,                        &
    gtol=gtol,stol=stol,ftol=merge(ftol,ftol_sce,mometh/=SCE_mometh),           &
    mometh=mometh,                                                              &
    gmeth=gmeth,itermax=itermax,maxFev=maxFev,stpmax=stepmax,trustRad0=trustRad,&
    imeth_qn=imeth_qn,hmeth_qn=hmeth_qn,himeth_qn=himeth_qn,& !qnewtonUnwise_qn=qnewtonUnwise_qn,&
    mMem_lbfgs=mMem_lbfgs,  &
    doepsf=.false., fcallstotal=fcalls,                                                   &
!          sceSettings_sce=sceSettings_sce,&
    xopt=xmin,fopt=fmin,gradopt=grad,epsFmax=epsFmax,                           &
    err=err,message=message, memadd=memadd)
  close(umult)
endselect
CALL cpu_time(t1)
if(err<0)err=0
!
! Report results
IF(err /= 0)THEN
  WRITE(unt,*) "ERROR: "//TRIM(message)
ELSE
  WRITE(unt,*) "SUCCESS: "//TRIM(message)
END IF
WRITE(unt,*)
WRITE(unt,*) 'Func. value at solution',fmin
WRITE(unt,*)
CALL write_matrix(unt=unt,header=" Solution vector:",v=xmin,nFig=12,&
  transp=.false.,vLabel="  var"//trimv(ns(arthsi(ndim)),3)//": ",err=err2,message=msg)
WRITE(unt,*)
WRITE(unt,*) "CPU time:",t1-t0,"secs"
WRITE(unt,*) "Func=",fmin
WRITE(unt,*)
WRITE(unt,*) 'Iterations:',iter
WRITE(unt,*) 'Reported Func. evals:',fcalls
WRITE(unt,*) 'Reported Grad. evals:',gcalls
WRITE(unt,*) 'Reported Hess. evals:',hcalls
WRITE(unt,*)
CALL write_matrix(unt=unt,header=" Gradient at reported optimum:",v=grad,nFig=12,&
  transp=.false.,vLabel="  var"//trimv(ns(arthsi(ndim)),3)//": ",err=err2,message=msg)
!
! * estimate Hessian at the mode
selectcase(gmeth)
case(0)       ! get Hessian from the gradient
  call getHessFromGrad(evalFunc=evalFunc,x=xmin,gradFx=grad,xscale=xscale,epsF=epsFmax,&
    useHxDef=.true.,&
    hmeth=2,hessfd=hess,gcalls=fcalls0,err=err2,message=msg)
case default  ! get Hessian from function
  call getHessFromFunc(evalFunc=evalFunc,x=xmin,fx=fmin,xscale=xscale,epsF=epsFmax,&
    useHxDef=.true.,&
    hmeth=1,hessfd=hess,fcalls=fcalls0,err=err2,message=msg)
!       call getHessDiagFromFunc(evalFunc,dataIN,dataOUT,x,fx,xscale,epsF,hx,useHxDef,&
!           hmeth,hessDiag,fcalls0,err2,msg)
endselect
if(err2/=0) then
    write(*,*) 'Hessian computation:Fatal error to be handled!'! The calculations for the present gauge skipped!'
    close (uOut); close(unt)
    return
!    call fatalError(trim(msg))
endif

IF (inverseHessian)THEN
  hessOpt=hess; dumM=hess ! get power specter of the hessian matrix
  call getEigValsVecsRealSym(a2evex=hessOpt,doEigVex=.true.,emeth=2,&
    sortEigs=+1,eigvals=eigvals,err=err2,message=msg)
  if(err2/=0)then
    err=err2;message='f-optimerWrap/BUG?/'//TRIM(msg)
    close (uOut); close(unt)
    return
  endif
  CALL write_matrix(unt=unt,header=" Eigenvalues:",v=eigvals,nFig=12,&
    transp=.false.,vLabel="  eig"//trimv(ns(arthsi(ndim)),3)//": ",err=err2,message=msg)
  if(err2/=0)then
    err=err2;message='f-optimerWrap/BUG?/'//TRIM(msg)
    close (uOut); close(unt)
    return
  endif
  condMax=1.e-5_mrk/epsRe; condEff=one/sqrt(epsRe); checkEigs=.true.
  if(checkEigs)then
    call invertEigs(eigvals=eigvals,DefStat=+1,condMax=condMax,&
      nmodE=nmodE,okE=okE,condOrig=cn1,cond=cn2,err=err2,message=msg)
  else
    eigvals=one/eigvals; nmodE=0; cn1=maxval(abs(eigvals))/minval(abs(eigvals)); cn2=cn1
  endif
  WRITE(unt,*) 'Condition number of original Hessian:', cn1
  WRITE(unt,*)
  if(nmodE==0)then
    write(unt,'(a)')"Original covariance safely positive definite"
  else
    write(unt,'(a)')"Warning - original covariance not positive definite [eigen-fixed]"
    WRITE(unt,*) 'Condition number of modified covariance:', cn2
    WRITE(unt,*)
  endif   ! construct safely positive definite modified covariance
  call eig_makemat(eigvals=eigvals,eigvecs=hessOpt,a=hess,err=err2,message=msg)
  hessOpt=matmul(hess,dumM) ! modification error
  call covar_to_correl(covar=hess,correl=hessOpt,keepSdev=.true.,&
                       iBad=iBad,err=err2,message=msg)
  if(iBad/=0)then
    err=iBad;message="BUG - covariance not positive definite slipped"//trim(msg);
    close (uOut); close(unt)
    return
  elseif(err2/=0)then
    err=err2;message='f-optimerWrap/BUG/'//TRIM(msg);
    close (uOut); close(unt)
    return
  endif
  CALL write_matrix(unt=unt,header=' Approximate std dev/correlation matrix:',&
    m=hessOpt,nFig=12,display=-1,err=err2,message=msg,&
    vLabel="  var"//trimv(ns(arthsi(ndim)),3)//": ")
  WRITE(unt,*)
ENDIF
CLOSE(uout)
CLOSE (unt)
! End procedure here
END SUBROUTINE optimiserWrap
!----------------------------------------------------

end module Optimisation_tools

!*******************************************************************
subroutine myApplication(runMode)
!-----------------------------------------------------------------
! Sample myApplication procedure that governs the execution of any
! DMSL-FORTWIN applications.
! Arguments:
! * runMode determines whether startUp,normalRun,or finish is proceeding.
! * message is used for any strongly abnormal finish that could not
!   be handled otherwise. Be prepared to crash/hangup in many such cases.
!   Report to author such cases for investigation.
! Comments
! It is STRONGLY recommended that the following template be
! used to ensure correct interfacing with DMSL
!-----------------------------------------------------------------
use iotools_dmsl_mod, only:startupMode,actualRun,& ! essentials
  dmsl_AppName,dmsl_exceptHand,&  ! other runtime OS features
  dmsl_windowX_nchars,dmsl_windowY_nchars,&  ! other runtime OS features
  fatalError   ! any other procedures
use utilities_dmsl_kit,only:ns=>number_string
use kinds_dmsl_kit

implicit none
! dummies
integer,intent(in)::runMode
! locals 

interface
endinterface
! Start procedure here
selectcase(runMode)
case(startupMode) ! opportunity to change any OS variables
  ! This is the first statement in DMSL that is executed.
  dmsl_AppName="DMSL Application" ! can change the window title
  dmsl_exceptHand=.true.          ! and floating point handling
  dmsl_windowX_nchars=100    ! Specify window width
  dmsl_windowY_nchars=40    ! and window height
case(actualRun)   ! run main Fortran code
!  call GEV_3_Likelihood(N,x,F,pp,PrOfNonExceedanceTreshold,polozenie,skala,ksztalt,kwantyl,lnL,AIC,MaxL,Dmax,xTreshold,NumberOfElementsAboveTreshold,pval,Delta,C_or_V)

case default
  call fatalError("BAD ARGUMENT TO MYAPPLICATION: "//&
                    trim(ns(runMode))//"&YOUR'RE ABOUT TO CRASH!!!")
endselect
! End procedure here
endsubroutine myApplication
!*******************************************************************
subroutine myFinish(finishMode,messageIN,statOUT,messageOUT)
! Purpose: shutdown routine.
! * message is used for any strongly abnormal finish that could not
!   be handled otherwise. Be prepared to crash/hangup in many such cases.
!   Report to author such cases for investigation.
use ioTools_dmsl_mod,only:alertDialog
implicit none
! dummies
integer(4),intent(in),optional::finishMode
character(*),intent(in),optional::messageIN
integer(4),intent(inout),optional::statOUT
character(*),intent(out),optional::messageOUT
! Start procedure here
if(present(messageIN))call alertdialog(trim(messageIN))
if(present(statOUT))statOUT=0
if(present(messageOUT))messageOUT="nada"
! * do anything else here - close data files, etc.
! * but do not stop (to allow freeing any resources by DMSL)
! End procedure here
endsubroutine myFinish
!*******************************************************************

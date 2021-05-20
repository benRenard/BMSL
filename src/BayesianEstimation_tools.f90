module BayesianEstimation_tools

!~**********************************************************************
!~* Purpose: Simple Bayesian-MCMC fitting of a dist. to some 1-d data
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified: 08/07/2010
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1.GetLogPrior
!~*    2.
!~*    3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use Distribution_tools,only:covariate, IsCovNeeded

implicit none
Private
public :: GetLogPrior, Bayes_EstimPar, Bayes_Predictive, LoadMCMCSample,&
          GetMCMCSummary,GeneratePred,GetMode,GetPredictiveEstimates,GetModalEstimates,&
          CreateFlatPriorList

! variables globally available to this module
real(mrk),allocatable:: Xdata(:)
character(100)::Xdist
type, public:: PriorListType
    character(250)::dist
    real(mrk), allocatable::par(:)
end type PriorListType
type(PriorListType), allocatable:: XpriorList(:)
type(covariate),allocatable::Xcov(:)

Contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GetLogPrior(teta,PriorList,lp, feas, isnull,err,mess)

!^**********************************************************************
!^* Purpose: Compute prior distribution of vector teta - assumes prior independence
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:08/07/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.teta, par vector
!^*    2.PriorList, list of priors & their par (see PriorListType above)
!^*    3.
!^* OUT
!^*    1.lp, log-prior
!^*    2.feas, are priors feasible?
!^*    3.isnull, is log-prior=0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************

use Distribution_tools, only: GetParNumber, GetParFeas, GetPdf

real(mrk), intent(in)::teta(:)
Type(PriorListType), intent(in)::PriorList(:)
real(mrk), intent(out)::lp
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
integer(mik)::n,i,npar
real(mrk), allocatable::pdf(:)

err=0;mess='';feas=.true.;isnull=.false.;lp=undefRN

! check size
if(size(teta)/=size(PriorList)) then
    feas=.false.;err=1;mess="GetLogPrior:FATAL:size mismatch [teta]/[PriorList]";
    return;
endif
n=size(teta)

! check PriorList
do i=1,n
    ! check prior par are given
    if(.not.(allocated(PriorList(i)%par))) then
        feas=.false.;err=1
        mess="GetLogPrior:FATAL:[PriorList(i)%par] is not allocated";
        return
    endif
    ! check size
    CALL GetParNumber(PriorList(i)%dist, npar, err, mess)
    If(err>0) then
        mess='GetLogPrior:'//trim(mess)
        feas=.false.;return
    endif
    if(size(PriorList(i)%par)/=npar) then
        feas=.false.;err=1;mess="GetLogPrior:FATAL:size mismatch [PriorList(i)%par]";
        return;
    endif
    ! check that prior par are feasible
    CALL GetParFeas(PriorList(i)%dist, PriorList(i)%par, feas, err, mess)
    If(err>0) then
        mess='GetLogPrior:'//trim(mess)
        feas=.false.;return
    endif
    if(.not. feas) then
        !err=1;mess="GetLogPrior:FATAL:unfeasible [PriorList(i)%par]"
        return
    endif
enddo

! Everything's OK, compute!
if(allocated(pdf)) deallocate(pdf)
allocate(pdf(n))
do i=1,n
    call GetPdf(distID=PriorList(i)%dist,x=teta(i),par=PriorList(i)%par,&
                loga=.true.,pdf=pdf(i),feas=feas,isnull=isnull,err=err,mess=mess)
    If(err>0) then
        mess='GetLogPrior:'//trim(mess)
        feas=.false.;return
    endif
    if(.not. feas) then 
        !err=1;mess="GetLogPrior:FATAL:unfeasible [PriorList(i)%par]"
        return
    endif
    if(isnull) return ! impossible teta value, prior=0
enddo

! sum log-prior
lp=sum(pdf)

end subroutine GetLogPrior

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GetLogPost(teta,X, distID, PriorList,cov,lp, feas, isnull,err,mess)

!^**********************************************************************
!^* Purpose: Compute unnormalized posterior
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:08/07/2010
!^**********************************************************************
!^* Comments: mv handled out of this sub
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.teta, par vector
!^*    2.X, data
!^*    3.distID, distribution of data
!^*    4.PriorList, list of priors & their par (see PriorListType above)
!^*    5.[cov], covariate for distId
!^* OUT
!^*    1.lp, log-post
!^*    2.feas,  feasible?
!^*    3.isnull, is log-post=0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************

use Distribution_tools, only: GetSamplePdf, GetParNumber,covariate

real(mrk), intent(in)::teta(:), X(:)
character(*), intent(in)::distID
Type(PriorListType), intent(in)::PriorList(:)
type(covariate), intent(in), optional::cov(:)
real(mrk), intent(out)::lp
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
integer(mik):: npar
real(mrk)::prior, lkh

err=0;mess='';feas=.true.;isnull=.false.;lp=undefRN

! check size
CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
    mess='GetLogPost:'//trim(mess)
    feas=.false.;return
endif
if(size(teta)/=npar) then
    feas=.false.;err=1;mess="GetLogPost:FATAL:size mismatch [teta]";
    return
endif

! get log-prior
call GetLogPrior(teta,PriorList,prior, feas, isnull,err,mess)
If(err>0) then
    mess='GetLogPost:'//trim(mess)
    feas=.false.;return
endif
if(.not. feas) then
    !err=1;mess="GetLogPost:FATAL:unfeasible [PriorList(i)%par]"
    return
endif
if(isnull) return

! get lkh
call GetSamplePdf(DistId=DistId,x=X,par=teta,&
            loga=.true.,cov=cov,pdf=lkh,feas=feas,isnull=isnull,err=err,mess=mess)
If(err>0) then
    mess='GetLogPost:'//trim(mess)
    feas=.false.;return
endif
if(.not. feas) then
    !err=1;mess="GetLogPost:FATAL:unfeasible [PriorList(i)%par]"
    return
endif
if(isnull) return

lp=prior+lkh

end subroutine GetLogPost

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Posterior_wrapper(x,feas,isnull,fx,fAux,err,mess)

!^**********************************************************************
!^* Purpose: wrapper to GetLogPost to comply with MCMC interface
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:08/07/2010
!^**********************************************************************
!^* Comments: mv handled out of here
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, teta
!^* OUT
!^*    1. see GetLogPost
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^* INOUT
!^*    1.
!^*    2.
!^*    3.
!^**********************************************************************
use Distribution_tools, only:IsCovNeeded

real(mrk),intent(in)::x(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess

if(IsCovNeeded(Xdist)) then
    if(.not.allocated(Xcov)) then
        err=1;mess="Posterior_wrapper:FATAL:Covariate needed for this distribution";return
    endif
    call GetLogPost(teta=x,X=XData, distID=Xdist, PriorList=XpriorList,cov=Xcov,lp=fx, &
                feas=feas, isnull=isnull,err=err,mess=mess)
else
    call GetLogPost(teta=x,X=XData, distID=Xdist, PriorList=XpriorList,lp=fx, &
                feas=feas, isnull=isnull,err=err,mess=mess)
endif
if(present(fAux)) fAux=UndefRN

end subroutine Posterior_wrapper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Bayes_EstimPar(X,distID,mv,PriorList,cov,& ! Data & Model
                ! Tuning of the MCMC sampler
                start,startStd,&
                nAdapt,nCycles,&
                MinMoveRate,MaxMoveRate,&
                DownMult,UpMult,&
                OutFile, headers,&
                ! error handling
                err,mess)

!^**********************************************************************
!^* Purpose: Parameter estimation using Bayes-MCMC
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:08/07/2010
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.X, data (1-d vector)
!^*    2.distID, distribution (see distribution_tools module)
!^*    3.mv, code for missing values
!^*    4.PriorList, list of priors & their par (see PriorListType above)
!^*    5.[cov], covariate for distId
!^*    5.start, starting point
!^*    6.startStd, starting std for each marginal jump distribution
!^*    7.nAdapt, number of iterations before adapting the jump std
!^*    8.nCycles, number of adaptations (total number of iteration is hence nAdapt*nCycles
!^*        If nCycles<1, user will be asked stop/continue
!^*    9.MinMoveRate, objective move rate (lower bound)
!^*    10.MaxMoveRate, objective move rate (upper bound)
!^*    11.DownMult, value multiplying the jump std whan move rate is too low
!^*    12.UpMult, value multiplying the jump std whan move rate is too high
!^*    13.OutFile, address of output file for MCMC samples
!^*    14.[heaaders], headers of OutFile
!^* OUT
!^*    1.
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^* INOUT
!^*    1.
!^*    2.
!^*    3.
!^**********************************************************************

use MCMCStrategy_tools, only:Adaptive_Metro_OAAT

real(mrk), intent(in)::X(:), mv, MinMoveRate,MaxMoveRate,DownMult,UpMult
real(mrk), intent(inout)::start(:), startStd(:)
Character(*), intent(in)::distID
Type(PriorListType), intent(in)::PriorList(:)
type(covariate), intent(in), optional::cov(:)
integer(mik), intent(in)::nAdapt,nCycles                
character(*), intent(in)::OutFile
character(*), intent(in),optional::headers(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::n0,n,npar,i
logical, allocatable::mvlist(:)
real(mrk)::fx

err=0;mess=''

! handle mv
n0=size(X)
if(allocated(mvlist)) deallocate(mvlist)
allocate(mvlist(n0))
mvlist=(X/=mv)
n=count(mvlist)
if(n==0) then
    err=1;mess='Bayes_EstimPar:FATAL: no non-missing values!';return
endif

! populate global variables 
if(allocated(Xdata)) deallocate(Xdata)
allocate(Xdata(n))
Xdata=pack(X,mvlist)
Xdist=distID
npar=size(PriorList)
if(allocated(XpriorList)) deallocate(XpriorList)
allocate(XpriorList(npar))
XpriorList=PriorList
! optional covariate
if(IsCovNeeded(Xdist)) then
    if(.not.present(cov)) then
        err=1;mess="Bayes_EstimPar:FATAL:Covariate needed for this distribution";return
    else
        if(allocated(Xcov)) deallocate(Xcov)
        allocate(Xcov(size(cov)))
        do i=1,size(cov)
            if(associated(Xcov(i)%val)) nullify(Xcov(i)%val)
            allocate(Xcov(i)%val(n))
            Xcov(i)%val=pack(cov(i)%val,MVlist)
        enddo
    endif
endif

! MCMC
call Adaptive_Metro_OAAT(f=Posterior_wrapper,x=start,&
                fx=fx,std=startStd,&
                nAdapt=nAdapt,nCycles=nCycles,&
                MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
                DownMult=DownMult,UpMult=UpMult,&
                OutFile=OutFile,headers=headers,err=err,mess=mess)

end subroutine Bayes_EstimPar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Bayes_Predictive(File,Nburn,Nread,Nslim,& ! IN: Read properties
                                Nrep, distID, & ! IN: number of reps for EACH MCMC sample + data dist
                                modal, pred,& ! OUT: modal par estimate & replicates from the predictive
                                err, mess)

!^**********************************************************************
!^* Purpose: Generate deviates from the predictive distribution
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:08/07/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. File, containing MCMC samples
!^*    2. Nburn, number of lines to discard
!^*    3. Nread, number of lines to read
!^*    4. Nslim, only one line every Nslim will be used
!^*    5. Nrep, for each line, Nrep deviates will be generated
!^*(advise: keep Nrep ~ Nslim - Nrep>>Nslim will yield a huge number of replicates!)
!^* OUT
!^*    1.modal, modal parameter estimate
!^*    2.pred, deviates from the predictive
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************

use Distribution_tools, only: GetParNumber, Generate

integer(mik), intent(in)::Nburn, Nread, Nslim, Nrep
character(*),intent(in)::File, distID
real(mrk), intent(out)::modal(:)
real(mrk), pointer::pred(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,k,errcode,npar,n
real(mrk), allocatable::Temp(:,:), mcmc(:,:)

err=0;mess=''

if(IsCovNeeded(Xdist)) then ! Warn but try to go anyway!
   err=-1;
   mess="Bayes_Predictive:WARNING:Dist needing covariates not yet &
        &propely handled by this subroutine" 
   write(*,*) trim(mess)
   !;return
endif

open(unit=1,file=File,status='old')
read(1,*) !headers
do i=1,Nburn ! burn
    read(1,*,iostat=errcode)
    if(errcode/=0) then ! problem reading file
        err=1
        mess='Bayes_Predictive:FATAL:problem reading file during burnin - cannot proceed!'
        return
    endif
enddo

! Get number of col and allocate Temp
CALL GetParNumber(distID, npar, err, mess)
If(err>0) then
    mess='Bayes_Predictive:'//trim(mess)
    return
endif
if(allocated(Temp)) deallocate(Temp)
allocate(Temp(Nread,npar+1))

do i=1,Nread ! read used lines
    k=i
    read(1,*,iostat=errcode) Temp(i,:)
    if(errcode/=0) then ! problem reading file
        err=-1;k=k-1
        mess='Bayes_Predictive:WARNING:problem reading file before [Nread] lines!'
        EXIT
    endif
enddo
close(1)

! Slim
n=k/Nslim
if(allocated(mcmc)) deallocate(mcmc)
allocate(mcmc(n,Npar+1))
mcmc=Temp(1:n:nSlim,:)

! get modal estimate
if(size(modal)/=npar) then
    err=1; mess='Bayes_Predictive:SEVERE:size mismath [modal]'
    modal=undefRN
else
    call GetMode(mcmc(:,1:npar),mcmc(:,npar+1),modal,err,mess)
    If(err>0) then
      mess='Bayes_Predictive:'//trim(mess);return
    endif
endif

! get deviates from predictive
call GeneratePred(mcmc(:,1:npar),distID,nrep,pred,err,mess)
If(err>0) then
  mess='Bayes_Predictive:'//trim(mess);return
endif

end subroutine Bayes_Predictive

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine LoadMCMCsample(File,BurnFactor,Nslim,distID,npar,sep,& ! IN: Read properties
                                mcmc,LogPost,err,mess)

!^**********************************************************************
!^* Purpose: Load existing MCMC sample
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
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
!^*    1.File
!^*    2.BurnFactor
!^*    3.Nslim
!^*    4.distID
!^*    5.[npar], optional: if present, passes the number of inferred parameter
!^*                        if not present, deduced from distID
!^*    6.[seo], optional: separator
!^* OUT
!^*    1.mcmc
!^*    2.[LogPost]
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:GetSpareUnit
use DataRW_tools,only:DatRead,ReadSeparatedFile
use Distribution_tools, only:GetParNumber

integer(mik), intent(in)::Nslim
character(*),intent(in)::File,distID
real(mrk), intent(in)::BurnFactor
integer(mik), intent(in),optional::npar
character(*),intent(in),optional::sep
real(mrk),pointer::mcmc(:,:)
real(mrk),pointer, optional::LogPost(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::i,Nburn,Nread,N,np
real(mrk), pointer::y(:,:)
integer(mik), allocatable::rows(:)
character(250), allocatable::headers(:)

err=0;mess=''

CALL GetParNumber(trim(distID), np, err, mess)
If(err>0) then
    mess='LoadMCMCsample:'//trim(mess)
    return
endif
if(present(npar)) np=npar
if(allocated(headers)) deallocate(headers)
allocate(headers(np+1))
! Read MCMC file
if(present(sep)) then
    call ReadSeparatedFile(file=trim(File),sep=sep,nhead=1,ncol=np+1,y=y,err=err,mess=mess)
else
    call DatRead(file=trim(File),ncol=np+1,y=y,headers=headers,err=err,mess=mess)
endif
If(err>0) then
  mess='LoadMCMCsample:'//trim(mess);return
endif
N=size(y,dim=1)
Nburn=NINT(N*BurnFactor)
Nread=size((/(i,i= Nburn+1,N,Nslim)/))
if(allocated(rows)) deallocate(rows)
allocate(rows(Nread))
rows=(/(i,i= Nburn+1,N,Nslim)/)
Nread=size(rows)
if(associated(mcmc)) nullify(mcmc)
allocate(mcmc(Nread,Np))
mcmc=y(rows,1:Np)
if(present(LogPost)) then
  if(associated(LogPost)) nullify(LogPost)
  allocate(LogPost(Nread))
  LogPost=y(rows,Np+1)
endif

end subroutine LoadMCMCsample

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GetMCMCSummary(mcmc,LogPost,OutFile,parnames,&
                          modal,median,mean,&
                          sdev,IC90,correl,&
                          err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea LYON
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
!^*    1.mcmc
!^*    2.LogPost
!^*    3.[OutFile]
!^*    4.[parnames]
!^* OUT
!^*    1.[modal]
!^*    2.[median]
!^*    3.[mean]
!^*    4.[sdev]
!^*    5.[IC90]
!^*    6.[correl]
!^*    7.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    8.mess, error message
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalStats,GetEmpiricalQuantile
use numerix_dmsl_kit, only: getvar
use utilities_dmsl_kit, only:GetSpareUnit,number_string

real(mrk),intent(in)::mcmc(:,:),LogPost(:)
character(*), intent(in), optional::OutFile,parnames(:)
real(mrk), intent(out), optional::modal(:),median(:),mean(:),&
                               sdev(:),IC90(:,:),correl(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk), parameter::lop=0.05_mrk,hip=0.95_mrk
integer(mik)::n,p,i,j,unt
real(mrk)::mode(size(mcmc,dim=2)),average(size(mcmc,dim=2)),&
           med(size(mcmc,dim=2)),sd(size(mcmc,dim=2)),&
           ic(size(mcmc,dim=2),2),&
           cov(size(mcmc,dim=2),size(mcmc,dim=2)),&
           cor(size(mcmc,dim=2),size(mcmc,dim=2))
character(250)::left(size(mcmc,dim=2)),fmt

n=size(mcmc,dim=1)
p=size(mcmc,dim=2)
! Get modal estimate
call GetMode(mcmc,LogPost,mode,err,mess)
If(err>0) then
  mess='GetMCMCSummary:'//trim(mess);return
endif
! Get standard stats
do i=1,p
  call GetEmpiricalStats(x=mcmc(:,i), & !data,
                mean=average(i),median=med(i),&
                std=sd(i),err=err,mess=mess)
  If(err>0) then
    mess='GetMCMCSummary:'//trim(mess);return
  endif
  call GetEmpiricalQuantile(p=lop,x=mcmc(:,i),q=ic(i,1),err=err,mess=mess)
  If(err>0) then
    mess='GetMCMCSummary:'//trim(mess);return
  endif
  call GetEmpiricalQuantile(p=hip,x=mcmc(:,i),q=ic(i,2),err=err,mess=mess)
  If(err>0) then
    mess='GetMCMCSummary:'//trim(mess);return
  endif
enddo
! Get Correlation matrix
cov=getvar(x=mcmc,package=2,method="c")
forall(i=1:p,j=1:p) cor(i,j)=cov(i,j)/(sd(i)*sd(j))

! populate optional outputs
if(present(modal)) modal=mode
if(present(median)) median=med
if(present(mean)) mean=average
if(present(sdev)) sdev=sd
if(present(IC90)) IC90=ic
if(present(correl)) correl=cor

if(present(OutFile)) then
  call GetSpareUnit(unt,err,mess)
  If(err>0) then
    mess='GetMCMCSummary'//trim(mess);return
  endif
  open(unit=unt,file=trim(OutFile),status="REPLACE",iostat=err)
  If(err/=0) then
    mess='GetMCMCSummary'//trim(mess);return
  endif
  ! get lefters
  if(present(parnames)) then
    left=parnames
  else
    do i=1,p
      left(i)="par"//trim(number_string(i))
    enddo
  endif
  ! Write parameters stats
  write(unt,'(7A12)') "parameter","mode","median","mean","stdev","q0.05","q0.95"
  do i=1,p
    write(unt,'(A12,6e12.3)') trim(left(i)),mode(i),&
                    med(i),average(i),&
                    sd(i),ic(i,1),ic(i,2)
  enddo
  ! Write correl matrix
  write(unt,*) "Posterior Correlation Matrix"
  fmt='(' // trim(number_string(p)) // 'F7.2)'
  do i=1,p
    write(unt,trim(fmt)) cor(i,:)
  enddo
  close(unt)
endif

end subroutine GetMCMCSummary
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GetModalEstimates(mcmc,LogPost,distID,&
                          mode,&
                          x,pdf,cdf,&
                          p,Qp,IC,&
                          level,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea LYON
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
!^*    1.mcmc
!^*    2.LogPost
!^*    3.x, points where pdf/cdf are computed
!^*    4.p, points where quantiles are computed
!^*    5.[level], confidence level for IC (default 0.9)
!^* OUT
!^*    1.mode
!^*    2.pdf
!^*    3.cdf
!^*    4.Qp
!^*    5.IC
!^*    6.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    7.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetQuantile, Getpdf, Getcdf
use EmpiricalStats_tools, only: GetEmpiricalQuantile
real(mrk),intent(in)::mcmc(:,:),LogPost(:)
real(mrk), intent(in)::x(:),p(:)
real(mrk), intent(in), optional::level
character(*), intent(in)::distID
real(mrk), intent(out)::mode(size(mcmc,dim=2)),pdf(size(x)),cdf(size(x)),&
                        Qp(size(p)),IC(size(p),2)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk)::ConfLevel,lop,hip
integer(mik)::n,i,j,nmcmc
logical::feas,isnull
real(mrk)::temp(size(mcmc,dim=1))

err=0;mess='';mode=undefRN;pdf=undefRN;cdf=undefRN;Qp=undefRN;IC=undefRN
if(present(level)) then
    ConfLevel=level
else
    ConfLevel=0.9_mrk
endif
lop=0.5_mrk*(1._mrk-ConfLevel)
hip=1._mrk-lop

nmcmc=size(mcmc,dim=1)
! mode
call GetMode(mcmc,LogPost,mode,err,mess)
If(err>0) then
  mess='GetModalEstimates:'//trim(mess);return
endif
! compute pdf/cdf
n=size(x)
do i=1,n
    call GetPdf(DistId=DistId,x=x(i),par=mode,&
                loga=.false.,pdf=pdf(i),&
                feas=feas,isnull=isnull,err=err,mess=mess)
    If(err>0) then
      mess='GetModalEstimates:'//trim(mess);return
    endif
    call GetCdf(DistId=DistId,x=x(i),par=mode,cdf=cdf(i),&
                feas=feas,err=err,mess=mess)
    If(err>0) then
      mess='GetModalEstimates:'//trim(mess);return
    endif
enddo
! Compute quantile
n=size(p)
do i=1,n
    ! modal quantiles
    call GetQuantile(DistId=DistId,p=p(i),par=mode,q=Qp(i),&
                feas=feas,err=err,mess=mess)
    If(err>0) then
      mess='GetModalEstimates:'//trim(mess);return
    endif
    ! IC computation
    do j=1,nmcmc
        call GetQuantile(DistId=DistId,p=p(i),par=mcmc(j,:),q=temp(j),&
                    feas=feas,err=err,mess=mess)
        If(err>0) then
          mess='GetModalEstimates:'//trim(mess);return
        endif
    enddo
    call GetEmpiricalQuantile(p=lop,x=temp,q=IC(i,1),err=err,mess=mess)
    If(err>0) then
      mess='GetModalEstimates:'//trim(mess);return
    endif
    call GetEmpiricalQuantile(p=hip,x=temp,q=IC(i,2),err=err,mess=mess)
    If(err>0) then
      mess='GetModalEstimates:'//trim(mess);return
    endif
enddo

end subroutine GetModalEstimates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GetPredictiveEstimates(pred,&
                          x,pdf,cdf,&
                          p,Qp,&
                          ppFormula, &! optional, formula for plotting position [default Hazen]
                          kernelID,h,& !optional, parameters of kde [defaults Gaussian + built_in h value]
                          err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea LYON
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
!^*    1.pred, replicates from the predictive distribution
!^*    2.x, points where pdf/cdf are computed
!^*    3.p, points where quantiles are computed
!^*    4.[ppFormula]
!^*    5.[kernelID]
!^*    6.[h]
!^* OUT
!^*    1.pdf
!^*    2.cdf
!^*    3.Qp
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use EmpiricalStats_tools, only: GetEmpiricalQuantile,KDE,GetEmpiricalPval

real(mrk), intent(in)::pred(:),x(:),p(:)
real(mrk), intent(out)::pdf(size(x)),cdf(size(x)),Qp(size(p))
real(mrk), intent(in), optional::h
character(*), intent(in), optional::ppFormula, kernelID
character(*), intent(out)::mess
integer(mik), intent(out)::err
!locals
integer(mik)::i,n
err=0;mess='';pdf=undefRN;cdf=undefRN;Qp=undefRN

! Get pdf
n=size(x)
do i=1,n
    call KDE(kernelID,x(i),pred,h,pdf(i),err,mess)
    if(err/=0) then
        mess='GetPredictiveEstimates:'//trim(mess);return
    endif
enddo
! Get cdf
do i=1,n
    call GetEmpiricalPval(obs=x(i),x=pred,ppFormula=ppFormula,&
                           pv=cdf(i),err=err,mess=mess)
    if(err>0) then
        mess='GetPredictiveEstimates:'//trim(mess);return
    endif
enddo
! Get quantile
n=size(p)
do i=1,n
    call GetEmpiricalQuantile(p=p(i),x=pred,ppFormula=ppFormula,&
                              q=Qp(i),err=err,mess=mess)
    if(err>0) then
        mess='GetPredictiveEstimates:'//trim(mess);return
    endif
enddo
end subroutine GetPredictiveEstimates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GetMode(mcmc,LogPost,mode,err,mess)
real(mrk),intent(in)::mcmc(:,:),LogPost(:)
real(mrk), intent(out)::mode(size(mcmc,dim=2))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::ml(1)

err=0;mess='';mode=UndefRN
if(size(LogPost)/=size(mcmc,dim=1)) then
    err=1;mess='GetMode::size mismatch[LogPost,mcmc]';return
endif
ml=Maxloc(LogPost)
mode=mcmc(ml(1),:)
end subroutine GetMode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GeneratePred(mcmc,distID,nrep,pred,err,mess)
! get deviates from predictive
use Distribution_tools, only:Generate
real(mrk),intent(in)::mcmc(:,:)
character(*), intent(in)::distID
integer(mik), intent(in)::nrep
real(mrk), pointer::pred(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::n,npred,i,j,k
logical::feas

err=0;mess=''
n=size(mcmc,dim=1)
npred=n*nrep
if(associated(pred)) nullify(pred)
allocate(pred(npred))
k=0
do i=1,n
    do j=1,nrep
        k=k+1
        call Generate(DistId=distID,par=mcmc(i,:),&
                    gen=pred(k),feas=feas,err=err,mess=mess)
        If(err>0) then
            mess='GeneratePred: '//trim(mess);return
        endif
        if(.not. feas) then ! shouldn't happen - MCMC should NOT generate infeasible par!
            err=1;
            mess='GeneratePred: unfeasible parameter!';return
        endif
    enddo
enddo
end subroutine GeneratePred

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine CreateFlatPriorList(npar,PriorList,err,mess)

!^**********************************************************************
!^* Purpose: Create a priorlist with improper flat priors for all npar parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:19/09/2013
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.npar
!^* OUT
!^*    1.PriorList
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************

integer(mik), intent(in)::npar
Type(PriorListType), intent(out)::PriorList(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::i
err=0;mess=''
! size check
if(size(PriorList)/=npar) then
    err=1;mess='CreateFlatPriorList:FATAL:size mismatch [PriorList,npar]';return
endif
do i=1,npar
    PriorList(i)%dist='FlatPrior'
    if(allocated(PriorList(i)%par)) deallocate(PriorList(i)%par);allocate(PriorList(i)%par(0))
enddo
end subroutine CreateFlatPriorList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module BayesianEstimation_tools

module SuniBay_tools

!~**********************************************************************
!~* Purpose: SuniBay-specific subroutines
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref
!~**********************************************************************
!~* Last modified:27/12/2010
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. SuniBay_PostProcess: post-process MCMC file
!~*		2. RedData: guess what? read data
!~*		3. GoodBye,Welcome: politness messages
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SuniBay_PostProcess,GoodBye,Welcome, ReadData

Contains

subroutine SuniBay_PostProcess(MCMCFile,Nburn,Nread,Nslim,Nrep,& ! Read properties
                        obs,indx,Ncov,cov,& ! data
                        obsHeaders, IndxHeaders, CovHeaders, & ! headers
                        DistID,nteta,&    ! Distribution ID 
                        DataFile, & ! File containing Data (indx and obs only)
                        DataEmpStatFile, & ! File containing empirical statistics of data
                        MCMCEmpStatFile, & ! File containing empirical statistics of MCMC samples
                        PCdfFile, & ! File containing pdf-cdf evaluations
                        SpaghettiFile, & ! File containing Spaghetti Quantile plots
                        QuantilePlotFile, & ! File containing Quantile plots
                        PredictiveFile, & ! File containing replicates from the predictive distribution
                        err,mess)

!^**********************************************************************
!^* Purpose: 
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
use EmpiricalStats_tools, only:GetEmpiricalStats,GetEmpiricalQuantile,GetEmpiricalPval
use utilities_dmsl_kit, only: getSpareUnit,cumsum,number_string
use Distribution_tools, only:covariate, GetParName, GetQuantile, GetPdf,GetCdf
use BayesianEstimation_tools, only:Bayes_Predictive

character(*), intent(in)::MCMCFile,DataEmpStatFile,MCMCEmpStatFile,&
                          PCdfFile,SpaghettiFile,QuantilePlotFile,&
                          DataFile,PredictiveFile
character(*), intent(in)::obsHeaders,IndxHeaders,CovHeaders(:)
integer(mik), intent(in)::Ncov,Nteta,Nburn,Nread,Nslim,Nrep
real(mrk), intent(in)::obs(:),indx(:)
type(covariate), intent(in)::cov(:)
character(*), intent(in)::DistID
integer(mik), intent(out)::err        
character(*),intent(out)::mess                   
!locals
real(mrk), parameter::pLo=0.05_mrk, pHi=0.95_mrk

character(250), dimension(16),parameter::lefters=(/"N              ",&
                                                   "Minimum        ",&
                                                   "Maximum        ",&
                                                   "Range          ",&
                                                   "Mean           ",&
                                                   "Median         ",&
                                                   "Q10%           ",&
                                                   "Q25%           ",&
                                                   "Q75%           ",&
                                                   "Q90%           ",&
                                                   "St.Dev.        ",&
                                                   "Variance       ",&
                                                   "CV             ",&
                                                   "Skewness       ",&
                                                   "Kurtosis       ",&
                                                   "MaxPost        "/)
character(15), dimension(5),parameter::PCdfHeaders=(/"x              ",&
                                                     "modal_pdf(x)   ",&
                                                     "modal_cdf(x)   ",&
                                                     "pred_pdf(x)    ",&
                                                     "pred_cdf(x)    "/)
character(15), dimension(5),parameter::QplotHeaders=(/"T              ",&
                                                      "Q(T)_maxpost   ",&
                                                      "Q(T)_CI_low    ",&
                                                      "Q(T)_CI_Hi     ",&
                                                      "Q(T)_pred      "/)

integer(mik), parameter::nprobs=63
real(mrk), dimension(nprobs)::probs=&
    (/0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,&
      0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,&
      0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,&
      0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,&
      0.991,0.992,0.993,0.994,0.995,0.996,0.997,0.998,0.999,&
      0.9991,0.9992,0.9993,0.9994,0.9995,0.9996,0.9997,0.9998,0.9999,&
      0.99991,0.99992,0.99993,0.99994,0.99995,0.99996,0.99997,0.99998,0.99999/)
integer(mik)::unt,Nmcmc,k,i,j,errcode,n,nH
real(mrk)::Ts(nprobs),Qs_maxpost(nprobs),Qs_pred(nprobs),&
        pdf_pred(nprobs),cdf_pred(nprobs),&
        pdf_maxpost(nprobs),cdf_maxpost(nprobs)
real(mrk), allocatable::Temp(:,:),mcmc(:,:),Summary(:,:),spag(:,:),Qs_CI(:,:)
character(250), allocatable::headers(:),SpagHeaders(:)
real(mrk)::ml(1), modal(nteta)
real(mrk), pointer::pred(:)
logical feas, isnull

!init
err=0;mess='';

! Size checks
if( size(obs)/=size(Indx) ) then
    err=1;mess='SuniBay_PostProcess: Fatal: size mismatch in input data';return
endif

! Rewrite data
call getSpareUnit(unt,err,mess)
If(err>0) then
    mess='SuniBay_PostProcess: '//trim(mess);return
endif
open(unit=unt,file=trim(DataFile),status='replace')
nH=2
write(unt,'(<nH>A15)') IndxHeaders,obsheaders
do i=1,size(obs)
    write(unt,'(<nH>e15.6)') indx(i), obs(i)
enddo
close(unt)

! Read MCMC file
call getSpareUnit(unt,err,mess)
If(err>0) then
    mess='SuniBay_PostProcess: '//trim(mess);return
endif
open(unit=unt,file=trim(MCMCFile),status='old')
read(unt,*) !headers
do i=1,Nburn ! burn
	read(unt,*,iostat=errcode)
	if(errcode/=0) then ! problem reading file
		err=1
		mess='SuniBay_PostProcess:FATAL:problem reading file during burnin - cannot proceed!'
		return
	endif
enddo
if(allocated(Temp)) deallocate(Temp)
allocate(Temp(Nread,nteta+1))
do i=1,Nread ! read used lines
	k=i
	read(unt,*,iostat=errcode) Temp(i,:)
	if(errcode/=0) then ! problem reading file
		err=-1;k=k-1
		mess='SuniBay_PostProcess:WARNING:problem reading file before [Nread] lines!'
		EXIT
	endif
enddo
close(unt)
! Slim
n=k/Nslim
if(allocated(mcmc)) deallocate(mcmc)
allocate(mcmc(n,nteta+1))
mcmc=Temp(1:n:nSlim,:)
nMCMC=size(mcmc,dim=1)

! Data Empirical Stats File
if(allocated(Summary)) deallocate(Summary)
allocate(Summary(15,2+nCov))
call GetEmpiricalStats(x=indx, all15=Summary(:,1),err=err,mess=mess)
if(err>0) then
    mess="SuniBay_PostProcess: "//trim(mess);return
endif
call GetEmpiricalStats(x=obs, all15=Summary(:,2),err=err,mess=mess)
if(err>0) then
    mess="SuniBay_PostProcess: "//trim(mess);return
endif
do i=1,nCov
    call GetEmpiricalStats(x=cov(i)%val, all15=Summary(:,2+i),err=err,mess=mess)
    if(err>0) then
        mess="SuniBay_PostProcess: "//trim(mess);return
    endif
enddo
!Write2File
call getSpareUnit(unt,err,mess)
If(err>0) then
    mess='SuniBay_PostProcess: '//trim(mess);return
endif
open(unit=unt,file=trim(DataEmpStatFile),status='replace')
nH=1+2+nCov
write(unt,'(<nH>A15)') "Stat           ",IndxHeaders,obsheaders,CovHeaders
do i=1,15
    write(unt,'(A15,<nH>e15.6)') lefters(i), Summary(i,:)
enddo
close(unt)

! MCMC empirical stats
if(allocated(Summary)) deallocate(Summary)
if(allocated(headers)) deallocate(headers)
allocate(Summary(15,nteta),headers(nteta))
do i=1,nteta
    call GetEmpiricalStats(x=mcmc(:,i), all15=Summary(:,i),err=err,mess=mess)
    if(err>0) then
        mess="SuniBay_PostProcess: "//trim(mess);return
    endif
enddo
ml=maxloc(mcmc(:,nteta+1))
! Create Headers
Call GetParName(DistID=DistID, name=headers, err=err, mess=mess)
If(err>0) then
    mess='SuniBay_PostProcess: '//trim(mess);return
endif
!Write2File
call getSpareUnit(unt,err,mess)
If(err>0) then
    mess='SuniBay_PostProcess: '//trim(mess);return
endif
open(unit=unt,file=trim(MCMCEmpStatFile),status='replace')
nh=nteta+1
write(unt,'(<nH>A15)') "Stat           ",headers
do i=1,15
    write(unt,'(A15,<nH>e15.6)') lefters(i), Summary(i,:)
enddo
write(unt,'(A15,<nH>e15.6)') lefters(16), mcmc(ml(1),1:(nh-1))
close(unt)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! IMPORTANT NOTE: most diagnostics below should not work with 
! distributions requiring covariates, except XXX_HISTO
! Specialized diagnostics will be implemented later on
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Get Predictive Distribution
Call Bayes_Predictive(File=trim(MCMCFile),Nburn=Nburn,Nread=Nread,Nslim=Nslim,& ! IN: Read properties
								Nrep=Nrep, distID=distID, & ! IN: number of reps for EACH MCMC sample + data dist
								modal=modal, pred=pred,& ! OUT: modal par estimate & replicates from the predictive
								err=err, mess=mess)

If(err>0) then
    mess='SuniBay_PostProcess:'//trim(mess);return
endif

!Write2File
call getSpareUnit(unt,err,mess)
If(err>0) then
    mess='SuniBay_PostProcess: '//trim(mess);return
endif
open(unit=unt,file=trim(PredictiveFile),status='replace')
write(unt,'(A15)') "Replicates"
do i=1,size(pred)
    write(unt,'(e15.6)') pred(i)
enddo
close(unt)

! Compute x's where pdf & cdf will be evaluated
do i=1,nprobs
    call GetQuantile(DistId=DistId,p=probs(i),par=modal,q=Qs_maxpost(i),&
            feas=feas,err=err,mess=mess)
    If(err>0) then
        mess='SuniBay_PostProcess: '//trim(mess);return
    endif
    If(.not.feas) then
        mess='SuniBay_PostProcess: GetQuantile: Fatal: Unfeasible parameters';return
    endif
enddo

! Compute pdfs
do i=1, nprobs
    ! From modal estimate
    call GetPdf(DistId=DistId,x=Qs_maxpost(i),par=modal,loga=.false.,&
                pdf=pdf_maxpost(i),feas=feas,isnull=isnull,err=err,mess=mess)
    If(err>0) then
        mess='SuniBay_PostProcess: '//trim(mess);return
    endif
    call GetCdf(DistId=DistId,x=Qs_maxpost(i),par=modal,&
                cdf=cdf_maxpost(i),feas=feas,err=err,mess=mess)
    If(err>0) then
        mess='SuniBay_PostProcess: '//trim(mess);return
    endif
    !from predictive
    Call GaussianKDE(x=Qs_maxpost(i),sample=pred,pdf=pdf_pred(i),err=err,mess=mess)
    If(err>0) then
        mess='SuniBay_PostProcess: '//trim(mess);return
    endif
    Call GetEmpiricalPval(obs=Qs_maxpost(i),x=pred,pv=cdf_pred(i),err=err,mess=mess)
    If(err>0) then
        mess='SuniBay_PostProcess: '//trim(mess);return
    endif
enddo
!Write2File
call getSpareUnit(unt,err,mess)
If(err>0) then
    mess='SuniBay_PostProcess: '//trim(mess);return
endif
open(unit=unt,file=trim(PCdfFile),status='replace')
write(unt,'(5A15)') PCDFheaders
do i=1,nprobs
    write(unt,'(5e15.6)') Qs_maxpost(i),pdf_maxpost(i),&
                        cdf_maxpost(i),pdf_pred(i),cdf_pred(i)
enddo
close(unt)

! Spaghetti quantile plot
Ts=1._mrk/real(1-probs,mrk)
if(allocated(spag)) deallocate(spag)
allocate(spag(nprobs,nmcmc))
if(allocated(spagHeaders)) deallocate(spagHeaders)
allocate(spagHeaders(nmcmc))

do i=1,Nmcmc
    spagHeaders(i)="QT_"//trim(number_string(i))
    do j=1,nprobs
        call GetQuantile(DistId=DistId,p=probs(j),par=mcmc(i,1:nteta),&
            q=spag(j,i),feas=feas,err=err,mess=mess)
    enddo
    If(err>0) then
        mess='SuniBay_PostProcess: '//trim(mess);return
    endif
    If(.not.feas) then
        mess='SuniBay_PostProcess: GetQuantile: Fatal: Unfeasible parameters';return
    endif
enddo
!Write Spaghetti files
call getSpareUnit(unt,err,mess)
If(err>0) then
    mess='SuniBay_PostProcess: '//trim(mess);return
endif
open(unit=unt,file=trim(SpaghettiFile),status='replace')
nH=1+nMCMC
write(unt,'(<nH>A15)') "T              ",spagHeaders
do i=1,nprobs
    write(unt,'(<nH>e15.6)') Ts(i), spag(i,:)
enddo
close(unt)

! Summary quantile plot
! get Confidence intervals from spag & predictive quantiles
if(allocated(Qs_CI)) deallocate(Qs_CI)
allocate(Qs_CI(nmcmc,2))
do i=1,nprobs
    call GetEmpiricalQuantile(p=pLo,x=Spag(i,:),q=Qs_CI(i,1),err=err,mess=mess)
    if(err>0) then
        mess="SuniBay_PostProcess: "//trim(mess);return
    endif
    call GetEmpiricalQuantile(p=pHi,x=Spag(i,:),q=Qs_CI(i,2),err=err,mess=mess)
    if(err>0) then
        mess="SuniBay_PostProcess: "//trim(mess);return
    endif
    call GetEmpiricalQuantile(p=probs(i),x=pred,q=Qs_pred(i),err=err,mess=mess)
    if(err>0) then
        mess="SuniBay_PostProcess: "//trim(mess);return
    endif

enddo
!Write To File
call getSpareUnit(unt,err,mess)
If(err>0) then
    mess='SuniBay_PostProcess: '//trim(mess);return
endif
open(unit=unt,file=trim(QuantilePlotFile),status='replace')
nH=5
write(unt,'(<nH>A15)') QplotHeaders
do i=1,nprobs
    write(unt,'(<nH>e15.6)') Ts(i), Qs_maxpost(i),Qs_CI(i,:),Qs_pred(i)
enddo
close(unt)

end subroutine SuniBay_PostProcess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ReadData(DataFile,nrow,ncol,nHeader,X,headers,err,mess)

!^**********************************************************************
!^* Purpose: Read matrix-like data file
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:04/10/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.DataFile
!^*		2.nrow (WITHOUT header lines)
!^*		3.ncol
!^*		4.nHeader, number of skipped header lines
!^* OUT
!^*		1.X
!^*		2.headers, headers
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:number_string

character(*), intent(in)::DataFile
integer(mik), intent(in)::nrow,ncol,nHeader
real(mrk), intent(out)::X(nrow,ncol)
character(250), intent(out)::headers(ncol)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i

err=0;mess='';X=undefRN

open(unit=1,file=DataFile, status='old',iostat=err)
if(err/=0) then
    mess="ReadData: Fatal: problem opening file";return
endif

do i=1,nHeader
    read(1,*,iostat=err) headers
enddo

if(err/=0 .or. nHeader/=1) then
    do i=1,ncol
        headers(i)="Var"//trim(number_string(i))
    enddo
endif

do i=1,nrow
    read(1,*,iostat=err) X(i,:)
    if(err/=0) then
        mess="ReadData: Fatal: problem reading file, row: "//trim(number_string(i));return
    endif
enddo
close(1)

end subroutine ReadData

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Welcome

write(*,*) " Welcome To...                        "
write(*,*) "  ______                       __        _______                       "
write(*,*) " /      \                     |  \      |       \                      "
write(*,*) "|  @@@@@@\ __    __  _______   \@@      | @@@@@@@\  ______   __    __  "
write(*,*) "| @@___\@@|  \  |  \|       \ |  \      | @@__/ @@ |      \ |  \  |  \ "
write(*,*) " \@@    \ | @@  | @@| @@@@@@@\| @@      | @@    @@  \@@@@@@\| @@  | @@ "
write(*,*) " _\@@@@@@\| @@  | @@| @@  | @@| @@      | @@@@@@@\ /      @@| @@  | @@ "
write(*,*) "|  \__| @@| @@__/ @@| @@  | @@| @@      | @@__/ @@|  @@@@@@@| @@__/ @@ "
write(*,*) " \@@    @@ \@@    @@| @@  | @@| @@      | @@    @@ \@@    @@ \@@    @@ "
write(*,*) "  \@@@@@@   \@@@@@@  \@@   \@@ \@@       \@@@@@@@   \@@@@@@@ _\@@@@@@@ "
write(*,*) "                                                            |  \__| @@ "
write(*,*) "                                                             \@@    @@ "
write(*,*) "                                                              \@@@@@@  "
write(*,*) "           |                       "
write(*,*) "         \ _ /                     "
write(*,*) "       -= (_) =-                   "
write(*,*) "         /   \         _\/_        "
write(*,*) "           |           //o\  _\/_  "
write(*,*) "    _____ _ __ __ ____ _ | __/o\\ _"
write(*,*) "  =-=-_-__=_-= _=_=-=_,-=|===,-|-,_"
write(*,*) "   =- _=-=- -_=-=_,-==         |   "
write(*,*) "     =- =- -=.--==                 "
write(*,*) "                                   "
write(*,*) "  Statistical estimation of UNIvariate "
write(*,*) " distributions using BAYesian inference"
write(*,*) "                                      "
write(*,*) "---------------------------------------"
write(*,*) "                                      "

end subroutine Welcome

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GoodBye

write(*,*) "---------------------------------------"
write(*,*) "                                      "
write(*,*) " Thanx For Using Me!"
write(*,*) " Press [Enter] and I'll go away"
write(*,*) "                                      "
write(*,*) "      _~     "
write(*,*) "   _ )_)_~   "
write(*,*) "  )_))_))_)  "
write(*,*) "  _!__!__!_  "
write(*,*) "  \_______/  "
write(*,*) "~~~~~~~~~~~~~"
write(*,*) "                                      "
read(*,*)

end subroutine GoodBye

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GaussianKDE(x,sample,pdf,err,mess)

!^**********************************************************************
!^* Purpose: Quick-and-dirty kde implementation - see with D to put this
!^* in DMSL. In particular bw is imposed here
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified: 27/12/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.x;value at which pdf is evaluated
!^*		2.sample, obs used for performing estimation
!^* OUT
!^*		1.pdf
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use numerix_dmsl_kit, only:getvar,normal_logp

real(mrk), intent(in)::x,sample(:)
real(mrk),intent(out)::pdf
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
integer(mik)::N,i
real(mrk)::sigma,h,toto(size(sample))

err=0;mess="";
N=size(sample)
sigma=sqrt(getvar(sample,"c"))
h=(1.06_mrk*sigma)/((real(N,mrk))**(0.2))
do i=1,N
    toto(i)=exp(normal_logp( ( x-sample(i) )/h,0._mrk,1._mrk))
enddo
pdf=(real(1,mrk)/(real(N*h,mrk)))*sum(toto)

end subroutine GaussianKDE


end module SuniBay_tools
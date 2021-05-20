module EmpiricalStats_tools

!~**********************************************************************
!~* Purpose: Compute empirical stats from a data vector x
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified:10/08/2010
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: enhance extrapolation for pval & quantile evaluators
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. GettPP: plotting position catalogue
!~*		2. GetEmpiricalPval
!~*		3. GetEmpiricalQuantile
!~*		4. GetEmpiricalStats, 15 statistics frequently used in statistical summaries
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: GetPP,GetEmpiricalPval,GetEmpiricalQuantile,GetEmpiricalStats,KDE,&
            GetEmpiricalDistSummary,P2T,T2p,Get4LMoments,GetRank,GetPearsonCorr

! aliases for Kernel names
character(250), parameter, public:: unif_K="Uniform_Kernel",uniform_K="Uniform_Kernel",&
                         triangle_K="Triangular_Kernel",triangular_K="Triangular_Kernel",&
                         epan_K="Epanechnikov_Kernel",epanechnikov_K="Epanechnikov_Kernel",&
                         quartic_K="Biweight_Kernel",biweight_K="Biweight_Kernel",&
                         triweight_K="Triweight_Kernel",&
                         tricube_K="Tricube_Kernel",&
                         gaussian_K="Gaussian_Kernel",gauss_K="Gaussian_Kernel",&
                         cosine_K="Cosine_Kernel",cos_K="Cosine_Kernel"

Contains

pure subroutine GetPP(i,n,formula,pp,err,mess)

!^**********************************************************************
!^* Purpose: plotting position evaluator
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 10/08/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. i, rank
!^*		2. n, total 
!^*		3. [formula], default is Hazen - see below for list 
!^* OUT
!^*		1. pp
!^*		2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3. mess, error message
!^**********************************************************************

! "Standard"     (i)/(n)
! "MinusOne"     (i-1)/(n)
! "Hazen"        (i-0.5)/(n)
! "Weibull"      (i)/(n+1) 
! "Benard"       (i-0.3)/(n+0.4)
! "Cunnane"      (i-0.4)/(n+0.2)
! "Beard"        (i-0.31)/(n+0.38)
! "Blom"         (i-0.375)/(n+0.25)
! "Gringorten"   (i-0.44)/(n+0.12)
! "Landwehr"     (i-0.35)/(n)
! "Tukey"        (3i-1)/(3n+1)

integer(mik), intent(in)::i,n
character(*), intent(in), optional::formula
real(mrk), intent(out)::pp
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250)::form

err=0;mess='';pp=undefRN

if ( (i<=0) .or. (i>n) ) then
    err=1; mess="GetPP:FATAL:i should be in [1:n]";return
endif
if(present(formula)) then
    form=trim(formula)
else
    form='Hazen'
endif

select case (trim(form))
case("Standard")
    pp=real(i,mrk)/real(n,mrk)
case("MinusOne")
    pp=real(i-1,mrk)/real(n,mrk)
case("Hazen")
    pp=(real(i,mrk)-0.5_mrk)/real(n,mrk)
case("Weibull")
    pp=real(i,mrk)/real(n+1,mrk)
case("Benard")
    pp=( real(i,mrk)-0.3_mrk )/( real(n,mrk)+0.4_mrk )
case("Cunnane")
    pp=( real(i,mrk)-0.4_mrk )/( real(n,mrk)+0.2_mrk )
case("Beard")
    pp=( real(i,mrk)-0.31_mrk )/( real(n,mrk)+0.38_mrk )
case("Blom")
    pp=( real(i,mrk)-0.375_mrk )/( real(n,mrk)+0.25_mrk )
case("Gringorten")
    pp=( real(i,mrk)-0.44_mrk )/( real(n,mrk)+0.12_mrk )
case("Landwehr")
    pp=( real(i,mrk)-0.35_mrk )/real(n,mrk)
case("Tukey")
    pp=real(3*i-1,mrk)/real(3*n+1,mrk)
case default
    err=1; mess="GetPP:FATAL:unknown formula";return
end select

end subroutine GetPP

subroutine GetRank(xIN,mv,xOut,err,mess)

!^**********************************************************************
!^* Purpose: Get rank
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 22/10/2012
!^**********************************************************************
!^* Comments: No randomization for ties 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Implement options for ties
!^**********************************************************************
!^* IN
!^*		1. xIN, original data
!^*		2. mv, value for missing data
!^* OUT
!^*		1. xOUT, ranks
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use numerix_dmsl_kit, only:indexx_qsort
real(mrk),intent(in)::xIN(:),mv
real(mrk), intent(out)::xOUT(size(xIN))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk), allocatable::dummy(:)
integer(mik), allocatable::indx(:),location(:)
logical::mask(size(xIN))
integer(mik)::i,n

err=0;mess='';xout=mv
mask=xIN/=mv
n=count(mask)
if(allocated(dummy)) deallocate(dummy);allocate(dummy(n))
if(allocated(indx)) deallocate(indx);allocate(indx(n))
if(allocated(location)) deallocate(location);allocate(location(n))
dummy=pack(xIN,mask)
location=pack((/(i,i=1,size(xIN))/),mask)
call indexx_qsort(arr=dummy,indx=indx,ascnd=.true.,err=err,message=mess)
!quicksort(arr=dummy,ascnd=.true.,err=err)
if(err>0) then
    mess="GetRank:FATAL: "//trim(mess);return
endif
do i=1,n
    xOUT(location(indx(i)))=real(i,mrk)
enddo
end subroutine GetRank

subroutine GetEmpiricalPval(obs,x,ppFormula,IsXSorted,pv,err,mess)

!^**********************************************************************
!^* Purpose: Compute the empirical pval of obs into the sample x
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 10/08/2010
!^**********************************************************************
!^* Comments: Dummy handling of tails for the time being
!^* if obs<min(x), pv=0
!^* if obs>max(x), pv=1
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Try smarter extrapolation schemes for tails
!^**********************************************************************
!^* IN
!^*		1. obs, value whose pval is sought
!^*		2. x, sample used to compute the empirical cdf
!^*		3. [ppFormula] (default is Hazen)
!^*		4. [IsXSorted] (default is false) - might save time to sort out of this sub
!^* OUT
!^*		1.pv
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:linearInterp
use numerix_dmsl_kit, only:quicksort

real(mrk), intent(in)::obs,x(:)
character(*), intent(in), optional::ppFormula
logical, intent(in),optional::IsXSorted
real(mrk), intent(out)::pv
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,n
real(mrk),allocatable::F(:),arr(:)
logical::IsXS

err=0;mess='';pv=undefRN
if(present(IsXSorted)) then
    IsXS=IsXSorted
else
    IsXS=.false.
endif
if(obs<minval(x)) then
    err=-1;mess="GetEmpiricalPval:WARNING:extrapolation not enabled"
    pv=0._mrk
elseif (obs>maxval(x)) then
    err=-1;mess="GetEmpiricalPval:WARNING:extrapolation not enabled"
    pv=1._mrk
else
    n=size(x)
    if(allocated(arr)) deallocate(arr)
    if(allocated(F)) deallocate(F)
    allocate(arr(n),F(n))
    arr=x
    if(.not.IsXS) call quicksort(arr=arr,ascnd=.true.,err=err)
    if(err>0) then
        mess="GetEmpiricalPval:FATAL: "//trim(mess);return
    endif
    do i=1,n
        call GetPP(i=i,n=n,formula=ppFormula,pp=F(i),err=err,mess=mess)
        if(err>0) then
            mess="GetEmpiricalPval:FATAL: "//trim(mess);return
        endif
    enddo
    call linearInterp(x=obs,xx=arr,yy=F,y=pv,err=err,message=mess)
    if(err>0) then
        mess="GetEmpiricalPval:FATAL: "//trim(mess);return
    endif
endif

end subroutine GetEmpiricalPval

subroutine GetEmpiricalQuantile(p,x,ppFormula,IsXSorted,q,err,mess)

!^**********************************************************************
!^* Purpose: Compute the empirical p-quantile from the dist.of sample x
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 10/08/2010
!^**********************************************************************
!^* Comments: Dummy handling of tails for the time being
!^* if p<pp(1,n), q=-HugeRe
!^* if p>pp(n,n), q=+HugeRe
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Try smarter extrapolation schemes for tails
!^**********************************************************************
!^* IN
!^*		1. p, p-value
!^*		2. x, sample used to compute the empirical cdf
!^*		3. [ppFormula] (default is Hazen)
!^*		4. [IsXSorted] (default is false) - might save time to sort out of this sub
!^* OUT
!^*		1. q, the empirical p-quantile
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:linearInterp
use numerix_dmsl_kit, only:quicksort

real(mrk), intent(in)::p,x(:)
character(*), intent(in), optional::ppFormula
logical, intent(in),optional::IsXSorted
real(mrk), intent(out)::q
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,n
real(mrk)::mini,maxi
real(mrk),allocatable::F(:),arr(:)
logical::IsXS

err=0;mess='';q=undefRN
if(present(IsXSorted)) then
    IsXS=IsXSorted
else
    IsXS=.false.
endif
n=size(x)
call GetPP(i=1,n=n,formula=ppFormula,pp=mini,err=err,mess=mess)
if(err>0) then
    mess="GetEmpiricalQuantile:FATAL: "//trim(mess);return
endif
call GetPP(i=n,n=n,formula=ppFormula,pp=maxi,err=err,mess=mess)
if(err>0) then
    mess="GetEmpiricalQuantile:FATAL: "//trim(mess);return
endif


if(p<mini) then
    err=-1;mess="GetEmpiricalQuantile:WARNING:extrapolation not enabled"
    q=-HugeRe
elseif (p>maxi) then
    err=-1;mess="GetEmpiricalQuantile:WARNING:extrapolation not enabled"
    q=HugeRe
else
    if(allocated(arr)) deallocate(arr)
    if(allocated(F)) deallocate(F)
    allocate(arr(n),F(n))
    arr=x
    if(.not.IsXS) call quicksort(arr=arr,ascnd=.true.,err=err)
    if(err>0) then
        mess="GetEmpiricalQuantile:FATAL: "//trim(mess);return
    endif
    do i=1,n
        call GetPP(i=i,n=n,formula=ppFormula,pp=F(i),err=err,mess=mess)
        if(err>0) then
            mess="GetEmpiricalQuantile:FATAL: "//trim(mess);return
        endif
    enddo
    call linearInterp(x=p,xx=F,yy=arr,y=q,err=err,message=mess)
    if(err>0) then
        mess="GetEmpiricalQuantile:FATAL: "//trim(mess);return
    endif
endif

end subroutine GetEmpiricalQuantile

subroutine GetEmpiricalStats(x, & !data,
                ppFormula, & ! Formula for plotting position - used for quantile stats
                n,mini,maxi,range, &! basics
                mean,median,LeftDecile,LeftQuartile,RightQuartile, RightDecile,& ! location stats
                std,var,CV,& ! dispersion stats
                skewness,kurtosis, & ! shape stats
                all15,& ! pack all stats into a single vector
                err,mess)

!^**********************************************************************
!^* Purpose: Get a standard statistical summary for data x
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:10/08/2010
!^**********************************************************************
!^* Comments: all output variables are optional (except err & mess)
!^* output "all15" can be used to pack all stats into a single vector of length 15
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.x, data vector
!^*		2.[ppFormula], Formula for plotting position (default is Hazen)
!^* OUT
!^*		1. [n], sample size
!^*		2. [mini], sample minimum
!^*		3. [maxi], sample maximum
!^*		4. [range], max-min
!^*		5. [mean]
!^*		6. [median]
!^*		7. [LeftDecile]
!^*		8. [LeftQuartile]
!^*		9. [RightDecile]
!^*		10. [RightQuartile]
!^*		11. [std]
!^*		12. [var]
!^*		13. [CV]
!^*		14. [skewness]
!^*		15. [kurtosis] EXCESS kurtosis
!^*		16. [all15]
!^*		17.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		18.mess, error message
!^**********************************************************************
use numerix_dmsl_kit, only:getmean, getvar, getCV,getmoments

real(mrk), intent(in)::x(:)
character(*), intent(in), optional::ppFormula
integer(mik), intent(out), optional::n
real(mrk), intent(out), optional::mini,maxi,range,mean,median,LeftDecile,&
LeftQuartile,RightQuartile, RightDecile,std,var,CV,skewness,kurtosis, all15(15)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::Xn
real(mrk)::Xmini,Xmaxi,Xrange,Xmean,Xmedian,XLeftDecile,&
XLeftQuartile,XRightQuartile, XRightDecile,Xstd,Xvar,XCV,Xskewness,Xkurtosis
real(mrk)::crap

err=0;mess=''

if( present(n) .or. present(all15) ) then
    Xn=size(x)
    if(present(n)) n=Xn
    if(present(all15)) all15(1)=Xn
endif

if( present(mini) .or. present(all15) ) then
    Xmini=minval(x)
    if(present(mini)) mini=Xmini
    if(present(all15)) all15(2)=Xmini
endif

if( present(maxi) .or. present(all15) ) then
    Xmaxi=maxval(x)
    if(present(maxi)) maxi=Xmaxi
    if(present(all15)) all15(3)=Xmaxi
endif

if( present(range) .or. present(all15) ) then
    Xrange=maxval(x)-minval(x)
    if(present(range)) range=Xrange
    if(present(all15)) all15(4)=Xrange
endif

if( present(mean) .or. present(all15) ) then
    Xmean=getmean(x)
    if(present(mean)) mean=Xmean
    if(present(all15)) all15(5)=Xmean
endif

if( present(median) .or. present(all15) ) then
    call GetEmpiricalQuantile(p=0.5_mrk,x=x,ppFormula=ppFormula,q=Xmedian,err=err,mess=mess)
    if(err>0) then
        mess="EmpiricalStats_tools: "//trim(mess);return
    endif    
    if(present(median)) median=Xmedian
    if(present(all15)) all15(6)=Xmedian
endif

if( present(LeftDecile) .or. present(all15) ) then
    call GetEmpiricalQuantile(p=0.1_mrk,x=x,ppFormula=ppFormula,q=XLeftDecile,err=err,mess=mess)
    if(err>0) then
        mess="EmpiricalStats_tools: "//trim(mess);return
    endif    
    if(present(LeftDecile)) LeftDecile=XLeftDecile
    if(present(all15)) all15(7)=XLeftDecile
endif

if( present(LeftQuartile) .or. present(all15) ) then
    call GetEmpiricalQuantile(p=0.25_mrk,x=x,ppFormula=ppFormula,q=XLeftQuartile,err=err,mess=mess)
    if(err>0) then
        mess="EmpiricalStats_tools: "//trim(mess);return
    endif    
    if(present(LeftQuartile)) LeftQuartile=XLeftQuartile
    if(present(all15)) all15(8)=XLeftQuartile
endif

if( present(RightQuartile) .or. present(all15) ) then
    call GetEmpiricalQuantile(p=0.75_mrk,x=x,ppFormula=ppFormula,q=XRightQuartile,err=err,mess=mess)
    if(err>0) then
        mess="EmpiricalStats_tools: "//trim(mess);return
    endif    
    if(present(RightQuartile)) RightQuartile=XRightQuartile
    if(present(all15)) all15(9)=XRightQuartile
endif

if( present(RightDecile) .or. present(all15) ) then
    call GetEmpiricalQuantile(p=0.9_mrk,x=x,ppFormula=ppFormula,q=XRightDecile,err=err,mess=mess)
    if(err>0) then
        mess="EmpiricalStats_tools: "//trim(mess);return
    endif    
    if(present(RightDecile)) RightDecile=XRightDecile
    if(present(all15)) all15(10)=XRightDecile
endif

if( present(std) .or. present(all15) ) then
    Xstd=sqrt(getvar(x,method="c"))
    if(present(std)) std=Xstd
    if(present(all15)) all15(11)=Xstd
endif

if( present(var) .or. present(all15) ) then
    Xvar=getvar(x,method="c")
    if(present(var)) var=Xvar
    if(present(all15)) all15(12)=Xvar
endif

if( present(CV) .or. present(all15) ) then
    XCV=getCV(x,method="c")
    if(present(CV)) CV=XCV
    if(present(all15)) all15(13)=XCV
endif

if( present(skewness) .or. present(all15) ) then
    call getmoments(x=x,mean=Xmean, adev=crap,sdev=Xstd,var=Xvar,skew=Xskewness,kurt=Xkurtosis,err=err,message=mess)
    !call getmoments(x=x,ave=Xmean, adev=crap,sdev=Xstd,var=Xvar,skew=Xskewness,kurt=Xkurtosis,err=err,message=mess)
    if(err>0) then
        mess="EmpiricalStats_tools: Warning: "//trim(mess) !;return Do Not return: put garbage in slew&kurt instead
        if(present(skewness)) skewness=undefRN
        if(present(all15)) all15(14)=undefRN
        err=-10 ! turn into a warning message
        write(*,*) trim(mess)
    else    
        if(present(skewness)) skewness=Xskewness
        if(present(all15)) all15(14)=Xskewness
    endif
endif

if( present(kurtosis) .or. present(all15) ) then
    call getmoments(x=x,mean=Xmean, adev=crap,sdev=Xstd,var=Xvar,skew=Xskewness,kurt=Xkurtosis,err=err,message=mess)
    !call getmoments(x=x,ave=Xmean, adev=crap,sdev=Xstd,var=Xvar,skew=Xskewness,kurt=Xkurtosis,err=err,message=mess)
    if(err>0) then
        mess="EmpiricalStats_tools: Warning: "//trim(mess) !;return
        if(present(kurtosis)) kurtosis=undefRN
        if(present(all15)) all15(15)=undefRN
        err=-10 ! turn into a warning message
        write(*,*) trim(mess)
    else    
        if(present(kurtosis)) kurtosis=Xkurtosis
        if(present(all15)) all15(15)=Xkurtosis
    endif
endif

end subroutine GetEmpiricalStats

subroutine KDE(kernelID,u,X,h,out,err,mess)

!^**********************************************************************
!^* Purpose: kernel density estimator
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:26/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.[kernelID], default Gaussian
!^*		2.u, value where estimate is sought
!^*		3.X, data
!^*		4.[h], smoothing parameter. Default value = Silverman's rule of
!^*		       thumb, h = 1.06*std(X)*n^(-1/5)
!^* OUT
!^*		1.out, estimated pdf(u)
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*),intent(in),optional::kernelID
real(mrk), intent(in)::u, X(:)
real(mrk), intent(in),optional::h
real(mrk), intent(out)::out
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
real(mrk)::smooth,std,k(size(X))
integer(mik)::n,i
character(250)::kid

err=0;mess=""
if(present(kernelID)) then
    kid=kernelID
else
    kid="Gaussian_Kernel"
endif
n=size(X)
! handle optional h
if(present(h)) then
  smooth=h
else
  call GetEmpiricalStats(x=X,std=std,err=err,mess=mess)
  If(err>0) then
    mess='KDE'//trim(mess);return
  endif
  smooth=(1.06_mrk*std)/(real(n,mrk)**0.2_mrk)
endif

if(smooth<=0._mrk) then
  err=1;mess="KDE:FATAL:smoothing parameter should be strictly positive";return
endif

do i=1,n
  call kernelfunk(kernelID=kid,u=(u-X(i))/smooth,k=k(i),err=err,mess=mess)
  If(err>0) then
    mess='KDE'//trim(mess);return
  endif
enddo
out=(1._mrk/(real(n,mrk)*smooth))*sum(k)

end subroutine KDE

subroutine GetEmpiricalDistSummary(X,& ! data
                        ppFormula, &! optional,, formula for plotting position [default Hazen]
                        kernelID,h,& !optional, parameters of kde [defaults Gaussian + built_in h value]
                        p2Tfactor,&! optional, factor for converting probabilities to return periods, default 1
                        invertT,&! optional, Revert p-T relation? (yes if large T for small values)
                        Xname,& ! optional, name given to X, will be written in header of output file
                        SortedX,Epdf,Ecdf,ET,&! outputs, all optionals
                        Outfile,&! optional, file to write results
                        err,mess)

!^**********************************************************************
!^* Purpose: kernels definition
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:26/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.kernel
!^*		2.u
!^* OUT
!^*		1.k (=kernel(u))
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:GetSpareUnit
use numerix_dmsl_kit, only:quicksort
real(mrk),intent(in)::X(:)
character(*), optional, intent(in)::OutFile,kernelID,ppFormula,Xname
real(mrk), optional, intent(in)::h,p2Tfactor
Logical, optional, intent(in)::invertT
real(mrk), optional, intent(out)::sortedX(:),Epdf(:),Ecdf(:),ET(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
real(mrk),dimension(size(X))::sX,pdf,cdf,T
integer(mik)::i,n,unt
character(250)::varname
err=0;mess='';
n=size(X)
sX=X
call quicksort(sX,.true.,err)
if(err/=0) then
    mess='GetEmpiricalDistSummary:SortingError';return
endif
do i=1,n
    ! Get Ecdf
    call GetPP(i,n,ppFormula,cdf(i),err,mess)
    if(err/=0) then
        mess='GetEmpiricalDistSummary:'//trim(mess);return
    endif
    ! Transform to T
    T(i)=p2T(cdf(i),p2Tfactor,invertT)
    ! Get KDE estimate of pdf
    call KDE(kernelID,sx(i),sX,h,pdf(i),err,mess)
    if(err/=0) then
        mess='GetEmpiricalDistSummary:'//trim(mess);return
    endif
enddo    
if(present(sortedX)) sortedX=sX
if(present(Epdf)) Epdf=pdf
if(present(Ecdf)) Ecdf=cdf
if(present(ET)) ET=T
if(present(OutFile)) then
  call GetSpareUnit(unt,err,mess)
  If(err>0) then
    mess='GetEmpiricalDistSummary'//trim(mess);return
  endif
  open(unit=unt,file=trim(OutFile),status="REPLACE",iostat=err)
  If(err/=0) then
    mess='GetEmpiricalDistSummary'//trim(mess);return
  endif
  ! Write parameters stats
  if(present(XName)) then
    varname=Xname
  else
    varname="Sorted data"
  endif
  write(unt,'(4A15)') trim(varname),"Epdf","Ecdf","ET"
  do i=1,n
    write(unt,'(4e15.6)') sx(i),pdf(i),cdf(i),T(i)
  enddo
  close(unt)
endif
end subroutine GetEmpiricalDistSummary

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Get4LMoments(x,LMOM,err,mess)

!^**********************************************************************
!^* Purpose: Compute the first 4 sample L-moments according to Hosking,
!^* ie l1,l2 (loc/scale) and ratios r3 (L-skew) and r4(L-kurt)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 07/04/2012
!^**********************************************************************
!^* Comments: Those 4 are sufficient to estimate the parameters of most 
!^* distributions, at least those useful in hydro-applications 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Generalize it to individual or first K L-moments
!^* + properly distinguish L-moments and L-moment ratios
!^**********************************************************************
!^* IN
!^*		1. x, sample used to compute L-MOM
!^* OUT
!^*		1.LMOM, real vector of size 4
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use numerix_dmsl_kit, only:quicksort
real(mrk), intent(in)::x(:)
real(mrk), intent(out)::LMOM(4)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::i,n
real(mrk)::b0,b1,b2,b3
real(mrk), dimension(size(X))::sX

err=0;mess='';LMOM=undefRN
n=size(X)
if(n<=3) then
    err=1;mess='Get4LMoments:Fatal:n<=3';return
endif
sX=X
call quicksort(sX,.true.,err)
if(err/=0) then
    mess='Get4LMoments:SortingError';return
endif
! L1
b0=sum(sX)/real(n,mrk)
LMOM(1)=b0
! L2
b1=0._mrk
do i=1,n
    b1=b1+(real(i-1,mrk)/real(n-1,mrk))*sX(i)
enddo
b1=b1/real(n,mrk)
LMOM(2)=b0*pstar(1,0)+b1*pstar(1,1)
! L3
b2=0._mrk
do i=1,n
    b2=b2+(real((i-1)*(i-2),mrk)/real((n-1)*(n-2),mrk))*sX(i)
enddo
b2=b2/real(n,mrk)
LMOM(3)=(b0*pstar(2,0)+b1*pstar(2,1)+b2*pstar(2,2))/LMOM(2)
! L4
b3=0._mrk
do i=1,n
    b3=b3+(real((i-1)*(i-2)*(i-3),mrk)/real((n-1)*(n-2)*(n-3),mrk))*sX(i)
enddo
b3=b3/real(n,mrk)
LMOM(4)=(b0*pstar(3,0)+b1*pstar(3,1)+b2*pstar(3,2)+b3*pstar(3,3))/LMOM(2)

end subroutine Get4LMoments
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function p2T(p,p2TFactor,invertT)
! probability to return period transformation
! note: p2Tfactor is optional,and corresponds to the average number of 
!       data per year (for a return period in years) - default is 1
! p=1-1/(T*p2Tfactor),T=1/(1-p)*p2Tfactor
! if invertT is true, then p=1/(T*p2Tfactor),T=1/(p*p2Tfactor)

real(mrk),intent(in)::p
real(mrk),intent(in),optional::p2Tfactor
logical,intent(in),optional::invertT
real(mrk)::p2T
real(mrk)::f
logical::iT

if(p>=1._mrk) then
    p2T=undefRN;return
elseif(p<0._mrk) then
    p2T=-undefRN;return
endif 
if(present(p2TFactor)) then
    f=p2TFactor
else
    f=1._mrk
endif
if(present(invertT)) then
    iT=invertT
else
    iT=.false.
endif
if(iT) then
    p2T=1._mrk/(p*f)
else
    p2T=1._mrk/((1._mrk-p)*f)
endif
end function p2T

elemental function T2p(T,p2TFactor,invertT)
! return period to probability transformation
! note: p2Tfactor is optional,and corresponds to the average number of 
!       data per year (for a return period in years) - default is 1
! p=1-1/(T*p2Tfactor),T=1/(1-p)*p2Tfactor
! if invertT is true, then p=1/(T*p2Tfactor),T=1/(p*p2Tfactor)
real(mrk),intent(in)::T
real(mrk),intent(in),optional::p2Tfactor
logical,intent(in),optional::invertT
real(mrk)::T2p
real(mrk)::f
logical:: iT

if(present(p2TFactor)) then
    f=p2TFactor
else
    f=1._mrk
endif
if(present(invertT)) then
    iT=invertT
else
    iT=.false.
endif
if(T*f<1._mrk) then
    T2p=-undefRN;return
endif 
if(iT) then
    T2p=1._mrk/(T*f)
else
    T2p=1._mrk-1._mrk/(T*f)
endif
end function T2p

Subroutine GetPearsonCorr(x,y,ro,err,mess)
!^**********************************************************************
!^* Purpose: Pearson correlation
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 25/11/2013
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. x,y, samples used to compute correlation
!^* OUT
!^*		1.ro, correlation
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
real(mrk),intent(in)::x(:),y(:)
real(mrk),intent(out)::ro
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
real(mrk)::mx,sx,my,sy

err=0;mess='';ro=undefRN
if(size(x)/=size(y)) then
    err=1;mess='GetPearsonCorr:FATAL:size mismatch[x,y]';return;
endif
call GetEmpiricalStats(x=x,mean=mx,std=sx,err=err,mess=mess)
if(err/=0) then;mess='GetPearsonCorr:'//trim(mess);return;endif
if(sx==0._mrk) then
    err=1;mess='GetPearsonCorr:FATAL:std[x]==0';return;
endif
call GetEmpiricalStats(x=y,mean=my,std=sy,err=err,mess=mess)
if(err/=0) then;mess='GetPearsonCorr:'//trim(mess);return;endif
if(sy==0._mrk) then
    err=1;mess='GetPearsonCorr:FATAL:std[y]==0';return;
endif
ro=dot_product(x-mx,y-my)/(real(size(x),mrk)*sx*sy)
               
end Subroutine GetPearsonCorr

!!!!!!!!!!!
! PRIVATE !
!!!!!!!!!!!

pure subroutine kernelfunk(kernelID,u,k,err,mess)

!^**********************************************************************
!^* Purpose: kernels definition
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:26/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.kernel
!^*		2.u
!^* OUT
!^*		1.k (=kernel(u))
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:pi

character(*), intent(in)::kernelID
real(mrk), intent(in)::u
real(mrk), intent(out)::k
integer(mik), intent(out)::err
character(*),intent(out)::mess

err=0;mess=""
select case(trim(kernelID))
case(unif_K)
  if(abs(u)<=1._mrk) then
    k=0.5_mrk
  else
    k=0._mrk
  endif
case(triangle_K)
  if(abs(u)<=1._mrk) then
    k=1._mrk-abs(u)
  else
    k=0._mrk
  endif
case(epan_K)
  if(abs(u)<=1._mrk) then
    k=0.75_mrk*(1._mrk-u**2)
  else
    k=0._mrk
  endif
case(biweight_K)
  if(abs(u)<=1._mrk) then
    k=(15._mrk/16._mrk)*((1._mrk-u**2)**2)
  else
    k=0._mrk
  endif
case(triweight_K)
  if(abs(u)<=1._mrk) then
    k=(35._mrk/32._mrk)*((1._mrk-u**2)**3)
  else
    k=0._mrk
  endif
case(tricube_K)
  if(abs(u)<=1._mrk) then
    k=(70._mrk/81._mrk)*((1._mrk-abs(u)**3)**3)
  else
    k=0._mrk
  endif
case(gaussian_K)
  k=(1._mrk/sqrt(2._mrk*pi))*exp(-0.5_mrk*u**2)
case(cosine_K)
  if(abs(u)<=1._mrk) then
    k=0.25_mrk*pi*cos(0.5_mrk*pi*u)
  else
    k=0._mrk
  endif
case default
  err=1; mess="kernelfunk:FATAL:unknown kernelID";return
end select
end subroutine kernelfunk

function pstar(r,k)
use utilities_dmsl_kit, only:bico
! Used for L-moments computation
integer(mik), intent(in)::r,k
real(mrk)::pstar
if(r==k) then
    pstar=bico(r,k)*bico(r+k,k)
else
    pstar=(-1._mrk)**(r-k)*bico(r,k)*bico(r+k,k)
endif
end function pstar
end module EmpiricalStats_tools
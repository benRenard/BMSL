module MultivariateDistribution_tools
!~**********************************************************************
!~* Purpose: Utilities for multivariate distributions
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon, Univerity of Adelaide
!~**********************************************************************
!~* Created: 09/07/2019
!~**********************************************************************
!~* Comments:
!~*     1. Parameterization is more tricky than for one-dimensional
!~*        distributions: it is not practical to stack all parameters
!~*        in a vector (e.g. think about the symetric covariance matrix
!~*        of a multi-Gaussian).
!~*     2. => use of a derived type mDistPar (see below) for parameters
!~*        of a p-dimensional distribution.
!~*     3. A mDistPar parameter for a p-dimensional distribution contains:
!~*        (i)   a vector of scalar parameters 'spar' (e.g. mean / sdev for an iid Gaussian)
!~*        (ii)  a p-vector vpar (typically, mean vector)
!~*        (iii) a p*p-matrix mpar (typically, covariance matrix)
!~*        (iv)  various other fields for additional parameters and/or computational tricks
!~*     4. In particular, for distributions involving covariance inversion,
!~*        isCov tells if mpar is a covariance (=> inversion needed),
!~*        saveInvert tells if inverted matrix needs to be saved (in which
!~*        case it overwrites covariance in mpar and isCov is turned to .false.)
!~*     5. So far mostly pdf evaluation, need to be completed later
!~*        (simulation, cdf, various types of quantiles, etc.)
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~*     1. Complete available utilities
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. initPar: initialise parameter object
!~*		2. getParFeas: assess parameter feasability
!~*		3. GetLogPdf: compute joint log-pdf
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use types_dmsl_kit,only:vectorA1_r_type,vectorA2_r_type,vectorA1_i_type ! non-rectangular arrays

implicit none
Private
Public::initPar,getParFeas,GetLogPdf,getLogPdf_singleCompUpdate,NNGP_getNN

! derived type used for storing the parameters of a multivariate distribution
type, public::mDistPar
    real(mrk), allocatable::spar(:) ! Scalar parameters; e.g. mean and sdev for 'Gaussian_IID'
    real(mrk), allocatable::vpar(:) ! Vector parameter; size p; e.g. mean vector for 'Gaussian'
    real(mrk), allocatable::mpar(:,:) ! Matrix parameter; size p*p; e.g. covariance (or precision) for 'Gaussian'
    type(vectorA1_i_type),allocatable::i1(:) ! list of integer vectors; e.g. i1(k)%v(:) may contain the list of neighboring parents for site k
    type(vectorA1_r_type),allocatable::r1(:) ! list of real vectors; e.g.  additional scalar of vector parameters
    type(vectorA2_r_type),allocatable::r2(:) ! list of real matrices; e.g.  r2(k)%v(:,:) contains the inverse covariance of neighboring parents for site k
    logical:: isCov=.true. ! does matrix parameter represent covariance or rather precision?
    logical:: saveInvert=.true. ! overwrite covariance with precision when inversion performed as part of pdf computation?
    real(mrk):: logDet=undefRN ! log-determinant
    real(mrk), allocatable::logps(:) ! components of the joint pdf for distribution that can be written as joint_logp = sum(logps) (e.g. iid, AR1 or NNGP)
end type mDistPar

Character(100), parameter:: & ! Shortcut for distribution names
            ! Uniform distributions
            FLAT='Flat',& ! Improper flat distribution (log-pdf=0 everywhere)
            UNIF='Uniform',& ! independant uniform distributions for each margin
            UNIF_IID='Uniform_IID',& ! iid uniform distribution
            ! Gaussian distributions
            GAUSS='Gaussian', & ! Multivariate Gaussian Distribution
            GAUSS_IID='Gaussian_IID', & ! iid Gaussian Distribution
            GAUSS_AR1='Gaussian_AR1', & ! AR(1) Gaussian process
            GAUSS_AR1_vmean='Gaussian_AR1_vmean', & ! AR(1) Gaussian process with a variable mean
            NNGP='NNGP' ! Nearest Neighbors Gaussian process

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine initPar(DistId,p,par,err,mess)
!^**********************************************************************
!^* Purpose: initialize parameter for a p-dimensional DistId
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 05/09/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, distribution ID
!^*    2.p, dimension
!^* OUT
!^*    1.par, initialized parameter
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************
character(*), intent(in)::DistID
integer(mik),intent(in)::p
type(mDistPar),intent(out)::par
integer(mik), intent(out)::err
character(*),intent(out)::mess
!Locals
character(250),parameter::procname='initPar'
integer(mik)::i

err=0;mess=''

if(p<=0) then
    mess=trim(procname)//':Fatal:p<=0'
    err=1;return
endif

select case(DistID)
case(FLAT)! no parameter
    ! nothing to do
case(UNIF) ! bounds in r1
    allocate(par%r1(p))
    do i=1,p
        allocate(par%r1(i)%v(2)) ! (lower bound, higher bound) for ith margin
    enddo
case(UNIF_IID)! 2 scalar pars: lower bound in spar(1), upper bound in spar(2)
    allocate(par%spar(2));par%spar=undefRN
    allocate(par%logps(p));par%logps=undefRN
case(GAUSS_IID)! 2 scalar pars: mean in spar(1), stdev in spar(2)
    allocate(par%spar(2));par%spar=undefRN
    allocate(par%logps(p));par%logps=undefRN
case(GAUSS_AR1)! 3 scalar pars: constant in spar(1), innovation stdev in spar(2), rho in spar(3)
    allocate(par%spar(3));par%spar=undefRN
    allocate(par%logps(p));par%logps=undefRN
case(GAUSS_AR1_vmean)! mean in vpar, then 2 scalar pars: innovation stdev in spar(1), rho in spar(2)
    allocate(par%vpar(p));par%vpar=undefRN
    allocate(par%spar(2));par%spar=undefRN
    allocate(par%logps(p));par%logps=undefRN
case(GAUSS) ! mean in vpar, covariance (or precision) in mpar
    allocate(par%vpar(p));par%vpar=undefRN
    allocate(par%mpar(p,p));par%mpar=undefRN
case(NNGP) ! mean in vpar, covariance in mpar, nearest parents in i1, mini-covs in r2
    allocate(par%spar(1));par%spar=undefRN ! not used in computations, but can store max number of neighbors if useful
    allocate(par%vpar(p));par%vpar=undefRN
    allocate(par%mpar(p,p));par%mpar=undefRN
    allocate(par%i1(p))
    allocate(par%r2(p))
    allocate(par%logps(p));par%logps=undefRN
case default
    mess=trim(procname)//':Fatal:Unavailable Dist'
    err=1;return
end select

end subroutine initPar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getParFeas(DistId,par,feas,err,mess)
!^**********************************************************************
!^* Purpose: parameter feasability & check for allocation and sizes
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 09/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, distribution ID
!^*    2.par, parameters
!^* OUT
!^*    1.pdf, result
!^*    2.feas, is par feasible?
!^*    3.isnull, pdf==0?
!^*    4.err, error code
!^*    5.mess, error message
!^**********************************************************************
character(*), intent(in)::DistID
type(mDistPar),intent(in)::par
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!Locals
character(250),parameter::procname='getParFeas'
integer(mik)::i

feas=.true.;err=0;mess=''

select case(DistID)
case(FLAT)! no parameter
    ! nothing to do
case(UNIF) ! bounds in r1
    if(.not.allocated(par%r1)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[r1] not allocated';return;endif
    do i=1,size(par%r1)
        if(size(par%r1(i)%v)/=2) then;err=1;mess=trim(procname)//':'//trim(DistID)//':wrong size [r1]';return;endif
        if(par%r1(i)%v(1)>=par%r1(i)%v(2)) then
            ! mess=trim(procname)//':'//trim(DistID)//':unfeasible bounds (lower>=higher)'
            ! err=1;return
            feas=.false.;return
        endif
    enddo
case(UNIF_IID)! 2 scalar pars: lower bound in spar(1), upper bound in spar(2)
    if(.not.allocated(par%spar)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[spar] not allocated';return;endif
    if(size(par%spar)/=2) then;err=1;mess=trim(procname)//':'//trim(DistID)//':wrong size [spar]';return;endif
    if(par%spar(2)<=par%spar(1)) then;feas=.false.;return;endif
case(GAUSS_IID)! 2 scalar pars: mean in spar(1), stdev in spar(2)
    if(.not.allocated(par%spar)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[spar] not allocated';return;endif
    if(size(par%spar)/=2) then;err=1;mess=trim(procname)//':'//trim(DistID)//':wrong size [spar]';return;endif
    if(par%spar(2)<=0._mrk) then;feas=.false.;return;endif
case(GAUSS_AR1)! 3 scalar pars: constant in spar(1), innovation stdev in spar(2), rho in spar(3)
    if(.not.allocated(par%spar)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[spar] not allocated';return;endif
    if(size(par%spar)/=3) then;err=1;mess=trim(procname)//':'//trim(DistID)//':wrong size [spar]';return;endif
    if(par%spar(2)<=0._mrk .or. abs(par%spar(3))>=1.0_mrk) then;feas=.false.;return;endif
case(GAUSS_AR1_vmean)! mean in vpar, then 2 scalar pars: innovation stdev in spar(1), rho in spar(2)
    if(.not.allocated(par%vpar)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[vpar] not allocated';return;endif
    if(.not.allocated(par%spar)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[spar] not allocated';return;endif
    if(size(par%spar)/=2) then;err=1;mess=trim(procname)//':'//trim(DistID)//':wrong size [spar]';return;endif
    if(par%spar(1)<=0._mrk .or. abs(par%spar(2))>=1.0_mrk) then;feas=.false.;return;endif
case(GAUSS) ! mean in vpar, covariance in mpar
    if(.not.allocated(par%vpar)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[vpar] not allocated';return;endif
    if(.not.allocated(par%mpar)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[mpar] not allocated';return;endif
    if(size(par%vpar)/=size(par%mpar,dim=1) .or. size(par%vpar)/=size(par%mpar,dim=2)) then
        err=1;mess=trim(procname)//':'//trim(DistID)//':wrong size [vpar and/or mpar]';return
    endif
case(NNGP) ! mean in vpar, covariance in mpar, nearest parents in i1, mini-covs in r2
    if(.not.allocated(par%vpar)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[vpar] not allocated';return;endif
    if(.not.allocated(par%mpar)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[i1] not allocated';return;endif
    if(.not.allocated(par%i1)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[i1] not allocated';return;endif
    if(.not.allocated(par%r2)) then;err=1;mess=trim(procname)//':'//trim(DistID)//':[r2] not allocated';return;endif
    if(size(par%vpar)/=size(par%mpar,dim=1) .or. size(par%vpar)/=size(par%mpar,dim=2) .or. &
       size(par%vpar)/=size(par%i1) .or. size(par%vpar)/=size(par%r2) ) then
        mess=trim(procname)//':'//trim(DistID)//':wrong size [vpar and/or mpar and/or i1 and/or r2]'
        err=1;return
    endif
case default
    mess=trim(procname)//':Fatal:Unavailable Dist'
    err=1;return
end select

end subroutine getParFeas
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getLogPdf(DistId,x,par,logp,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: compute logpdf(x|par) for DistID
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 09/07/2019
!^**********************************************************************
!^* Comments: par is INOUT because it might be modified during the
!^*           pdf computation (e.g. covariance matrix is inverted)
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, distribution ID
!^*    2.x, vector of values
!^* INOUT
!^*    3.par, parameters
!^* OUT
!^*    1.logp, result
!^*    2.feas, is par feasible?
!^*    3.isnull, is pdf=0?
!^*    3.err, error code
!^*    4.mess, error message
!^**********************************************************************
character(*), intent(in)::DistID
real(mrk),intent(in)::x(:)
type(mDistPar),intent(inout)::par
real(mrk), intent(out)::logp
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!Locals
character(250),parameter::procname='getLogPdf'
real(mrk)::sInv(size(x),size(x)),ldet,a,b
integer(mik)::i

!Init
logp=UndefRN;feas=.true.;isnull=.false.;err=0;mess=''

!Feasability
call GetParFeas(DistId,par,feas,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not.feas) return

!Compute
select case(DistID)
case(FLAT)
    logp=0._mrk
case(UNIF) ! bounds in r1
    do i=1,size(par%r1)
        a=par%r1(i)%v(1);b=par%r1(i)%v(2)
        if(x(i)<a .or. x(i)>b) then ! out of bounds
            isnull=.true.; return
        else
            par%logps(i)=-1.0_mrk*log(b-a)
        endif
    enddo
    logp=sum(par%logps)
case(UNIF_IID)
    if(any(x<par%spar(1)) .or. any(x>par%spar(2))) then ! out of bounds
        isnull=.true.; return
    else
        logp=-1.0_mrk*real(size(x),mrk)*log(par%spar(2)-par%spar(1))
    endif
case(GAUSS)
    if(par%isCov) then ! slow version where covariance needs to be inverted
        call gauss_logpdf(x=x,mu=par%vpar,sigma=par%mpar,logp=logp,&
                          sigmaInv=sInv,logDet=ldet,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(par%saveInvert) then ! save logdet and overwrite covariance in par%mpar
            par%logDet=ldet
            par%mpar=sInv;par%isCov=.false.
        endif
    else ! fast version where cov. already inverted
        call gauss_logpdf_fast(x=x,mu=par%vpar,sigmainv=par%mpar,logdet=par%logDet,logp=logp,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    endif
case(GAUSS_IID)
    call gaussiid_logpdf(x=x,mu=par%spar(1),sdev=par%spar(2),&
                         logps=par%logps,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    logp=sum(par%logps)
case(GAUSS_AR1)
    call gaussAR1_logpdf(x=x,c=par%spar(1)+0._mrk*x,sdev=par%spar(2),rho=par%spar(3),&
                         logps=par%logps,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    logp=sum(par%logps)
case(GAUSS_AR1_vmean)
    call gaussAR1_logpdf(x=x,c=par%vpar,sdev=par%spar(1),rho=par%spar(2),&
                         logps=par%logps,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    logp=sum(par%logps)
case(NNGP)
    call NNGP_logpdf(x=x,mu=par%vpar,sigma=par%mpar,Pa=par%i1,invCovs=par%r2,&
                     doInvert=par%isCov,logps=par%logps,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    logp=sum(par%logps)
    par%isCov=.false. ! mini-covariances have been inverted
case default
    mess=trim(procname)//':Fatal:Unavailable Dist'
    err=1;return
end select

end subroutine getLogPdf
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getLogPdf_singleCompUpdate(DistId,x,par,oldLogp,iComp,inc,&
                                      newLogp,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: update logpdf(x|par) when a single component 'iComp' of
!^*          vector 'x' is incremented by 'inc'.
!^*          This is useful for one-at-time MCMC samplers because a
!^*          smart logpdf update is often way more efficient than
!^*          recomputing the logpdf from scratch as if x was all new.
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 05/09/2019
!^**********************************************************************
!^* Comments: Assumes any inversion of covariance has already been done.
!^*           Will exit with an error if not the case.
!^*           Note that par is INOUT to allow updating par%logps
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, distribution ID
!^*    2.x, vector of values
!^*    3.oldLogp, log-pdf value before update
!^*    4.iComp, component of 'x' that is updated
!^*    5.inc, increment ( newx(iComp)=x(iComp)+inc )
!^* INOUT
!^*    1.par, parameters
!^* OUT
!^*    1.newLogp, updated log-pdf
!^*    2.feas, is par feasible?
!^*    3.isnull, newLogp==0?
!^*    4.err, error code
!^*    5.mess, error message
!^**********************************************************************
character(*), intent(in)::DistID
real(mrk),intent(in)::x(:),oldLogp,inc
type(mDistPar),intent(inout)::par
integer(mik),intent(in)::iComp
real(mrk), intent(out)::newLogp
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!Locals
character(250),parameter::procname='getLogPdf_singleCompUpdate'
real(mrk)::newx(size(x)),beforeInc,m,s,a,b

!Init
newLogp=oldLogp;feas=.true.;isnull=.false.;err=0;mess=''

!Compute
select case(DistID)
case(FLAT)
    newLogp=0._mrk
case(UNIF)
    a=par%r1(iComp)%v(1);b=par%r1(iComp)%v(2)
    if(x(iComp)+inc<a .or. x(iComp)+inc>b) then ! out of bounds
        isnull=.true.; return
    else ! still within bounds, log-pdf does not change
        newLogp=oldLogp
    endif
case(UNIF_IID)
    a=par%spar(1);b=par%spar(2)
    if(x(iComp)+inc<a .or. x(iComp)+inc>b) then ! out of bounds
        isnull=.true.; return
    else ! still within bounds, log-pdf does not change
        newLogp=oldLogp
    endif
case(GAUSS)
    if(par%isCov) then ! return with error, cov should already be inverted at this stage
        mess=trim(procname)//':covariance matrix should already be inverted'
        err=1;return
    endif
    newLogp=oldLogp - 0.5_mrk * (2._mrk*inc*sum((x-par%vpar)*par%mpar(:,iComp))+par%mpar(iComp,iComp)*inc**2)
case(GAUSS_IID)
    m=par%spar(1)
    s=par%spar(2)
    beforeInc=par%logps(iComp)
    par%logps(iComp)=dnorm(x(iComp)+inc,m,s)
    newLogp=oldLogp+par%logps(iComp)-beforeInc
case(GAUSS_AR1)
    call gaussAR1_singleCompUpdate(x=x,c=par%spar(1)+0._mrk*x,&
                                   sdev=par%spar(2),rho=par%spar(3),&
                                   oldLogp=oldLogp,iComp=iComp,inc=inc,&
                                   logps=par%logps,newLogp=newLogp,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
case(GAUSS_AR1_vmean)
    call gaussAR1_singleCompUpdate(x=x,c=par%vpar,&
                                   sdev=par%spar(2),rho=par%spar(3),&
                                   oldLogp=oldLogp,iComp=iComp,inc=inc,&
                                   logps=par%logps,newLogp=newLogp,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
case(NNGP)
    call NNGP_singleCompUpdate(x=x,mu=par%vpar,sigma=par%mpar,&
                               Pa=par%i1,invCovs=par%r2,&
                               oldLogp=oldLogp,iComp=iComp,inc=inc,&
                               logps=par%logps,newLogp=newLogp,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
case default ! no smart computer available, compute logpdf from scratch
    newx=x;newx(iComp)=x(iComp)+inc
    call getLogPdf(DistId=DistId,x=newx,par=par,&
                   logp=newLogp,feas=feas,isnull=isnull,err=err,mess=mess)
end select

end subroutine getLogPdf_singleCompUpdate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!*********!
! PRIVATE !
!*********!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine gauss_logpdf(x,mu,sigma,logp,sigmaInv,logDet,err,mess)
!^**********************************************************************
!^* Purpose: gaussian log pdf, slow version (sigma needs to be inverted)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 09/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, vector of values
!^*    2.mu, mean vector
!^*    3.sigma, covariance
!^* OUT
!^*    1.logp, result
!^*    2.sigmaInv, inverse covariance
!^*    3.logdet, log-determinant
!^*    4.err, error code
!^*    5.mess, error message
!^**********************************************************************
use linalg_dmsl_kit,only:choles_invrt
real(mrk), intent(in)::x(:),mu(:),sigma(:,:)
real(mrk), intent(out)::logp,sigmaInv(:,:),logDet
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='gauss_logpdf'
logical::feas

logp=UndefRN;sigmaInv=undefRN;logDet=undefRN;err=0;mess=''

sigmaInv=sigma
call choles_invrt(ainv=sigmaInv,logDet=logDet,posDefinite=feas,err=err,message=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not.feas) return

call gauss_logpdf_fast(x,mu,sigmaInv,logDet,logp,err,mess)

end subroutine gauss_logpdf
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine gauss_logpdf_fast(x,mu,sigmainv,logdet,logp,err,mess)
!^**********************************************************************
!^* Purpose: gaussian log pdf, fast version (sigma already inverted)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 09/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, vector of values
!^*    2.mu, mean vector
!^*    3.sigmainv, inverse covariance
!^*    4.logdet, log-determinant
!^* OUT
!^*    1.logp, result
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:lnTwoPi
real(mrk), intent(in)::x(:),mu(:),sigmainv(:,:),logdet
real(mrk), intent(out)::logp
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='gauss_logpdf_fast'
integer(mik)::n
real(mrk)::centred(size(x)),z

logp=UndefRN;err=0;mess=''
n=size(x)
centred=x-mu
z=dot_product(matmul(centred,sigmainv),centred)
logp=-0.5_mrk*(real(n,mrk)*lnTwoPi+logdet+z)

end subroutine gauss_logpdf_fast
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine gaussiid_logpdf(x,mu,sdev,logps,err,mess)
!^**********************************************************************
!^* Purpose: iid gaussian log pdf COMPONENTS
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 09/07/2019
!^**********************************************************************
!^* Comments: WARNING: the full log-pdf is sum(logps)
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, vector of values
!^*    2.mu, mean
!^*    3.sdev, standard deviation
!^* OUT
!^*    1.logps, result
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:lnTwoPi
real(mrk), intent(in)::x(:),mu,sdev
real(mrk), intent(out)::logps(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='gaussiid_logpdf'
integer(mik)::n
real(mrk)::centred(size(x)),z

logps=UndefRN;err=0;mess=''
logps=dnorm(x,mu,sdev)

end subroutine gaussiid_logpdf
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine gaussAR1_logpdf(x,c,sdev,rho,logps,err,mess)
!^**********************************************************************
!^* Purpose: gaussian AR(1) log pdf COMPONENTS
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 26/07/2019
!^**********************************************************************
!^* Comments: WARNING: the full log-pdf is sum(logps)
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, vector of values
!^*    2.c, vector of constants
!^*    3.sdev, stdev of innovations
!^*    4.rho, lag-1 autocorr
!^* OUT
!^*    1.logps, result
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************
real(mrk), intent(in)::x(:),c(:),sdev,rho
real(mrk), intent(out)::logps(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='gaussAR1_logpdf'
integer(mik)::i
real(mrk)::lp

logps=UndefRN;err=0;mess=''

! first data, use marginal distribution
logps(1)=dnorm(x=x(1),mu=c(1)/(1.0_mrk-rho),sdev=sdev/sqrt(1.0_mrk-rho**2))
! Conditional distributions for subsequent data
do i=2,size(x)
    logps(i)=dnorm(x=x(i),mu=c(i)+rho*x(i-1),sdev=sdev)
enddo

end subroutine gaussAR1_logpdf
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NNGP_logpdf(x,mu,sigma,Pa,invCovs,doInvert,logps,err,mess)
!^**********************************************************************
!^* Purpose: Nearest Neighbors Gaussian Process log pdf COMPONENTS
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 04/09/2019
!^**********************************************************************
!^* Comments: WARNING: the full log-pdf is sum(logps)
!^**********************************************************************
!^* References: Datta et al. 2016 and Banerjee 2017
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, vector of values
!^*    2.mu, mean vector of the underlying Gaussian process
!^*    3.sigma, covariance matrix of the underlying Gaussian process
!^*    4.Pa, nearest parents. Pa(i)%v(:) contains the indices of the nearest parents of site i
!^*    5.doInvert. if false, assume invCovs contains already-inverted matrices, otherwise invert and store in invCovs
!^* INOUT
!^*    1.invCov. invCovs(i)%v(:,:) contains the (small!) covariance matrix for the nearest parents of site i
!^* OUT
!^*    1.logps, vector of logpdf components (full log-pdf is sum(logps))
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************
use linalg_dmsl_kit, only:choles_invrt
use numerix_dmsl_kit,only:normal_logp
use types_dmsl_kit,only:vectorA2_r_type,vectorA1_i_type
real(mrk),intent(in)::x(:),mu(:),sigma(:,:)
type(vectorA1_i_type),intent(in)::Pa(:)
type(vectorA2_r_type),intent(inout)::invCovs(:)
logical,intent(in)::doInvert
real(mrk),intent(out)::logps(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(*),parameter::procnam="NNGP_logpdf"
integer(mik)::i,n,k
real(mrk)::mi,vi,muCond,varCond,logDet
logical::posdef

logps=undefRN;err=0;mess=''
! init
n=size(x)
do i=1,n
    k=size(Pa(i)%v) ! number of parents
    if(k>0) then
        if(doInvert) then
            if(allocated(invCovs(i)%v)) deallocate(invCovs(i)%v);allocate(invCovs(i)%v(k,k))
            invCovs(i)%v=sigma(Pa(i)%v,Pa(i)%v)
            call choles_invrt(ainv=invCovs(i)%v,posDefinite=posdef,logDet=logDet,&
                              err=err,message=mess)
            if(err/=0) return
            if(.not. posdef) then;err=1;mess=trim(procnam)//':notPosDef';return;endif
        endif
        ! conditioning terms modifying the marginal mean/variance
        muCond=dot_product(matmul(sigma(i,Pa(i)%v),invCovs(i)%v),x(Pa(i)%v)-mu(Pa(i)%v))
        varCond=dot_product(matmul(sigma(i,Pa(i)%v),invCovs(i)%v),sigma(Pa(i)%v,i))
    else
        muCond=0._mrk
        varCond=0._mrk
    endif
    ! conditional mean and variance
    mi=mu(i)+muCond
    vi=sigma(i,i)-varCond
    logps(i)=normal_logp(x=x(i), mean=mi, var=vi)
enddo

end subroutine NNGP_logpdf
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NNGP_singleCompUpdate(x,mu,sigma,Pa,invCovs,&
                                 oldLogp,iComp,inc,&
                                 logps,newLogp,err,mess)
!^**********************************************************************
!^* Purpose: Nearest Neighbors Gaussian Process log pdf COMPONENTS
!^*          Implements the smart update of logps' when a single
!^*          component icomp is incremented.
!^*          invCovs contains already-inverted matrices
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 05/09/2019
!^**********************************************************************
!^* Comments: WARNING: the full log-pdf is sum(logps)
!^**********************************************************************
!^* References: Datta et al. 2016 and Banerjee 2017
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, vector of values
!^*    2.mu, mean vector of the underlying Gaussian process
!^*    3.sigma, covariance matrix of the underlying Gaussian process
!^*    4.Pa, nearest parents. Pa(i)%v(:) contains the indices of the nearest parents of site i
!^*    5.invCov. invCovs(i)%v(:,:) contains the (small!) covariance matrix for the nearest parents of site i
!^*    6.oldLogp, joint log-pdf before update
!^*    7.iComp, component of 'x' that is updated
!^*    8.inc, increment ( newx(iComp)=x(iComp)+inc )
!^* INOUT
!^*    1.logps, vector of logpdf components (full log-pdf is sum(logps))
!^* OUT
!^*    1.newLogp, joint log-pdf after update
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************
use types_dmsl_kit,only:vectorA2_r_type,vectorA1_i_type
real(mrk),intent(in)::x(:),mu(:),sigma(:,:),inc,oldLogp
integer(mik),intent(in)::iComp
type(vectorA1_i_type),intent(in)::Pa(:)
type(vectorA2_r_type),intent(in)::invCovs(:)
real(mrk),intent(inout)::logps(:)
real(mrk),intent(out)::newLogp
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(*),parameter::procnam="NNGP_singleCompUpdate"
integer(mik)::i,n,k
real(mrk)::mi,vi,muCond,sigmaCond,logDet,newx(size(x)),beforeInc
logical::posdef

newLogp=oldLogP;err=0;mess=''
! init
n=size(x)
newx=x;newx(iComp)=x(iComp)+inc
do i=1,n
    if(i==iComp .or. any(Pa(i)%v==iComp)) then ! needs update
        k=size(Pa(i)%v) ! number of parents
        if(k>0) then
            ! conditioning terms modifying the marginal mean/variance
            muCond=dot_product(matmul(sigma(i,Pa(i)%v),invCovs(i)%v),newx(Pa(i)%v)-mu(Pa(i)%v))
            sigmaCond=dot_product(matmul(sigma(i,Pa(i)%v),invCovs(i)%v),sigma(Pa(i)%v,i))
        else
            muCond=0._mrk
            sigmaCond=0._mrk
        endif
        ! conditional mean and variance
        mi=mu(i)+muCond
        vi=sigma(i,i)-sigmaCond
        beforeInc=logps(i)
        logps(i)=dnorm(x=newx(i),mu=mi,sdev=sqrt(vi))
        newLogp=newLogp+logps(i)-beforeInc
    endif
enddo

end subroutine NNGP_singleCompUpdate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine gaussAR1_singleCompUpdate(x,c,sdev,rho,&
                                     oldLogp,iComp,inc,&
                                     logps,newLogp,err,mess)
!^**********************************************************************
!^* Purpose: Gaussian AR(1) Process log pdf COMPONENTS
!^*          Implements the smart update of logps' when a single
!^*          component icomp is incremented.
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 09/02/2021
!^**********************************************************************
!^* Comments: WARNING: the full log-pdf is sum(logps)
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, vector of values
!^*    2.c, vector of constants
!^*    3.sdev, stdev of innovations
!^*    4.rho, lag-1 autocorr
!^*    5.oldLogp, joint log-pdf before update
!^*    6.iComp, component of 'x' that is updated
!^*    7.inc, increment ( newx(iComp)=x(iComp)+inc )
!^* INOUT
!^*    1.logps, vector of logpdf components (full log-pdf is sum(logps))
!^* OUT
!^*    1.newLogp, joint log-pdf after update
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************
real(mrk),intent(in)::x(:),c(:),sdev,rho,inc,oldLogp
integer(mik),intent(in)::iComp
real(mrk),intent(inout)::logps(:)
real(mrk),intent(out)::newLogp
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(*),parameter::procnam="gaussAR1_singleCompUpdate"
real(mrk)::beforeInc,m,s

! init

newLogp=oldLogP;err=0;mess=''
! update component iComp of AR1 pdf
if(iComp==1) then ! marginal gaussian
    m=c(1)/(1._mrk-rho)
    s=sdev/sqrt(1._mrk-rho**2)
else ! conditional Gaussian
    m=c(iComp)+rho*x(iComp-1)
    s=sdev
endif
beforeInc=logps(iComp)
logps(iComp)=dnorm(x(iComp)+inc,m,s)
newLogp=oldLogp+logps(iComp)-beforeInc
! update component iComp+1 of AR1 pdf, since it is conditioned on x(iComp)
if(iComp<size(x)) then
    m=c(iComp+1)+rho*(x(iComp)+inc)
    s=sdev
    beforeInc=logps(iComp+1)
    logps(iComp+1)=dnorm(x(iComp+1),m,s)
    newLogp=oldLogp+logps(iComp+1)-beforeInc
endif
end subroutine gaussAR1_singleCompUpdate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NNGP_getNN(k,sigma,Pa,err,mess)
!^**********************************************************************
!^* Purpose: get Nearest Neighbors from the covariance matrix sigma
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 12/09/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References: Datta et al. 2016 and Banerjee 2017
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.k, max number of neighbors
!^*    2.sigma, covariance matrix of the underlying Gaussian process
!^* OUT
!^*    1.Pa, nearest parents. Pa(i)%v(:) contains the indices of the nearest parents of site i
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************
use types_dmsl_kit,only:vectorA1_i_type
use numerix_dmsl_kit,only:indexx_qsort
real(mrk),allocatable::foo(:)
integer(mik),allocatable::ifoo(:)

integer(mik),intent(in)::k
real(mrk),intent(in)::sigma(:,:)
type(vectorA1_i_type),intent(out)::Pa(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(*),parameter::procname="NNGP_getNN"
integer(mik)::i,n,kk

err=0;mess=''
! init
n=size(Pa)
! first site has no parents
if(allocated(Pa(1)%v)) deallocate(Pa(1)%v);allocate(Pa(1)%v(0))
! loop on subsequent sites
do i=2,n
    ! covariances to parents
    if(allocated(foo)) deallocate(foo);allocate(foo(i-1))
    if(allocated(ifoo)) deallocate(ifoo);allocate(ifoo(i-1))
    foo=sigma(i,1:(i-1))
    ! get k most-correlated (i.e. nearest) parents
    call indexx_qsort(arr=foo,indx=ifoo,ascnd=.false.,err=err,message=mess) ! could be more efficient - only need the first k indices.
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    kk=min(i-1,k)
    allocate(Pa(i)%v(kk))
    Pa(i)%v=ifoo(1:kk)
  enddo

end subroutine NNGP_getNN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function dnorm(x,mu,sdev)
! function version for a 1-D gaussian log-pdf
use utilities_dmsl_kit,only:lnTwoPi
real(mrk),intent(in)::x,mu,sdev
real(mrk)::dnorm
! locals
real(mrk)::z

z= ((x-mu)**2)/(sdev**2)
dnorm= -1._mrk*log(sdev)-0.5_mrk*(lnTwoPi+z)

end function dnorm
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module MultivariateDistribution_tools

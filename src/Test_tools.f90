module Test_tools

!~**********************************************************************
!~* Purpose: Test utilities 
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 20/02/2012
!~**********************************************************************
!~* Comments: 
!~* 
!~* 
!~* 
!~* 
!~* 
!~* 
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
public :: KendallTest, RegionalTest

! Shortcut for test names
Character(100), parameter, PUBLIC:: MK='Mann-Kendall',Pettitt='Pettitt'

! Type "TestResultType" encapasulating all usefull characteristics & results of the test
type,public::TestResultType
	real(mrk)::level=undefRN ! error level (in interval (0;1)) at which test is applied
	logical::H0=.true. ! result; true = "H0 is true" ie do not reject it, false = the contrary
	real(mrk)::pval=undefRN ! p-value of the test
	real(mrk)::stat=undefRN !test statistic
	real(mrk)::critical=undefRN ! critical value of the test for error level "level"
	real(mrk)::Xtra_real=undefRN ! used to pass any useful by-product of the test (eg location of a step change, Kendall correlation, etc)
	real(mrk),pointer::Xtra_rv(:) => NULL() ! Same as above for passing a vector
	integer(mik)::Xtra_integer=undefIN ! Same as above for passing an integer
	integer(mik),pointer::Xtra_iv(:) => NULL()  ! Same as above for passing an integer vector
	character(250)::Xtra_string='' ! Same as above for passing a comment
	character(250),pointer::Xtra_sv(:) => NULL()! Same as above for passing many comments 
end type TestResultType

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine KendallTest(x,y,level,out,err,mess)

!^**********************************************************************
!^* Purpose: Kendall's tau test of association between x & y
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:23/02/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. x
!^*		2. y
!^*		3. level
!^* OUT
!^*		1.out, test result object
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetCdf,GetQuantile

real(mrk), intent(in)::x(:),y(:),level
type(TestResultType), intent(out)::out
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::i,j,n,S,xN,yN,xExt(size(x)),yExt(size(y))
real(mrk)::var,tau,v0,vt,vu,v1,v2,ti,s1,s2,&
            s3,s4,cdf,foo,n0,n1,n2
logical::Xties, Yties, feas

err=0;mess='';out%level=level
! Check Size
if(size(x)/=size(y)) then
    err=1;mess='KendallTest:Fatal:size mismatch [x,y]';return
endif
n=size(x)
if(n<2) then
    err=1;mess='KendallTest:Fatal:size[x,y]<2';return
endif

S=0;Xties=.false.;Yties=.false.
do i=1,n-1
    do j=i+1,n
        foo=sign(1._mrk,x(j)-x(i))*sign(1._mrk,y(j)-y(i)) ! discordant or concordant?
        if(x(i)==x(j)) then
            Xties=.true.;cycle !tied values - skip
        endif
        if(y(i)==y(j)) then
            Yties=.true.;cycle !tied values - skip
        endif
        S=S+nint(foo)
    enddo
enddo

! Start variance computation
! Basic variance component, no ties correction
v0=real(n*(n-1)*(2*n+5),mrk)
! ties-in-x-component
if(.not. Xties) then ! don't waste time 
    vt=0._mrk
    n1=0._mrk
else ! compute component
    call GetTies(x,xN,xExt)
    vt=0._mrk;n1=0._mrk
    do i=1,xN
        ti=real(xExt(i),mrk)
        vt=vt + ti * (ti-1._mrk) * (2._mrk*ti+5._mrk)
        n1=n1+ 0.5_mrk * ti * (ti-1._mrk)
    enddo
endif
! ties-in-y-component
if(.not. Yties) then ! don't waste time 
    vu=0._mrk
    n2=0._mrk
else ! compute component
    call GetTies(y,yN,yExt)
    vu=0._mrk;n2=0._mrk
    do i=1,yN
        ti=real(yExt(i),mrk)
        vu=vu + ti * (ti-1._mrk) * (2._mrk*ti+5._mrk)
        n2=n2+ 0.5_mrk * ti * (ti-1._mrk)
    enddo
endif
! ties-in-both-x-and-y-component
if( (.not. Xties) .or. (.not. Yties) )then ! don't waste time 
    v1=0._mrk;v2=0._mrk
else
    s1=0._mrk;s2=0._mrk;s3=0._mrk;s4=0._mrk
    do i=1,xN
        ti=real(xExt(i),mrk)
        s1=s1 + ti * (ti-1._mrk)
        s3=s3 + ti * (ti-1._mrk) * (ti-2._mrk)
    enddo
    do i=1,yN
        ti=real(yExt(i),mrk)
        s2=s2 + ti * (ti-1._mrk)
        s4=s4 + ti * (ti-1._mrk) * (ti-2._mrk)
    enddo
    v1=(s1*s2)/real(2*n*(n-1),mrk)
    v2=(s3*s4)/real(9*n*(n-1)*(n-2),mrk)
endif

! Get Kendall's tau
n0=0.5_mrk*n*(n-1)
tau=real(S,mrk)/sqrt((n0-n1)*(n0-n2))
out%Xtra_real=tau

! Final variance value
var=(v0-vt-vu)/18._mrk+v1+v2
if(var<=0._mrk) then
    err=1;mess='KendallTest:Fatal:Var[S]<=0!!';return
endif
!continuity correction
if(S>0) then
    out%stat=real(S-1,mrk)/sqrt(var)
elseif(S<0) then
    out%stat=real(S+1,mrk)/sqrt(var)
else
    out%stat=0._mrk
endif

! Test results
! p-value
call GetCdf(DistId='Gaussian',x=-1._mrk*abs(out%stat),par=(/0._mrk,1._mrk/),&
            cdf=cdf,feas=feas,err=err,mess=mess)
if(.not. feas) then
    mess='KendallTest:Fatal:Cdf computation unfeasible';err=1;return
endif
if(err>0) then
    mess='KendallTest:'//trim(mess);return
endif
out%pval=2._mrk*cdf
! Decision
out%H0=(out%pval>out%level)
! critical value
call GetQuantile(DistId='Gaussian',p=1._mrk-0.5_mrk*out%level,par=(/0._mrk,1._mrk/),&
            q=out%critical,feas=feas,err=err,mess=mess)
if(.not. feas) then
    mess='KendallTest:Fatal:Quantile computation unfeasible';err=1;return
endif
if(err>0) then
    mess='KendallTest:'//trim(mess);return
endif
! Messages
if(n<=10) out%Xtra_string='Warning: small sample: Asymptotic Gaussian approximation might be inaccurate'

end subroutine KendallTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine RegionalTest(y,x,level,mv,MVoption,out,err,mess)

!^**********************************************************************
!^* Purpose: Regional test, see Renard et al WRR 2008
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:22/10/2012
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Checks, especially for option 2
!^**********************************************************************
!^* IN
!^*		1. y (eg time)
!^*		2. x (data)
!^*		3. level
!^*		4. mv, value used for missing data
!^*		5. [MVoption], 1=accept no MV (safer), 2=handle MV (needs more check)
!^* OUT
!^*		1.out, test result object
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
use EmpiricalStats_tools, only:GetRank
use Distribution_tools, only:GetQuantile
use linalg_dmsl_kit, only:choles_invrt
use numerix_dmsl_kit,only:chi2p_cdf,chi2_icdf,mnormal_logp_sub
use utilities_dmsl_kit, only:number_string
real(mrk), intent(in)::y(:),x(:,:),level,mv
integer(mik), intent(in), optional::MVoption
type(TestResultType), intent(out)::out
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik),parameter::MinCommonVal=5, MVOdefault=1
integer(mik)::n,p,i,j,k,ncommon,indxlist(size(x,dim=2)),MVO
integer(mik),allocatable::indx(:)
real(mrk)::ytilde(size(y)),xtilde(size(x,dim=1),size(x,dim=2)),z(size(y)),&
            Sigma(size(x,dim=2),size(x,dim=2)),SigmaInv(size(x,dim=2),size(x,dim=2)),&
            s1bar,s2bar,g0(size(x,dim=2)),g1,g2,&
            z0(size(x,dim=2)),z1(size(x,dim=1)),z2(size(x,dim=1)),z3(size(x,dim=1)),&
            beta,LogDet,f0,f1,LL0,LL1,num,den,ones(size(x,dim=2)),LD(size(x,dim=1))
logical::feas,mask(size(x,dim=1)),mask2(size(x,dim=2))
real(mrk),allocatable::s1(:),s2(:),Sig(:,:),one(:),g4(:)
type SingInvStorage
    real(mrk),allocatable::mat(:,:)
end type SingInvStorage
type(SingInvStorage)::SigInv(size(x,dim=1))

err=0;mess='';out%level=level
if(present(MVOption)) then
    MVO=MVoption
else
    MVO=MVOdefault
endif
n=size(x,dim=1);p=size(x,dim=2);ones=1._mrk
if(any(y==mv)) then
    mess='RegionalTest:FATAL:MV not allowed in predictor y';err=1;return
endif
! center y
ytilde=y-sum(y)/real(n,mrk)

if(MVO==1) then 
    ! normal-score-transform x
    xtilde=x
    do i=1,p
        if(any(x(:,i)==mv)) then
            mess='RegionalTest:FATAL:MV not allowed within this option';err=1;return
        endif
        call GetRank(x(:,i),mv,z,err,mess)
        if(err>0) then
            mess='RegionalTest:'//trim(mess);return
        endif
        do j=1,n
            call GetQuantile(distID='Gaussian',p=(z(j)-0.5_mrk)/real(n,mrk),&
                    par=(/0._mrk,1._mrk/),q=xtilde(j,i),feas=feas,err=err,mess=mess)
            if(err>0) then
                mess='RegionalTest:'//trim(mess);return
            endif
        enddo
    enddo

    ! Get correlation matrix
    Sigma=(1._mrk/real(n,mrk))*matmul(transpose(xtilde),xtilde)
    ! Invert correlation matrix
    call choles_invrt(a=Sigma,ainv=SigmaInv,LogDet=LogDet,posDefinite=feas,err=err,message=mess)
    if(.not.feas) then
        err=1;mess='RegionalTest:Corr Matrix not Def-pos-np='//trim(number_string(n))//'-'//trim(number_string(p));return
    endif
    if(err>0) then
        mess='RegionalTest:'//trim(mess);return
    endif
    ! Compute trend estimate
    g0=matmul(ones,SigmaInv)
    z1=matmul(g0,transpose(xtilde))
    num=dot_product(z1,ytilde)
    g1=dot_product(ytilde,ytilde)
    g0=matmul(ones,SigmaInv)
    den=g1*dot_product(g0,ones)
    beta=num/den
    out%Xtra_real=beta
    ! Get likelihoods
    LL0=0._mrk;LL1=0._mrk
    do i=1,n
        ! Model M0
        call mnormal_logp_sub(x=xtilde(i,:),mean=0._mrk*ones,do_invrs=.false.,var_invrs=SigmaInv,logDet=LogDet,&
                        posdef=feas,logp=f0,err=err,message=mess)
        if(.not.feas) then
            err=1;mess='RegionalTest:Corr Matrix not Def-pos';return
        endif
        if(err>0) then
            mess='RegionalTest:'//trim(mess);return
        endif
        ! Model M1
        call mnormal_logp_sub(x=xtilde(i,:),mean=beta*ytilde(i)*ones,do_invrs=.false.,var_invrs=SigmaInv,logDet=LogDet,&
                        posdef=feas,logp=f1,err=err,message=mess)
        if(.not.feas) then
            err=1;mess='RegionalTest:Corr Matrix not Def-pos';return
        endif
        if(err>0) then
            mess='RegionalTest:'//trim(mess);return
        endif
        LL0=LL0+f0
        LL1=LL1+f1
    enddo
    out%stat=-2._mrk*(LL0-LL1)
    out%pval=1._mrk-chi2p_cdf(max(0._mrk,out%stat),1._mrk)
    out%H0=(out%pval>1._mrk-out%level)
    out%critical=chi2p_cdf(1._mrk-out%level,1._mrk)
else
    ! Option 2: handles MV, but poorly-checked so far...
    ! normal-score-transform x
    xtilde=x
    do i=1,p
        call GetRank(x(:,i),mv,z,err,mess)
        if(err>0) then
            mess='RegionalTest:'//trim(mess);return
        endif
        do j=1,n
            if(z(j)/=mv) then
                call GetQuantile(distID='Gaussian',p=(z(j)-0.5_mrk)/real(n,mrk),&
                        par=(/0._mrk,1._mrk/),q=xtilde(j,i),feas=feas,err=err,mess=mess)
                if(err>0) then
                    mess='RegionalTest:'//trim(mess);return
                endif
            else
                xtilde(j,i)=mv
            endif
        enddo
    enddo
    ! Get correlation matrix
    Sigma=1._mrk
    do i=2,p ! Note: pairwise computation makes the best use of available data but doesn't ensure definite-positiveness...
        do j=1,i-1
            mask=((xtilde(:,i)/=mv).and.(xtilde(:,j)/=mv))
            ncommon=count(mask)
            if(ncommon<MinCommonVal) then
                mess='RegionalTest:FATAL:Not enough common years to get a proper correlation';err=1;return
            endif
            if(allocated(s1)) deallocate(s1);allocate(s1(ncommon))
            if(allocated(s2)) deallocate(s2);allocate(s2(ncommon))
            s1=pack(xtilde(:,i),mask);s2=pack(xtilde(:,j),mask)
            s1bar=sum(s1)/real(ncommon,mrk);s2bar=sum(s2)/real(ncommon,mrk)
            Sigma(i,j)=(1._mrk/real(ncommon,mrk))*dot_product(s1-s1bar,s2-s2bar)
            Sigma(j,i)=Sigma(i,j)
        enddo
    enddo
    ! Invert correlation matrix
    ! Needs to be repeated at each time step => computationaly expansive!
    ! 2DO: implement Xun's approach to avoid repeating the same inversion
    indxlist=(/(i,i=1,p)/)
    do i=1,n
        ! Get data available at this time step
        mask2=xtilde(i,:)/=mv;k=count(mask2);if(k==0) cycle
        if(allocated(indx)) deallocate(indx);allocate(indx(k))
        indx=pack(indxlist,mask2) ! index of available data
        ! get corr sub-matrix, invert it and store it for subsequent use in likelihood computation
        if(allocated(Sig)) deallocate(Sig);allocate(Sig(k,k))
        Sig=Sigma(indx,indx)
        if(allocated(SigInv(i)%mat)) deallocate(SigInv(i)%mat);allocate(SigInv(i)%mat(k,k))
        call choles_invrt(a=Sig,ainv=SigInv(i)%mat,LogDet=LD(i),posDefinite=feas,err=err,message=mess)
        if(.not.feas) then
            err=1;mess='RegionalTest:Corr Matrix not Def-pos';return
        endif
        if(err>0) then
            mess='RegionalTest:'//trim(mess);return
        endif
        ! Compute trend estimate
        ! Pack available data
        if(allocated(s1)) deallocate(s1);allocate(s1(k))
        s1=pack(xtilde(i,:),mask2)
        ! columnwise and entire sums of inverse matrix
        if(allocated(g4)) deallocate(g4);allocate(g4(k))
        g4=sum(SigInv(i)%mat,dim=2)
        z2(i)=sum(SigInv(i)%mat)
        ! dot product of g0 with packed data
        z1(i)=dot_product(s1,g4)
        ! square predictor
        z3(i)=ytilde(i)**2
    enddo
    ! compute trend
    num=dot_product(ytilde,z1)
    den= dot_product(z3,z2)
    beta=num/den
    out%Xtra_real=beta
    ! Get likelihoods
    LL0=0._mrk;LL1=0._mrk
    do i=1,n
        ! Get available data 
        mask2=xtilde(i,:)/=mv;k=count(mask2);if(k==0) cycle
        if(allocated(s1)) deallocate(s1);allocate(s1(k))
        s1=pack(xtilde(i,:),mask2)
        if(allocated(one)) deallocate(one);allocate(one(k))
        one=1._mrk
        ! Model M0
        call mnormal_logp_sub(x=s1,mean=0._mrk*one,do_invrs=.false.,var_invrs=SigInv(i)%mat,logDet=LD(i),&
                        posdef=feas,logp=f0,err=err,message=mess)
        if(.not.feas) then
            err=1;mess='RegionalTest:Corr Matrix not Def-pos';return
        endif
        if(err>0) then
            mess='RegionalTest:'//trim(mess);return
        endif
        ! Model M1
        call mnormal_logp_sub(x=s1,mean=beta*ytilde(i)*one,do_invrs=.false.,var_invrs=SigInv(i)%mat,logDet=LD(i),&
                        posdef=feas,logp=f1,err=err,message=mess)
        if(.not.feas) then
            err=1;mess='RegionalTest:Corr Matrix not Def-pos';return
        endif
        if(err>0) then
            mess='RegionalTest:'//trim(mess);return
        endif
        LL0=LL0+f0
        LL1=LL1+f1
    enddo
    out%stat=-2._mrk*(LL0-LL1)
    out%pval=1._mrk-chi2p_cdf(max(0._mrk,out%stat),1._mrk)
    out%H0=(out%pval>1._mrk-out%level)
    out%critical=chi2p_cdf(1._mrk-out%level,1._mrk)
endif
end subroutine RegionalTest
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Private

subroutine GetTies(x,Ngroups,Extent)
! Get tied values from a series x
! Ngroups = number of groups with identical values
! Extent = extent of each group
! Example: (1,2,3,1,1,3,4,0) => we have 3 times the value "1" and twice the value "3"
! => Ngroups=2, Extent = (3,2,undefIN,undefIN,undefIN,undefIN,undefIN,undefIN)
! Note that Extent has same size than x, but only the first Ngroups values are meaningful
! This is admitedly a bit ugly, but is made as a quick-fix to avoid using pointers that could be associated
! with zero-sized arrays when there is no ties (I'm not sure how the compiler handles this...) 
use numerix_dmsl_kit, only:quicksort

real(mrk), intent(in)::x(:)
integer(mik), intent(out)::Ngroups,Extent(size(x))
! locals
integer(mik)::i,n,err,p
real(mrk)::arr(size(x))
logical::mask(size(x))

Extent=undefIN
n=size(x)
arr=x
call quicksort(arr=arr,ascnd=.true.,err=err)
i=1;Ngroups=0
do while(i<=n)
    mask=(arr==arr(i))
    p=count(mask) ! number of ties for this value
    if(p==0) then
        write(*,*) "GetTies: Fatal: this shouldn't happen..."
        write(*,*) "You're likely to enter an infinite while loop..."
        read(*,*)
    endif
    if(p>1) then ! ties exist
        Ngroups=Ngroups+1
        Extent(Ngroups)=p
    endif
    i=i+p !increment
enddo

end subroutine GetTies

end module Test_tools
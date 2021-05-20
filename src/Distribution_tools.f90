module Distribution_tools

!~**********************************************************************
!~* Purpose: Distribution utilities
!~**********************************************************************
!~* Programmer: Ben Renard, University of Newcastle
!~**********************************************************************
!~* Last modified:24/11/2011
!~**********************************************************************
!~* Comments: subs GetSamplePdf and GenerateSample are vectorial version
!~* of GetPdf and Generate - usefull for non-idependent samples (eg AR1).
!~* BE CAREFULL WITH PARAMETERIZATION: in particular, standard deviation
!~* and scale are used in this library, rather than variance or squared
!~* scale. This allows checking positivity of standard deviation/scale
!~* parameters. DMSL uses the alternative parameterization (ie uses
!~* variances and squared scales)
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: GetMoments? GetCharacteristicFunction? GetSupport? GetQuantile?
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1. GetParNumber(DistID, npar, err, mess)
!~*    2. GetParName(DistID, name, err, mess)
!~*    3. GetParFeas(DistID, par, feas, err, mess)
!~*    4. GetPdf(DistId,x,par,loga,cov,pdf,feas,isnull,err,mess)
!~*    5. GetCdf(DistId,x,par,cov,cdf,feas,err,mess)
!~*    6. GetQuantile(DistId,p,par,cov,q,feas,err,mess)
!~*    7. GetMean(DistId,par,cov,m,feas,isfinite,err,mess)
!~*    8. GetStdev(DistId,par,cov,s,feas,isfinite,err,mess)
!~*    9. Generate(DistId,par,cov,gen,feas,err,mess)
!~*    10. GetSamplePdf(DistId,x,par,loga,cov,pdf,feas,isnull,err,mess)
!~*    11. GenerateSample(DistId,par,cov,gen,feas,err,mess)
!~*     12. IsCovNeeded(DistID) - logical function
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: GetParNumber, GetParName, GetParFeas, GetPdf, GetCdf, GetQuantile,&
          Generate, GetSamplePdf, GenerateSample, IsCovNeeded,GetRoughEstimate,&
          GetMode, GetMean, GetStdev
Character(100), parameter, PUBLIC:: & ! Shortcut for distribution names
            ! Continuous distributions
            GAUSS='Gaussian', GAUSSIAN='Gaussian', NORMAL='Gaussian',&
            Gauss_reparameterized='Gaussian_reparameterized',Gaussian_reparameterized='Gaussian_reparameterized',&
            UNIF='Uniform', UNIFORM='Uniform', &
            InvChi2='Inverse_Chi2', &
            GEORGE='George distribution', & ! N(0,sig^2+((a+b*cov))^2). sum of 2 errors
            VAGUEPRIOR='FlatPrior', FLATPRIOR='FlatPrior',& ! flat prior (=1)
            AR1='Gaussian AR(1)',& !stationary AR(1) gaussian model - (Xt-mu)=ro(X(t-1)-mu)+et, et~N(0,sigma^2)
            GEV='GEV',& !GEV(location,scale,shape)
            GEV_reparameterized='GEV_reparameterized',& !reparameterization for regionalisation purposes: GEV(location,lambda=(scale/loc),shape)
            GEV_min='GEV_min',& ! GEV for minima
            GEV_min_pos='GEV_min_pos',& ! GEV for minima constrained to be positive
            Gumbel='Gumbel',GUM='Gumbel',&
            GUM_reparameterized='Gumbel_reparameterized',Gumbel_reparameterized='Gumbel_reparameterized',&!reparameterization for regionalisation purposes: GUM(location,lambda=(scale/loc))
            Gumbel_min='Gumbel_min',GUM_min='Gumbel_min',& ! Gumbel for minima
            EXPO='Exponential',&!Exp(location,scale)
            GPD='GPD',&!GPD(location,scale,shape)
            GEV_Histo='GEV_Histo',& !GEV with historical info as covariates.
                                  ! cov(1)=type of data
                                  !         = 0 if exact data
                                  !         = 1 if data "smaller than" value in cov(2)
                                  !         = 2 if data "larger than" value in cov(2)
                                  !         = 3 if data "between" values in cov(2) and cov(3)
            TRIANGLE='Triangle', TRIANGULAR='Triangle',&
            FlatPriorPositive='FlatPrior+',&
            FlatPriorNegative='FlatPrior-',&
            PearsonIII='PearsonIII',&
            PearsonIII_reparameterized='PearsonIII_reparameterized',&
            LogN='LogNormal',LogNormal='LogNormal',&
            LogN3='LogNormal3',LogNormal3='LogNormal3',&
            Beta='Beta',&
            Kumaraswamy='Kumaraswamy',Kuma='Kumaraswamy',&
            MSPrior='MSPrior',& ! Martins & Stedinger's prior for GEV shape (some beta dist. translated between -0.5 and 0.5)
            ! Discrete distributions
            Bernoulli='Bernoulli',&
            GEOM='Geometric',&
            Poisson='Poisson',&
            Binomial='Binomial',&
            NegBinomial='NegBinomial',&
            NegBinomial2='NegBinomial2'

! type needed for allowing 1-D or 2-D covariates in subs GetSamplePdf and GenerateSample
! this is basically a 'vector of vectors'
type, public::covariate
    real(mrk), pointer::val(:) =>NULL()
end type covariate

Contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetParNumber(DistID, npar, err, mess)

!^**********************************************************************
!^* Purpose: returns the number of parameters of a distribution
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:30/05/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution ID
!^* OUT
!^*    1.nPar, the numbr of parameters
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************

character(*),intent(in)::DistID
integer(mik),intent(out)::nPar
integer(mik), intent(out)::err
character(*),intent(out)::mess

err=0;mess='';nPar=UndefIN
select case(DistID)
case(FLATPRIOR,FlatPriorPositive,FlatPriorNegative,MSPrior)
    npar=0
case(Bernoulli, GEOM, Poisson)
    npar=1
case(GAUSS, UNIF, InvChi2, Gumbel, Binomial, NegBinomial, NegBinomial2,&
     GUM_reparameterized,Gauss_reparameterized,LogN,GUM_min,GEV_min_pos,&
     EXPO,Beta,Kuma)
    npar=2
case(AR1, GEORGE, GEV,GEV_Histo,GPD,TRIANGLE,GEV_reparameterized,PearsonIII,&
    PearsonIII_reparameterized,GEV_min,LogN3)
    npar=3
case default
    err=1;mess='GetParNumber:Fatal:Unavailable Dist'
end select

end subroutine GetParNumber

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetParName(DistID, name, err, mess)

!^**********************************************************************
!^* Purpose: returns the names of parameters of a distribution
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:30/05/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution ID
!^* OUT
!^*    1.name, parameters names
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************

character(*),intent(in)::DistID
character(*),intent(out)::name(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::npar

err=0;mess='';name=undefCH

! check size
call GetParNumber(DistID, npar, err, mess)
if(err>0) then
    mess='GetParName: '//trim(mess);return
endif
if(size(name)/=npar) then
    err=2;mess='GetParName: dimension mismatch';return
endif

select case(DistID)
case(GAUSS)
    name(1)='mean'
    name(2)='standard_deviation'
case(Gauss_reparameterized)
    name(1)='mean'
    name(2)='cv'
case(UNIF)
    name(1)='lower_bound'
    name(2)='higher_bound'
case(InvChi2)
    name(1)='dof'
    name(2)='scale'
case(GEORGE)
    name(1)='remanent_std'
    name(2)='measurement_std_intercept'
    name(3)='measurement_std_slope'
case(FLATPRIOR,FlatPriorPositive,FlatPriorNegative)
    !name(1)='starting_point'
    !name(2)='starting_std'
case(AR1)
    name(1)='mean'
    name(2)='standard_deviation'
    name(3)='lag1_autocorrelation'
case(GEV,GEV_Histo,GEV_min)
    name(1)='location'
    name(2)='scale'
    name(3)='shape'
case(GEV_reparameterized)
    name(1)='location'
    name(2)='scale-to-location ratio'
    name(3)='shape'
case(GPD)
    name(1)='threshold'
    name(2)='scale'
    name(3)='shape'
case(Gumbel,Gumbel_min,GEV_min_pos)
    name(1)='location'
    name(2)='scale'
case(Gumbel_reparameterized)
    name(1)='location'
    name(2)='scale-to-location ratio'
case(EXPO)
    name(1)='threshold'
    name(2)='scale'
case(TRIANGLE)
    name(1)='peak'
    name(2)='lower_bound'
    name(3)='higher_bound'
case(Bernoulli)
    name(1)='success_prob'
case(GEOM)
    name(1)='success_prob'
case(Poisson)
    name(1)='lambda'
case(Binomial)
    name(1)='success_prob'
    name(2)='n'
case(NegBinomial)
    name(1)='success_prob'
    name(2)='pre-defined_fail_nb'
case(NegBinomial2)
    name(1)='mean'
    name(2)='pre-defined_fail_nb'
case(PearsonIII)
    name(1)='location'
    name(2)='scale'
    name(3)='shape'
case(PearsonIII_reparameterized)
    name(1)='location'
    name(2)='scale-to-location ratio'
    name(3)='shape'
case(LogN)
    name(1)='mean_log'
    name(2)='standard_deviation_log'
case(LogN3)
    name(1)='threshold'
    name(2)='mean_log_excess'
    name(3)='standard_deviation_log_excess'
case(Beta,Kuma)
    name(1)='shape1'
    name(2)='shape2'
case default
    err=1;mess='GetParName:Fatal:Unavailable Dist'
end select

end subroutine GetParName

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetParFeas(DistID, par, feas, err, mess)

!^**********************************************************************
!^* Purpose: check parameters feasability
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:30/05/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution ID
!^*    2.par, parameter vector
!^* OUT
!^*    1.feas, feasability
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************

character(*),intent(in)::DistID
real(mrk), intent(in)::par(:)
logical,intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! Locals
logical::ok

err=0;mess='';feas=.true.

! check size
call CheckParSize(DistID, par, ok, err, mess)
if(.not.ok) then
    err=2;mess='GetParFeas: dimension mismatch';return
endif

select case(DistID)
case(GAUSS,LogN)
    if (par(2)<=0.0_mrk) feas=.false.
case(LogN3)
    if (par(3)<=0.0_mrk) feas=.false.
case(Gauss_reparameterized)
    if(par(1)==0._mrk .or. par(1)*par(2)<=0._mrk) feas=.false.
case(UNIF)
    if (par(2)<=par(1)) feas=.false.
case(InvChi2)
    if (par(1)<=0.0_mrk .OR. par(2)<=0.0_mrk) feas=.false.
case(GEORGE)
    if (par(1)<=0.0_mrk) feas=.false.
case(FLATPRIOR,FlatPriorPositive,FlatPriorNegative,MSPrior)
    !if (par(2)<=0.0_mrk) feas=.false.
case(AR1)
    if (par(2)<=0.0_mrk .OR. abs(par(3))>=1.0_mrk) feas=.false.
case(GEV,Gumbel,GEV_Histo, GPD,GEV_min,Gumbel_min,EXPO)
    if (par(2)<=0.0_mrk) feas=.false.
case(GEV_min_pos)
    if (par(1)<=0.0_mrk .or. par(2)<=0.0_mrk) feas=.false.
case(GEV_reparameterized,Gumbel_reparameterized)
    if (par(2)*par(1) <= 0.0_mrk) feas=.false.
case(TRIANGLE)
    if (par(3)<=par(2) .or. par(1)<=par(2) .or. par(1)>=par(3)) feas=.false.
case(Bernoulli)
    if (par(1)<0 .or. par(1)>1) feas=.false.
case(GEOM)
    if (par(1)<0 .or. par(1)>1) feas=.false.
case(Poisson)
    if (par(1)<=0) feas=.false.
case(Binomial)
    if (par(1)<0 .or. par(1)>1) feas=.false.
    if (par(2)<=0) feas=.false.
case(NegBinomial)
    if (par(1)<0 .or. par(1)>1) feas=.false.
    if (par(2)<=0) feas=.false.
case(NegBinomial2)
    !if (par(1)<0 .or. par(1)>1) feas=.false.
    if (par(2)<=0 .or. par(1)<=0) feas=.false.
case(PearsonIII)
    if (par(2)==0._mrk .or. par(3)<=0._mrk) feas=.false.
case(PearsonIII_reparameterized)
    if(par(1)==0._mrk .or. par(2)==0._mrk .or. par(3)<=0._mrk) feas=.false.
case(Beta,Kuma)
    if (par(1)<=0 .or. par(2)<=0) feas=.false.
case default
    err=1;mess='GetParFeas:Fatal:Unavailable Dist'
end select

end subroutine GetParFeas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GetRoughEstimate(DistID, X, par, err, mess)

!^**********************************************************************
!^* Purpose: Get a rough estimation of parameters - mostly used for
!^* automatically setting the starting point in MCMC
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:25/01/2012
!^**********************************************************************
!^* Comments: data assumed without any missing value
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution ID
!^*    2.X, data
!^* OUT
!^*    1.par, parameter rough estimate
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalStats
character(*),intent(in)::DistID
real(mrk), intent(in)::X(:)
real(mrk), intent(out)::par(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! Locals
real(mrk),parameter::eps=0.000001_mrk,xi0=0.001_mrk ! 0.1_mrk
logical::ok
real(mrk)::mean,std,scale,skew,mini

err=0;mess='';

! check size
call CheckParSize(DistID, par, ok, err, mess)
if(.not.ok) then
    err=2;mess='GetRoughEstimate: dimension mismatch';return
endif

select case(DistID)
case(FlatPrior,FlatPriorPositive,FlatPriorNegative)
    ! do nothing
case(GAUSS)
    call GetEmpiricalStats(x=X,mean=par(1),std=par(2),err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
case(GAUSS_reparameterized)
    call GetEmpiricalStats(x=X,mean=par(1),std=scale,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    if(par(1)==0._mrk) then
        err=1;mess='GetRoughEstimate:Fatal:zero-mean';return
    else
        par(2)=scale/par(1)
    endif
case(UNIF)
    call GetEmpiricalStats(x=X,mini=par(1),maxi=par(2),err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
case(GEV)
    call GetEmpiricalStats(x=X,mean=mean,std=std,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(2)=0.7797_mrk*std
    par(1)=mean-0.5772_mrk*par(2)
    par(3)=xi0
case(GEV_min)
    call GetEmpiricalStats(x=-1._mrk*X,mean=mean,std=std,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(2)=0.7797_mrk*std
    par(1)=-1._mrk*(mean-0.5772_mrk*par(2))
    par(3)=xi0
case(GEV_min_pos)
    call GetEmpiricalStats(x=-1._mrk*X,mean=mean,std=std,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(2)=0.7797_mrk*std
    par(1)=-1._mrk*(mean-0.5772_mrk*par(2))
case(GEV_reparameterized)
    call GetEmpiricalStats(x=X,mean=mean,std=std,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    scale=0.7797_mrk*std
    par(1)=mean-0.5772_mrk*scale
    par(2)=scale/par(1)
    par(3)=xi0
case(Gumbel)
    call GetEmpiricalStats(x=X,mean=mean,std=std,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(2)=0.7797_mrk*std
    par(1)=mean-0.5772_mrk*par(2)
case(Gumbel_min)
    call GetEmpiricalStats(x=-1._mrk*X,mean=mean,std=std,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(2)=0.7797_mrk*std
    par(1)=-1._mrk*(mean-0.5772_mrk*par(2))
case(Gumbel_reparameterized)
    call GetEmpiricalStats(x=X,mean=mean,std=std,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    scale=0.7797_mrk*std
    par(1)=mean-0.5772_mrk*scale
    par(2)=scale/par(1)
case(EXPO)
    call GetEmpiricalStats(x=X,mean=mean,mini=par(1),err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(2)=mean-par(1)
case(GPD)
    call GetEmpiricalStats(x=X,mean=mean,mini=par(1),err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(2)=mean-par(1)
    par(3)=xi0
case(TRIANGLE)
    call GetEmpiricalStats(x=X,mean=par(1),mini=par(2),maxi=par(3),err=err,mess=mess)
    par(2)=par(2)-eps !*abs(par(2))
    par(3)=par(3)+eps !*abs(par(3))
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
case(PearsonIII)
    call GetEmpiricalStats(x=X,mini=mini,std=std,skewness=skew,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(1)=mini-0.01_mrk*std
    if(skew==0._mrk) then
        par(3)=4._mrk
    else
        par(3)=4._mrk/(skew**2)
    endif
    par(2)=std/sqrt(par(3))
case(PearsonIII_reparameterized)
    call GetEmpiricalStats(x=X,mini=mini,std=std,skewness=skew,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(1)=mini-0.01_mrk*std
    if(skew==0._mrk) then
        par(3)=4._mrk
    else
        par(3)=4._mrk/(skew**2)
    endif
    par(2)=(std/sqrt(par(3)))/par(1)
case(LogN)
    call GetEmpiricalStats(x=log(X),mean=par(1),std=par(2),err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
case(LogN3)
    call GetEmpiricalStats(x=X,mini=mini,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
    par(1)=mini-eps
    call GetEmpiricalStats(x=log(X-par(1)),mean=par(2),std=par(3),err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
case(Poisson)
    call GetEmpiricalStats(x=X,mean=par(1),std=scale,err=err,mess=mess)
    if(err>0) then
        mess="GetRoughEstimate: "//trim(mess);return
    endif
case(GEORGE,NegBinomial2,NegBinomial,Binomial,GEOM,Bernoulli,AR1,GEV_Histo,InvChi2)
    err=3;mess='GetRoughEstimate:Fatal:Not yet coded for this Dist';return
case(Beta,Kuma) ! corresponds to a U[0,1]. 2DO: could be improved
    par(1)=1._mrk
    par(2)=1._mrk
case default
    err=1;mess='GetRoughEstimate:Fatal:Unavailable Dist';return
end select
end subroutine GetRoughEstimate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetPdf(DistId,x,par,loga,cov,pdf,feas,isnull,err,mess)

!^**********************************************************************
!^* Purpose: compute pdf(x|par) for DistID
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:30/05/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    2.x, value
!^*    3.par, parameters
!^*    4.loga, log-pdf or natural pdf?
!^*    5.[cov], covariate
!^* OUT
!^*    1.pdf, result
!^*    2.feas, is par feasible?
!^*    3.isnull, pdf==0? (usefull for log-pdf of bounded-support distribution)
!^*    4.err, error code
!^*    5.mess, error message
!^**********************************************************************

use numerix_dmsl_kit, only:normal_logp, invChi2_logp,poisson_logPmf,binomial_pmf,&
                      negBinomial_pmf,beta_logP,exp_logp

character(*), intent(in)::DistID
real(mrk), intent(in)::x, par(:)
logical, intent(in)::loga
real(mrk), intent(in), optional::cov(:)
real(mrk), intent(out)::pdf
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!Locals
real(mrk)::measSTD,cdf, cdf1,cdf2,bet,pTemp
integer(mik)::xi

!Init
pdf=UndefRN;feas=.true.;isnull=.false.;err=0;mess=''
xi=real2int(x)

!Feasability
call GetParFeas(DistID, par, feas, err, mess)
if(err>0) then
    mess='GetPdf: '//trim(mess);return
endif
if(.not.feas) return

!Compute
select case(DistID)
case(GAUSS)
    pdf=normal_logp(x=x, mean=par(1), var=(par(2))**2)
    if(.not.loga) pdf=exp(pdf)
case(GAUSS_reparameterized)
    pdf=normal_logp(x=x, mean=par(1), var=(par(1)*par(2))**2)
    if(.not.loga) pdf=exp(pdf)
case(UNIF)
    call Uni_logP(x=x, a=par(1), b=par(2), logp=pdf, feas=feas, isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(InvChi2)
    if(x<=0.0_mrk) then
        isnull=.true.
     else
        pdf=invChi2_logp(x=x, v=par(1), s2=(par(2))**2)
    endif
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(GEORGE)
    if(.not.present(cov)) then
        err=10;mess='GetPdf:Fatal:argument [cov] missing'
        return
    endif
    measSTD=par(2)+par(3)*cov(1)
    if(measSTD<=0.0_mrk) then
        feas=.false.;return
    endif
    pdf=normal_logp(x=x, mean=0.0_mrk, var=(par(1))**2+(measSTD)**2)
    if(.not.loga) pdf=exp(pdf)
case(FLATPRIOR)
    err=-1;mess='GetPdf:Warning: improper [FLATPRIOR]';pdf=0.0_mrk
    if(.not.loga) pdf=1.0_mrk
case(FlatPriorPositive)
    err=-1;mess='GetPdf:Warning: improper [FLATPRIOR]';
    if(x>=0._mrk) then
        if(loga) then
            pdf=0.0_mrk
        else
            pdf=1._mrk
        endif
    else
        isnull=.true.
        if(loga) then
            pdf=-HugeRe
        else
            pdf=0._mrk
        endif
    endif
case(FlatPriorNegative)
    err=-1;mess='GetPdf:Warning: improper [FLATPRIOR]';
    if(x<=0._mrk) then
        if(loga) then
            pdf=0.0_mrk
        else
            pdf=1._mrk
        endif
    else
        isnull=.true.
        if(loga) then
            pdf=-HugeRe
        else
            pdf=0._mrk
        endif
    endif
case(AR1) ! Marginal distribution is used - GetSamplePdf should be preferred for autoregressive models!
    err=-1;mess='GetPdf:Warning: sub [GetPdf] inadequate for [AR1] - use [GetSamplePdf]';
    pdf=normal_logp(x=x, mean=par(1), var=(par(2))**2/(1.0_mrk-par(3)**2))
    if(.not.loga) pdf=exp(pdf)
case(GEV)
    call GEV_logP(x=x,loc=par(1),scale=par(2),shape=par(3),logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(GEV_min)
    call GEV_logP(x=-1._mrk*x,loc=-1._mrk*par(1),scale=par(2),shape=par(3),logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(GEV_min_pos)
    call GEV_logP(x=-1._mrk*x,loc=-1._mrk*par(1),scale=par(2),shape=par(2)/par(1),logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(GEV_reparameterized)
    call GEV_logP(x=x,loc=par(1),scale=par(2)*par(1),shape=par(3),logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(EXPO)
    if(x<par(1)) isnull=.true.
    pdf=exp_logp(x=x-par(1), beta=1._mrk/par(2))
    if(.not.loga) then
        if(isnull) then
           pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(GPD)
    call GPD_logP(x=x,loc=par(1),scale=par(2),shape=par(3),logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(GEV_Histo)
    if(.not.present(cov)) then
        err=-10;mess='GetPdf:WARNING:argument [cov] missing - standard GEV computed'
        call GEV_logP(x=x,loc=par(1),scale=par(2),shape=par(3),logp=pdf,feas=feas,isnull=isnull)
    else
        select case(int(cov(1)))
        case(0)
            call GEV_logP(x=x,loc=par(1),scale=par(2),shape=par(3),logp=pdf,feas=feas,isnull=isnull)
        case(1)
            call GEV_cdf(x=cov(2),loc=par(1),scale=par(2),shape=par(3),cdf=cdf,feas=feas)
            if(cdf==0.) then
                isnull=.true.
                pdf=undefRN
            else
                pdf=log(cdf)
            endif
        case(2)
            call GEV_cdf(x=cov(2),loc=par(1),scale=par(2),shape=par(3),cdf=cdf,feas=feas)
            if(cdf==1.) then
                isnull=.true.
                pdf=undefRN
            else
                pdf=log(1.-cdf)
            endif
        case(3)
            call GEV_cdf(x=cov(2),loc=par(1),scale=par(2),shape=par(3),cdf=cdf1,feas=feas)
            call GEV_cdf(x=cov(3),loc=par(1),scale=par(2),shape=par(3),cdf=cdf2,feas=feas)
            if(cdf2<=cdf1) then
                isnull=.true.
                pdf=undefRN
            else
                pdf=log(cdf2-cdf1)
            endif
        endselect
    endif
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(Gumbel)
    call GEV_logP(x=x,loc=par(1),scale=par(2),shape=0._mrk,logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(Gumbel_min)
    call GEV_logP(x=-1._mrk*x,loc=-1._mrk*par(1),scale=par(2),shape=0._mrk,logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(Gumbel_reparameterized)
    call GEV_logP(x=x,loc=par(1),scale=par(1)*par(2),shape=0._mrk,logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(TRIANGLE)
    call Triangle_logP(x=x, peak=par(1),a=par(2), b=par(3), logp=pdf, feas=feas, isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(Bernoulli)
    call Ber_pdf(x=xi,p=par(1),pdf=pdf,isnull=isnull)
    if(loga) then
        if(isnull) then
            pdf=-HugeRe
        else
            pdf=log(pdf)
        endif
    endif
case(GEOM)
    ! Probability of total number of try until success
    call Geo_pdf(x=xi,p=par(1),pdf=pdf,isnull=isnull)
    if(loga) then
        if(isnull) then
            pdf=-HugeRe
        else
            pdf=log(pdf)
        endif
    endif
case(Poisson)
    if (xi<0_mik) then
        isnull=.true.
    else
        pdf=poisson_logPmf(theta=xi,rate=par(1))
    endif
    if(.not. loga) then
        if(isnull) then
            pdf=0._mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(Binomial)
    if (xi<0_mik) isnull=.true.
    if (isnull) then
        pdf=0.0_mrk
    else
        pdf=binomial_pmf(theta=xi,n=real2int(par(2)),p=par(1))
    endif
    if(loga) then
        if(isnull) then
            pdf=-HugeRe
        else
            pdf=log(pdf)
        endif
    endif
case(NegBinomial)
    ! Probability of total number of try until par(2) failures occur
    if (xi<0_mik) isnull=.true.
    if (isnull) then
        pdf=0.0_mrk
    else
        bet=1._mrk/par(1)-1._mrk
        pdf=negBinomial_pmf(theta=xi,alpha=real2int(par(2)),beta=bet)
    endif
    if(loga) then
        if(isnull) then
            pdf=-HugeRe
        else
            pdf=log(pdf)
        endif
    endif
case(NegBinomial2)
    ! Probability of total number of try until par(2) failures occur
    if (xi<0_mik) isnull=.true.
    if (isnull) then
        pdf=0.0_mrk
    else
        pTemp=par(1)/(par(1)+par(2))
        bet=1._mrk/pTemp-1._mrk
        pdf=negBinomial_pmf(theta=xi,alpha=real2int(par(2)),beta=bet)
    endif
    if(loga) then
        if(isnull) then
            pdf=-HugeRe
        else
            pdf=log(pdf)
        endif
    endif
case(PearsonIII)
    call PearsonIII_logP(x=x,loc=par(1),scale=par(2),shape=par(3),logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(PearsonIII_reparameterized)
    call PearsonIII_logP(x=x,loc=par(1),scale=par(2)*par(1),shape=par(3),logp=pdf,feas=feas,isnull=isnull)
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(LogN)
    if(x<=0._mrk) then
        isnull=.true.
    else
        pdf=normal_logp(x=log(x), mean=par(1), var=(par(2))**2)-log(x)
    endif
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(LogN3)
    if(x-par(1)<=0._mrk) then
        isnull=.true.
    else
        pdf=normal_logp(x=log(x-par(1)), mean=par(2), var=(par(3))**2)-log(x-par(1))
    endif
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(Beta)
    if(x<=0.0_mrk .or. x>=1.0_mrk) then
        isnull=.true.
     else
        pdf=beta_logP(x=x, alpha=par(1), beta=par(2))
    endif
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(Kuma)
    if(x<=0.0_mrk .or. x>=1.0_mrk) then
        isnull=.true.
     else
        pdf=log(par(1)) + log(par(2)) + (par(1)-1._mrk)*log(x) + (par(2)-1._mrk)*log(1._mrk-x**par(1))
    endif
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case(MSprior)
    if(x<=-0.5_mrk .or. x>=0.5_mrk) then
        isnull=.true.
    else
        pdf=beta_logP(x=x+0.5_mrk,alpha=6._mrk,beta=9._mrk)
    endif
    if(.not.loga) then
        if(isnull) then
            pdf=0.0_mrk
        else
            pdf=exp(pdf)
        endif
    endif
case default
    err=1;mess='GetPdf:Fatal:Unavailable Dist'
end select

end subroutine GetPdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetCdf(DistId,x,par,cov,cdf,feas,err,mess)

!^**********************************************************************
!^* Purpose: compute Cdf(x|par) for DistID
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:30/05/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    2.x, value
!^*    3.par, parameters
!^*    4.[cov], covariates
!^* OUT
!^*    1.cdf, result
!^*    2.feas, is par feasible?
!^*    3.err, error code
!^*    4.mess, error message
!^**********************************************************************

use numerix_dmsl_kit, only:normal_cdf,poisson_cdf,binomial_cdf,exp_cdf
use utilities_dmsl_kit, only:betai

character(*), intent(in)::DistID
real(mrk), intent(in)::x, par(:)
real(mrk), intent(in), optional::cov(:)
real(mrk), intent(out)::cdf
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
real(mrk)::measSTD,tempcdf
integer(mik)::xi,i
logical::isnull

!Init
cdf=UndefRN;feas=.false.;err=0;mess=''
xi=real2int(x)

!Feasability
call GetParFeas(DistID, par, feas, err, mess)
if(err>0) then
    mess='GetCdf: '//trim(mess);return
endif
if(.not.feas) return

!Compute
select case(DistID)
case(GAUSS)
    cdf=normal_cdf(x=x, mean=par(1), sdev=par(2))
case(GAUSS_reparameterized)
    cdf=normal_cdf(x=x, mean=par(1), sdev=par(2)*par(1))
case(UNIF)
    call Uni_cdf(x=x, a=par(1), b=par(2), cdf=cdf, feas=feas)
case(InvChi2)
    call InvChi2_cdf(x=x, DoF=par(1), scale=par(2), cdf=cdf, feas=feas)
case(GEORGE)
    if(.not.present(cov)) then
        err=10;mess='GetCdf:Fatal:argument [cov] missing'
        return
    endif
    measSTD=par(2)+par(3)*cov(1)
    if(measSTD<=0.0_mrk) then
        feas=.false.;return
    endif
    cdf=normal_cdf(x=x, mean=0.0_mrk, sdev=sqrt((par(1))**2+(measSTD)**2))
case(FLATPRIOR)
    err=20;mess='GetCdf:Fatal: No CDF for improper [FLATPRIOR]';return
case(FLATPRIORPositive)
    err=20;mess='GetCdf:Fatal: No CDF for improper [FLATPRIOR+]';return
case(FLATPRIORNegative)
    err=20;mess='GetCdf:Fatal: No CDF for improper [FLATPRIOR-]';return
case(AR1) ! Marginal distribution is used
    err=-1;mess='GetCdf:Warning: Marginal CDF computed for [AR1]';
    cdf=normal_cdf(x=x, mean=par(1), sdev=par(2)/sqrt(1.0_mrk-par(3)**2))
case(GEV,GEV_Histo)
    call GEV_cdf(x=x,loc=par(1),scale=par(2),shape=par(3),cdf=cdf,feas=feas)
case(GEV_min)
    call GEV_cdf(x=-1._mrk*x,loc=-1._mrk*par(1),scale=par(2),shape=par(3),cdf=cdf,feas=feas)
    cdf=1._mrk-cdf
case(GEV_min_pos)
    call GEV_cdf(x=-1._mrk*x,loc=-1._mrk*par(1),scale=par(2),shape=par(2)/par(1),cdf=cdf,feas=feas)
    cdf=1._mrk-cdf
case(GEV_reparameterized)
    call GEV_cdf(x=x,loc=par(1),scale=par(1)*par(2),shape=par(3),cdf=cdf,feas=feas)
case(EXPO)
    if(x<par(1)) then
        cdf=0._mrk
    else
        cdf=exp_cdf(x=x-par(1),beta=1._mrk/par(2))
    endif
case(GPD)
    call GPD_cdf(x=x,loc=par(1),scale=par(2),shape=par(3),cdf=cdf,feas=feas)
case(Gumbel)
    call GEV_cdf(x=x,loc=par(1),scale=par(2),shape=0._mrk,cdf=cdf,feas=feas)
case(Gumbel_min)
    call GEV_cdf(x=-1._mrk,loc=-1._mrk*par(1),scale=par(2),shape=0._mrk,cdf=cdf,feas=feas)
    cdf=1._mrk-cdf
case(Gumbel_reparameterized)
    call GEV_cdf(x=x,loc=par(1),scale=par(1)*par(2),shape=0._mrk,cdf=cdf,feas=feas)
case(TRIANGLE)
    call Triangle_cdf(x=x, peak=par(1),a=par(2), b=par(3), cdf=cdf, feas=feas)
case(Bernoulli)
    if (x<0._mrk) cdf=0._mrk
    if (x>=1._mrk) cdf=1._mrk
    if (x>=0._mrk .and. x<1._mrk) cdf=1-par(1)
case(GEOM)
    cdf=1._mrk-(1._mrk-par(1))**xi
case(Poisson)
    cdf=poisson_cdf(x=xi,Ex=par(1))
case(Binomial)
    !cdf=binomial_cdf(x=xi,n=real2int(par(2)),p=par(1)) !not works, why?

    cdf=0._mrk
    do i=0,xi
        call GetPdf(DistId=Binomial,x=real(i,mrk),par=par,loga=.false.,pdf=tempcdf,&
                feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then
            mess='GetCdf: '//trim(mess);return
        endif
        cdf=cdf+tempcdf
    enddo
case(NegBinomial)
    cdf=0._mrk
    do i=0,xi
        call GetPdf(DistId=NegBinomial,x=real(i,mrk),par=par,loga=.false.,pdf=tempcdf,&
                feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then
            mess='GetCdf: '//trim(mess);return
        endif
        cdf=cdf+tempcdf
    enddo
case(NegBinomial2)
    cdf=0._mrk
    do i=0,xi
        call GetPdf(DistId=NegBinomial2,x=real(i,mrk),par=par,loga=.false.,pdf=tempcdf,&
                feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then
            mess='GetCdf: '//trim(mess);return
        endif
        cdf=cdf+tempcdf
    enddo
case(PearsonIII)
    call PearsonIII_cdf(x=x,loc=par(1),scale=par(2),shape=par(3),cdf=cdf,feas=feas)
case(PearsonIII_reparameterized)
    call PearsonIII_cdf(x=x,loc=par(1),scale=par(2)*par(1),shape=par(3),cdf=cdf,feas=feas)
case(LogN)
    if(x<=0._mrk) then
        cdf=0._mrk
    else
        cdf=normal_cdf(x=log(x), mean=par(1), sdev=par(2))
    endif
case(LogN3)
    if(x-par(1)<=0._mrk) then
        cdf=0._mrk
    else
        cdf=normal_cdf(x=log(x-par(1)), mean=par(2), sdev=par(3))
    endif
case(Beta)
    if(x<=0._mrk) then
        cdf=0._mrk
    elseif(x>=1._mrk) then
        cdf=1._mrk
    else
        cdf=betai(a=par(1),b=par(2),x=x)
    endif
case(Kuma)
    if(x<=0._mrk) then
        cdf=0._mrk
    elseif(x>=1._mrk) then
        cdf=1._mrk
    else
        cdf=1._mrk-(1._mrk-x**par(1))**par(2)
    endif
case default
    err=1;mess='GetCdf:Fatal:Unavailable Dist'
end select

end subroutine GetCdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GetQuantile(DistId,p,par,cov,q,feas,err,mess)

!^**********************************************************************
!^* Purpose: compute InvCdf(p|par) for DistID
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:08/11/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    2.p, nonexceedance probability (0 and 1 not allowed)
!^*    3.par, parameters
!^*    4.[cov], covariates
!^* OUT
!^*    1.q, result
!^*    2.feas, is par feasible?
!^*    3.err, error code
!^*    4.mess, error message
!^**********************************************************************

use numerix_dmsl_kit, only:normal_icdf,exp_icdf

character(*), intent(in)::DistID
real(mrk), intent(in)::p, par(:)
real(mrk), intent(in), optional::cov(:)
real(mrk), intent(out)::q
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals

!Init
q=UndefRN;feas=.false.;err=0;mess=''

!Feasability
if((p<=0._mrk).or.(p>=1._mrk)) then
    err=1;mess='GetQuantile:Fatal:invalid [p]';return
endif
call GetParFeas(DistID, par, feas, err, mess)
if(err>0) then
    mess='GetQuantile: '//trim(mess);return
endif
if(.not.feas) return

!Compute
select case(DistID)
case(UNIF)
    q=par(1)+p*(par(2)-par(1))
case(GAUSS)
    q=normal_icdf(p=p,mean=par(1), sdev=par(2))
case(GAUSS_reparameterized)
    q=normal_icdf(p=p,mean=par(1), sdev=par(2)*par(1))
case(GEV,GEV_Histo)
    call GEV_Quant(p=p, loc=par(1),scale=par(2),shape=par(3),q=q, feas=feas)
case(GEV_min)
    call GEV_Quant(p=1._mrk-p, loc=-1._mrk*par(1),scale=par(2),shape=par(3),q=q, feas=feas)
    q=-1._mrk*q
case(GEV_min_pos)
    call GEV_Quant(p=1._mrk-p, loc=-1._mrk*par(1),scale=par(2),shape=par(2)/par(1),q=q, feas=feas)
    q=-1._mrk*q
case(GEV_reparameterized)
    call GEV_Quant(p=p, loc=par(1),scale=par(2)*par(1),shape=par(3),q=q, feas=feas)
case(EXPO)
    q=par(1)+exp_icdf(p=p,beta=1._mrk/par(2))
case(GPD)
    call GPD_Quant(p=p, loc=par(1),scale=par(2),shape=par(3),q=q, feas=feas)
case(Gumbel)
    call GEV_Quant(p=p, loc=par(1),scale=par(2),shape=0._mrk,q=q, feas=feas)
case(Gumbel_min)
    call GEV_Quant(p=1._mrk-p, loc=-1._mrk*par(1),scale=par(2),shape=0._mrk,q=q, feas=feas)
    q=-1._mrk*q
case(Gumbel_reparameterized)
    call GEV_Quant(p=p, loc=par(1),scale=par(1)*par(2),shape=0._mrk,q=q, feas=feas)
case(TRIANGLE)
    call Triangle_Quant(p=p, peak=par(1),a=par(2),b=par(3),q=q,feas=feas)
case(FLATPRIOR)
    err=20;mess='GetQuantile:Fatal: No quant for improper [FLATPRIOR]';return
case(FLATPRIORPositive)
    err=20;mess='GetQuantile:Fatal: No quant for improper [FLATPRIOR+]';return
case(FLATPRIORNegative)
    err=20;mess='GetQuantile:Fatal: No quant for improper [FLATPRIOR-]';return
case(PearsonIII)
    call PearsonIII_Quant(p=p, loc=par(1),scale=par(2),shape=par(3),q=q,feas=feas)
case(PearsonIII_reparameterized)
    call PearsonIII_Quant(p=p, loc=par(1),scale=par(2)*par(1),shape=par(3),q=q,feas=feas)
case(LogN)
    q=normal_icdf(p=p,mean=par(1), sdev=par(2));q=exp(q)
case(LogN3)
    q=normal_icdf(p=p,mean=par(2), sdev=par(3));q=exp(q)+par(1)
case(Beta)
    call Beta_Quant(p=p,a=par(1),b=par(2),q=q,feas=feas)
case(Kuma)
    q=(1._mrk-(1._mrk-p)**(1._mrk/par(2)))**(1._mrk/par(1))
case(Binomial)
    q=real(Binomial_Quant(p=p,success_prob=par(1), n=nint(par(2))),mrk)
case default
    err=1;mess='GetQuantile:Fatal:Unavailable Dist'
end select

end subroutine GetQuantile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetMean(DistId,par,cov,m,feas,isfinite,err,mess)

!^**********************************************************************
!^* Purpose: compute the mean of DistID
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:31/05/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    2.par, parameters
!^*    3.[cov], covariates
!^* OUT
!^*    1.m, result
!^*    2.feas, is par feasible?
!^*    3.isfinite, finite mean?
!^*    4.err, error code
!^*    5.mess, error message
!^**********************************************************************

use numerix_dmsl_kit, only:normal_logp
use utilities_dmsl_kit, only:gammaln

character(*), intent(in)::DistID
real(mrk), intent(in)::par(:)
real(mrk), intent(in), optional::cov(:)
real(mrk), intent(out)::m
logical, intent(out)::feas, isfinite
integer(mik), intent(out)::err
character(*),intent(out)::mess

!Init
m=UndefRN;feas=.false.;isfinite=.true.;err=0;mess=''

!Feasability
call GetParFeas(DistID, par, feas, err, mess)
if(err>0) then
    mess='GetMean: '//trim(mess);return
endif
if(.not.feas) return

!Compute
select case(DistID)
case(GAUSS,GAUSS_reparameterized)
    m=par(1)
case(InvChi2)
    if(par(1)<=2.0_mrk) then
        isfinite=.false.
    else
        m=(par(1)*par(2)**2)/(par(1)-2.0_mrk)
    endif
case(UNIF)
    m=0.5_mrk*(par(1)+par(2))
case(EXPO)
    m=par(1)+par(2)
case(GEORGE)
    m=0.0_mrk
case(AR1)
    m=par(1)
case(TRIANGLE)
    err=1;mess='GetMean:Fatal:Not yet coded for dist=TRIANGLE'
case(GEV,Gumbel,GEV_Histo,GPD,GEV_reparameterized,Gumbel_reparameterized,GEV_min,Gumbel_min,GEV_min_pos)
    err=1;mess='GetMean:Fatal:Not yet coded for this distribution'
case(Bernoulli)
    m=par(1)
case(GEOM)
    m=1._mrk/par(1)
case(Poisson)
    m=par(1)
case(Binomial)
    m=par(1)*par(2)
case(NegBinomial)
    m=par(1)*par(2)/(1._mrk-par(1))
case(NegBinomial2)
    m=par(1)
case(FLATPRIOR)
    err=20;isfinite=.false.
    mess='GetMean:Fatal: No Mean for improper [FLATPRIOR]';return
case(FLATPRIORPositive)
    err=20;isfinite=.false.
    mess='GetMean:Fatal: No Mean for improper [FLATPRIOR+]';return
case(FLATPRIORNegative)
    err=20;isfinite=.false.
    mess='GetMean:Fatal: No Mean for improper [FLATPRIOR-]';return
case(PearsonIII,PearsonIII_reparameterized)
    err=1;mess='GetMean:Fatal:Not yet coded for dist=PearsonIII'
case(LogN)
    m=exp(par(1) + 0.5_mrk * (par(2)**2) )
case(LogN3)
    m=par(1) + exp(par(2) + 0.5_mrk * (par(3)**2) )
case(Beta)
    m=par(1)/(par(1)+par(2))
case(Kuma)
    m= exp(log(par(2)) + gammaln(1._mrk+1._mrk/par(1)) + gammaln(par(2)) - gammaln(1._mrk+1._mrk/par(1)+par(2)))
case default
    err=1;mess='GetMean:Fatal:Unavailable Dist'
end select

end subroutine GetMean

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetMode(DistId,par,cov,m,feas,err,mess)

!^**********************************************************************
!^* Purpose: compute the mode of DistID
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:15/07/2014
!^**********************************************************************
!^* Comments: When mode is not uniquely defined, one possible mode is returned
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    2.par, parameters
!^*    3.[cov], covariates
!^* OUT
!^*    1.m, result
!^*    2.feas, is par feasible?
!^*    3.err, error code
!^*    4.mess, error message
!^**********************************************************************

character(*), intent(in)::DistID
real(mrk), intent(in)::par(:)
real(mrk), intent(in), optional::cov(:)
real(mrk), intent(out)::m
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess

!Init
m=UndefRN;feas=.false.;err=0;mess=''

!Feasability
call GetParFeas(DistID, par, feas, err, mess)
if(err>0) then;mess='GetMode: '//trim(mess);return;endif
if(.not.feas) return

!Compute
select case(DistID)
case(GAUSS,GAUSS_reparameterized)
    m=par(1)
case(InvChi2)
    m=(par(1)*par(2)**2)/(par(1)+2.0_mrk)
case(UNIF)
    m=0.5_mrk*(par(1)+par(2))
case(GEORGE)
    m=0.0_mrk
case(AR1)
    m=par(1)
case(TRIANGLE)
    m=par(1)
case(GEV,GEV_Histo,GEV_reparameterized,GEV_min,GEV_min_pos)
    err=1;mess='GetMode:Fatal:Not Yet Coded for GEV'
case(Gumbel,Gumbel_reparameterized,Gumbel_min)
    m=par(1)
case(EXPO)
    m=par(1)
case(GPD)
    m=par(1)
case(Bernoulli)
    if(par(1)>=1._mrk-par(1)) then;m=1._mrk;else;m=0._mrk;endif
case(GEOM)
    m=1._mrk
case(Poisson)
    m=floor(par(1))
case(Binomial)
    err=1;mess='GetMode:Fatal:Not Yet Coded for Binomial'
case(NegBinomial,NegBinomial2)
    err=1;mess='GetMode:Fatal:Not Yet Coded for Negative Binomial'
case(FLATPRIOR)
    feas=.false.;mess='GetMode:Fatal: No mode for improper [FLATPRIOR]';return
case(FLATPRIORPositive)
    feas=.false.;mess='GetMode:Fatal: No mode for improper [FLATPRIOR+]';return
case(FLATPRIORNegative)
    feas=.false.;mess='GetMode:Fatal: No mode for improper [FLATPRIOR-]';return
case(PearsonIII,PearsonIII_reparameterized)
    err=1;mess='GetMode:Fatal:Not Yet Coded for Pearson III'
case(LogN)
    m=exp(par(1)-par(2)**2)
case(LogN3)
    m=par(1)+exp(par(2)-par(3)**2)
case(Beta)
    if(par(1)>1._mrk .and. par(2)>1._mrk) then
        m=(par(1)-1._mrk)/(par(1)+par(2)-2._mrk)
    elseif (par(1)==1._mrk .and. par(2)==1._mrk) then
        m=0.5_mrk
    else
        m=undefRN
    endif
case(Kuma)
   if (par(1)==1._mrk .and. par(2)==1._mrk) then
        m=0.5_mrk
   elseif(par(1)>=1._mrk .and. par(2)>=1._mrk) then
        m=( (par(1)-1._mrk)/(par(1)*par(2)-1._mrk) )**(1._mrk/par(1))
    else
        m=undefRN
    endif
case default
    err=1;mess='GetMode:Fatal:Unavailable Dist'
end select

end subroutine GetMode

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetStdev(DistId,par,cov,s,feas,isfinite,err,mess)

!^**********************************************************************
!^* Purpose: compute the Std dev of DistID
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:31/05/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    2.par, parameters
!^*    3.[cov], covariates
!^* OUT
!^*    1.s, result
!^*    2.feas, is par feasible?
!^*    3.isfinite, finite stdev?
!^*    4.err, error code
!^*    5.mess, error message
!^**********************************************************************

use numerix_dmsl_kit, only:normal_logp
use utilities_dmsl_kit, only:gammaln

character(*), intent(in)::DistID
real(mrk), intent(in)::par(:)
real(mrk), intent(in), optional::cov(:)
real(mrk), intent(out)::s
logical, intent(out)::feas, isfinite
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk)::m1,m2

!Init
s=UndefRN;feas=.false.;isfinite=.true.;err=0;mess=''

!Feasability
call GetParFeas(DistID, par, feas, err, mess)
if(err>0) then
    mess='GetStdev: '//trim(mess);return
endif
if(.not.feas) return

!Compute
select case(DistID)
case(GAUSS)
    s=par(2)
case(GAUSS_reparameterized)
    s=par(2)*par(1)
case(EXPO)
    s=par(2)
case(UNIF)
    s=(par(2)-par(1))/sqrt(12.0_mrk)
case(InvChi2)
    if(par(1)<=4.0_mrk) then
        isfinite=.false.
    else
        s=(sqrt(2.0_mrk)*par(1)*par(2)**2)/((par(1)-2.0_mrk)*sqrt(par(1)-4.0_mrk))
    endif
case(GEORGE)
    if(.not.present(cov)) then
        err=10;mess='GetStdev:Fatal:optional cov needed for dist GEORGE';return
    endif
    s=sqrt(par(1)**2+(par(2)+par(3)*cov(1))**2)
case(AR1) ! Marginal Std!!!!
    s=par(2)/sqrt(1.0_mrk-par(3)**2)
case(GEV,Gumbel,GEV_Histo,GPD,GEV_reparameterized,Gumbel_reparameterized,GEV_min,Gumbel_min,GEV_min_pos)
    err=1;mess='GetStdev:Fatal:Not yet coded for this distribution'
case(TRIANGLE)
    err=1;mess='GetStdev:Fatal:Not yet coded for dist=TRIANGLE'
case(Bernoulli)
    s=sqrt(par(1)*(1._mrk-par(1)))
case(GEOM)
    s=sqrt(1._mrk-par(1))/par(1)
case(Poisson)
    s=par(1)
case(Binomial)
    s=sqrt(par(1)*par(2)*(1._mrk-par(1)))
case(NegBinomial)
    s=sqrt(par(1)*par(2))/(1._mrk-par(1))
case(NegBinomial2)
    s=sqrt(par(1)+par(1)*par(1)/par(2))
case(FLATPRIOR)
    err=20;isfinite=.false.
    mess='GetStdev:Fatal: No Stdev for improper [FLATPRIOR]';return
case(FLATPRIORPositive)
    err=20;isfinite=.false.
    mess='GetStdev:Fatal: No Stdev for improper [FLATPRIOR+]';return
case(FLATPRIORNegative)
    err=20;isfinite=.false.
    mess='GetStdev:Fatal: No Stdev for improper [FLATPRIOR-]';return
case(PearsonIII,PearsonIII_reparameterized)
    err=1;mess='GetStdev:Fatal:Unavailable PearsonIII'
case(LogN)
    s=exp(par(1) + 0.5_mrk * (par(2)**2) )*sqrt(exp(par(2)**2)-1._mrk)
case(LogN3)
    s=exp(par(2) + 0.5_mrk * (par(3)**2) )*sqrt(exp(par(3)**2)-1._mrk)
case(Beta)
    s=sqrt((par(1)*par(2))/(par(1)+par(2)+1._mrk)) / (par(1)+par(2))
case(Kuma)
    m1= exp(log(par(2)) + gammaln(1._mrk+1._mrk/par(1)) + gammaln(par(2)) - gammaln(1._mrk+1._mrk/par(1)+par(2)))
    m2= exp(log(par(2)) + gammaln(1._mrk+2._mrk/par(1)) + gammaln(par(2)) - gammaln(1._mrk+2._mrk/par(1)+par(2)))
    s=sqrt(m2-m1**2)
case default
    err=1;mess='GetStdev:Fatal:Unavailable Dist'
end select

end subroutine GetStdev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Generate(DistId,par,cov,gen,feas,err,mess)

!^**********************************************************************
!^* Purpose: Generate a random number from DistID
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:04/07/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    2.par, parameters
!^*    3.[cov], covariate
!^* OUT
!^*    1.gen, result
!^*    2.feas, is par feasible?
!^*    3.err, error code
!^*    4.mess, error message
!^**********************************************************************

use numerix_dmsl_kit, only:normaldev, uniran, invChi2dev, poidev, &
                 bnldev, negbindev,expdev

character(*), intent(in)::DistID
real(mrk), intent(in)::par(:)
real(mrk), intent(in), optional::cov(:)
real(mrk), intent(out)::gen
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
real(mrk)::dev, measSTD, bet,pTemp,prob

!Init
gen=UndefRN;feas=.true.;err=0;mess=''

!Feasability
call GetParFeas(DistID, par, feas, err, mess)
if(err>0) then
    mess='Generate: '//trim(mess);return
endif
if(.not.feas) return

!Compute
select case(DistID)
case(GAUSS)
    call normaldev(mean=par(1),sdev=par(2),gdev=gen,err=err,message=mess)
    if(err>0) then
        mess='Generate: '//trim(mess);return
    endif
case(GAUSS_reparameterized)
    call normaldev(mean=par(1),sdev=par(2)*par(1),gdev=gen,err=err,message=mess)
    if(err>0) then
        mess='Generate: '//trim(mess);return
    endif
case(UNIF)
    call uniran(dev)
    gen=(par(2)-par(1))*dev+par(1)
case(InvChi2)
    !call invChi2dev(ndeg=par(1),scaler=par(2),invchi2=gen,err=err)
    !if(err>0) then
    !    mess='Generate:error using [invChi2dev]';return
    !endif
    err=30;mess='Generate:Fatal: [InvChi2] disabled for the moment to maintain compatibility with several DMSL versions';return
case(FLATPRIOR)
    err=20;mess='Generate:Fatal: No Generation for improper [FLATPRIOR]';return
case(FLATPRIORPositive)
    err=20;mess='Generate:Fatal: No Generation for improper [FLATPRIOR+]';return
case(FLATPRIORNegative)
    err=20;mess='Generate:Fatal: No Generation for improper [FLATPRIOR-]';return
case(GEORGE)
    if(.not.present(cov)) then
        err=10;mess='Generate:Fatal:argument [cov] missing'
        return
    endif
    measSTD=par(2)+par(3)*cov(1)
    if(measSTD<=0.0_mrk) then
        feas=.false.;return
    endif
    call normaldev(mean=0.0_mrk,sdev=sqrt((par(1))**2+(measSTD)**2),gdev=gen,err=err,message=mess)
    if(err>0) then
        mess='Generate: '//trim(mess);return
    endif
case(AR1) ! Sampling from the marginal - GenerateSample should be preferred for autoregressive models!
    call normaldev(mean=par(1),sdev=par(2)/sqrt(1.0_mrk-par(3)**2),gdev=gen,err=err,message=mess)
    if(err>0) then
        mess='Generate: '//trim(mess);return
    endif
    err=-1;mess='Generate:Warning: sub [Generate] inadequate for [AR1] - use [GenerateSample]';
case(GEV,GEV_Histo)
    call GEV_dev(loc=par(1),scale=par(2),shape=par(3),dev=gen,feas=feas)
case(GEV_min)
    call GEV_dev(loc=-1._mrk*par(1),scale=par(2),shape=par(3),dev=gen,feas=feas)
    gen=-1._mrk*gen
case(GEV_min_pos)
    call GEV_dev(loc=-1._mrk*par(1),scale=par(2),shape=par(2)/par(1),dev=gen,feas=feas)
    gen=-1._mrk*gen
case(GEV_reparameterized)
    call GEV_dev(loc=par(1),scale=par(2)*par(1),shape=par(3),dev=gen,feas=feas)
case(EXPO)
    call expdev(beta=par(2),harvest=gen,err=err)
    err=0 ! error seems to be uninitialized in DMSL - put it to zero, will generate -HugeRe if par(2)<=0
    gen=gen+par(1)
case(GPD)
    call GPD_dev(loc=par(1),scale=par(2),shape=par(3),dev=gen,feas=feas)
case(Gumbel)
    call GEV_dev(loc=par(1),scale=par(2),shape=0._mrk,dev=gen,feas=feas)
case(Gumbel_min)
    call GEV_dev(loc=-1._mrk*par(1),scale=par(2),shape=0._mrk,dev=gen,feas=feas)
    gen=-1._mrk*gen
case(Gumbel_reparameterized)
    call GEV_dev(loc=par(1),scale=par(2)*par(1),shape=0._mrk,dev=gen,feas=feas)
case(TRIANGLE)
    call uniran(prob)
    call Triangle_Quant(prob,par(1),par(2),par(3),gen,feas)
case(Bernoulli)
    gen=Ber_dev(par(1))
case(GEOM)
    gen=Geo_dev(par(1))
case(Poisson)
    gen=poidev(lambda=par(1))
case(Binomial)
    gen=bnldev(n=real2int(par(2)),pp=par(1))
case(NegBinomial)
    bet=1._mrk/par(1)-1._mrk
    call negbindev(alpha=par(2),beta=bet,nbdev=gen,err=err,message=mess)
case(NegBinomial2)
    pTemp=par(1)/(par(1)+par(2))
    bet=1._mrk/pTemp-1._mrk
    call negbindev(alpha=par(2),beta=bet,nbdev=gen,err=err,message=mess)
case(PearsonIII)
    call PearsonIII_dev(loc=par(1),scale=par(2),shape=par(3),dev=gen,feas=feas)
case(PearsonIII_reparameterized)
    call PearsonIII_dev(loc=par(1),scale=par(2)*par(1),shape=par(3),dev=gen,feas=feas)
case(LogN)
    call normaldev(mean=par(1),sdev=par(2),gdev=gen,err=err,message=mess)
    if(err>0) then
        mess='Generate: '//trim(mess);return
    endif
    gen=exp(gen)
case(LogN3)
    call normaldev(mean=par(2),sdev=par(3),gdev=gen,err=err,message=mess)
    if(err>0) then
        mess='Generate: '//trim(mess);return
    endif
    gen=par(1)+exp(gen)
case(Beta,Kuma)
    call uniran(prob)
    call GetQuantile(DistId=DistID,p=prob,par=par,q=gen,feas=feas,err=err,mess=mess)
    if(err>0) then;mess='Generate: '//trim(mess);return;endif
case default
    err=1;mess='Generate:Fatal:Unavailable Dist'
end select

end subroutine Generate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetSamplePdf(DistId,x,par,loga,cov,pdf,feas,isnull,err,mess)

!^**********************************************************************
!^* Purpose: compute pdf(x|par) of a sample x(:) for DistID
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:04/07/2008
!^**********************************************************************
!^* Comments: this is where derived type covariate is useful, since cov
!^* can be 1D (only one covariate) or 2D (several covariates)
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    2.x, sample values
!^*    3.par, parameters
!^*    4.loga, log-pdf or natural pdf?
!^*    5.[cov], covariate
!^* OUT
!^*    1.pdf, result
!^*    2.feas, is par feasible?
!^*    3.isnull, pdf==0? (usefull for log-pdf of bounded-support distribution)
!^*    4.err, error code
!^*    5.mess, error message
!^**********************************************************************

use numerix_dmsl_kit, only:normal_logp

character(*), intent(in)::DistID
real(mrk), intent(in)::x(:), par(:)
logical, intent(in)::loga
type(covariate), intent(in), optional::cov(:)
real(mrk), intent(out)::pdf
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::n,i,j,ncov
real(mrk), allocatable:: density(:), coco(:)

!Init
pdf=UndefRN;feas=.false.;isnull=.false.;err=0;mess=''
n=size(x)
if(allocated(density)) deallocate(density)
allocate(density(n))
if(present(cov)) then
If(allocated(coco)) deallocate(coco)
    ncov=size(cov)
    allocate(coco(ncov))
endif

!Feasability
call GetParFeas(DistID, par, feas, err, mess)
if(err>0) then
    mess='GetSamplePdf: '//trim(mess);return
endif
if(.not.feas) return

!Compute
select case(DistID)
case(GAUSS, UNIF,InvChi2,GEORGE,FLATPRIOR, GEV,Gumbel,GEV_Histo,GPD,TRIANGLE,&
        GEV_reparameterized,GUM_reparameterized,FlatPriorPositive,FlatPriorNegative,&
        GAUSS_reparameterized,PearsonIII,PearsonIII_reparameterized,LogN,LogN3,&
        GEV_min,Gumbel_min,GEV_min_pos,EXPO,Beta,Kuma)
    SampleLoop :do i=1,n
        if(present(cov)) then
            forall(j=1:ncov) coco(j)=cov(j)%val(i)
            call GetPdf(DistID=DistId,x=x(i),par=par,loga=loga,cov=coco,&
                    pdf=density(i),feas=feas,isnull=isnull,err=err,mess=mess)
        else
            call GetPdf(DistID=DistId,x=x(i),par=par,loga=loga,&
                    pdf=density(i),feas=feas,isnull=isnull,err=err,mess=mess)
        endif
        if(err>0) then
            mess='GetSamplePdf: '//trim(mess);return
        endif
        if(.not.feas) return
        if(isnull) exit SampleLoop
    enddo SampleLoop
    if(isnull) then
        if(loga) then
            pdf=-HugeRe
        else
            pdf=0.0_mrk
        endif
    else
        if(loga) then
            pdf=sum(density)
        else
            pdf=product(density)
        endif
    endif
    deallocate(density)
    if(present(cov)) deallocate(coco)
case(AR1) ! likelihood - x1 from the marginal distribution (x0 unknown)
    ! First value - from the marginal
    call GetPdf(DistID=GAUSS,x=x(1),par=(/par(1), par(2)/sqrt(1.0_mrk-par(3)**2)/),loga=loga,&
                    pdf=density(1),feas=feas,isnull=isnull,err=err,mess=mess)
    if(err>0) then
        mess='GetSamplePdf: '//trim(mess);return
    endif
    if(.not.feas) return
    if(isnull) then
        if(loga) then
            pdf=-HugeRe
        else
            pdf=0.0_mrk
        endif
        return
    endif
    do i=2,n ! Subsequent values - conditioned on t-1
        call GetPdf(DistID=GAUSS,x=x(i),par=(/(1.0_mrk-par(3))*par(1)+par(3)*x(i-1), par(2)/),loga=loga,&
                    pdf=density(i),feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then
            mess='GetSamplePdf: '//trim(mess);return
        endif
        if(.not.feas) return
        if(isnull) then
            if(loga) then
                pdf=-HugeRe
            else
                pdf=0.0_mrk
            endif
            return
        endif
    enddo
    if(loga) then
        pdf=sum(density)
    else
        pdf=product(density)
    endif
case default
    err=1;mess='GetSamplePdf:Fatal:Unavailable Dist'
end select

end subroutine GetSamplePdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GenerateSample(DistId,par,cov,gen,feas,err,mess)

!^**********************************************************************
!^* Purpose: Generate a random sample from DistID
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:04/07/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    2.par, parameters
!^*    3.[cov], covariate
!^* OUT
!^*    1.gen, result
!^*    2.feas, is par feasible?
!^*    3.err, error code
!^*    4.mess, error message
!^**********************************************************************

use numerix_dmsl_kit, only:normal_logp

character(*), intent(in)::DistID
real(mrk), intent(in):: par(:)
type(covariate), intent(in), optional::cov(:)
real(mrk), intent(out)::gen(:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::n,i,j,ncov
real(mrk), allocatable:: coco(:)

!Init
gen=UndefRN;feas=.false.;err=0;mess=''
n=size(gen)
if(present(cov)) then
    if(allocated(coco)) deallocate(coco)
    ncov=size(cov)
    allocate(coco(ncov))
endif

!Feasability
call GetParFeas(DistID, par, feas, err, mess)
if(err>0) then
    mess='GenerateSample: '//trim(mess);return
endif
if(.not.feas) return

!Compute
select case(DistID)
case(GAUSS, UNIF,InvChi2,GEORGE,FLATPRIOR, GEV,Gumbel,GEV_Histo,GPD,TRIANGLE,&
        GEV_reparameterized,GUM_reparameterized,FlatPriorPositive,FlatPriorNegative,&
        GAUSS_reparameterized,PearsonIII,PearsonIII_reparameterized,&
        Bernoulli,GEOM,Poisson,Binomial,NegBinomial,NegBinomial2,LogN,LogN3,&
        GEV_min,Gumbel_min,GEV_min_pos,EXPO,Beta,Kuma)
    SampleLoop :do i=1,n
        if(present(cov)) then
            forall(j=1:ncov) coco(j)=cov(j)%val(i)
            call Generate(DistID=DistId,par=par,cov=coco,&
                    gen=gen(i),feas=feas,err=err,mess=mess)
        else
            call Generate(DistID=DistId,par=par,&
                    gen=gen(i),feas=feas,err=err,mess=mess)
        endif
        if(err>0) then
            mess='GenerateSample: '//trim(mess);return
        endif
        if(.not.feas) return
    enddo SampleLoop
    if(present(cov)) deallocate(coco)
case(AR1)
    ! First from the marginal
    call Generate(DistID=GAUSS,par=(/par(1), par(2)/sqrt(1.0_mrk-par(3)**2)/),&
                    gen=gen(1),feas=feas,err=err,mess=mess)
    if(err>0) then
        mess='GenerateSample: '//trim(mess);return
    endif
    if(.not.feas) return
    do i=2,n ! subsequent - conditional on t-1
        call Generate(DistID=GAUSS,par=(/(1.0_mrk-par(3))*par(1)+par(3)*gen(i-1), par(2)/),&
                        gen=gen(i),feas=feas,err=err,mess=mess)
        if(err>0) then
            mess='GenerateSample: '//trim(mess);return
        endif
        if(.not.feas) return
    enddo
case default
    err=1;mess='GenerateSample:Fatal:Unavailable Dist'
end select

end subroutine GenerateSample

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elemental function IsCovNeeded(distID)

!#**********************************************************************
!#* Purpose: Is Covariate needed for using this distribution?
!#**********************************************************************
!#* Programmer: Ben Renard, University of Newcastle
!#**********************************************************************
!#* Last modified: 17/07/2008
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*        1.distID
!#* OUT
!#*        1.IsCovNeeded, logical
!#**********************************************************************

character(*), intent(in)::DistID
logical:: IsCovNeeded

if(DistID==GEORGE .or. DistID==GEV_Histo) then
    IsCovNeeded=.true.
else
    IsCovNeeded=.false.
endif

end function IsCovNeeded



!*********!
! PRIVATE !
!*********!

pure subroutine CheckParSize(DistID, par, ok, err, mess)

!^**********************************************************************
!^* Purpose: check size(par)=ParNumber(dist)
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:30/05/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.DistID, the distribution
!^*    1.par, parameter vector
!^* OUT
!^*    1.ok
!^*    2.err, error code
!^*    3.mess, error message
!^**********************************************************************

character(*), intent(in)::DistID
real(mrk), intent(in)::par(:)
logical, intent(out)::ok
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::npar

err=0;mess='';ok=.false.
call GetParNumber(DistID, npar, err, mess)
if(err>0) then
    mess='CheckParSize: '//trim(mess);return
endif
if(size(par)==npar) ok=.true.

end subroutine CheckParSize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine Uni_logP(x, a, b, logp, feas, isnull)

!^**********************************************************************
!^* Purpose: Log-density of a U[a,b]
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:16/07/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. a, lower bound
!^*    3. b, upper bound
!^* OUT
!^*    1.logp, log-density
!^*    2.feas, feasible?
!^*    3.isnull, is density=0?
!^**********************************************************************

real(mrk), intent(in)::x, a, b
real(mrk), intent(out):: logp
logical, intent(out)::feas, isnull

logp=UndefRN;feas=.true.;isnull=.false.

if(b<=a) then
    feas=.false.
    return
endif

if(x<a .OR. x>b) then
    isnull=.true.
    return
else
    logp=-1.0_mrk*log(b-a)
endif

end subroutine Uni_Logp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine Uni_cdf(x, a, b, cdf, feas)

!^**********************************************************************
!^* Purpose: cdf of a U[a,b]
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:16/07/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. a, lower bound
!^*    3. b, upper bound
!^* OUT
!^*    1.cdf, log-density
!^*    2.feas, feasible?
!^**********************************************************************

real(mrk), intent(in)::x, a, b
real(mrk), intent(out):: cdf
logical, intent(out)::feas

cdf=UndefRN;feas=.true.

if(b<=a) then
    feas=.false.
    return
endif

if(x<a) then
    cdf=0.0_mrk
elseif(x>b) then
    cdf=1.0_mrk
else
    cdf=(x-a)/(b-a)
endif

end subroutine Uni_cdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine Triangle_logP(x, peak, a, b, logp, feas, isnull)

!^**********************************************************************
!^* Purpose: Log-density of a Tiangle[peak,a,b]
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 02/08/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. peak, peak
!^*    3. a, lower bound
!^*    4. b, upper bound
!^* OUT
!^*    1.logp, log-density
!^*    2.feas, feasible?
!^*    3.isnull, is density=0?
!^**********************************************************************

real(mrk), intent(in)::x, peak, a, b
real(mrk), intent(out):: logp
logical, intent(out)::feas, isnull

logp=UndefRN;feas=.true.;isnull=.false.

if(b<=a .or. peak<a .or. peak>b) then
    feas=.false.
    return
endif

if(x<=a .OR. x>=b) then
    isnull=.true.
    return
else
    if(x<=peak) then
        logp=log(2._mrk)+log(x-a)-log(b-a)-log(peak-a)
    else
        logp=log(2._mrk)+log(b-x)-log(b-a)-log(b-peak)
    endif
endif

end subroutine Triangle_Logp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine Triangle_cdf(x, peak, a, b, cdf, feas)

!^**********************************************************************
!^* Purpose: cdf of a Triangle[peak,a,b]
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 02/08/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. peak, peak
!^*    3. a, lower bound
!^*    4. b, upper bound
!^* OUT
!^*    1.cdf
!^*    2.feas, feasible?
!^**********************************************************************

real(mrk), intent(in)::x, peak,a, b
real(mrk), intent(out):: cdf
logical, intent(out)::feas

cdf=UndefRN;feas=.true.

if(b<=a .or. peak<a .or. peak>b) then
    feas=.false.
    return
endif

if(x<a) then
    cdf=0.0_mrk
elseif(x>b) then
    cdf=1.0_mrk
elseif(x<=peak) then
    cdf=((x-a)**2)/( (b-a)*(peak-a) )
else
    cdf= (peak-a)/(b-a) + ( ( (x-peak)*(2*b-x-peak) )/( (b-a)*(b-peak) ) )
endif

end subroutine Triangle_cdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine Triangle_Quant(p, peak, a, b, q, feas)

!^**********************************************************************
!^* Purpose: quantile of a Triangle[peak,a,b]
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 12/01/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. p, probability
!^*    2. peak, peak
!^*    3. a, lower bound
!^*    4. b, upper bound
!^* OUT
!^*    1.q, quantile
!^*    2.feas, feasible?
!^**********************************************************************

real(mrk), intent(in)::p, peak,a, b
real(mrk), intent(out):: q
logical, intent(out)::feas
real(mrk)::A1

q=UndefRN;feas=.true.

if(b<=a .or. peak<a .or. peak>b) then
    feas=.false.
    return
endif

A1=(peak-a)/(b-a)
if(p<=A1) then
    q= a+sqrt(p*(b-a)*(peak-a))
else
    q= b-sqrt((1-p)*(b-a)*(b-peak))
endif

end subroutine Triangle_Quant

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GEV_logP(x, loc, scale, shape,logp, feas, isnull)

!^**********************************************************************
!^* Purpose: Log-density of a GEV(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:23/10/2008
!^**********************************************************************
!^* Comments: if abs(shape)<0.001, Gumbel is computed
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.logp, log-density
!^*    2.feas, feasible?
!^*    3.isnull, is density=0?
!^**********************************************************************

real(mrk), intent(in)::x, loc, scale, shape
real(mrk), intent(out):: logp
logical, intent(out)::feas, isnull
!locals
real(mrk)::redx, redx2

logp=UndefRN;feas=.true.;isnull=.false.

if(scale<=0._mrk) then
    feas=.false.;return
endif

redx=(x-loc)/scale;redx2=1._mrk-shape*redx
if(abs(shape)<0.001_mrk) then !Gumbel
    if(redx<=-log(HugeRe)) then ! overflow ahead!
        isnull=.true.;return
    endif
    logp=-log(scale)-redx-exp(-redx)
else !GEV
    if((redx2<=0._mrk).or.(log(redx2)/shape>=log(HugeRe))) then !overflow
        isnull=.true.;return
    else
        logp=-log(scale)+((1._mrk/shape)-1._mrk)*log(redx2)-exp(log(redx2)/shape)
    endif
endif

end subroutine GEV_Logp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine PearsonIII_Logp(x, loc, scale, shape,logp, feas, isnull)

!^**********************************************************************
!^* Purpose: Log-density of a PearsonIII(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:03/04/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.logp, log-density
!^*    2.feas, feasible?
!^*    3.isnull, is density=0?
!^**********************************************************************
use utilities_dmsl_kit, only:gammaln
real(mrk), intent(in)::x, loc, scale, shape
real(mrk), intent(out):: logp
logical, intent(out)::feas, isnull
!locals
real(mrk)::redx

logp=UndefRN;feas=.true.;isnull=.false.

if(scale==0._mrk .or. shape<=0._mrk) then
    feas=.false.;return
endif
redx=(x-loc)/scale
if(redx<0._mrk) then
    isnull=.true.;return
endif

logp=-log(abs(scale))-gammaln(shape)+(shape-1._mrk)*log(redx)-redx

end subroutine PearsonIII_Logp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GPD_logP(x, loc, scale, shape,logp, feas, isnull)

!^**********************************************************************
!^* Purpose: Log-density of a GPD(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:26/11/2010
!^**********************************************************************
!^* Comments: if abs(shape)<0.001, expo is computed
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.logp, log-density
!^*    2.feas, feasible?
!^*    3.isnull, is density=0?
!^**********************************************************************

real(mrk), intent(in)::x, loc, scale, shape
real(mrk), intent(out):: logp
logical, intent(out)::feas, isnull
!locals
real(mrk)::redx, redx2

logp=UndefRN;feas=.true.;isnull=.false.

if(scale<=0._mrk) then
    feas=.false.;return
endif

redx=(x-loc)/scale;redx2=1._mrk-shape*redx
if(abs(shape)<0.001_mrk) then !expo
    if(redx<=-log(HugeRe)) then ! overflow ahead!
        isnull=.true.;return
    endif
    if(redx>=0._mrk) then
        logp=-log(scale)-redx
    else
        isnull=.true.;return
    endif
else !GPD
    if((redx2<=0._mrk).or.(redx<0._mrk)) then !overflow
        isnull=.true.;return
    else
        logp=-log(scale)+((1._mrk/shape)-1._mrk)*log(redx2)
    endif
endif

end subroutine GPD_Logp


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GEV_cdf(x, loc,scale,shape,cdf, feas)

!^**********************************************************************
!^* Purpose: cdf of a GEV(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:23/10/2008
!^**********************************************************************
!^* Comments: if abs(shape)<0.001, Gumbel is computed
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.cdf, log-density
!^*    2.feas, feasible?
!^**********************************************************************

real(mrk), intent(in)::x, loc,scale,shape
real(mrk), intent(out):: cdf
logical, intent(out)::feas
!locals
real(mrk)::redx, redx2

cdf=UndefRN;feas=.true.

if(scale<=0._mrk) then
    feas=.false.;return
endif

redx=(x-loc)/scale;redx2=1._mrk-shape*redx
if(abs(shape)<0.001_mrk) then !Gumbel
    if(redx<=-log(HugeRe)) then ! overflow ahead!
        cdf=0._mrk;return
    endif
    cdf=exp(-exp(-redx))
else !GEV
    if(redx2<=0._mrk) then ! out of support
        if(shape>0._mrk) then
            cdf=1._mrk;return
        else
            cdf=0._mrk;return
        endif
    else
        if(log(redx2)/shape>=log(HugeRe)) then !overflow
            cdf=0._mrk
        else
            cdf=exp(-exp(log(redx2)/shape))
        endif
    endif
endif

end subroutine GEV_cdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine PearsonIII_cdf(x, loc,scale,shape,cdf, feas)

!^**********************************************************************
!^* Purpose: cdf of a PearsonIII(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:03/04/2012
!^**********************************************************************
!^* Comments: Stedinger's approximation - not sur how robust this is...
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.cdf, cdf
!^*    2.feas, feasible?
!^**********************************************************************
use numerix_dmsl_kit, only:normal_cdf

real(mrk), intent(in)::x, loc,scale,shape
real(mrk), intent(out):: cdf
logical, intent(out)::feas
!locals
real(mrk)::moy,sd,skew,foo,foobar,up

cdf=UndefRN;feas=.true.

if((scale==0._mrk).or.(shape<=0._mrk)) then
    feas=.false.;return
endif

moy = loc + shape*scale
sd = sqrt(shape)*abs(scale)
skew=(2._mrk*(scale/abs(scale)))/sqrt(shape)
foo=((x-moy)/sd+(2._mrk/skew))*0.5_mrk*skew
if(foo<0._mrk) then
    if(scale>0._mrk) then
        cdf=0._mrk;return
    else
        cdf=1._mrk;return
    endif
endif
foobar=foo**(1._mrk/3._mrk)
up=(foobar-1._mrk+(skew**2)/36._mrk)*(6._mrk/skew)
cdf=normal_cdf(x=up, mean=0._mrk, sdev=1._mrk)
end subroutine PearsonIII_cdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine PearsonIII_cdf_bugged(x, loc,scale,shape,cdf, feas)

!^**********************************************************************
!^* Purpose: cdf of a PearsonIII(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:03/04/2012
!^**********************************************************************
!^* Comments: Apparently only works for strictly-positively-skewed PIII
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.cdf, log-density
!^*    2.feas, feasible?
!^**********************************************************************
use utilities_dmsl_kit, only:gammaln,gammq,gammp

real(mrk), intent(in)::x, loc,scale,shape
real(mrk), intent(out):: cdf
logical, intent(out)::feas
!locals
real(mrk)::redx

cdf=UndefRN;feas=.true.

if(scale==0._mrk .or. shape<=0._mrk) then
    feas=.false.;return
endif
redx=(x-loc)/scale
if(redx<=0._mrk) then
    if(scale>0._mrk) then
        cdf=0._mrk;return
    else
        cdf=1._mrk;return
    endif
else
    !if(scale>0._mrk) then
        cdf=gammp(shape,redx)/exp(gammaln(shape))
    !else
    !    cdf=gammq(shape,redx)/exp(gammaln(shape))
    !endif
endif
end subroutine PearsonIII_cdf_bugged

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GPD_cdf(x, loc,scale,shape,cdf, feas)

!^**********************************************************************
!^* Purpose: cdf of a GPD(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:26/11/2010
!^**********************************************************************
!^* Comments: if abs(shape)<0.001, Expo is computed
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. x, value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.cdf, log-density
!^*    2.feas, feasible?
!^**********************************************************************

real(mrk), intent(in)::x, loc,scale,shape
real(mrk), intent(out):: cdf
logical, intent(out)::feas
!locals
real(mrk)::redx, redx2

cdf=UndefRN;feas=.true.

if(scale<=0._mrk) then
    feas=.false.;return
endif

redx=(x-loc)/scale;redx2=1._mrk-shape*redx
if(redx<=0._mrk) then
    cdf=0._mrk;return
endif
if(abs(shape)<0.001_mrk) then !Expo
    if(redx<=-log(HugeRe)) then ! overflow ahead!
        cdf=0._mrk;return
    endif
    cdf=1._mrk-exp(-redx)
else !GPD
    if(redx2<=0._mrk) then ! out of support
        if(shape>0._mrk) then
            cdf=1._mrk;return
        else
            cdf=0._mrk;return
        endif
    else
        if(log(redx2)/shape>=log(HugeRe)) then !overflow
            cdf=0._mrk
        else
            cdf=1._mrk-exp(log(redx2)/shape)
        endif
    endif
endif
end subroutine GPD_cdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GEV_Quant(p, loc,scale,shape,q, feas)

!^**********************************************************************
!^* Purpose: quantile of a GEV(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:08/11/2008
!^**********************************************************************
!^* Comments: if abs(shape)<0.001, Gumbel is computed
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. p, p-value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.q, quantile
!^*    2.feas, feasible?
!^**********************************************************************

real(mrk), intent(in)::p, loc,scale,shape
real(mrk), intent(out):: q
logical, intent(out)::feas
!locals

q=UndefRN;feas=.true.

if((scale<=0._mrk).or.(p<=0._mrk).or.(p>=1._mrk)) then
    feas=.false.;return
endif

if(abs(shape)<0.001_mrk) then !Gumbel
    q=loc-scale*log(-log(p))
else !GEV
    q=loc+(scale/shape)*(1._mrk-(-log(p))**shape)
endif

end subroutine GEV_Quant

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine PearsonIII_Quant(p, loc,scale,shape,q, feas)

!^**********************************************************************
!^* Purpose: quantile of a PearsonIII(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:03/04/2012
!^**********************************************************************
!^* Comments: Stedinger's approximation - not sur how robust this is...
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. p, p-value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.q, quantile
!^*    2.feas, feasible?
!^**********************************************************************
use numerix_dmsl_kit, only:normal_icdf

real(mrk), intent(in)::p, loc,scale,shape
real(mrk), intent(out):: q
logical, intent(out)::feas
!locals
real(mrk)::moy,sd,skew,up,kp

q=UndefRN;feas=.true.

if((scale==0._mrk).or.(shape<=0._mrk).or.(p<=0._mrk).or.(p>=1._mrk)) then
    feas=.false.;return
endif

moy = loc + shape*scale
sd = sqrt(shape)*abs(scale)
skew=(2._mrk*(scale/abs(scale)))/sqrt(shape)
up=normal_icdf(p=p,mean=0._mrk, sdev=1._mrk)
kp=(2._mrk/skew)*(1._mrk+(skew*up)/6._mrk-(skew**2)/36._mrk)**3-(2._mrk/skew)

q=moy+sd*kp

end subroutine PearsonIII_Quant
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GPD_Quant(p, loc,scale,shape,q, feas)

!^**********************************************************************
!^* Purpose: quantile of a GPD(location, scale, shape)
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:26/11/2010
!^**********************************************************************
!^* Comments: if abs(shape)<0.001, Gumbel is computed
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. p, p-value
!^*    2. loc, location
!^*    3. scale
!^*    3. shape
!^* OUT
!^*    1.q, quantile
!^*    2.feas, feasible?
!^**********************************************************************

real(mrk), intent(in)::p, loc,scale,shape
real(mrk), intent(out):: q
logical, intent(out)::feas

q=UndefRN;feas=.true.

if((scale<=0._mrk).or.(p<=0._mrk).or.(p>=1._mrk)) then
    feas=.false.;return
endif

if(abs(shape)<0.001_mrk) then !Expo
    q=loc-scale*log(1._mrk-p)
else !GPD
    q=loc+(scale/shape)*(1._mrk-(1._mrk-p)**shape)
endif

end subroutine GPD_Quant

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GEV_dev(loc,scale,shape,dev,feas)

!^**********************************************************************
!^* Purpose: deviate from a GEV
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified: 08/07/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. loc, location
!^*    2. scale
!^*    3. shape
!^* OUT
!^*    1.dev, deviate
!^*    2.feas, feasible?
!^**********************************************************************

use numerix_dmsl_kit, only: uniran

real(mrk), intent(in)::loc,scale,shape
real(mrk), intent(out):: dev
logical, intent(out)::feas
!locals
real(mrk)::u

feas=.true.;dev=undefRN

call uniran(u)

call GEV_Quant(u, loc,scale,shape,dev, feas)

end subroutine GEV_dev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine PearsonIII_dev(loc,scale,shape,dev,feas)

!^**********************************************************************
!^* Purpose: deviate from a PearsonIII
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified: 03/04/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. loc, location
!^*    2. scale
!^*    3. shape
!^* OUT
!^*    1.dev, deviate
!^*    2.feas, feasible?
!^**********************************************************************

use numerix_dmsl_kit, only: uniran

real(mrk), intent(in)::loc,scale,shape
real(mrk), intent(out):: dev
logical, intent(out)::feas
!locals
real(mrk)::u

feas=.true.;dev=undefRN

call uniran(u)

call PearsonIII_Quant(u, loc,scale,shape,dev, feas)

end subroutine PearsonIII_dev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine GPD_dev(loc,scale,shape,dev,feas)

!^**********************************************************************
!^* Purpose: deviate from a GPD
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified: 26/11/2010
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. loc, location
!^*    2. scale
!^*    3. shape
!^* OUT
!^*    1.dev, deviate
!^*    2.feas, feasible?
!^**********************************************************************

use numerix_dmsl_kit, only: uniran

real(mrk), intent(in)::loc,scale,shape
real(mrk), intent(out):: dev
logical, intent(out)::feas
!locals
real(mrk)::u

feas=.true.;dev=undefRN

call uniran(u)

call GPD_Quant(u, loc,scale,shape,dev, feas)

end subroutine GPD_dev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine InvChi2_cdf(x, DoF, scale, cdf, feas)

!^**********************************************************************
!^* Purpose: cdf of InvChi2(Dof, scale)
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified: 16/07/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x
!^*    2.DoF
!^*    3.scale
!^* OUT
!^*    1.cdf
!^*    2.feas
!^**********************************************************************
use utilities_dmsl_kit, only: gammq
real(mrk), intent(in) :: x, DoF, scale
real(mrk), intent(out) :: cdf
logical, intent(out) :: feas
!locals
real(mrk)::halfdof

cdf=UndefRN;feas=.true.
if(DoF<=0.0_mrk .OR. scale<=0.0_mrk) then
    feas=.false.;return
endif

if(x<=0.0_mrk) then
    cdf=0.0_mrk
else
    halfdof=0.5_mrk*DoF
    cdf=gammq(halfdof,halfdof*((scale)**2/x))
endif
end subroutine InvChi2_cdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure function Real2Int(r)
!^**********************************************************************
!^* Purpose: Transforme a real to the nearest integer.
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:24/11/2011
!^**********************************************************************
!^* Comments: For a value 0.5->??? lol I don't know.
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.real
!^* OUT
!^*     1.integer
!^**********************************************************************
real(mrk),intent(in)::r
integer(mik)::Real2Int
Real2Int=int((r+0.5_mrk),mik)

end function Real2Int

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine Ber_pdf(x,p,pdf,isnull)
!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:24/11/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.x
!^*     2.p
!^* OUT
!^*     1.pdf
!^*     2.isnull
!^**********************************************************************
integer(mik),intent(in)::x
real(mrk),intent(in)::p
real(mrk),intent(out)::pdf
logical,intent(out)::isnull

if (x/=1_mik .and. x/=0_mik) isnull=.true.
if (isnull) pdf=0.0_mrk
if (x==1_mik) pdf=p
if (x==0_mik) pdf=1._mrk-p

end subroutine Ber_pdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine Geo_pdf(x,p,pdf,isnull)
!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:24/11/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.x
!^*     2.p
!^* OUT
!^*     1.pdf
!^*     2.isnull
!^**********************************************************************
integer(mik),intent(in)::x
real(mrk),intent(in)::p
real(mrk),intent(out)::pdf
logical,intent(out)::isnull

if (x<=0_mik) isnull=.true.
if (isnull) then
    pdf=0.0_mrk
else
    pdf=((1._mrk-p)**(x-1_mik))*p
endif

end subroutine Geo_pdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ber_dev(p)
!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:24/11/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.x
!^*     2.p
!^* OUT
!^*     1.pdf
!^*     2.isnull
!^**********************************************************************
use numerix_dmsl_kit, only:uniran

real(mrk),intent(in)::p
integer(mik)::Ber_dev
real(mrk)::dev

call uniran(dev)
if (dev<p) then
    Ber_dev=1_mik
else
    Ber_dev=0_mik
endif

end function Ber_dev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Geo_dev(p)
!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:24/11/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.x
!^*     2.p
!^* OUT
!^*     1.pdf
!^*     2.isnull
!^**********************************************************************
use numerix_dmsl_kit, only:uniran

real(mrk),intent(in)::p
integer(mik)::geo_dev
real(mrk)::dev
integer(mik)::temp

temp=1
call uniran(dev)
do while (dev>p)
    temp=temp+1
    call uniran(dev)
enddo
Geo_dev=temp

end function Geo_dev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Binomial_Quant(p,success_prob, n)
use numerix_dmsl_kit,only:binomial_pmf
real(mrk),intent(in)::p,success_prob
integer(mik),intent(in)::n
integer(mik)::Binomial_Quant
! locals
integer(mik)::k
real(mrk)::cumul

k=0
cumul=binomial_pmf(theta=k,n=n,p=success_prob)
do while(p>cumul)
    k=k+1;
    cumul=cumul+binomial_pmf(theta=k,n=n,p=success_prob)
enddo
Binomial_Quant=k
end function Binomial_Quant

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine Beta_Quant(p,a,b,q,feas)

!^**********************************************************************
!^* Purpose: quantile of a beta(a,b)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:02/01/2021
!^**********************************************************************
!^* Comments: Simple bisection - easy but certainly not the fastest.
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: look in the literature for faster algorithms
!^**********************************************************************
!^* IN
!^*    1. p, nonexceedance probability
!^*    2. a, shape1
!^*    3. b, shape2
!^* OUT
!^*    1.q, quantile
!^*    2.feas, feasible?
!^**********************************************************************

use utilities_dmsl_kit, only:betai

real(mrk), intent(in)::p,a,b
real(mrk), intent(out):: q
logical, intent(out)::feas
!locals
real(mrk)::xl,xu,fl,fu,x,fx
integer(mik)::i
real(mrk), parameter::tolF=0.000000001_mrk,tolX=0.0000000001_mrk
integer(mik), parameter::itmax=100000

q=UndefRN;feas=.true.
xl=0._mrk;xu=1._mrk;fl=0._mrk-p;fu=1._mrk-p

do i=1,itmax
    x=xl-(fl/(fu-fl))*(xu-xl)
    fx=betai(a,b,x)-p
    if(abs(fx)<tolF) then;q=x;return;endif
    if(fx>0._mrk) then
        xu=x;fu=fx
    else
        xl=x;fl=fx
   endif
   if(0.5_mrk*(xu-xl)<tolX) then;q=x;return;endif
enddo
! not yet returned => itmax exceeded
feas=.false.
end subroutine Beta_Quant

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Distribution_tools

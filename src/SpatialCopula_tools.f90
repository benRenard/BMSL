module SpatialCopula_tools

!~**********************************************************************
!~* Purpose: copulas for spatial analysis
!~**********************************************************************
!~* Programmer: Ben Renard, University of Newcastle
!~**********************************************************************
!~* Last modified: 08/06/2008
!~**********************************************************************
!~* Comments: the term 'spatial copula' might be a poor choice - the idea
!~* is to consider copulas parametrized by a dependence matrix - leading 
!~* to useful analogy with geostatistics... Mainly elliptical copulas 
!~* for the time being, maybe more sophisticated stuff in the future...
!~* For computational efficiency purpose, there is no pre-check for 
!~* definite positiveness - directly incorporated in pdf and cdf computation
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
public :: ApplyCop, ApplyCopPrime

character(100), parameter, PUBLIC:: GAUSSIAN_COP='Gaussian', &
									STUDENT_COP='Student'

Contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine ApplyCop(CopId,u,par,isInv,logdet,C,feas,err,mess)

!^**********************************************************************
!^* Purpose: compute C(u) - allow deriving multivariate cdf
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:08/06/2008
!^**********************************************************************
!^* Comments: A big issue is that multivariate cdf rarely have 
!^* a closed-form expression... left for future work; see here:
!^* http://www.math.wsu.edu/faculty/genz/homepage for numerical methods
!^* Optional arguments aimed at increasing computational efficiency
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.CopID, the distribution
!^*		2.u, values
!^*		3.par, parameters (dependence matrix)
!^*		4.[IsInv], is the dependence matrix already inverted? assume no by default
!^*		5.[LogDet], log-determinant of the dependence matrix
!^* OUT
!^*		1.C, result
!^*		2.feas, is par/x feasible?
!^*		3.err, error code
!^*		4.mess, error message
!^**********************************************************************


character(*), intent(in)::CopID
real(mrk), intent(in)::u(:), par(:,:)
logical, intent(in), optional::IsInv
real(mrk),intent(in),optional::logDet
real(mrk), intent(out)::C
logical, intent(out)::feas
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
logical:: IsI
!Init
C=UndefRN;feas=.false.;err=0;mess=''

if(present(isInv)) then
	IsI=IsInv
else
	IsI=.false.
endif
!Feasability - u in [0,1]?
if((any(u<0.0_mrk)).or.(any(u>1.0_mrk))) return
!size OK?
if(size(u)/=size(par,dim=1)) then
	err=2;mess='ApplyCop:Fatal: Dimension mismatch';return
endif

!Compute
select case(CopID)
case(GAUSSIAN_COP)
	err=3;mess='ApplyCop:Fatal:Gaussian Copula:No closed-form expression'
case(STUDENT_COP)
	err=3;mess='ApplyCop:Fatal:STUDENT Copula:No closed-form expression'
case default
	err=1;mess='ApplyCop:Fatal:Unavailable Copula'
end select

end subroutine ApplyCop


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ApplyCopPrime(CopId,u,par,AdditionalPar,isInv,logdet,loga,C,feas,isnull,err,mess)

!^**********************************************************************
!^* Purpose: compute C'(u) - this is the one used for deriving the pdf
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
!^**********************************************************************
!^* Last modified:08/06/2008
!^**********************************************************************
!^* Comments: Strongly advised to handle matrix inversion and determinant
!^* computation outside of this sub for computational eficiency
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.CopID, the distribution
!^*		2.u, values
!^*		3.par, parameters (dependence matrix)
!^*		4.[AdditionalPar], additional parameters (eg tail dependence)
!^*		5.[IsInv], is the dependence matrix already inverted? assume no by default
!^*		6.[LogDet], log-determinant of the dependence matrix
!^*		7.Loga, is the log of c' required?
!^* OUT
!^*		1.C, result
!^*		2.feas, is par/x feasible?
!^*		3.isNull, c'==0?
!^*		4.err, error code
!^*		5.mess, error message
!^**********************************************************************

use linalg_dmsl_kit, only:choles_invrt
use numerix_dmsl_kit, only:normal_icdf,tdist_icdf,tdist_icdf_rough,tdist_logp

character(*), intent(in)::CopID
real(mrk), intent(in)::u(:), par(:,:)
logical, intent(in), optional::IsInv
real(mrk),intent(in),optional::logDet,AdditionalPar(:)
logical, intent(in)::loga
real(mrk), intent(out)::C
logical, intent(out)::feas, isNull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
logical:: IsI, pD
real(mrk):: logc,logd,logp
real(mrk), allocatable::sigInv(:,:), Id(:,:),vect(:), vect2(:)
integer(mik)::p, i
real(mrk), parameter::BotchThresh=0.000001_mrk

!Init
C=UndefRN;feas=.true.;isNull=.false.;err=0;mess=''
p=size(par,dim=1)

if(present(isInv)) then
	IsI=IsInv
else
	IsI=.false.
endif
!Feasability - u in [0,1]?
if((any(u<0.0_mrk)).or.(any(u>1.0_mrk))) then
	feas=.false.;return
endif
! Any uk=0 or 1 means an impossible value - zero density
if((any(u==0.0_mrk)).or.(any(u==1.0_mrk))) then
	isnull=.true.;return
endif

!size OK?
if(size(u)/=p) then
	err=2;mess='ApplyCopPrime:Fatal: Dimension mismatch';return
endif

!Compute
select case(CopID)
case(GAUSSIAN_COP)
	!Allocate matrix for sigma inverse
	if(allocated(sigInv)) deallocate(sigInv)
	allocate(sigInv(p,p))
	if(allocated(vect)) deallocate(vect)
	allocate(vect(p))
	if(allocated(Id)) deallocate(Id)
	allocate(Id(p,p))
	! Inverse variance matrix / get logDet if needed
	if((.not.IsI).or.(.not.(present(logdet)))) then
		call choles_invrt(a=par,ainv=sigInv,posDefinite=pD,&
				logDet=logd,err=err,message=mess)
		If(err>0) then
			mess='ApplyCopPrime:'//trim(mess);return
		endif
	endif
	! preliminaries
	if(IsI) sigInv=par
	if(present(logDet)) logd=logDet
	Id=0.0_mrk
	do i=1,p
		Id(i,i)=1.0_mrk
		vect(i)=normal_icdf(u(i),0._mrk,1._mrk)
	enddo
	! Start computation
	logc=-0.5_mrk*logd-0.5_mrk*dot_product(matmul(vect,sigInv-Id),vect)
	if(loga) then
		C=logc
	else
		C=exp(logc)
	endif
case(Student_COP)
	!!!!!!!!!!!!!MEGA-BOTCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(any(u<BotchThresh).or. any(1._mrk-u<BotchThresh)) then
		isnull=.true.;return
	endif
	!!!!!!!!!!!!!END MEGA-BOTCH!!!!!!!!!!!!!!!!!!!!!!!!

	! Check AdditionalPar (tail dependence dof) is present
	if(.not.present(AdditionalPar)) then
		err=3;mess='ApplyCopPrime:Fatal:STUDENT_COP:[AdditionalPar] needed for the student copula';return
	endif
	!Allocate matrix for sigma inverse
	if(allocated(sigInv)) deallocate(sigInv)
	allocate(sigInv(p,p))
	if(allocated(vect)) deallocate(vect)
	allocate(vect(p))
	if(allocated(vect2)) deallocate(vect2)
	allocate(vect2(p))
	! Inverse variance matrix / get logDet if needed
	if((.not.IsI).or.(.not.(present(logdet)))) then
		call choles_invrt(a=par,ainv=sigInv,posDefinite=pD,&
				logDet=logd,err=err,message=mess)
		If(err>0) then
			mess='ApplyCopPrime:'//trim(mess);return
		endif
	endif
	! preliminaries
	if(IsI) sigInv=par
	if(present(logDet)) logd=logDet
	do i=1,p
		vect(i)=tdist_icdf(u(i),AdditionalPar(1))
		vect2(i)=tdist_logp(x=vect(i),v=AdditionalPar(1))
	enddo
	! Start computation
	CALL MultivariateStudent_LogPdf(x=vect,dof=AdditionalPar(1),&
					mu=(/(0._mrk,i=1,p)/),SigmaInv=sigInv,&
					SigmaLogDet=logd,&
					logP=logp, feas=feas, IsNull=IsNull,&
					err=err,mess=mess)
	If(err>0) then
		mess='ApplyCopPrime:'//trim(mess);return
	endif
	If(IsNull) then
		return
	endif
	logc=logp - sum(vect2)
	if(loga) then
		C=logc
	else
		C=exp(logc)
	endif
case default
	err=1;mess='ApplyCopPrime:Fatal:Unavailable Copula'
end select

end subroutine ApplyCopPrime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!
! LOCAL FUNCTIONS !
!!!!!!!!!!!!!!!!!!!

pure subroutine MultivariateStudent_LogPdf(x,dof,mu,SigmaInv,SigmaLogDet,logP, feas, IsNull, err,mess)

!^**********************************************************************
!^* Purpose: log-pdf of the multivariate student
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:20/01/2010
!^**********************************************************************
!^* Comments: 1/ SigmaInv and SigmaDet already computed elsewhere
!^*		2/ Sigma is not equal to the covariance; cov=(dof/(dof-2))*sigma
!^**********************************************************************
!^* References: Wikipedia
!^**********************************************************************
!^* 2Do List: a "MultivariateDistributionLib" module....
!^**********************************************************************
!^* IN
!^*		1. x, input vector
!^*		2. dof, degrees of freedom
!^*		3. mu, mean vector
!^*		4. SigmaInv, Inverse of the scatter matrix
!^*		5. SigmaLogDet, log-Determinant of the scatter matrix
!^* OUT
!^*		1.logP, log-pdf
!^*		2.feas, is par/x feasible?
!^*		3.isNull, LogP==0?
!^*		4.err, error code
!^*		5.mess, error message
!^**********************************************************************

use utilities_dmsl_kit, only: gammaln, pi

real(mrk), intent(in)::x(:),dof,mu(:),SigmaInv(:,:), SigmaLogDet
real(mrk),intent(out)::LogP
logical, intent(out)::feas,IsNull
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
integer(mik)::d
real(mrk)::toto

!Init
LogP=UndefRN;feas=.true.;isNull=.false.;err=0;mess=''
d=size(x)

! Check Dim
if(size(mu)/=d) then
	err=1;mess='MultivariateStudent_LogPdf:Fatal: Dimension mismatch [mu]';return
endif
if((size(SigmaInv,dim=1)/=d) .or. (size(SigmaInv,dim=2)/=d)) then
	err=1;mess='MultivariateStudent_LogPdf:Fatal: Dimension mismatch [SigmaInv]';return
endif

! Check Feasability
if(dof<=0) then
	feas=.false.;return
endif

toto=dot_product(matmul(x-mu,SigmaInv),x-mu)/dof
LogP=gammaln((dof+d)*0.5_mrk)-gammaln(dof*0.5_mrk) &
    -0.5_mrk*d*log(dof)-0.5_mrk*d*log(pi)-0.5_mrk*SigmaLogDet &
	-((dof+d)*0.5_mrk)*log(1._mrk+toto)

end subroutine MultivariateStudent_LogPdf


end module SpatialCopula_tools
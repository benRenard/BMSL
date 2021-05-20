module Dependogram_tools

!~**********************************************************************
!~* Purpose: 
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref
!~**********************************************************************
!~* Last modified: 11/10/2008
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

implicit none
Private
public :: GetDependParNumber, GetDependParName, &
		  GetDependParFeas, GetDepMatrix

Character(100), parameter, PUBLIC:: &
            Expo_dep='Exponential',Power='Power'

Contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetDependParNumber(dependID, npar, err, mess)

!^**********************************************************************
!^* Purpose: returns the number of parameters of the dependogram
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:12/10/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.dependID, the dependogram ID
!^* OUT
!^*		1.nPar, the numbr of parameters
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************

character(*),intent(in)::dependID
integer(mik),intent(out)::nPar
integer(mik), intent(out)::err
character(100),intent(out)::mess

err=0;mess='';nPar=UndefIN
select case(DependID)
case(Expo_dep)
    npar=2
case(Power)
    npar=1
case("Gauss")
    npar=2
case default
    err=1;mess='GetDependParNumber:Fatal:Unavailable [dependID]'
end select

end subroutine GetDependParNumber

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetDependParName(dependID, name, err, mess)

!^**********************************************************************
!^* Purpose: returns the names of parameters of the dependogram
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:12/10/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.DependID, the dependogram ID
!^* OUT
!^*		1.name, parameters names
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************

character(*),intent(in)::DependID
character(100),intent(out)::name(:)
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
integer(mik)::npar

err=0;mess='';name=undefCH

! check size
call GetDependParNumber(DependID, npar, err, mess)
if(err>0) then
    mess='GetDependParName: '//trim(mess);return
endif
if(size(name)/=npar) then
    err=2;mess='GetDependParName: dimension mismatch';return
endif

select case(dependID)
case(EXPO_DEP)
    name(1)='range'
    name(2)='nugget'
case("Gauss")
    name(1)='range'
    name(2)='nugget'
case(Power)
    name(1)='par'
case default
    err=1;mess='GetDependParName:Fatal:Unavailable [dependID]'
end select

end subroutine GetDependParName

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetDependParFeas(dependID, par, feas, err, mess)

!^**********************************************************************
!^* Purpose: check parameters feasability
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:12/10/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.dependID, the dependogram ID
!^*		2.par, parameter vector
!^* OUT
!^*		1.feas, feasability
!^*		2.err, error code
!^*		3.mess, error message
!^**********************************************************************

character(*),intent(in)::dependID
real(mrk), intent(in)::par(:)
logical,intent(out)::feas
integer(mik), intent(out)::err
character(100),intent(out)::mess
! Locals
logical::ok
integer(mik)::npar

err=0;mess='';feas=.true.

! check size
call GetDependParNumber(DependID, npar, err, mess)
if(err>0) then
    mess='GetDependParFeas: '//trim(mess);feas=.false.;return
endif
if(size(par)/=npar) then
    err=2;mess='GetDependParFeas:Fatal:size[par]';feas=.false.;return
endif

select case(DependID)
case(EXPO_DEP)
    if ((par(1)<0.0_mrk).or.(par(2)<=0.0_mrk).or.(par(2)>1.0_mrk)) feas=.false.
case("Gauss")
    if ((par(1)<0.0_mrk).or.(par(2)<=0.0_mrk).or.(par(2)>1.0_mrk)) feas=.false.
case(Power)
    if (par(1)<0.0_mrk .or. par(1)>1._mrk) feas=.false.
case default
    err=1;mess='GetDependParFeas:Fatal:Unavailable [dependID]'
end select

end subroutine GetDependParFeas

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine GetDepMatrix(dependID,D,par,V,feas,err,mess)

!^**********************************************************************
!^* Purpose: compute dependence matrix from distance matrix
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified: 11/10/2008
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.dependID, the dependogram function
!^*		2.D, distance matrix
!^*		3.par, parameters of dependogram function
!^* OUT
!^*		1.V, the dependence matrix
!^*		2.feas, feasibility
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*), intent(in)::dependID
real(mrk), intent(in)::D(:,:), par(:)
real(mrk), intent(out)::V(:,:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(100),intent(out)::mess
!locals
integer(mik)::i,j,p

!init
V=undefRN;feas=.true.;err=0;mess=''
!Check feasibility
call GetDependParFeas(dependID, par, feas, err, mess)
If(err>0) then
	mess='GetDepMatrix'//trim(mess);feas=.false.;return
endif
if(.not.feas) return
if(any(D<0.0_mrk)) then
	feas=.false.;return
endif
!Check sizes
p=size(D,dim=1)
if((size(D,dim=2)/=p).or.(size(V,dim=2)/=p).or.(size(V,dim=1)/=p)) then
	err=1;mess='GetDepMatrix:Fatal:size [D or V]';return
endif

select case(DependID)
case(EXPO_DEP)
    forall(i=1:p,j=1:p) V(i,j)=Expo_dependogram(D(i,j),par)
case("Gauss")
    forall(i=1:p,j=1:p) V(i,j)=Gauss_dependogram(D(i,j),par)
case(Power)
    forall(i=1:p,j=1:p) V(i,j)=par(1)**D(i,j)
case default
    err=1;mess='GetDepMatrix:Fatal:Unavailable [dependID]'
end select

end subroutine GetDepMatrix

!*********!
! PRIVATE !
!*********!
pure function Expo_dependogram(dist, par)

!#**********************************************************************
!#* Purpose: returns the dependence coeff using exponential dependogram
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref
!#**********************************************************************
!#* Last modified: 12/10/2008
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*		1.dist, distance
!#*		2.par, range and nugget
!#* OUT
!#*		1.Expo_dependogram
!#**********************************************************************
real(mrk), intent(in)::dist, par(2)
real(mrk)::Expo_dependogram

if(dist<0.0_mrk) then
	Expo_dependogram=undefRN;return
endif

if(dist==0.0_mrk) then
	Expo_dependogram=1.0_mrk
else
	Expo_dependogram=par(2)*exp(-par(1)*dist)
endif

end function Expo_dependogram

pure function Gauss_dependogram(dist, par)

!#**********************************************************************
!#* Purpose: returns the dependence coeff using exponential dependogram
!#**********************************************************************
!#* Programmer: Ben Renard, Cemagref
!#**********************************************************************
!#* Last modified: 12/10/2008
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*		1.dist, distance
!#*		2.par, range and nugget
!#* OUT
!#*		1.Expo_dependogram
!#**********************************************************************
real(mrk), intent(in)::dist, par(2)
real(mrk)::Gauss_dependogram

if(dist<0.0_mrk) then
	Gauss_dependogram=undefRN;return
endif

if(dist==0.0_mrk) then
	Gauss_dependogram=1.0_mrk
else
	Gauss_dependogram=par(2)*exp(-par(1)*dist*dist)
endif

end function Gauss_dependogram

end module Dependogram_tools
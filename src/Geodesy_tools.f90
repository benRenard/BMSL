module Geodesy_tools

!~**********************************************************************
!~* Purpose: tools for geodesy, mostly distance computations 
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Created: 08/02/2017
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: GetDistance

interface GetDistance
  module procedure GetDistance_2,GetDistance_n
endinterface GetDistance

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetDistance_2(formula,pt1,pt2,par,cov1,cov2,d,err,mess)
!^**********************************************************************
!^* Purpose: Compute distance between two points
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Created: 08/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.formula
!^*		2.pt1
!^*		3.pt2
!^*		4.[par] optional parameter vector 
!^*		5.[cov1] optional covariate vector at site 1
!^*		6.[cov2] optional covariate vector at site 2
!^* OUT
!^*		1.d, distance
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
character(*),intent(in)::formula
real(mrk),intent(in)::pt1(:),pt2(:)
real(mrk),intent(in),optional::par(:),cov1(:),cov2(:)
real(mrk),intent(out)::d
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='GetDistance_2'
real(mrk),parameter::R = 6371._mrk ! earth's average radius in km
real(mrk),parameter:: toRad = 0.017453292519943295 ! pi/180
real(mrk)::diff(size(pt1))

err=0;mess='';d=undefRN

select case(trim(formula))
case('euclidean','Euclidean')
    diff=pt1-pt2
    d=sqrt(dot_product(diff,diff))
case('haversine','Haversine')
    ! Using Salvador's script in this page:
    ! https://stackoverflow.com/questions/27928/calculate-distance-between-two-latitude-longitude-points-haversine-formula
    diff=(pt1-pt2)*toRad
    d=0.5_mrk-cos(diff(2))/2._mrk+&
      cos(pt1(2)*toRad)*cos(pt2(2)*toRad)*(1._mrk-cos(diff(1)))/2._mrk
    d=2._mrk*R*asin(sqrt(d))
case default
    err=1;mess=trim(procname)//':Fatal:Unavailable formula';return
end select

end subroutine GetDistance_2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetDistance_n(formula,pts,par,cov,D,err,mess)
!^**********************************************************************
!^* Purpose: Compute distance matrix between n points
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Created: 08/02/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.formula
!^*		2.pts
!^*		3.[par] optional parameter vector 
!^*		4.[cov] optional covariate matrix (row: sites, columns: covariates)
!^* OUT
!^*		1.D, distance matrix
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
character(*),intent(in)::formula
real(mrk),intent(in)::pts(:,:)
real(mrk),intent(in),optional::par(:),cov(:,:)
real(mrk),intent(out)::D(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='GetDistance_n'
integer(mik)::n,i,j

err=0;mess='';D=0._mrk

n=size(pts,dim=1)
if(size(d,dim=1)/=n .or. size(d,dim=2)/=n) then
    err=1;mess=trim(procname)//':Fatal:incorrect size [D]';return
endif
if(n==1) return
do i=2,n
    do j=1,i-1
        if(present(cov) .and. present(par)) then
            call GetDistance_2(formula=formula,pt1=pts(i,:),pt2=pts(j,:),&
                               par=par,cov1=cov(i,:),cov2=cov(j,:),&
                               d=D(i,j),err=err,mess=mess)
        else
            call GetDistance_2(formula=formula,pt1=pts(i,:),pt2=pts(j,:),&
                               d=D(i,j),err=err,mess=mess)
        endif
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        D(j,i)=D(i,j)
    enddo
enddo

end subroutine GetDistance_n
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module Geodesy_tools

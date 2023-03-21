module CircularStats_tools

!~**********************************************************************
!~* Purpose: Compute circular stats from angle data x in RADIANS
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified:25/07/2017
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References: https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Analysis.pdf
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1. GetCircularStats, 5 statistics frequently used in circular statistics
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: GetCircularStats

Contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetCircularStats(x,& ! input data
                            n,& !sample size
                            meandir,meanlength,& ! dir & length of mean vector
                            std,var,& ! std and variance
                            package,& ! pack all stats into a single vector
                            err,mess)

!^**********************************************************************
!^* Purpose: Get a summary for circular data x
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:25/07/2017
!^**********************************************************************
!^* Comments: all output variables are optional (except err & mess)
!^* output "package" can be used to pack all stats into a single vector
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, data vector (angles in radians
!^* OUT
!^*    1. [n], sample size
!^*    2. [meandir], direction of mean vector
!^*    3. [meanlength], length of mean vector
!^*    4. [std], circular standard deviation
!^*    5. [var], circular variance
!^*    6. [package], all stats above packed into a single vector
!^*    7.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    8.mess, error message
!^**********************************************************************
use numerix_dmsl_kit, only:getmean, getvar, getCV,getmoments
use utilities_dmsl_kit,only:pi,twopi

real(mrk), intent(in)::x(:)
integer(mik), intent(out), optional::n
real(mrk), intent(out), optional::meandir,meanlength,std,var,package(5)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::p,Xn
real(mrk)::C,S,R,Xmeandir,Xmeanlength,Xstd,Xvar


err=0;mess=''

p=size(x)
C=sum(cos(x))/real(p,mrk)
S=sum(sin(x))/real(p,mrk)
R=sqrt(C**2+S**2)

if( present(n) .or. present(package) ) then
    Xn=p
    if(present(n)) n=Xn
    if(present(package)) package(1)=Xn
endif

if( present(meandir) .or. present(package) ) then
    if(C>0._mrk .and. S>0._mrk) then
        Xmeandir=atan(S/C)
    else if(C>0._mrk .and. S<0._mrk) then
        Xmeandir=atan(S/C)+twopi
    else if(C<0._mrk) then
        Xmeandir=atan(S/C)+pi
    else
        Xmeandir=undefRN
    endif
    if(present(meandir)) meandir=Xmeandir
    if(present(package)) package(2)=Xmeandir
endif

if( present(meanlength) .or. present(package) ) then
    Xmeanlength=R
    if(present(meanlength)) meanlength=Xmeanlength
    if(present(package)) package(3)=Xmeanlength
endif

if( present(std) .or. present(package) ) then
    Xstd=sqrt(-2._mrk*log(R))
    if(present(std)) std=Xstd
    if(present(package)) package(4)=Xstd
endif

if( present(var) .or. present(package) ) then
    Xvar=1._mrk-R
    if(present(var)) var=Xvar
    if(present(package)) package(5)=Xvar
endif

end subroutine GetCircularStats

end module CircularStats_tools

module Geodesy_tools

!~**********************************************************************
!~* Purpose: tools for geodesy, e.g. distance computations,
!~*          raster grids, interpolation, etc.
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Created: 08/02/2017 - Last Update: 03/08/2022
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: GetDistance, validateRasterGrid, rasterToPoints, interpolateOnGrid, &
          ReadRasterASC, WriteRasterASC

type, public:: rasterGridType
    integer(mik) :: ncols=-99,nrows=-99
    character(250)::name='AintGotNoName'
    real(mrk)::xllcorner=undefRN, yllcorner=undefRN, cellsize=undefRN
    real(mrk) :: mv=undefRN ! missing values code
end type rasterGridType

interface GetDistance
  module procedure GetDistance_2,GetDistance_1n,GetDistance_1grid,GetDistance_n
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
!^*     1.formula
!^*     2.pt1
!^*     3.pt2
!^*     4.[par] optional parameter vector
!^*     5.[cov1] optional covariate vector at site 1
!^*     6.[cov2] optional covariate vector at site 2
!^* OUT
!^*     1.d, distance
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
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
subroutine GetDistance_1n(formula,pt1,pts,par,cov1,cov,D,err,mess)
!^**********************************************************************
!^* Purpose: Compute distance vector between one and n points
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Created: 03/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.formula
!^*     2.pt1, single point
!^*     3.pts, n points
!^*     4.[par] optional parameter vector
!^*     5.[cov1] optional covariate vector for single site
!^*     6.[cov] optional covariate matrix (row: sites, columns: covariates)
!^* OUT
!^*     1.D, distance vector
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
character(*),intent(in)::formula
real(mrk),intent(in)::pt1(:),pts(:,:)
real(mrk),intent(in),optional::par(:),cov1(:),cov(:,:)
real(mrk),intent(out)::D(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='GetDistance_1n'
integer(mik)::n,i

err=0;mess='';D=0._mrk

n=size(pts,dim=1)
if(size(D)/=n) then
    err=1;mess=trim(procname)//':Fatal:incorrect size [D]';return
endif
do i=1,n
    if(present(cov1) .and. present(cov) .and. present(par)) then
        call GetDistance_2(formula=formula,pt1=pt1,pt2=pts(i,:),&
                           par=par,cov1=cov1,cov2=cov(i,:),&
                           d=D(i),err=err,mess=mess)
    else
        call GetDistance_2(formula=formula,pt1=pt1,pt2=pts(i,:),&
                           d=D(i),err=err,mess=mess)
    endif
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo

end subroutine GetDistance_1n

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine GetDistance_1grid(formula,pt1,grid,D,err,mess)
!^**********************************************************************
!^* Purpose: Compute distance vector between one point and a grid
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Created: 03/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.formula
!^*     2.pt1, single point
!^*     3.grid, grid on which interpolated values should be computed (rasterGrid object)
!^* OUT
!^*     1.D, distance matrix
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
character(*),intent(in)::formula
real(mrk),intent(in)::pt1(:)
type(rasterGridType), intent(in)::grid
real(mrk),intent(out)::D(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='GetDistance_1grid'
integer(mik)::i,j
real(mrk)::pt2(2)

err=0;mess='';D=undefRN

! Validate grid
call validateRasterGrid(grid,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
! Validate D
if(size(D,dim=1)/=grid%nrows .or. size(D,dim=2)/=grid%ncols) then
    err=1;mess=trim(procname)//':Fatal:size mismatch [D,grid]';return
endif

do i=1,grid%nrows
    do j=1,grid%ncols
        pt2(1)=grid%xllcorner+real(j-0.5,mrk)*grid%cellsize
        pt2(2)=grid%yllcorner+real(i-0.5,mrk)*grid%cellsize
        call GetDistance_2(formula=formula,pt1=pt1,pt2=pt2,&
                           D=D(grid%nrows-i+1,j),err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    enddo
enddo

end subroutine GetDistance_1grid

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
!^*     1.formula
!^*     2.pts
!^*     3.[par] optional parameter vector
!^*     4.[cov] optional covariate matrix (row: sites, columns: covariates)
!^* OUT
!^*     1.D, distance matrix
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
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
subroutine validateRasterGrid(grid,err,mess)
!^**********************************************************************
!^* Purpose: Validate raster grid object.
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Created: 03/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. grid, grid definition (rasterGrid object)
!^* OUT
!^*    1.err, error code; ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
type(rasterGridType), intent(in)::grid
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='validateRasterGrid'
logical::ncolOK,nrowOK,sizeOK

err=0;mess=''

ncolOK=(grid%ncols>0)
nrowOK=(grid%nrows>0)
sizeOK=(grid%cellsize>0)

if(ncolOK .and. nrowOK .and. sizeOK) then
    mess='ok'
else
    err=1;mess=trim(procname)//':FATAL'
    if(.not. ncolOK) mess=trim(mess)//' - invalid ncols'
    if(.not. nrowOK) mess=trim(mess)//' - invalid nrows'
    if(.not. sizeOK) mess=trim(mess)//' - invalid cellsize'
endif

end subroutine validateRasterGrid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rasterToPoints(grid,M,pts,err,mess)
!^**********************************************************************
!^* Purpose: Transform a raster dataset in long matrix format with
!^*          columns (x,y,val).
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Created: 03/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. grid, grid definition (rasterGrid object)
!^*    2. [M], data matrix. If not present, third column of "pts" is filled with junk.
!^* OUT
!^*     1.points, (nrow*ncol) x 3 array with columns (x,y,val)
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
type(rasterGridType), intent(in)::grid
real(mrk),intent(in),optional::M(:,:)
real(mrk),allocatable,intent(out)::pts(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='rasterToPoints'
integer(mik)::i,j,k

err=0;mess=''

call validateRasterGrid(grid,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif

allocate(pts(grid%ncols*grid%nrows,3))
pts=undefRN;k=0
do i=1,grid%nrows
    do j=1,grid%ncols
        k=k+1
        pts(k,1)=grid%xllcorner+real(j-0.5,mrk)*grid%cellsize
        pts(k,2)=grid%yllcorner+real(i-0.5,mrk)*grid%cellsize
        if(present(M)) pts(k,3)=M(grid%nrows-i+1,j)
    enddo
enddo

end subroutine rasterToPoints

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine interpolateOnGrid(pts,grid,method,par,M,err,mess)
!^**********************************************************************
!^* Purpose: Interpolate a set of (x,y,value) points on a regular grid
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Created: 03/08/2022
!^**********************************************************************
!^* Comments: For inverse-distance weighting, zero-distances are replaced
!^*           by a small number (1/1000000 of grid resolution).
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. pts, n*3 array with columns (x,y,value) to be interpolated
!^*     2. grid, grid on which interpolated values should be computed (rasterGrid object)
!^*     3. method, interpolation method
!^*     4.[par] optional parameter vector
!^* OUT
!^*     1.M, array containing interpolated values
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:number_string,changeCase
real(mrk),intent(in)::pts(:,:)
type(rasterGridType), intent(in)::grid
character(*),intent(in)::method
real(mrk),intent(in),optional::par(:)
real(mrk), allocatable, intent(out)::M(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='interpolateOnGrid'
real(mrk),parameter::IDWpar_def=1._mrk
integer(mik)::n,p,i
real(mrk),allocatable::distances(:,:,:),sumDist(:,:),weights(:,:,:)
real(mrk)::smallDist,pow

err=0;mess=''

! Validate grid
call validateRasterGrid(grid,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
! Validate points to be interpolated
n=size(pts,dim=1);p=size(pts,dim=2)
if(n<1) then
    mess=trim(procname)//':[points] should have at least 1 row: n='//trim(number_string(n))
    err=1;return
endif
if(p<3) then
    mess=trim(procname)//':[points] should have at least 3 columns: p='//trim(number_string(p))
    err=1;return
endif

allocate(M(grid%nrows,grid%ncols))
M=0._mrk

select case(changeCase(trim(method),'lower'))
case('inversedistance','inversedistanceweighting','idw')
    ! handle parameter
    if(present(par)) then;pow=par(1);else;pow=IDWpar_def;endif
    ! Compute distance between each point and grid
    allocate(distances(grid%nrows,grid%ncols,n),weights(grid%nrows,grid%ncols,n))
    allocate(sumDist(grid%nrows,grid%ncols))
    smallDist=grid%cellsize/1000000._mrk
    do i=1,n
        call GetDistance(formula='euclidean',pt1=pts(i,1:2),grid=grid,&
                         D=distances(:,:,i),err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        where(distances(:,:,i)<=0._mrk) distances(:,:,i)=smallDist
    enddo
    ! Deduce weights
    sumDist=sum(1._mrk/distances**pow,dim=3)
    do i=1,n
        weights(:,:,i)=(1._mrk/distances(:,:,i)**pow)/sumDist
        M=M+weights(:,:,i)*pts(i,3)
    enddo
case default
    err=1;mess=trim(procname)//':Fatal:Unavailable method';return
end select

end subroutine interpolateOnGrid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ReadRasterASC(file,gridOnly,grid,M,err,mess)

!^**********************************************************************
!^* Purpose: read Raster grid in ASC format
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 03/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, address of data file
!^*    2.[gridOnly], only read grid properties? default .false.
!^* OUT
!^*    1.grid, grid definition (rasterGrid object)
!^*    2.[M], data matrix
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:number_string
character(*), intent(in)::file
logical, optional, intent(in)::gridOnly
type(rasterGridType), intent(out)::grid
real(mrk), allocatable, intent(out), optional::M(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='ReadRasterASC'
integer(mik)::unt,j
character(250)::txt
logical::gOnly

err=0;mess=''

if(present(gridOnly)) then;gOnly=gridOnly;else;gOnly=.false.;endif

call getSpareUnit(unt,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file: '//trim(file);return;endif

! Read grid properties
read(unt,*,iostat=err) txt,grid%ncols
if(err>0) then;mess=trim(procname)//':problem reading ncols';close(unt);return;endif
read(unt,*,iostat=err) txt,grid%nrows
if(err>0) then;mess=trim(procname)//':problem reading nrows';close(unt);return;endif
read(unt,*,iostat=err) txt,grid%xllcorner
if(err>0) then;mess=trim(procname)//':problem reading xllcorner';close(unt);return;endif
read(unt,*,iostat=err) txt,grid%yllcorner
if(err>0) then;mess=trim(procname)//':problem reading yllcorner';close(unt);return;endif
read(unt,*,iostat=err) txt,grid%cellsize
if(err>0) then;mess=trim(procname)//':problem reading cellsize';close(unt);return;endif
read(unt,*,iostat=err) txt,grid%mv
if(err>0) then;mess=trim(procname)//':problem reading NODATA_value';close(unt);return;endif

if(.not.gOnly) then
    if(.not.present(M)) then
        mess=trim(procname)//':M should be present unless gridOnly=.true.'
        err=1;return
    endif
    allocate(M(grid%nrows,grid%ncols))
    do j=1, grid%nrows
        read(unt,*,iostat=err) M(j,:)
        if(err>0) then
            mess=trim(procname)//':problem reading data at row:'//trim(number_string(j))
            close(unt);return
        endif
    enddo
endif

close(unt)

end subroutine ReadRasterASC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine WriteRasterASC(grid,M,file,err,mess)

!^**********************************************************************
!^* Purpose: write Raster grid in ASC format
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 03/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. grid, grid definition (rasterGrid object)
!^*    2. M, data matrix
!^*    3. file, address of data file
!^* OUT
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:number_string
type(rasterGridType), intent(in)::grid
real(mrk), intent(in)::M(:,:)
character(*), intent(in)::file
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='WriteRasterASC'
integer(mik)::unt,j
character(250)::fmat

err=0;mess=''

! Check M format
if( size(M,dim=1)/=grid%nrows .or. size(M,dim=2)/=grid%ncols ) then
    mess=trim(procname)//':M dimension incompatible with grid'
    err=1;return
endif

call getSpareUnit(unt,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file: '//trim(file);return;endif

! write grid properties
write(unt,'(G0,4X,I0)',iostat=err) adjustl('ncols'),grid%ncols
if(err>0) then;mess=trim(procname)//':problem writing ncols';close(unt);return;endif
write(unt,'(G0,4X,I0)',iostat=err) 'nrows',grid%nrows
if(err>0) then;mess=trim(procname)//':problem writing nrows';close(unt);return;endif
write(unt,'(G0,4X,e16.8)',iostat=err) 'xllcorner',grid%xllcorner
if(err>0) then;mess=trim(procname)//':problem writing xllcorner';close(unt);return;endif
write(unt,'(G0,4X,e16.8)',iostat=err) 'yllcorner',grid%yllcorner
if(err>0) then;mess=trim(procname)//':problem writing yllcorner';close(unt);return;endif
write(unt,'(G0,4X,e16.8)',iostat=err) 'cellsize',grid%cellsize
if(err>0) then;mess=trim(procname)//':problem writing cellsize';close(unt);return;endif
write(unt,'(G0,4X,e16.8)',iostat=err) 'NODATA_value',grid%mv
if(err>0) then;mess=trim(procname)//':problem writing NODATA_value';close(unt);return;endif

fmat="("//trim(number_string(grid%ncols))//"e16.8)"
do j=1, grid%nrows
    write(unt,fmat,iostat=err) M(j,:)
    if(err>0) then
        mess=trim(procname)//':problem writing data at row:'//trim(number_string(j))
        close(unt);return
    endif
enddo

close(unt)

end subroutine WriteRasterASC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Geodesy_tools

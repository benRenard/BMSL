module Netcdf_tools

!~**********************************************************************
!~* Purpose: Tools for handling netcdf file format
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!~**********************************************************************
!~* Last modified: 12/11/2012
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
public :: NETCDF_load,NETCDF2TimeSeries,NETCDF2TimeSeries_XYPT,&
          NETCDF2TimeSeries_XYT,&
          NETCDF_SubSample

type, public::NETCDFdim
    character(250)::nickname=''
    integer(mik)::n=undefRN ! number of values
    real(mrk),allocatable::val(:) ! values
end type NETCDFdim

type, public::NETCDFType
    integer(mik):: Ndim=undefIN ! number of dimensions
    type(NETCDFdim),allocatable::d(:) ! descriptors for each dimension
    ! botchy use of a 7-shaped array due to the inability to deal with unknown-shape arrays in fortran...
    ! Also note the use of single-precision for the time being, to decrease memory requirements
    real(4), allocatable::field(:,:,:,:,:,:,:) 
end type NETCDFType
Contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NETCDF_load(file,y,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:12/11/2012
!^**********************************************************************
!^* Comments: Only accept singled-variable netcdf for the moment
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. file, address of the netcdf file
!^*		2. 
!^*		3.
!^* OUT
!^*		1.y, object of type NCDFType
!^*		2.err
!^*		3.mess
!^*		4.
!^*		5.
!^**********************************************************************
use netcdf90
character(*), intent(in)::file
type(NETCDFType),intent(out)::y
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::ncid,nVariables,nAttributes,i,FieldDim(7)
real(4), allocatable::S1(:),S2(:,:),S3(:,:,:),S4(:,:,:,:),S5(:,:,:,:,:),S6(:,:,:,:,:,:),S7(:,:,:,:,:,:,:)
integer(mik)::dimids(NF90_MAX_VAR_DIMS)

err=0;mess=''
err=nf90_open(trim(file), 0, ncid)
if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_open';return;endif
err=nf90_Inquire(ncid, y%ndim, nVariables, nAttributes)
if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_Inquire';return;endif
if(nVariables/=y%ndim+1) then;err=1;mess='NETCDF_load:Fatal:multi-variable not allowed';return;endif
if(allocated(y%d)) deallocate(y%d);allocate(y%d(y%ndim))
err = nf90_Inquire_Variable(ncid, y%ndim+1, dimids=dimids)
do i=1,y%ndim
    err=nf90_Inquire_Dimension(ncid, dimids(i), y%d(i)%nickname, y%d(i)%n)
    if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_Inquire_Dimension';return;endif
    if(allocated(y%d(i)%val)) deallocate(y%d(i)%val);allocate(y%d(i)%val(y%d(i)%n))
    err=nf90_get_var(ncid=ncid,varid=dimids(i),values=y%d(i)%val)
    if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_get_var';return;endif
enddo
FieldDim=1
FieldDim(1:y%ndim)=y%d(1:y%ndim)%n
if(allocated(y%field)) deallocate(y%field);allocate(y%field(FieldDim(1),FieldDim(2),FieldDim(3),FieldDim(4),FieldDim(5),FieldDim(6),FieldDim(7)))
select case(y%ndim)
    case(7)
        if(allocated(S7)) deallocate(S7);allocate(S7(FieldDim(1),FieldDim(2),FieldDim(3),FieldDim(4),FieldDim(5),FieldDim(6),FieldDim(7)))
        err=nf90_get_var(ncid=ncid,varid=y%ndim+1,values=S7)
        if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_get_var';return;endif
        y%field(:,:,:,:,:,:,:)=S7
    case(6)
        if(allocated(S6)) deallocate(S6);allocate(S6(FieldDim(1),FieldDim(2),FieldDim(3),FieldDim(4),FieldDim(5),FieldDim(6)))
        err=nf90_get_var(ncid=ncid,varid=y%ndim+1,values=S6)
        if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_get_var';return;endif
        y%field(:,:,:,:,:,:,1)=S6
    case(5)
        if(allocated(S5)) deallocate(S5);allocate(S5(FieldDim(1),FieldDim(2),FieldDim(3),FieldDim(4),FieldDim(5)))
        err=nf90_get_var(ncid=ncid,varid=y%ndim+1,values=S5)
        if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_get_var';return;endif
        y%field(:,:,:,:,:,1,1)=S5
    case(4)
        if(allocated(S4)) deallocate(S4);allocate(S4(FieldDim(1),FieldDim(2),FieldDim(3),FieldDim(4)))
        err=nf90_get_var(ncid=ncid,varid=y%ndim+1,values=S4)
        if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_get_var';return;endif
        y%field(:,:,:,:,1,1,1)=S4
    case(3)
        if(allocated(S3)) deallocate(S3);allocate(S3(FieldDim(1),FieldDim(2),FieldDim(3)))
        err=nf90_get_var(ncid=ncid,varid=y%ndim+1,values=S3)
        if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_get_var';return;endif
        y%field(:,:,:,1,1,1,1)=S3
    case(2)
        if(allocated(S2)) deallocate(S2);allocate(S2(FieldDim(1),FieldDim(2)))
        err=nf90_get_var(ncid=ncid,varid=y%ndim+1,values=S2)
        if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_get_var';return;endif
        y%field(:,:,1,1,1,1,1)=S2
    case(1)
        if(allocated(S1)) deallocate(S1);allocate(S1(FieldDim(1)))
        err=nf90_get_var(ncid=ncid,varid=y%ndim+1,values=S1)
        if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_get_var';return;endif
        y%field(:,1,1,1,1,1,1)=S1
    case default
        err=2;mess='NETCDF_load:Fatal:y%ndim not in 1:7';return
end select
err=nf90_close(ncid)
if(err/=0) then;mess='NETCDF_load:Fatal error:nf90_close';return;endif
end subroutine NETCDF_load
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NETCDF_SubSample(y,step,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:12/11/2012
!^**********************************************************************
!^* Comments: Only accept singled-variable netcdf for the moment
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. file, address of the netcdf file
!^*		2. 
!^*		3.
!^* OUT
!^*		1.y, object of type NCDFType
!^*		2.err
!^*		3.mess
!^*		4.
!^*		5.
!^**********************************************************************
type(NETCDFType),intent(INOUT)::y
integer(mik),intent(in)::step(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
type(NETCDFType)::foo
integer(mik)::n,i,k,step7(7),k7(7),j

err=0;mess=''
n=size(step)
foo%ndim=y%ndim
if(allocated(foo%d)) deallocate(foo%d);allocate(foo%d(foo%ndim))
! sub-sample each dimension
k7=1;k7(1:y%ndim)=y%d(:)%n
do i=1,min(n,y%Ndim)
    k=size((/(j,j=1,y%d(i)%n,step(i))/));k7(i)=k
    if(allocated(foo%d(i)%val)) deallocate(foo%d(i)%val);allocate(foo%d(i)%val(k))
    foo%d(i)%val=y%d(i)%val(1::step(i))
    deallocate(y%d(i)%val);allocate(y%d(i)%val(k));y%d(i)%val=foo%d(i)%val
    y%d(i)%n=k
enddo
! sub-sample field
! first fill-in step if necessary
step7=1;step7(1:min(n,y%Ndim))=step(1:min(n,y%Ndim))
! create field
if(allocated(foo%field)) deallocate(foo%field)
allocate(foo%field(k7(1),k7(2),k7(3),k7(4),k7(5),k7(6),k7(7)))
foo%field=y%field(1::step7(1),1::step7(2),1::step7(3),1::step7(4),1::step7(5),1::step7(6),1::step7(7))
deallocate(y%field);allocate(y%field(k7(1),k7(2),k7(3),k7(4),k7(5),k7(6),k7(7)))
y%field=foo%field
end subroutine NETCDF_SubSample
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NETCDF2TimeSeries(y,TimeStep,FirstDate,Tname,dates,X,err,mess)
!^**********************************************************************
!^* Purpose: transform a NETCDF object to a big matrix with rows corresponding to time
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:12/11/2012
!^**********************************************************************
!^* Comments: Not properly and thoroughly checked, unsafe... 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. y, object of type NCDFType
!^*		2. TimeStep in the original netcdf file ('daily' or 'monthly')
!^*		3. FirstDate in the original netcdf file
!^*		4. [Tname], character identifying the time dimension - default 'T'
!^* OUT
!^*		1. dates, time series of dates corresponding to each row of X 
!^*		2. X, big matrix
!^*		3. err
!^*		4. mess
!^**********************************************************************
use Dates_tools
type(NETCDFType),intent(in)::y
character(*), intent(in)::TimeStep
type(DateType),intent(in)::FirstDate
character(*), intent(in), optional::Tname
type(DateType), pointer::dates(:)
real(mrk), pointer::X(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
character(250), parameter::TN_def='T'
character(250)::TN
integer(mik)::i,Tindx,n,p,order(y%ndim)

err=0;mess=''
! Handle optional variables
if(present(Tname)) then;TN=Tname;else;TN=TN_def;endif
! Get time dimension
Tindx=UndefIN
do i=1,y%ndim
    if(trim(y%d(i)%nickname)==trim(TN)) then
        Tindx=i;exit
    endif
enddo
if(Tindx==UndefIN) then;err=1;mess='NETCDF2TimeSeries:Fatal:Time dimension not found in y.';return;endif
! Get sizes
n=y%d(Tindx)%n
p=product(shape(y%field))/n
if(associated(X)) deallocate(X);allocate(X(n,p))
if(associated(dates)) deallocate(dates);allocate(dates(n))
!order=(/(i,i=1,y%ndim)/)
!order(y%ndim)=Tindx;order(Tindx)=y%ndim
X=reshape(y%field,(/n,p/),order=(/2,1/))
select case(trim(TimeStep))
    case('daily')
        dates=DailyList(FirstDate,InNDays(FirstDate,n-1))
    case('monthly')
        dates(1)=FirstDate
        do i=2,n
            dates(i)=InOneMonth(dates(i-1))
        enddo
    case default
        err=2;mess='NETCDF2TimeSeries:Fatal:Unknown TimeStep';return
end select
end subroutine NETCDF2TimeSeries
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NETCDF2TimeSeries_XYPT(y,TimeStep,FirstDate,dates,lon,lat,X,err,mess)
!^**********************************************************************
!^* Purpose: transform a NETCDF object to a big matrix with rows corresponding to time
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:12/11/2012
!^**********************************************************************
!^* Comments: Assumes dimensions are in the order: X Y P T (standard for NCEP export of geop. height) 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. y, object of type NCDFType
!^*		2. TimeStep in the original netcdf file ('daily' or 'monthly')
!^*		3. FirstDate in the original netcdf file
!^* OUT
!^*		1. dates, time series of dates corresponding to each row of X 
!^*		1. lon, longitude variation along columns of X 
!^*		1. lat, latitude variation along columns of X  
!^*		2. X, big matrix
!^*		3. err
!^*		4. mess
!^**********************************************************************
use Dates_tools
type(NETCDFType),intent(in)::y
character(*), intent(in)::TimeStep
type(DateType),intent(in)::FirstDate
type(DateType), pointer::dates(:)
real(mrk), pointer::lat(:),lon(:)
real(mrk), pointer::X(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::i,j,n,p,nlat,nlon

err=0;mess=''
! Get sizes
n=y%d(4)%n
p=product(shape(y%field))/n
if(associated(X)) deallocate(X);allocate(X(n,p))
if(associated(lat)) deallocate(lat);allocate(lat(p))
if(associated(lon)) deallocate(lon);allocate(lon(p))
if(associated(dates)) deallocate(dates);allocate(dates(n))
nlat=y%d(2)%n;nlon=y%d(1)%n
do i=1,n
    do j=1,nlat
        X(i,((j-1)*nlon+1):(j*nlon))=y%field(:,j,1,i,1,1,1)
        lat(((j-1)*nlon+1):(j*nlon))=y%d(2)%val(j)
        lon(((j-1)*nlon+1):(j*nlon))=y%d(1)%val(:)
    enddo
enddo
select case(trim(TimeStep))
    case('daily')
        dates=DailyList(FirstDate,InNDays(FirstDate,n-1))
    case('monthly')
        dates(1)=FirstDate
        do i=2,n
            dates(i)=InOneMonth(dates(i-1))
        enddo
    case default
        err=2;mess='NETCDF2TimeSeries_XYPT:Fatal:Unknown TimeStep';return
end select
end subroutine NETCDF2TimeSeries_XYPT
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine NETCDF2TimeSeries_XYT(y,TimeStep,FirstDate,dates,lon,lat,X,err,mess)
!^**********************************************************************
!^* Purpose: transform a NETCDF object to a big matrix with rows corresponding to time
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:14/11/2012
!^**********************************************************************
!^* Comments: Assumes dimensions are in the order: X Y T (standard for NCEP export of eg surface temp) 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. y, object of type NCDFType
!^*		2. TimeStep in the original netcdf file ('daily' or 'monthly')
!^*		3. FirstDate in the original netcdf file
!^* OUT
!^*		1. dates, time series of dates corresponding to each row of X 
!^*		1. lon, longitude variation along columns of X 
!^*		1. lat, latitude variation along columns of X  
!^*		2. X, big matrix
!^*		3. err
!^*		4. mess
!^**********************************************************************
use Dates_tools
type(NETCDFType),intent(in)::y
character(*), intent(in)::TimeStep
type(DateType),intent(in)::FirstDate
type(DateType), pointer::dates(:)
real(mrk), pointer::lat(:),lon(:)
real(mrk), pointer::X(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::i,j,n,p,nlat,nlon

err=0;mess=''
! Get sizes
n=y%d(3)%n
p=product(shape(y%field))/n
if(associated(X)) deallocate(X);allocate(X(n,p))
if(associated(lat)) deallocate(lat);allocate(lat(p))
if(associated(lon)) deallocate(lon);allocate(lon(p))
if(associated(dates)) deallocate(dates);allocate(dates(n))
nlat=y%d(2)%n;nlon=y%d(1)%n
do i=1,n
    do j=1,nlat
        X(i,((j-1)*nlon+1):(j*nlon))=y%field(:,j,i,1,1,1,1)
        lat(((j-1)*nlon+1):(j*nlon))=y%d(2)%val(j)
        lon(((j-1)*nlon+1):(j*nlon))=y%d(1)%val(:)
    enddo
enddo
select case(trim(TimeStep))
    case('daily')
        dates=DailyList(FirstDate,InNDays(FirstDate,n-1))
    case('monthly')
        dates(1)=FirstDate
        do i=2,n
            dates(i)=InOneMonth(dates(i-1))
        enddo
    case default
        err=2;mess='NETCDF2TimeSeries_XYT:Fatal:Unknown TimeStep';return
end select
end subroutine NETCDF2TimeSeries_XYT


end module Netcdf_tools

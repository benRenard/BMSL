module ExtraFlo_utilities

!~**********************************************************************
!~* Purpose: 
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 25/03/2012
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

public :: ReadCVlist,ReadCovFile,ReadXYFile,ReadParFile,ReadAmax,&
          PackAmaxByRegion,PackAmax,PackAmaxByDecomposition,ReNumberRegion
! derived types
type,public:: partype
  real(mrk)::pmean,pmed,pstd,pmode
end type partype
type,public:: AtSitePar_type
  integer(mik)::n
  real(mrk)::mean,med
  type(partype)::loc,scale,shape
end type AtSitePar_type
type, public::VoVtype
    real(mrk), pointer::val(:) =>NULL()
end type VoVtype

Contains
subroutine ReadCVlist(file,n,liste,err,mess)
  character(*), intent(in)::file
  integer(mik), intent(out)::n
  character(*), pointer::liste(:)
  integer(mik), intent(out)::err
  character(250), intent(out)::mess
  !locals
  integer(mik)::errcode,i

  err=0;mess='';n=undefIN;
  
  open(unit=1, file=trim(file), status='OLD', iostat=errcode)
  If(errcode/=0) then
    err=1;mess="problem reading file "//trim(file);return
  endif
  errcode=0;n=0
  do while(errcode==0)
	  n=n+1
	  read(1,*,iostat=errcode)
  enddo
  n=n-1
  if(associated(liste)) deallocate(liste)
  allocate(liste(n))
  rewind(1)
  do i=1,n
    read(1,*,iostat=errcode) liste(i)
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(file);return
    endif
  enddo
  close(1)
end subroutine ReadCVlist
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ReadCovFile(liste,folder,region,cov,err,mess)
  character(*), intent(in)::folder,liste(:),region
  real(mrk), intent(out)::cov(:,:)
  integer(mik), intent(out)::err
  character(250), intent(out)::mess
  !locals
  integer(mik)::errcode,i,n
  
  n=size(liste)
  do i=1,n
    open(unit=1, file=trim(folder)//trim(liste(i))//'_'//trim(region)//'_cov.txt', status='OLD', iostat=errcode)
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    read(1,*,iostat=errcode) ! headers
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    read(1,*,iostat=errcode) cov(i,:)
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    close(1)
  enddo
end subroutine ReadCovFile
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ReadXYFile(liste,folder,xy,err,mess)
  character(*), intent(in)::folder,liste(:)
  real(mrk), intent(out)::xy(:,:)
  integer(mik), intent(out)::err
  character(250), intent(out)::mess
  !locals
  integer(mik)::errcode,i,n
  
  n=size(liste)
  do i=1,n
    open(unit=1, file=trim(folder)//trim(liste(i))//'_xy.txt', status='OLD', iostat=errcode)
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    read(1,*,iostat=errcode) ! headers
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    read(1,*,iostat=errcode) xy(i,:)
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    close(1)
  enddo
end subroutine ReadXYFile
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ReadParFile(liste,folder,distID,AtSitePar,err,mess)
  character(*), intent(in)::folder,liste(:),distID
  type(AtSitePar_type)::AtSitePar(:)
  integer(mik), intent(out)::err
  character(250), intent(out)::mess
  !locals
  integer(mik)::errcode,i,n
  real(mrk)::crapR1,crapR2,crapR3,crapR4
  character(250)::crapC1, crapC2, crapC3, crapC4,dist
  
  n=size(liste)
  do i=1,n
    open(unit=1, file=trim(folder)//trim(liste(i))//'_'//trim(distID)//'_par.txt', status='OLD', iostat=errcode)
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    read(1,*,iostat=errcode) ! headers
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    read(1,*,iostat=errcode) crapC1, AtSitePar(i)%n,AtSitePar(i)%mean,&
                             AtSitePar(i)%med,AtSitePar(i)%loc%pmean,&
                             AtSitePar(i)%scale%pmean,AtSitePar(i)%shape%pmean
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    read(1,*,iostat=errcode) crapC1, AtSitePar(i)%n,AtSitePar(i)%mean,&
                             AtSitePar(i)%med,AtSitePar(i)%loc%pstd,&
                             AtSitePar(i)%scale%pstd,AtSitePar(i)%shape%pstd
    If(errcode/=0) then
      err=1;mess="problem reading file "//trim(liste(i));return
    endif
    close(1)
  enddo
end subroutine ReadParFile
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ReadAmax(file,Amax,pmv,year,err,mess)
use utilities_dmsl_kit, only:number_string,GetSpareUnit
character(*), intent(in)::file
real(mrk), pointer::Amax(:),pmv(:)
integer(mik), pointer,optional::year(:)
character(*), intent(out)::mess
integer(mik), intent(out)::err
!locals
integer(mik)::errcode,n,k,y,unt,i
character(250)::date

err=0;mess=""
call GetSpareUnit(unt,err,mess)
open(unit=unt, file=file, status='OLD',iostat=err)
if(err/=0) return
errcode=0;n=0
read(unt,*) !header
do while(errcode==0)
    n=n+1
    read(unt,*,iostat=errcode)
enddo
n=n-1
if(associated(Amax)) deallocate(Amax)
allocate(Amax(n))
if(associated(pmv)) deallocate(pmv)
allocate(pmv(n))
if(present(year)) then
    if(associated(year)) deallocate(year)
    allocate(year(n))
endif

rewind(unt)
read(unt,*) !header
do i=1,n
    read(unt,*,iostat=errcode) k,Amax(i),pmv(i),y,date
    if(present(year)) year(i)=y
    if(errcode/=0) then
        err=-1;mess="ReadAmax:FATAL:line"//trim(number_string(i));return
    endif
enddo
close(unt)
end subroutine ReadAmax
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine PackAmaxByRegion(Y,cov,list,CatchmentSize,regioncol,pmvMax,&
                    path2Amax,LogSpace,packedAmax,err,mess)
! IN
real(mrk),intent(in)::Y(:),cov(:,:),CatchmentSize(:),pmvMax
integer(mik),intent(in)::regioncol
character(*), intent(in)::path2Amax,list(:)
logical, intent(in)::LogSpace
! OUT
integer(mik),intent(out)::err
character(*),intent(out)::mess
type(VoVtype),pointer::packedAmax(:)
!locals
integer(mik), parameter::BigDummyN=100000
integer(mik)::reg,nregion,p,foo(1),NRegionalDist,site,Ndata,N
logical::mask(size(cov,dim=1))
real(mrk), allocatable::Ymasked(:),CatchmentSize_masked(:)
character(250),allocatable::listmasked(:)
real(mrk), pointer::Amax(:),Amax_raw(:),pmv(:)
logical,allocatable::mask2(:)
real(mrk)::NormalizedData(BigDummyN)

foo=maxval(nint(cov(:,regioncol),mik))
nregion=foo(1)
if(associated(packedAmax)) deallocate(packedAmax);allocate(packedAmax(nregion))
do reg=1,nregion
    ! Select stations in current region
    mask=cov(:,regioncol)==reg
    p=count(mask)
    if(p==0) then
        write(*,*) "This region is EMPTY!!!!"
        write(*,*), "region:",reg
        write(*,*) "It will be skipped but subsequent computations may become problematic..."
        CYCLE
    endif
    if(allocated(Ymasked)) deallocate(Ymasked);allocate(Ymasked(p))
    if(allocated(listmasked)) deallocate(listmasked);allocate(listmasked(p))
    if(allocated(CatchmentSize_masked)) deallocate(CatchmentSize_masked);allocate(CatchmentSize_masked(p))
    Ymasked=pack(Y,mask)
    listmasked=pack(list,mask)
    CatchmentSize_masked=pack(CatchmentSize,mask)
    ! Aggregate all non-missing data from this region
    NRegionalDist=0
    do site=1,p
        call ReadAmax(file=trim(path2Amax)//trim(listmasked(site))//"_ann.txt",&
                      Amax=Amax_raw,pmv=pmv,err=err,mess=mess)
        if(err>0) then
            mess='PackAmax'//trim(mess);return
        endif
        ! Handle MV
        Ndata=size(pmv)
        if(allocated(mask2)) deallocate(mask2)
        allocate(mask2(Ndata))
        mask2 = (pmv<=pmvMax)
        N=count(mask2)
        if(associated(Amax)) deallocate(Amax)
        allocate(Amax(N))
        Amax=pack(Amax_raw,mask2)
        ! Convert in mm
        Amax=(Amax*86.4_mrk)/CatchmentSize_masked(site)
        if(LogSpace) then
            NormalizedData((NRegionalDist+1):(NRegionalDist+N))=log(Amax)/Ymasked(site) !-Ymasked(site)
        else
            NormalizedData((NRegionalDist+1):(NRegionalDist+N))=Amax/Ymasked(site)
        endif
        NRegionalDist=NRegionalDist+N
    enddo
    if(associated(packedAmax(reg)%val)) deallocate(packedAmax(reg)%val)
    allocate(packedAmax(reg)%val(NRegionalDist))
    packedAmax(reg)%val=NormalizedData(1:NRegionalDist)
enddo
end subroutine PackAmaxByRegion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine PackAmax(list,CatchmentSize,pmvMax,&
                    path2Amax,LogSpace,packedAmax,err,mess)
! IN
real(mrk),intent(in)::CatchmentSize(:),pmvMax
character(*), intent(in)::path2Amax,list(:)
logical, intent(in)::LogSpace
! OUT
integer(mik),intent(out)::err
character(*),intent(out)::mess
type(VoVtype),pointer::packedAmax(:)
!locals
integer(mik)::site,Nsite,Ndata,N
logical,allocatable::mask(:)
real(mrk), pointer::Amax_raw(:),pmv(:)

Nsite=size(CatchmentSize)
if(associated(packedAmax)) deallocate(packedAmax);allocate(packedAmax(Nsite))
do site=1,Nsite
     call ReadAmax(file=trim(path2Amax)//trim(list(site))//"_ann.txt",&
                   Amax=Amax_raw,pmv=pmv,err=err,mess=mess)
     if(err>0) then
         mess='PackAmax'//trim(mess);return
     endif
     ! Handle MV
     Ndata=size(pmv)
     if(allocated(mask)) deallocate(mask)
     allocate(mask(Ndata))
     mask = (pmv<=pmvMax)
     N=count(mask)
     if(associated(packedAmax(site)%val)) deallocate(packedAmax(site)%val)
     allocate(packedAmax(site)%val(N))
     packedAmax(site)%val=pack(Amax_raw,mask)
     ! Convert in mm
     packedAmax(site)%val=(packedAmax(site)%val*86.4_mrk)/CatchmentSize(site)
     if(LogSpace) packedAmax(site)%val=log( packedAmax(site)%val)
enddo
end subroutine PackAmax
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine PackAmaxByDecomposition(list,CatchmentSize,pmvMax,mvcode,&
                    path2Decomposition,path2Amax,LogSpace,packedAmax1,packedAmax2,&
                    suffix1,suffix2,err,mess)
! IN
real(mrk),intent(in)::CatchmentSize(:),pmvMax,mvcode
character(*), intent(in)::path2Amax,path2Decomposition,list(:),suffix1,suffix2
logical, intent(in)::LogSpace
! OUT
integer(mik),intent(out)::err
character(*),intent(out)::mess
type(VoVtype),pointer::packedAmax1(:),packedAmax2(:)
!locals
integer(mik)::site,Nsite,Ndata
logical,allocatable::mask(:)
real(mrk), pointer::Amax_raw(:),pmv(:)
real(mrk), allocatable::Amax_raw1(:),Amax_raw2(:),temp1(:),temp2(:)
logical, allocatable::mask1(:),mask2(:)
integer(mik), pointer::year(:),year1(:),year2(:)
real(mrk)::dummy(2050)

Nsite=size(CatchmentSize)
if(associated(packedAmax1)) deallocate(packedAmax1);allocate(packedAmax1(Nsite))
if(associated(packedAmax2)) deallocate(packedAmax2);allocate(packedAmax2(Nsite))
do site=1,Nsite
     ! Read At-site Amax
     call ReadAmax(trim(path2Amax)//trim(list(site))//"_ann.txt",&
                   Amax_raw,pmv,year,err,mess)
     if(err>0) then
         mess='PackAmaxByDecomposition'//trim(mess);return
     endif
     ! handle MV, part 1
     Ndata=size(pmv)
     if(allocated(mask)) deallocate(mask)
     allocate(mask(Ndata))
     mask = (pmv>pmvMax)
     where(mask) Amax_raw=mvcode
     ! Read decomposition
     call ReadDecomposition(file=trim(path2Decomposition)//trim(list(site)),&
                            suffix1=suffix1,suffix2=suffix2,&
                            year1=year1,year2=year2,err=err,mess=mess)
     if(err>0) then
         mess='PackAmaxByDecomposition'//trim(mess);return
     endif
     ! Finalize
     if(allocated(temp1)) deallocate(temp1);allocate(temp1(size(year1)))
     if(allocated(temp2)) deallocate(temp2);allocate(temp2(size(year2)))
     if(allocated(mask1)) deallocate(mask1);allocate(mask1(size(year1)))
     if(allocated(mask2)) deallocate(mask2);allocate(mask2(size(year2)))
     dummy(year)=Amax_raw
     temp1=dummy(year1)
     temp2=dummy(year2)  
     ! Remove MV
     mask1=(temp1/=mvcode)
     mask2=(temp2/=mvcode)
     if(associated(packedAmax1(site)%val)) deallocate(packedAmax1(site)%val)
     allocate(packedAmax1(site)%val(count(mask1)))
     if(associated(packedAmax2(site)%val)) deallocate(packedAmax2(site)%val)
     allocate(packedAmax2(site)%val(count(mask2)))
     packedAmax1(site)%val=(pack(temp1,mask1)*86.4_mrk)/CatchmentSize(site)  
     packedAmax2(site)%val=(pack(temp2,mask2)*86.4_mrk)/CatchmentSize(site)  
     if(LogSpace) then
        packedAmax1(site)%val=log(packedAmax1(site)%val) 
        packedAmax2(site)%val=log(packedAmax2(site)%val)
     endif
enddo
end subroutine PackAmaxByDecomposition 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ReadDecomposition(file,suffix1,suffix2,year1,year2,err,mess)
use utilities_dmsl_kit, only:GetSpareUnit,number_string
! IN
character(*), intent(in)::file,suffix1,suffix2
integer(mik), pointer::year1(:),year2(:)
! OUT
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::unt1,unt2,errcode,n1,n2,i

call GetSpareUnit(unt1,err,mess)
open(unit=unt1, file=trim(file)//trim(suffix1)//'.txt', status='OLD',iostat=err)
if(err/=0) return
call GetSpareUnit(unt2,err,mess)
open(unit=unt2, file=trim(file)//trim(suffix2)//'.txt', status='OLD',iostat=err)
if(err/=0) return
errcode=0;n1=0
do while(errcode==0)
    n1=n1+1
    read(unt1,*,iostat=errcode)
enddo
n1=n1-1
errcode=0;n2=0
do while(errcode==0)
    n2=n2+1
    read(unt2,*,iostat=errcode)
enddo
n2=n2-1
if(associated(year1)) deallocate(year1);allocate(year1(n1))
if(associated(year2)) deallocate(year2);allocate(year2(n2))
rewind(unt1);rewind(unt2)
do i=1,n1
    read(unt1,*,iostat=errcode) year1(i)
    if(errcode/=0) then
        err=-1;mess="ReadDecomposition:FATAL:line"//trim(number_string(i));return
    endif
enddo
do i=1,n2
    read(unt2,*,iostat=errcode) year2(i)
    if(errcode/=0) then
        err=-1;mess="ReadDecomposition:FATAL:line"//trim(number_string(i));return
    endif
enddo
end subroutine ReadDecomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ReNumberRegion(region)
real(mrk),intent(inout)::region(:)
!locals
integer(mik)::nreg,k,i,m
real(mrk)::foo(1),dummy(size(region))

foo=maxval(region)
nreg=nint(foo(1))
k=0;dummy=undefRN
do i=1,nreg
    m=count(nint(region)==i)
    if(m==0) then
        cycle
    else
        k=k+1
        where(nint(region)==i) dummy=real(k,mrk)
    endif
enddo
region=dummy
end subroutine ReNumberRegion
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end module ExtraFlo_utilities

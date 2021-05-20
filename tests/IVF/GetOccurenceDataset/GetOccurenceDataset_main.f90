program GetOccurenceDataset_main

use kinds_dmsl_kit
use Dates_tools
use DataRW_tools, only:ExtraRead,WriteSeparatedFile,ReadRHN
use TimeSeries_tools,only:OneDTimeSeriesType,FillTheGap_daily
use HydroIndices_tools,only:AnnualIndices_daily
use utilities_dmsl_kit, only:getNumItemsInFile
use EMpiricalStats_tools, only:GetEmpiricalQuantile
implicit none

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
character(250),parameter::forma="EXTRA" ! "RHN" or "EXTRA"
character(250),parameter::dir_data="C:\BEN\DATA\Q_France209\"! 'C:\BEN\DATA\Q_RHN3\Output_files\Floods\' !,
character(250),parameter::dir_results=""
character(250),parameter::xtractWhat="Amax",threshold_type="quantile" !"absolute" ! 
type(DateType), parameter:: YearStart=DateType(1,10,1,0,0,0),YearEnd=DateType(1,12,31,0,0,0)
real(mrk),parameter::threshold=0.5,pmvmax=15.
integer(mik),parameter::nmin=30 !40 
character(250),parameter::file_list="FileList.txt",file_xy="XY.txt"
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
character(250),parameter::result_xy="result_XY.txt",result_occ="result_data.txt",result_years="result_years.txt"
character(1),parameter::sep=";"
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! locals
type(OneDTimeSeriesType)::PJ,PJ0,Qcode
integer(mik)::err,i,j,nsite,neff,k,y0,yn,ney,compt
character(250)::mess,jchar
character(250),allocatable::files(:)
real(mrk), pointer::xtract(:),Pmv(:)
real(mrk)::th
integer(mik), pointer::YearList(:)
integer(mik), allocatable::dataset(:,:),years(:),finaldataset(:,:),finalyears(:)
character(2500)::foo
integer(mik)::foo2(1)
logical,allocatable::mask(:)
type::OccType    
    integer(mik)::year,val
end type OccType
type::OccSeries
    character(250)::name 
    integer(mik)::n   
    type(OccType),allocatable::series(:)
end type OccSeries
type(OccSeries),allocatable::occy(:)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! preliminaries
open(unit=1,file=trim(dir_data)//trim(file_list),status='old')
call getNumItemsInFile(unt=1,preRwnd=.true.,nskip=0,nitems=nsite,postPos=0,jchar=jchar,err=err,message=mess)
allocate(files(nsite),occy(nsite))
do i=1,nsite
    ! Read files
    read(1,*) files(i)
enddo
close(1)
if(trim(file_xy)/="") then
    open(unit=666,file=trim(dir_data)//trim(file_xy),status='old')
    read(666,'(A2500)') ! header
    open(unit=667,file=trim(result_xy),status='replace')
endif

k=0
do i=1,nsite
    write(*,*) "site:",i,"/",nsite
	! Read XY
	if(trim(file_xy)/="") then; read(666,'(A2500)') foo;endif
    ! Read PJ
    if(trim(forma)=="RHN") then
        Call ReadRHN(trim(dir_data)//trim(files(i)),PJ0,Qcode,err,mess)
    else
        call ExtraRead(trim(dir_data)//trim(files(i)),PJ0,Qcode,err,mess)
    endif
    deallocate(Qcode%ts)
    call FillTheGap_daily(xin=PJ0,xout=PJ,err=err,mess=mess)
    deallocate(PJ0%ts)
    ! extract desired variable
    select case(trim(xtractWhat))
    case('Amax')
        call AnnualIndices_daily(xin=PJ,&
				YearStart=YearStart,YearEnd=YearEnd, & 
				AMax=xtract,Pmv=Pmv,YearList=YearList,&
				err=err,mess=mess)
	end select
	! kill years with too many MV
	where(Pmv>pmvmax) xtract=PJ%mv
	neff=count(xtract/=PJ%mv)
	if(neff<nmin) then
	    write(*,*) "rejected: neff=",neff;cycle
	endif
	! Save XY
	if(trim(file_xy)/="") then; write(667,'(A2500)') foo;endif
	! get occurrences
	k=k+1
	occy(k)%n=size(xtract)
	occy(k)%name=files(i)
	allocate(occy(k)%series(size(YearList)))
	occy(k)%series(:)%year=YearList
	if(threshold_type=="quantile") then
	    call GetEmpiricalQuantile(p=threshold,x=pack(xtract,xtract/=PJ%mv),q=th,err=err,mess=mess)
	else
	    th=threshold
	endif
	
	where(xtract>th) occy(k)%series(:)%val=1
	where(xtract<=th) occy(k)%series(:)%val=0
	where(xtract==PJ%mv) occy(k)%series(:)%val=PJ%mv	
	deallocate(PJ%ts)
enddo
if(trim(file_xy)/="") then;close(666);close(667);endif
write(*,*) "n usable site:",k,"/",nsite

! Squarify dataset
! get 1st and last year
y0=6666;yn=-6666
do i=1,k
    if(occy(i)%series(1)%year<y0) y0=occy(i)%series(1)%year
    if(occy(i)%series(occy(i)%n)%year>yn) yn=occy(i)%series(occy(i)%n)%year    
enddo
! Squarify
allocate(dataset(yn-y0+1,k),years(yn-y0+1))
dataset=PJ%mv
years=(/(y0+i-1,i=1,yn-y0+1)/)
do i=1,k
    if(allocated(mask)) deallocate(mask)
    allocate(mask(occy(i)%n))
    do j=1,yn-y0+1
        mask=occy(i)%series(:)%year==years(j)
        if(count(mask)>0) then
            foo2=pack(occy(i)%series(:)%val,mask)
            dataset(j,i)=foo2(1)
        endif
    enddo
enddo
! remove empty years
ney=0
do j=1,yn-y0+1
    if(count(dataset(j,:)/=PJ%mv)>0) ney=ney+1
enddo
allocate(finaldataset(ney,k),finalyears(ney))
compt=0
do j=1,yn-y0+1
    if(count(dataset(j,:)/=PJ%mv)>0) then
         compt=compt+1
         finaldataset(compt,:)=dataset(j,:)
         finalyears(compt)=years(j)
    endif
enddo
! Write 2 file
call WriteSeparatedFile(file=trim(result_occ),sep=sep,y=transpose(finaldataset),err=err,mess=mess)
call WriteSeparatedFile(file=trim(result_years),sep=sep,y=reshape(finalyears,(/ney,1/)),err=err,mess=mess)

end program GetOccurenceDataset_main
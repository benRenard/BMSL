program GetSeasonalStats_main

use kinds_dmsl_kit
use DataRW_tools, only:LoadFromRepository,ExtraRead
use TimeSeries_tools,only:OneDTimeSeriesType,GetMonthlyInterannualStat
implicit none

character(250),parameter::dir_data="C:\BEN\DATA\P_Chine\data_copy\",dir_results=""
character(250),parameter::file_list="FileList.txt"
! locals
type(OneDTimeSeriesType), pointer::y(:),Qcode(:)
integer(mik)::err,i,nsite
character(250)::mess
real(mrk),allocatable::WhereAmax(:,:),mean(:,:),sdev(:,:)

call LoadFromRepository(dir_data,trim(dir_data)//trim(file_list),"",ExtraRead,y,Qcode,err,mess)
nsite=size(y)
allocate(WhereAmax(nsite,12),mean(nsite,12),sdev(nsite,12))
do i=1,nsite
    call GetMonthlyInterannualStat(y(i),"mean",mean(i,:),err,mess)
    call GetMonthlyInterannualStat(y(i),"sdev",sdev(i,:),err,mess)
    call GetMonthlyInterannualStat(y(i),"WhereAmax",WhereAmax(i,:),err,mess)
enddo

open(unit=1,file="Seasonal_mean.txt",status="replace")
open(unit=2,file="Seasonal_sdev.txt",status="replace")
open(unit=3,file="Seasonal_whereAmax.txt",status="replace")
do i=1,nsite
    write(1,'(12e14.6)') mean(i,:)
    write(2,'(12e14.6)') sdev(i,:)
    write(3,'(12e14.6)') WhereAmax(i,:)
enddo
close(1);close(2);close(3)

end program GetSeasonalStats_main
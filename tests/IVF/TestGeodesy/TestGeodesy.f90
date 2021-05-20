program TestGeodesy
use kinds_dmsl_kit
use Geodesy_tools

implicit none

integer(mik),parameter::n=7
real(mrk)::cities(n,2)
real(mrk)::D1(n,n), D2(n,n) ! distance matrices
integer(mik)::err,i
character(250)::mess

cities(1,:)=(/138.600739,-34.928497/) ! Adelaide
cities(2,:)=(/151.781677,-32.928272/) ! Newy
cities(3,:)=(/-74.005974,40.712776/) ! NYC
cities(4,:)=(/360-74.005974,40.712776/) ! NYC, different longitude convention
cities(5,:)=(/4.835659,45.764042/) ! Lyon
cities(6,:)=(/5.374150,43.296010/) ! Marseille Vieux Port
cities(7,:)=(/5.418680,43.211390/) ! Marseille Sormiou

call GetDistance(formula='Euclidean',pts=cities,D=D1,err=err,mess=mess)
do i=1,n
    write(*,'(<i>F10.3)') D1(i,1:i)
enddo

call GetDistance(formula='Haversine',pts=cities,D=D2,err=err,mess=mess)
do i=1,n
    write(*,'(<i>F10.3)') D2(i,1:i)
enddo
end program TestGeodesy
program TestTextFileModel

use kinds_dmsl_kit
use TextFile_model

implicit none 

character(250),parameter:: file='MyModel.txt'
integer(mik),parameter::nOBS=200
integer(mik)::err,i
character(250)::mess
logical::feas(nOBS)
real(mrk)::IN(nOBS,3),OUT(nOBS,2),theta(4)

call TxtMdl_load(file,err,mess)
write(*,*) '----------------------------------'
write(*,*) '1. LOADING'
write(*,*) 'Error: ',err
write(*,*) 'Message: ',trim(mess)

IN(:,1)=(/(i,i=1,nOBS)/)
IN(:,2)=(/(18.+4.*sin(2.*3.14*i/12),i=1,nOBS)/)
IN(:,3)=(/(12.+1.*cos(2.*3.14*i/12),i=1,nOBS)/)
theta=(/100.,0.005,1000.,2000./)
call TxtMdl_Apply(IN,theta,OUT,feas,err,mess)
write(*,*) '----------------------------------'
write(*,*) '1. COMPUTING'
write(*,*) 'Error: ',err
write(*,*) 'Message: ',trim(mess)
do i=1,size(OUT,dim=2)
    write(*,*) 'Result: ',OUT(:,i)
enddo
write(*,*) 'Feasability: ',feas

read(*,*)

end program TestTextFileModel
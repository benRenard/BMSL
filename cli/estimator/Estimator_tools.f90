module Estimator_tools

use kinds_dmsl_kit
implicit none
private
public:: CRead_data, CRead_Inference, CRead_MCMC,CRead_ResultFiles,&
        CRead_ResultOptions,RemoveMV,ApplyOption2,CRead_SystematicErrors,&
        CRead_ResultFiles_HBay,CRead_GetQT,&
        ! deprecated
        CRead_data_OLD,CRead_ResultOptions_OLD

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine CRead_data(file,DataFile,DataType,err)
                        !,mvOption,O1_value,O2_cond,O2_value,err)

character(*), intent(in)::file
character(*), intent(out)::DataFile,DataType
!character(*), intent(out)::mvOption,O2_cond
!real(mrk), intent(out)::O1_value,O2_value
integer(mik), intent(out)::err

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return
read(1,'(A250)',iostat=err) DataFile
if(err/=0) return
read(1,'(A250)',iostat=err) DataType
if(err/=0) return
!read(1,*,iostat=err) mvOption
!if(err/=0) return
!read(1,*,iostat=err) O1_value
!if(err/=0) return
!read(1,*,iostat=err) O2_cond
!if(err/=0) return
!read(1,*,iostat=err) O2_value
!if(err/=0) return
close(1)

end subroutine CRead_data

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! deprecated
subroutine CRead_data_OLD(file,DataFile,DataType,&
                          mvOption,O1_value,O2_cond,O2_value,err)

character(*), intent(in)::file
character(*), intent(out)::DataFile,DataType
character(*), intent(out)::mvOption,O2_cond
real(mrk), intent(out)::O1_value,O2_value
integer(mik), intent(out)::err

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return
read(1,'(A250)',iostat=err) DataFile
if(err/=0) return
read(1,'(A250)',iostat=err) DataType
if(err/=0) return
read(1,*,iostat=err) mvOption
if(err/=0) return
read(1,*,iostat=err) O1_value
if(err/=0) return
read(1,*,iostat=err) O2_cond
if(err/=0) return
read(1,*,iostat=err) O2_value
if(err/=0) return
close(1)

end subroutine CRead_data_OLD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine CRead_Inference(file,dist,parnames,priors,err)

use Distribution_tools, only:GetParNumber
use BayesianEstimation_tools, only:PriorListType

character(*), intent(in)::file
character(*), intent(out)::dist
character(*), pointer:: parnames(:)
type(PriorListType), pointer:: priors(:)
integer(mik), intent(out)::err
!locals
integer(mik)::npar,i,npriorpar
character(250)::mess

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return
read(1,'(A250)',iostat=err) dist
if(err/=0) return
! allocate pointers
call GetParNumber(DistID=trim(dist), npar=npar, err=err, mess=mess)
if(err/=0) return
if(associated(parnames)) deallocate(parnames)
allocate(parnames(npar))
if(associated(priors)) deallocate(priors)
allocate(priors(npar))
! keep reading
do i=1,npar
  read(1,'(A250)',iostat=err) parnames(i)
  if(err/=0) return
  read(1,'(A250)',iostat=err) priors(i)%dist
  if(err/=0) return
  call GetParNumber(DistID=trim(priors(i)%dist), npar=npriorpar, err=err, mess=mess)
  if(err/=0) return
  if(allocated(priors(i)%par)) deallocate(priors(i)%par)
  allocate(priors(i)%par(npriorpar))
  if(npriorpar>0) then
    read(1,*,iostat=err) priors(i)%par
    if(err/=0) return
  else
    read(1,*,iostat=err)
    if(err/=0) return
  endif
enddo
close(1)

end subroutine CRead_Inference

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CRead_SystematicErrors(file,npar,parnames,priors,err)

use Distribution_tools, only:GetParNumber
use BayesianEstimation_tools, only:PriorListType

character(*), intent(in)::file
integer(mik), intent(in)::npar
character(*), pointer:: parnames(:)
type(PriorListType), pointer:: priors(:)
integer(mik), intent(out)::err
!locals
integer(mik)::i,npriorpar
character(250)::mess

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return
! allocate pointers
if(associated(parnames)) deallocate(parnames)
allocate(parnames(npar))
if(associated(priors)) deallocate(priors)
allocate(priors(npar))
! keep reading
do i=1,npar
  read(1,'(A250)',iostat=err) parnames(i)
  if(err/=0) return
  read(1,'(A250)',iostat=err) priors(i)%dist
  if(err/=0) return
  call GetParNumber(DistID=trim(priors(i)%dist), npar=npriorpar, err=err, mess=mess)
  if(err/=0) return
  if(allocated(priors(i)%par)) deallocate(priors(i)%par)
  allocate(priors(i)%par(npriorpar))
  if(npriorpar>0) then
    read(1,*,iostat=err) priors(i)%par
    if(err/=0) return
  endif
enddo
close(1)

end subroutine CRead_SystematicErrors

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CRead_MCMC(file,&
                Nadapt,Ncycles,BurnFactor,Nslim,Nrep,MinMoveRate,&
                MaxMoveRate,DownMult,UpMult,MultFactor,Eps,err)

character(*), intent(in)::file
integer(mik), intent(out):: Nadapt,Ncycles,Nslim,Nrep,err
real(mrk),intent(out)::BurnFactor,MinMoveRate,MaxMoveRate,DownMult,&
                       UpMult,MultFactor,Eps

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return
read(1,*,iostat=err) Nadapt
if(err/=0) return
read(1,*,iostat=err) Ncycles
if(err/=0) return
read(1,*,iostat=err) BurnFactor
if(err/=0) return
read(1,*,iostat=err) Nslim
if(err/=0) return
read(1,*,iostat=err) Nrep
if(err/=0) return
read(1,*,iostat=err) MinMoveRate
if(err/=0) return
read(1,*,iostat=err) MaxMoveRate
if(err/=0) return
read(1,*,iostat=err) DownMult
if(err/=0) return
read(1,*,iostat=err) UpMult
if(err/=0) return
read(1,*,iostat=err) MultFactor
if(err/=0) return
read(1,*,iostat=err) Eps
if(err/=0) return
close(1)
end subroutine CRead_MCMC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine CRead_ResultFiles(file,&
                mc2_File,& !mc2_Save,&
                par_File,& !par_Save,&
                dat_File,& !dat_Save,&
                pdf_File,& !pdf_Save,&
                cdf_File,& !cdf_Save,&
                prd_File,& !prd_Save,&
                qtl_File,& !qtl_Save,&
                emp_File,& !emp_Save,&
                err)

character(*), intent(in)::file
character(*), intent(out)::mc2_File,par_File,dat_File,pdf_File,&
                           cdf_File,prd_File,qtl_File,emp_File
!logical,intent(out)::mc2_Save,par_Save,dat_Save,pdf_Save,cdf_Save,&
!                          prd_Save,qtl_Save,emp_Save
integer(mik), intent(out)::err
! locals
integer(mik)::i

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return

!read(1,*,iostat=err) i
!if(err/=0) return
!mc2_Save=(i==1)
read(1,'(A250)',iostat=err) mc2_File
if(err/=0) return

!read(1,*,iostat=err) i
!if(err/=0) return
!par_Save=(i==1)
read(1,'(A250)',iostat=err) par_File
if(err/=0) return

!read(1,*,iostat=err) i
!if(err/=0) return
!dat_Save=(i==1)
read(1,'(A250)',iostat=err) dat_File
if(err/=0) return

!read(1,*,iostat=err) i
!if(err/=0) return
!pdf_Save=(i==1)
read(1,'(A250)',iostat=err) pdf_File
if(err/=0) return

!read(1,*,iostat=err) i
!if(err/=0) return
!cdf_Save=(i==1)
read(1,'(A250)',iostat=err) cdf_File
if(err/=0) return

!read(1,*,iostat=err) i
!if(err/=0) return
!prd_Save=(i==1)
read(1,'(A250)',iostat=err) prd_File
if(err/=0) return

!read(1,*,iostat=err) i
!if(err/=0) return
!qtl_Save=(i==1)
read(1,'(A250)',iostat=err) qtl_File
if(err/=0) return

!read(1,*,iostat=err) i
!if(err/=0) return
!emp_Save=(i==1)
read(1,'(A250)',iostat=err) emp_File
if(err/=0) return

close(1)

end subroutine CRead_ResultFiles

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine CRead_ResultFiles_HBay(file,&
                mc2_File,mc2_Save,&
                par_File,par_Save,&
                pdf_File,pdf_Save,&
                cdf_File,cdf_Save,&
                prd_File,prd_Save,&
                qtl_File,qtl_Save,&
                err)

character(*), intent(in)::file
character(*), intent(out)::mc2_File,par_File,pdf_File,&
                           cdf_File,prd_File,qtl_File
logical,intent(out)::mc2_Save,par_Save,pdf_Save,cdf_Save,&
                          prd_Save,qtl_Save
integer(mik), intent(out)::err
! locals
integer(mik)::i

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return

read(1,*,iostat=err) i
if(err/=0) return
mc2_Save=(i==1)
read(1,'(A250)',iostat=err) mc2_File
if(err/=0) return

read(1,*,iostat=err) i
if(err/=0) return
par_Save=(i==1)
read(1,'(A250)',iostat=err) par_File
if(err/=0) return

read(1,*,iostat=err) i
if(err/=0) return
pdf_Save=(i==1)
read(1,'(A250)',iostat=err) pdf_File
if(err/=0) return

read(1,*,iostat=err) i
if(err/=0) return
cdf_Save=(i==1)
read(1,'(A250)',iostat=err) cdf_File
if(err/=0) return

read(1,*,iostat=err) i
if(err/=0) return
prd_Save=(i==1)
read(1,'(A250)',iostat=err) prd_File
if(err/=0) return

read(1,*,iostat=err) i
if(err/=0) return
qtl_Save=(i==1)
read(1,'(A250)',iostat=err) qtl_File
if(err/=0) return

close(1)

end subroutine CRead_ResultFiles_HBay

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine CRead_ResultOptions(file,&
                       xmin,xmax,xn,& ! options for pdf/cdf
                       pp,kernel,& ! options for empirical stats
                       pmin,pmax,pn,level,PtoT,invertT,& ! options for quantiles
                       err)

character(*), intent(in)::file
real(mrk),intent(out)::xmin,xmax,xn,pmin,pmax,pn,level,PtoT
logical, intent(out)::invertT
character(250),intent(out)::pp,kernel
integer(mik), intent(out)::err

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return
read(1,*,iostat=err) xmin
if(err/=0) return
read(1,*,iostat=err) xmax
if(err/=0) return
read(1,*,iostat=err) xn
if(err/=0) return
read(1,'(A250)',iostat=err) pp
if(err/=0) return
read(1,'(A250)',iostat=err) kernel
if(err/=0) return
read(1,*,iostat=err) pmin
if(err/=0) return
read(1,*,iostat=err) pmax
if(err/=0) return
read(1,*,iostat=err) pn
if(err/=0) return
read(1,*,iostat=err) level
if(err/=0) return
read(1,*,iostat=err) PtoT
if(err/=0) return
read(1,*,iostat=err) invertT
if(err/=0) return
close(1)

end subroutine CRead_ResultOptions

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! deprecated
subroutine CRead_ResultOptions_OLD(file,&
                       xmin,xmax,xn,& ! options for pdf/cdf
                       pp,kernel,& ! options for empirical stats
                       pmin,pmax,pn,level,PtoT,& ! options for quantiles
                       err)

character(*), intent(in)::file
real(mrk),intent(out)::xmin,xmax,xn,pmin,pmax,pn,level,PtoT
character(250),intent(out)::pp,kernel
integer(mik), intent(out)::err

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return
read(1,*,iostat=err) xmin
if(err/=0) return
read(1,*,iostat=err) xmax
if(err/=0) return
read(1,*,iostat=err) xn
if(err/=0) return
read(1,'(A250)',iostat=err) pp
if(err/=0) return
read(1,'(A250)',iostat=err) kernel
if(err/=0) return
read(1,*,iostat=err) pmin
if(err/=0) return
read(1,*,iostat=err) pmax
if(err/=0) return
read(1,*,iostat=err) pn
if(err/=0) return
read(1,*,iostat=err) level
if(err/=0) return
read(1,*,iostat=err) PtoT
if(err/=0) return
close(1)

end subroutine CRead_ResultOptions_OLD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine CRead_GetQT(file,dist,QorT,val,PtoT,invertT,level,mc2_File,QT_File,err)

character(*), intent(in)::file
character(*), intent(out)::dist,mc2_File,QT_File
character(1), intent(out)::QorT
real(mrk), intent(out)::val,PtoT,level
logical, intent(out)::invertT
integer(mik), intent(out)::err
!locals
character(250)::mess

open(unit=1, file=file,status="OLD",iostat=err)
if(err/=0) return
read(1,'(A250)',iostat=err) dist
if(err/=0) return
read(1,'(A1)',iostat=err) QorT
if(err/=0) return
read(1,*,iostat=err) val
if(err/=0) return
read(1,*,iostat=err) PtoT
if(err/=0) return
read(1,*,iostat=err) invertT
if(err/=0) return
read(1,*,iostat=err) level
if(err/=0) return
read(1,'(A250)',iostat=err) mc2_File
if(err/=0) return
read(1,'(A250)',iostat=err) QT_File
if(err/=0) return
close(1)

end subroutine CRead_GetQT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine RemoveMV(yin,mvOption,O1_value,&
              O2_cond,O2_value,&
              yout,err)

real(mrk), intent(in)::yin(:,:)
character(*), intent(in)::mvOption,O2_cond
real(mrk), intent(in)::O1_value,O2_value
real(mrk),pointer::yout(:,:)
integer(mik), intent(out)::err
!locals
integer(mik),parameter::O1col=2,O2col=3
integer(mik)::i,n,k
logical::binme(size(yin,dim=1))

err=0
n=size(yin,dim=1)
k=0
do i=1,n
  select case (trim(mvOption))
  case("Option1")
    binme(i)= (yin(i,O1col)==O1_value)
  case("Option2")
    binme(i)=ApplyOption2(yin(i,O2col),O2_cond,O2_value)
  case("Both")
    binme(i)=(yin(i,O1col)==O1_value).or. ApplyOption2(yin(i,O2col),O2_cond,O2_value)
  case default
    err=1;return
  end select
  if(.not.binme(i)) k=k+1
enddo
if(associated(yout)) deallocate(yout)
allocate(yout(k,size(yin,dim=2)))
do i=1,size(yin,dim=2)
  yout(:,i)=pack(yin(:,i),.not.binme)
enddo
end subroutine RemoveMV
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ApplyOption2(y,O2_cond,O2_value)
real(mrk),intent(in)::y,O2_value
character(*), intent(in)::O2_cond
logical::ApplyOption2

select case(trim(O2_cond))
case("==")
  ApplyOption2=(y==O2_value)
case("<")
  ApplyOption2=(y<O2_value)
case(">")
  ApplyOption2=(y>O2_value)
case("<=")
  ApplyOption2=(y<=O2_value)
case(">=")
  ApplyOption2=(y>=O2_value)
case("/=")
  ApplyOption2=(y/=O2_value)
end select
end function ApplyOption2

end module Estimator_tools

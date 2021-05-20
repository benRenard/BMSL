module ExtraFlo_tools

!~**********************************************************************
!~* Purpose: 
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 07/03/2012
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
public :: GetFF_modal,GetFF_pred,GetPval_modal,GetPval_pred,&
          WriteReliabilityFile_modal,WriteReliabilityFile_pred,&
          WritePvalFile_modal,WritePvalFile_pred,&
          WriteQuantFile_modal,WriteQuantFile_pred,&
          GetSpanT_modal,GetSpanT_pred,&
          WriteSpanFile_pred,WriteSpanFile_modal,&
          GetNT_modal,GetNT_pred

Contains

subroutine GetFF_modal(X,distID,modal,FF,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalStats
use Distribution_tools, only:GetCdf
real(mrk),intent(in)::X(:),modal(:)
character(*),intent(in):: distID
real(mrk),intent(out)::FF
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::N
real(mrk)::maxi,cdf
logical::feas

err=0;mess='';FF=undefRN
if(modal(1)==undefRN) return 
call GetEmpiricalStats(x=X,n=N,maxi=maxi,err=err,mess=mess)
If(err>0) then
  mess='GetFF_modal:'//trim(mess);return
endif
call GetCdf(DistId=DistId,x=maxi,par=modal,cdf=cdf,feas=feas,err=err,mess=mess)
If(err>0) then
  mess='GetFF_modal:'//trim(mess);return
endif
If(.not.feas) then
  mess='GetFF_modal:Fatal:Unfeasible Cdf computation';return
endif
FF=cdf**N
end subroutine GetFF_modal
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetFF_pred(X,pred,FF,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalStats,GetEmpiricalPval

real(mrk),intent(in)::X(:),pred(:)
real(mrk),intent(out)::FF
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::N
real(mrk)::maxi,cdf
logical::feas

err=0;mess='';FF=undefRN
if(pred(1)==undefRN) return 
call GetEmpiricalStats(x=X,n=N,maxi=maxi,err=err,mess=mess)
If(err>0) then
  mess='GetFF_pred:'//trim(mess);return
endif
call GetEmpiricalPval(obs=maxi,x=pred,pv=cdf,err=err,mess=mess)
If(err>0) then
  mess='GetFF_pred:'//trim(mess);return
endif
FF=cdf**N
end subroutine GetFF_pred
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetPval_modal(X,distID,modal,Pval,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalStats
use Distribution_tools, only:GetCdf
use numerix_dmsl_kit, only:quicksort

real(mrk),intent(in)::X(:),modal(:)
character(*),intent(in):: distID
real(mrk),intent(out)::Pval(size(X))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::N,i
logical::feas

err=0;mess='';Pval=undefRN
if(modal(1)==undefRN) return 
N=size(X)
do i=1,N
    call GetCdf(DistId=DistId,x=X(i),par=modal,cdf=Pval(i),feas=feas,err=err,mess=mess)
    If(err>0) then
      mess='GetPval_modal:'//trim(mess);return
    endif
    If(.not.feas) then
      mess='GetPval_modal:Fatal:Unfeasible Cdf computation';return
    endif
enddo
call quicksort(arr=Pval,ascnd=.true.,err=err)

end subroutine GetPval_modal
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetPval_pred(X,pred,Pval,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalPval
use numerix_dmsl_kit, only:quicksort

real(mrk),intent(in)::X(:),pred(:)
real(mrk),intent(out)::Pval(size(X))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::N,i
logical::feas
real(mrk)::predi(size(pred))

err=0;mess='';Pval=undefRN
if(pred(1)==undefRN) return 
predi=pred
call quicksort(arr=predi,ascnd=.true.,err=err)
N=size(X)
do i=1,N
    call GetEmpiricalPval(obs=X(i),x=predi,IsXSorted=.true.,pv=Pval(i),err=err,mess=mess)
    If(err>0) then
      mess='GetPval_pred:'//trim(mess);return
    endif
enddo
call quicksort(arr=Pval,ascnd=.true.,err=err)

end subroutine GetPval_pred
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetNT_modal(T,X,distID,modal,NT,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use Distribution_tools, only:GetQuantile

real(mrk),intent(in)::T(:),X(:),modal(:)
character(*),intent(in):: distID
real(mrk),intent(out)::NT(size(T))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::N,i
real(mrk)::QT
logical::feas

err=0;mess='';NT=undefRN
if(modal(1)==undefRN) return 
N=size(T)
do i=1,N
    call GetQuantile(DistId=DistId,p=(1._mrk-1._mrk/T(i)),par=modal,q=QT,feas=feas,err=err,mess=mess)
    If(err>0) then
      mess='GetNT_modal:'//trim(mess);return
    endif
    If(.not.feas) then
      mess='GetNT_modal:Fatal:Unfeasible Cdf computation';return
    endif
    NT(i)=count(X>QT)
enddo

end subroutine GetNT_modal
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetNT_pred(T,X,pred,NT,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalQuantile
use numerix_dmsl_kit, only:quicksort

real(mrk),intent(in)::T(:),X(:),pred(:)
real(mrk),intent(out)::NT(size(T))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::N,i
real(mrk)::QT
logical::feas
real(mrk)::predi(size(pred))

err=0;mess='';NT=undefRN
if(pred(1)==undefRN) return 
predi=pred
call quicksort(arr=predi,ascnd=.true.,err=err)
N=size(T)
do i=1,N
    call GetEmpiricalQuantile(p=(1._mrk-1._mrk/T(i)),x=predi,IsXSorted=.true.,q=QT,err=err,mess=mess)
    If(err>0) then
      mess='GetNT_pred:'//trim(mess);return
    endif
    NT(i)=count(X>QT)
enddo

end subroutine GetNT_pred
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine WriteReliabilityFile_modal(file,X,distID,modal,ntrim,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use ExtraFlo_utilities
use utilities_dmsl_kit, only:GetSpareUnit

character(*), intent(in)::file,distID
type(VoVtype),intent(in)::X(:)
real(mrk), intent(in)::modal(:,:)
integer(mik), intent(in), optional::ntrim
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::Nsite,site,unt
real(mrk)::FF,NT(2)
real(mrk), allocatable::Y(:)

Nsite=size(modal,dim=1)
call GetSpareUnit(unt,err,mess)
If(err>0) then
  mess='WriteReliabilityFile_modal:'//trim(mess);return
endif
open(unit=unt,file=trim(file),status='replace')
write(unt,'(5A14)') 'site','n','FF','N10','N100'
do site=1,Nsite
    if(present(ntrim)) then
        if(allocated(Y)) deallocate(Y);allocate(Y(ntrim))
        if(size(X(site)%val)<ntrim) then
            write(unt,'(2I14,e14.6,2I14)') site,UndefIN,undefRN,UndefIN,UndefIN
            cycle
        else
            Y=X(site)%val(1:ntrim)
        endif
    else
        if(allocated(Y)) deallocate(Y);allocate(Y(size(X(site)%val)))
        Y=X(site)%val
    endif
    ! FF
    call GetFF_modal(X=Y,distID=distID,modal=modal(site,:),&
                     FF=FF,err=err,mess=mess)
    If(err>0) then
      mess='WriteReliabilityFile_modal:'//trim(mess);return
    endif
    ! N10
    call GetNT_modal(T=(/10._mrk,100._mrk/),X=Y,distID=distID,&
                    modal=modal(site,:),NT=NT,err=err,mess=mess)
    If(err>0) then
      mess='WriteReliabilityFile_modal:'//trim(mess);return
    endif
    write(unt,'(2I14,e14.6,2I14)') site,size(Y),FF,nint(NT(1)),nint(NT(2))
enddo
close(unt)
end subroutine WriteReliabilityFile_modal
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine WriteReliabilityFile_pred(file,X,pred,ntrim,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use ExtraFlo_utilities
use utilities_dmsl_kit, only:GetSpareUnit

character(*), intent(in)::file
type(VoVtype),intent(in)::X(:)
real(mrk), intent(in)::pred(:,:)
integer(mik), intent(in), optional::ntrim
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::Nsite,site,unt
real(mrk)::FF,NT(2)
real(mrk), allocatable::Y(:)

Nsite=size(pred,dim=1)
call GetSpareUnit(unt,err,mess)
If(err>0) then
  mess='WriteReliabilityFile_pred:'//trim(mess);return
endif
open(unit=unt,file=trim(file),status='replace')
write(unt,'(5A14)') 'site','n','FF','N10','N100'
do site=1,Nsite
    if(present(ntrim)) then
        if(allocated(Y)) deallocate(Y);allocate(Y(ntrim))
        if(size(X(site)%val)<ntrim) then
            write(unt,'(2I14,e14.6,2I14)') site,UndefIN,undefRN,UndefIN,UndefIN
            cycle
        else
            Y=X(site)%val(1:ntrim)
        endif
    else
        if(allocated(Y)) deallocate(Y);allocate(Y(size(X(site)%val)))
        Y=X(site)%val
    endif
    ! FF
    call GetFF_pred(X=Y,pred=pred(site,:),&
                     FF=FF,err=err,mess=mess)
    If(err>0) then
      mess='WriteReliabilityFile_pred:'//trim(mess);return
    endif
    ! N10
    call GetNT_pred(T=(/10._mrk,100._mrk/),X=Y,&
                    pred=pred(site,:),NT=NT,err=err,mess=mess)
    If(err>0) then
      mess='WriteReliabilityFile_pred:'//trim(mess);return
    endif
    write(unt,'(2I14,e14.6,2I14)') site,size(Y),FF,nint(NT(1)),nint(NT(2))
enddo
close(unt)
end subroutine WriteReliabilityFile_pred
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine WritePvalFile_modal(file,X,distID,modal,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use ExtraFlo_utilities
use utilities_dmsl_kit, only:GetSpareUnit

character(*), intent(in)::file,distID
type(VoVtype),intent(in)::X(:)
real(mrk), intent(in)::modal(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::Nsite,site,unt,N
real(mrk),allocatable::pval(:)

Nsite=size(modal,dim=1)
call GetSpareUnit(unt,err,mess)
If(err>0) then
  mess='WritePvalFile_modal:'//trim(mess);return
endif
open(unit=unt,file=trim(file),status='replace')
write(unt,'(3A14)') 'site','n','pvals'
do site=1,Nsite
    N=size(X(site)%val)
    if(allocated(pval)) deallocate(pval);allocate(pval(N))
    ! pvals
    call GetPval_modal(X=X(site)%val,distID=distID,modal=modal(site,:),&
                      Pval=pval,err=err,mess=mess)
    If(err>0) then
      mess='WritePvalFile_modal:'//trim(mess);return
    endif
    write(unt,'(2I14,<N>e14.6)') site,N,pval
enddo
close(unt)
end subroutine WritePvalFile_modal
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine WritePvalFile_pred(file,X,pred,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use ExtraFlo_utilities
use utilities_dmsl_kit, only:GetSpareUnit

character(*), intent(in)::file
type(VoVtype),intent(in)::X(:)
real(mrk), intent(in)::pred(:,:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::Nsite,site,unt,N
real(mrk),allocatable::pval(:)

Nsite=size(pred,dim=1)
call GetSpareUnit(unt,err,mess)
If(err>0) then
  mess='WritePvalFile_pred:'//trim(mess);return
endif
open(unit=unt,file=trim(file),status='replace')
write(unt,'(3A14)') 'site','n','pvals'
do site=1,Nsite
    N=size(X(site)%val)
    if(allocated(pval)) deallocate(pval);allocate(pval(N))
    ! pvals
    call GetPval_pred(X=X(site)%val,pred=pred(site,:),&
                      Pval=pval,err=err,mess=mess)
    If(err>0) then
      mess='WritePvalFile_pred:'//trim(mess);return
    endif
    write(unt,'(2I14,<N>e14.6)') site,N,pval
enddo
close(unt)
end subroutine WritePvalFile_pred
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine WriteQuantFile_pred(file,Listp,pred,LogSpace,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalQuantile
use utilities_dmsl_kit, only:GetSpareUnit
use numerix_dmsl_kit,only:quicksort

character(*), intent(in)::file
real(mrk), intent(in)::pred(:,:),Listp(:)
logical, intent(in)::LogSpace
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::i,Nsite,NlistP,site,unt
real(mrk),allocatable::Qp(:)
real(mrk)::predi(size(pred,dim=2))

Nsite=size(pred,dim=1)
NlistP=size(Listp)
if(allocated(Qp)) deallocate(Qp);allocate(Qp(NlistP))
call GetSpareUnit(unt,err,mess)
If(err>0) then
  mess='WriteQuantFile_pred:'//trim(mess);return
endif
open(unit=unt,file=trim(file),status='replace')
write(unt,'(A14,<NlistP>e14.6)') 'site',Listp
do site=1,Nsite
    predi=pred(site,:)
    call quicksort(arr=predi,ascnd=.true.,err=err)
    do i=1,NlistP
       ! quantiles
        call GetEmpiricalQuantile(p=Listp(i),x=predi,isXSorted=.true.,q=Qp(i),err=err,mess=mess)
        If(err>0) then
          mess='WriteQuantFile_pred:'//trim(mess);return
        endif
    enddo
    if(LogSpace) Qp=exp(Qp)
    write(unt,'(I14,<NlistP>e14.6)') site,Qp
enddo
close(unt)
end subroutine WriteQuantFile_pred 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine WriteQuantFile_modal(file,Listp,distID,modal,LogSpace,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use utilities_dmsl_kit, only:GetSpareUnit
use Distribution_tools, only:GetQuantile

character(*), intent(in)::file,distID
real(mrk), intent(in)::modal(:,:),Listp(:)
logical, intent(in)::LogSpace
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::i,Nsite,NlistP,site,unt
real(mrk),allocatable::Qp(:)
logical::feas

Nsite=size(modal,dim=1)
NlistP=size(Listp)
if(allocated(Qp)) deallocate(Qp);allocate(Qp(NlistP))
call GetSpareUnit(unt,err,mess)
If(err>0) then
  mess='WriteQuantFile_modal:'//trim(mess);return
endif
open(unit=unt,file=trim(file),status='replace')
write(unt,'(A14,<NlistP>e14.6)') 'site',Listp
do site=1,Nsite
    do i=1,NlistP
       ! quantiles
        call GetQuantile(DistId=distID,p=Listp(i),par=modal(site,:),q=Qp(i),&
                feas=feas,err=err,mess=mess)
        If(err>0) then
          mess='WriteQuantFile_modal:'//trim(mess);return
        endif
    enddo
    if(LogSpace) Qp=exp(Qp)
    write(unt,'(I14,<NlistP>e14.6)') site,Qp
enddo
close(unt)
end subroutine WriteQuantFile_modal  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetSpanT_modal(T,distID,modal1,modal2,LogSpace,SpanT,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use Distribution_tools, only:GetQuantile

real(mrk),intent(in)::T(:),modal1(:),modal2(:)
character(*),intent(in):: distID
logical, intent(in)::LogSpace
real(mrk),intent(out)::SpanT(size(T))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::N,i
real(mrk)::QT1,QT2
logical::feas

err=0;mess='';SpanT=undefRN
if(modal1(1)==undefRN .or. modal2(1)==undefRN) return 
N=size(T)
do i=1,N
    call GetQuantile(DistId=DistId,p=(1._mrk-1._mrk/T(i)),par=modal1,q=QT1,feas=feas,err=err,mess=mess)
    If(err>0) then
      mess='GetSpanT_modal:'//trim(mess);return
    endif
    If(.not.feas) then
      mess='GetSpanT_modal:Fatal:Unfeasible Cdf computation';return
    endif
    if(LogSpace) QT1=exp(QT1)
    call GetQuantile(DistId=DistId,p=(1._mrk-1._mrk/T(i)),par=modal2,q=QT2,feas=feas,err=err,mess=mess)
    If(err>0) then
      mess='GetSpanT_modal:'//trim(mess);return
    endif
    If(.not.feas) then
      mess='GetSpanT_modal:Fatal:Unfeasible Cdf computation';return
    endif
    if(LogSpace) QT2=exp(QT2)
    SpanT(i)=abs(QT1-QT2)/(0.5_mrk*(QT1+QT2))
enddo

end subroutine GetSpanT_modal
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine GetSpanT_pred(T,pred1,pred2,LogSpace,SpanT,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalQuantile
use numerix_dmsl_kit, only:quicksort

real(mrk),intent(in)::T(:),pred1(:),pred2(:)
logical, intent(in)::LogSpace
real(mrk),intent(out)::SpanT(size(T))
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::N,i
real(mrk)::QT1,QT2,predi1(size(pred1)),predi2(size(pred2))
logical::feas

err=0;mess='';SpanT=undefRN
if(pred1(1)==undefRN .or. pred2(1)==undefRN) return 
predi1=pred1
call quicksort(arr=predi1,ascnd=.true.,err=err)
predi2=pred2
call quicksort(arr=predi2,ascnd=.true.,err=err)
N=size(T)
do i=1,N
    call GetEmpiricalQuantile(p=(1._mrk-1._mrk/T(i)),x=predi1,IsXSorted=.true.,q=QT1,err=err,mess=mess)
    If(err>0) then
      mess='GetSpanT_pred:'//trim(mess);return
    endif
    if(LogSpace) QT1=exp(QT1)
    call GetEmpiricalQuantile(p=(1._mrk-1._mrk/T(i)),x=predi2,ISXSorted=.true.,q=QT2,err=err,mess=mess)
    If(err>0) then
      mess='GetSpanT_pred:'//trim(mess);return
    endif
    if(LogSpace) QT2=exp(QT2)
    SpanT(i)=abs(QT1-QT2)/(0.5_mrk*(QT1+QT2))
enddo

end subroutine GetSpanT_pred
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine WriteSpanFile_pred(file,pred1,pred2,LogSpace,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use ExtraFlo_utilities
use utilities_dmsl_kit, only:GetSpareUnit

character(*), intent(in)::file
real(mrk), intent(in)::pred1(:,:),pred2(:,:)
logical, intent(in)::LogSpace
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::Nsite,site,unt
real(mrk)::SpanT(3)

Nsite=size(pred1,dim=1)
call GetSpareUnit(unt,err,mess)
If(err>0) then
  mess='WriteSpanFile_pred:'//trim(mess);return
endif
open(unit=unt,file=trim(file),status='replace')
write(unt,'(4A14)') 'site','Span10','Span100','Span1000'
do site=1,Nsite
    call GetSpanT_pred(T=(/10._mrk,100._mrk,1000._mrk/),pred1=pred1(site,:),&
                    pred2=pred2(site,:),LogSpace=LogSpace,SpanT=SpanT,err=err,mess=mess)
    If(err>0) then
      mess='WriteSpanFile_pred:'//trim(mess);return
    endif
    write(unt,'(I14,3e14.6)') site,SpanT
enddo
close(unt)
end subroutine WriteSpanFile_pred
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine WriteSpanFile_modal(file,distID,modal1,modal2,LogSpace,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.
!^*		3.
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use ExtraFlo_utilities
use utilities_dmsl_kit, only:GetSpareUnit

character(*), intent(in)::file,distID
real(mrk), intent(in)::modal1(:,:),modal2(:,:)
logical, intent(in)::LogSpace
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::Nsite,site,unt
real(mrk)::SpanT(3)

Nsite=size(modal1,dim=1)
call GetSpareUnit(unt,err,mess)
If(err>0) then
  mess='WriteSpanFile_modal:'//trim(mess);return
endif
open(unit=unt,file=trim(file),status='replace')
write(unt,'(4A14)') 'site','Span10','Span100','Span1000'
do site=1,Nsite
    call GetSpanT_modal(T=(/10._mrk,100._mrk,1000._mrk/),distID=distID,modal1=modal1(site,:),&
                    modal2=modal2(site,:),LogSpace=LogSpace,SpanT=SpanT,err=err,mess=mess)
    If(err>0) then
      mess='WriteSpanFile_modal:'//trim(mess);return
    endif
    write(unt,'(I14,3e14.6)') site,SpanT
enddo
close(unt)
end subroutine WriteSpanFile_modal
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end module ExtraFlo_tools
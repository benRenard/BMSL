module aggregation_tools
! Purpose: tools for aggregating functions sampled regularly or not
! Programmer: Ben Renard, Irstea & the Big-D. Kavetski, UoA
! Last modified: 09/05/2015
!*********************************************************************
! Comments:  
!*********************************************************************
! References:
!*********************************************************************
! 2Do List: 
!*********************************************************************
! Quick description of public procedures:
! 1. 
!*********************************************************************
use kinds_dmsl_kit ! numeric kind definitions from DMSL
implicit none
private
public::aggregator,basic_op
! Continuity flags: a 'y' flagged as discontinuous should not be interpolated with the PREVIOUS 'y'
integer(mik),parameter,public::CC_CONT=QC_PRESENT,CC_DISCONT=QC_MISSVAL
! Indices defining what's in opID:
! (/1. interpolation option, 2. function, 3. qc option, 4. cc option/)
integer(mik),parameter,public::ix_interpol=1,ix_funk=2,ix_qc=3,ix_cc=4
! Interpolation options. 
! NONE   = no interpolation at all, just apply the operator function to the (x,y) verifying xmin<=x<xmax
! LINEAR = linear interpolation. Includes interpolation one step before xmin and one step after xmax
integer(mik),parameter,public::NONE_interpol=0,LINEAR_interpol=1
! Operator functions
integer(mik),parameter,public::SUM_funk=0,AVE_funk=1,MIN_funk=2,MAX_funk=3,&
                               LEFTVAL_funk=4,RIGHTVAL_funk=5 ! useful for regridding
! Quality code options
! IGNORE = if no interpolation, just apply the operator function to non-missing (x,y) in [xmin xmax]
!          if interpolation, interpolate above the missing (x,y)
! INFECT = if any missing value is encountered, the result is also missing
integer(mik),parameter,public::INFECT_qc=0, IGNORE_qc=1
! Continuity code options - only used when interpolation requested
! CONSERV = if any discontinuity flag is encoutered, result is missing
! LIBERAL = if a discontinuity flag is encoutered, interpolate anyway but flag the result as discontinuous
integer(mik),parameter,public::CONSERV_cc=0,LIBERAL_cc=1


interface aggregator
  module procedure aggregator_intervals,aggregator_grid
endinterface aggregator


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine basic_op(opID,y,x,xmin,xmax,qc,cc,aux,& !IN
                         yOut,xOut,qcOut,ccOut,err,message) !OUT
! Purpose: Basic aggregation operators
! Programmer: Ben Renard, Irstea & the BigD Kavetski
! Last modified: 08/05/2015
!*********************************************************************
! Comments: Designed for functions sampled at regular or irregular intervals
!           Current options do not handle  y's represent an average over a time step
!           Assumes sampling points x are ordered lo->hi
!*********************************************************************
! References:
!*********************************************************************
! IN
!   1. opID, in the form (/interpolation option, function, qc option, cc option/)
!   2. y, the function sampled at discrete points
!   3. x, the sampling points
!   4. xmin, lower bound of the interval on which the operator is applied
!   5. xmax, upper bound of the interval on which the operator is applied
!   6. [qc], quality codes
!   7. [cc], continuity codes
!   8. [aux], auxiliaries
! OUT
!   1. yout, value(s) after applying the operator
!   2. xout, x-value associated with yout (defaults: center of [xmin xmax])
!   3. [qcOut], quality code(s) associated with yout
!   4. [ccOut], continuity code(s) associated with yout
!   5. err, error code
!   6. message
!*********************************************************************
! 2Do List:
!*********************************************************************
use utilities_dmsl_kit,only:half,assertEq,ns=>number_string,quickLinInterp
implicit none
! dummies
integer(mik),intent(in)::opID(:)
real(mrk),intent(in)::y(:),x(:),xmin,xmax
integer(mik),intent(in),optional::qc(:),cc(:)
real(mrk),intent(in),optional::aux(:)
real(mrk),intent(out)::yout(:),xout
integer(mik),intent(out),optional::qcOut(:),ccOut(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam='basic_op'
integer(mik)::n,ix_n
integer(mik),pointer::ix(:)
real(mrk),allocatable::xs(:),ys(:)
logical(mlk)::ok,haveCbreak,abortFlag
! start procedure here
err=EXIT_SUCCESS; message=procnam//'/ok'; yout=undefRN; xout=undefRN
if(present(qcOut))qcOut=QC_MISSVAL
if(present(ccOut))ccOut=CC_CONT
call assertEq(size(x),size(y),ok,n)
if(.not.ok)then;err=EXIT_FAILURE;message='f-'//procnam//'/dimError[x,y]';return;endif
if(xmin>=xmax)then;err=EXIT_FAILURE;message='f-'//procnam//'/error[xmin>=xmax]';return;endif
xout=half*(xmin+xmax) ! assign x for interval [xmin xmax]
call getActiveIndices(x,xmin,xmax,opID,qc,ix,abortFlag)
if(abortFlag)then
    if(associated(ix)) nullify(ix);return
endif
ix_n=size(ix)
allocate(xs(ix_n),ys(ix_n)); xs=x(ix);ys=y(ix) ! get relevant x's and y's
selectcase(opID(ix_interpol))   ! interpolation method
case(NONE_interpol)
  if(present(qcOut))qcOut=QC_PRESENT
  selectcase(opID(ix_funk))     ! aggregation method
  case(SUM_funk); yout(1)=sum(ys)
  case(AVE_funk); yout(1)=sum(ys)/real(ix_n,mrk)
  case(MAX_funk); yout(1)=maxval(ys)
  case(MIN_funk); yout(1)=minval(ys)
  case(LEFTVAL_funk); yout(1)=ys(1);xout=xmin
  case(RIGHTVAL_funk); yout(1)=ys(ix_n);xout=xmax
  case default
    err=EXIT_FAILURE;message='f-'//procnam//'/unknown[funk='//trim(ns(opID(ix_funk)))//']';return
  endselect
case(LINEAR_interpol)    ! perform linear interpolation at both ends
  ys(1)   =quickLinInterp(xs(1),ys(1),xs(2),ys(2),xmin);                 xs(1)=xmin
  ys(ix_n)=quickLinInterp(xs(ix_n-1),ys(ix_n-1),xs(ix_n),ys(ix_n),xmax); xs(ix_n)=xmax
  if(present(qcOut))qcOut=QC_PRESENT
  if(present(cc))then
    haveCbreak=.not.all(cc(ix(2:))==CC_CONT)
    if(haveCbreak)then
      selectcase(opID(ix_cc))
      case(LIBERAL_cc)      ! do interpolation anyway, just flag the discontinuity
        if(present(ccOut))ccOut=CC_DISCONT            
      case(CONSERV_cc)      ! abandon and flag quality code as missing
        qcOut=QC_MISSVAL;
        if(allocated(xs)) deallocate(xs)
        if(allocated(ys)) deallocate(ys)
        if(associated(ix)) nullify(ix)
        return
      case default
        err=EXIT_FAILURE;message='f-'//procnam//'/unknown[cc option='//trim(ns(opID(ix_cc)))//']';return
      endselect
    endif
  endif
  selectcase(opID(ix_funk)) ! aggregation method
  case(SUM_funk); yout(1)=quickTrap(xs,ys)
  case(AVE_funk); yout(1)=quickTrap(xs,ys)/(xmax-xmin)
  case(MAX_funk); yout(1)=maxval(ys)
  case(MIN_funk); yout(1)=minval(ys)
  case(LEFTVAL_funk); yout(1)=ys(1);xout=xmin
  case(RIGHTVAL_funk); yout(1)=ys(ix_n);xout=xmax
  case default
    err=EXIT_FAILURE;message='f-'//procnam//'/unknown[funk='//trim(ns(opID(ix_funk)))//']';return
  endselect
case default
  err=EXIT_FAILURE;message='f-'//procnam//'/unknown[interpol='//trim(ns(opID(ix_interpol)))//']';return
endselect
if(allocated(xs)) deallocate(xs)
if(allocated(ys)) deallocate(ys)
if(associated(ix)) nullify(ix)
! end procedure here
endsubroutine basic_op

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine aggregator_intervals(operateur,opID,y,x,intervals,qc,cc,aux,& !IN
                                     yOut,xOut,qcOut,ccOut,err,message) !OUT
! Purpose: Aggregator engine
! Programmer: Ben Renard, Irstea & the BigD Kavetski
! Last modified: 09/05/2015
!*********************************************************************
! References:
!*********************************************************************
! IN
!   1. operateur, subroutine implementing the operator
!   2. opID, the properties of the operator to apply:
!      (/interpolation option, aggregation function, qc option, cc option/)
!   3. y, the function sampled at discrete points
!   4. x, the sampling points
!   5. intervals, the list of intervals where the operator is applied
!   6. [qc], quality codes
!   7. [cc], continuity codes
!   8. [aux], auxiliaries
! OUT
!   1. yout, values after applying the operator
!   2. xout, x-values associated with yout
!   3. [qcOut], quality codes associated with yout
!   4. [ccOut], continuity codes associated with yout
!   5. err, error code
!   6. message
!*********************************************************************
! 2Do List:
!*********************************************************************
use utilities_dmsl_kit,only:assertEq
implicit none
interface
  pure subroutine operateur(opID,y,x,xmin,xmax,qc,cc,aux,&
                            yOut,xOut,qcOut,ccOut,err,message)
  use kinds_dmsl_kit
  implicit none
  integer(mik),intent(in)::opID(:)
  real(mrk),intent(in)::y(:),x(:),xmin,xmax
  integer(mik),intent(in),optional::qc(:),cc(:)
  real(mrk),intent(in),optional::aux(:)
  real(mrk),intent(out)::yout(:),xout
  integer(mik),intent(out),optional::qcOut(:),ccOut(:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine operateur
endinterface
! dummies
integer(mik),intent(in)::opID(:)
real(mrk),intent(in)::y(:),x(:),intervals(:,:)
integer(mik),intent(in),optional::qc(:),cc(:)
real(mrk),intent(in),optional::aux(:)
real(mrk),intent(out)::yout(:,:),xout(:)
integer(mik),intent(out),optional::qcOut(:,:),ccOut(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam='aggregator_intervals'
real(mrk)::xmin,xmax
integer(mik)::i,int_n,int_p
logical(mlk)::ok
! start procedure here
err=EXIT_SUCCESS; message=procnam//'/ok'; yout=undefRN; xout=undefRN
if(present(qcOut))qcOut=QC_MISSVAL
if(present(ccOut))ccOut=CC_CONT
! check sizes
call assertEq(size(intervals,2),2,ok,int_p)
if(.not.ok)then;err=EXIT_FAILURE;message='f-'//procnam//'/dimError[intervals]';return;endif
call assertEq(size(intervals,1),size(yout,1),ok,int_n)
if(.not.ok)then;err=EXIT_FAILURE;message='f-'//procnam//'/dimError[intervals,yout]';return;endif
call assertEq(size(intervals,1),size(xout,1),ok,int_n)
if(.not.ok)then;err=EXIT_FAILURE;message='f-'//procnam//'/dimError[intervals,xout]';return;endif
if(present(qcOut))then
  call assertEq(size(intervals,1),size(qcout,1),ok,int_n)
  if(.not.ok)then;err=EXIT_FAILURE;message='f-'//procnam//'/dimError[intervals,qcOut]';return;endif
endif
if(present(ccOut))then
  call assertEq(size(intervals,1),size(ccout,1),ok,int_n)
  if(.not.ok)then;err=EXIT_FAILURE;message='f-'//procnam//'/dimError[intervals,ccOut]';return;endif
endif
do i=1,size(intervals,1)
  xmin=intervals(i,1); xmax=intervals(i,2)
  if(xmin>=xmax)then;err=EXIT_FAILURE;message='f-'//procnam//'/error[xmin>=xmax]';return;endif
! Horrible if-elses to handle optional outputs agreed upon with D after a long debate
  if(present(qcOut).and.present(ccOut))then ! both optional outputs are present
    call operateur(opID=opID,y=y,x=x,xmin=xmin,xmax=xmax,qc=qc,cc=cc,aux=aux,&
      yOut=yOut(i,:),xOut=xOut(i),qcOut=qcOut(i,:),ccOut=ccOut(i,:),err=err,message=message)
    if(err/=EXIT_SUCCESS)then;message='f-'//procnam//'&'//message;return;endif
  elseif(present(qcOut))then                ! only qcOut is present
    call operateur(opID=opID,y=y,x=x,xmin=xmin,xmax=xmax,qc=qc,cc=cc,aux=aux,&
      yOut=yOut(i,:),xOut=xOut(i),qcOut=qcOut(i,:),err=err,message=message)
    if(err/=EXIT_SUCCESS)then;message='f-'//procnam//'&'//message;return;endif
  elseif(present(ccOut))then                ! only ccOut is present
    call operateur(opID=opID,y=y,x=x,xmin=xmin,xmax=xmax,qc=qc,cc=cc,aux=aux,&
      yOut=yOut(i,:),xOut=xOut(i),ccOut=ccOut(i,:),err=err,message=message)
    if(err/=EXIT_SUCCESS)then;message='f-'//procnam//'&'//message;return;endif
  else                                      ! none are present
    call operateur(opID=opID,y=y,x=x,xmin=xmin,xmax=xmax,qc=qc,cc=cc,aux=aux,&
      yOut=yOut(i,:),xOut=xOut(i),err=err,message=message)
    if(err/=EXIT_SUCCESS)then;message='f-'//procnam//'&'//message;return;endif
  endif
enddo
! end procedure here
endsubroutine aggregator_intervals

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure subroutine aggregator_grid(operateur,opID,y,x,grid,qc,cc,aux,& !IN
                                     yOut,xOut,qcOut,ccOut,err,message) !OUT
! Purpose: Aggregator engine
! Programmer: Ben Renard, Irstea & the BigD Kavetski
! Last modified: 10/05/2015
!*********************************************************************
! References:
!*********************************************************************
! IN
!   1. operateur, subroutine implementing the operator
!   2. opID, the properties of the operator to apply:
!      (/interpolation option, aggregation function, qc option, cc option/)
!   3. y, the function sampled at discrete points
!   4. x, the sampling points
!   5. grid, the grid defining the intervals where the operator is applied
!   6. [qc], quality codes
!   7. [cc], continuity codes
!   8. [aux], auxiliaries
! OUT
!   1. yout, values after applying the operator
!   2. xout, x-values associated with yout
!   3. [qcOut], quality codes associated with yout
!   4. [ccOut], continuity codes associated with yout
!   5. err, error code
!   6. message
!*********************************************************************
! 2Do List:
!*********************************************************************
use utilities_dmsl_kit,only:assertEq
implicit none
interface
  pure subroutine operateur(opID,y,x,xmin,xmax,qc,cc,aux,&
                            yOut,xOut,qcOut,ccOut,err,message)
  use kinds_dmsl_kit
  implicit none
  integer(mik),intent(in)::opID(:)
  real(mrk),intent(in)::y(:),x(:),xmin,xmax
  integer(mik),intent(in),optional::qc(:),cc(:)
  real(mrk),intent(in),optional::aux(:)
  real(mrk),intent(out)::yout(:),xout
  integer(mik),intent(out),optional::qcOut(:),ccOut(:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine operateur
endinterface
! dummies
integer(mik),intent(in)::opID(:)
real(mrk),intent(in)::y(:),x(:),grid(:)
integer(mik),intent(in),optional::qc(:),cc(:)
real(mrk),intent(in),optional::aux(:)
real(mrk),intent(out)::yout(:,:),xout(:)
integer(mik),intent(out),optional::qcOut(:,:),ccOut(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam='aggregator_grid'
real(mrk)::intervals(size(grid)-1,2)
integer(mik)::n

n=size(grid)
intervals(:,1)=grid(1:n-1);intervals(:,2)=grid(2:n)
call aggregator_intervals(operateur,opID,y,x,intervals,qc,cc,aux,& !IN
                      yOut,xOut,qcOut,ccOut,err,message)
if(err/=EXIT_SUCCESS)then;message='f-'//procnam//'&'//message;return;endif

end subroutine aggregator_grid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! private subs

pure subroutine getActiveIndices(x,xmin,xmax,opID,qc,ix,abortFlag)
! Purpose: Determine active indices ix, such that
! [xmin xmax] is inside [x(ix(1)) x(ix(2))] and ix(1) and ix(2) are non-missing.
use utilities_dmsl_kit,only:huntloc,arthsi
implicit none
! dummies
real(mrk),intent(in)::x(:),xmin,xmax
integer(mik),intent(in),optional::qc(:)
integer(mik),intent(in)::opID(:)
integer(mik),pointer::ix(:)
logical(mlk),intent(out)::abortFlag
! locals
logical(mlk),allocatable::qcMask(:)
integer(mik)::ix_first,ix_last,ix_nAll,ix_n
integer(mik)::i,n
! start procedure here
n=size(x); if(associated(ix))nullify(ix); abortFlag=.false.
ix_first=huntloc(x,xmin)              ! index of last x smaller than xmin
ix_last =huntloc(x,xmax,ix_first)     ! index of first x larger than xmax
if(ix_last<=0)then ! can't do anything
  abortFlag=.true.; return
endif
selectcase(opID(ix_interpol))
case(NONE_interpol)
  if(x(ix_last)==xmax) ix_last=ix_last-1 ! step back if x on the upper bound
case(LINEAR_interpol)
  if(x(ix_last)<xmax) ix_last=ix_last+1 ! unless x is on the upper bound, step forward
  if(present(qc))then   ! extend either side to hit non-missing values
    do i=ix_first,1,-1; if(qc(i)==QC_PRESENT)exit; enddo; ix_first=i
    do i=ix_last ,n,+1; if(qc(i)==QC_PRESENT)exit; enddo; ix_last=i
  endif
endselect
ix_nAll=ix_last-ix_first+1
if(ix_first<=0.or.ix_last>n.or.ix_nAll==0)then ! can't do anything
  if(associated(ix)) nullify(ix)
  if(allocated(qcMask)) deallocate(qcMask)
  abortFlag=.true.; return
endif
if(present(qc))then ! pack the indices ix to exclude missing values
  allocate(qcMask(ix_first:ix_last))
  qcMask=(qc(ix_first:ix_last)==QC_PRESENT);ix_n=count(qcmask)
  selectcase(opID(ix_qc))
  case(INFECT_qc)   ! abort if a missing value encountered in [ix_first,ix_last]
    if(.not.all(qcmask))then
        if(associated(ix)) nullify(ix)
        if(allocated(qcMask)) deallocate(qcMask)
        abortFlag=.true.; return
    endif
  case(IGNORE_qc)   ! abort only if all values are missing
    if(ix_n==0)then
      if(associated(ix)) nullify(ix)
      if(allocated(qcMask)) deallocate(qcMask)
      abortFlag=.true.; return
    endif
  endselect
  allocate(ix(ix_n))
  ix=pack(arthsi(ix_nAll),qcmask)
  deallocate(qcMask)
else                ! simpler case: no missing values
  ix_n=ix_nAll; allocate(ix(ix_n)); ix=arthsi(ix_n)
endif
ix=ix_first-1+ix
! end procedure here
endsubroutine getActiveIndices
   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pure function quickTrap(x,y)result(res)
! Purpose: Trapezoidal area
use utilities_dmsl_kit,only:zero,half
implicit none
! dummies
real(mrk),intent(in)::x(:),y(:)
real(mrk)::res
! locals
integer(mik)::i,n
! start procedure here
n=size(x); res=zero
do i=1,n-1
  res=res+half*(x(i+1)-x(i))*(y(i+1)+y(i))
enddo
! end procedure here
endfunction quickTrap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endmodule aggregation_tools

module HistoricBay_tools

!~**********************************************************************
!~* Purpose: Tools to perform flood frequency analyses with historical data
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified:01/08/2011
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1.
!~*    2.
!~*    3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use BayesianEstimation_tools, only:PriorListType

implicit none
Private
public :: ReadHistoricData,HistoricBay_EstimPar,HistoricBay_Predictive,FatalExit

type, public:: HistoDataType ! object for ONE historical value
    integer(mik)::DataType=undefIN ! Type of the historical data, defined as follows:
                           ! 0 = data equal to ...  (up to multiplicative systematic RC error, see below)
                           ! 1 = data smaller than ...
                           ! 2 = data greater than ...
                           ! 3 = data between ...
    real(mik)::Indx=undefIN ! indexing value (eg year or i)
    real(mrk):: IsEqualTo=undefRN ! Value for type-0 data; ignored if data is of any other type
    real(mrk):: IsSmallerThan=undefRN ! Value for type-1 data; ignored if data is of any other type
    real(mrk):: IsLargerThan=undefRN ! Value for type-2 data; ignored if data is of any other type
    real(mrk), dimension(2):: IsBetween=(/undefRN,undefRN/) ! Values for type-3 data; ignored if data is of any other type
    integer(mik):: SystErr=undefIN ! handling of multiplicative systematic error. Convention as follows:
                                   ! SystErr=0 means no error
                                   ! SystErr=k means the observed discharge is affected by the kth systematic error
                                   ! (in general, because this observation is computed with the kth rating curve)
    real(mrk)::QD=undefRN ! Highly experimental, be careful... QD is a daily discharge
                          ! Basically, the idea is the following: historical data are maximum stages, and therefore
                          ! correspond to peak (instantaneous) discharges. However, if some daily discharges are available
                          ! concomitently with peak discharge, one may try to jointly estimate the distribution of
                          ! (Qp,QD) by assuming a multiplicative "peak coefficient", assumed to be a realization of a
                          ! Gaussian distribution N(mu,sigma). mu and sigma are then estimated as part of the inferrence.
                          ! Note that if no QD are provided, only the distribution of Qp is estimated
end type HistoDataType

type, public:: HistoSampleType ! object for a sample of several historical values
    integer(mik)::nTot=undefIN ! total number of values in sample
    integer(mik)::n0=undefIN ! number of values of type 0 (see HistoDataType above)
    integer(mik)::n1=undefIN ! number of values of type 1 (see HistoDataType above)
    integer(mik)::n2=undefIN ! number of values of type 2 (see HistoDataType above)
    integer(mik)::n3=undefIN ! number of values of type 3 (see HistoDataType above)
    integer(mik)::nSystErr=undefIN ! number of systematic errors (in general, corresponding to the number of rating curves)
    integer(mik)::nQD=undefIN ! number of available QD values (see HistoDataType above)
    character(250)::nickname='AintGotNoName'
    type(HistoDataType), allocatable::Hsample(:) ! sample of historical values
    real(mrk)::mv=-9999._mrk ! value for missing data
end type HistoSampleType

type,public::HistoricBayType ! type containing all needed stuff for inference (data, model, prior, etc.)
    type(HistoSampleType)::X ! data, stored in an object HistoSampleType
    character(250)::distID='UnknownDist' ! Distribution of Data.
                               ! Any distribution (available in distribution_tools) can be used if there are no
                               ! systematic multiplicative errors; however, so far, only Gumbel/GEV are accepted if systematic
                               ! errors are included (because some distribution-specific computations are needed).
    integer(mik)::npar=undefRN ! number of parameters of the distribution above
    type(PriorListType), pointer::DistPrior(:),SystErrPrior(:) ! priors
    type(PriorListType)::DToPeakPrior(2) ! priors
    character(250)::nickname='AintGotNoName' ! name for this inference
end type HistoricBayType

! Object "HistoricBayType" globally available within this module
type(HistoricBayType)::Jbay

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadHistoricData(file,& ! path to file
                                 nskip,& ! number of headers row to be skipped
                                 nrow,ncol,& ! number of row/columns (excluding headers)
                                 IndxCol,& ! column containing the indexing variable (eg years)
                                 DataTypeCol,& ! column containing the data type (0 to 3)
                                 IsEqualToCol,& ! column containing value for type-0 data
                                 IsSmallerThanCol,& ! column containing value for type-1 data
                                 IsLargerThanCol,& ! column containing value for type-2 data
                                 IsBetweenCol,& ! 2 columns containing value for type-3 data
                                 QDcol,& ! column containing QD values (only works if corresponding Qpeak is of type 0)
                                 SystErrCol,& ! column containing the systematic error number (0 for no error)
                                 mv,& ! value denoting missing data
                                 nickname,& ! name for this dataset
                                 X,& ! data stored in an object HistoSampleType
                                 err,mess)
!^**********************************************************************
!^* Purpose: Read file containing historic data
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:01/08/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, path to file
!^*    2.nskip, number of headers row to be skipped
!^*    3.nrow, number of rows (excluding headers)
!^*    4.ncol, number of columns
!^*    5.IndxCol, column containing the indexing variable (eg years)
!^*    6.DataTypeCol, column containing the data type (0 to 3)
!^*    7.IsEqualToCol, column containing value for type-0 data
!^*    8.IsSmallerThanCol, column containing value for type-1 data
!^*    9.IsLargerThanCol, column containing value for type-2 data
!^*    10.IsBetweenCol, 2 columns containing value for type-3 data
!^*    11.<QDcol>, column containing QD values (only works if corresponding Qpeak is of type 0)
!^*    12.<SystErrCol>, column containing the systematic error number (0 for no error)
!^*    13.mv, value denoting missing data
!^*    14.<nickname>, name for this dataset
!^* OUT
!^*    1.X, data stored in an object HistoSampleType
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:GetSpareUnit

character(*), intent(in)::file
integer(mik), intent(in)::nskip,nrow,ncol,IndxCol,DataTypeCol,IsEqualToCol,&
                          IsSmallerThanCol,IsLargerThanCol,IsBetweenCol(2)
integer(mik), intent(in), optional::QDcol,SystErrCol
character(*), intent(in), optional::nickname
real(mrk), intent(in)::mv
type(HistoSampleType), intent(out)::X
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::unt,i,foo(1)
real(mrk)::dummyV(ncol)
err=0;mess=''

call getspareunit(unt,err=err,message=mess)
if(err>0) then
    mess="ReadHistoricData: "//trim(mess);return
endif

open(unit=unt,file=file,status='old')
do i=1,Nskip
    read(unt,*)
enddo
! Start populating X
X%ntot=nrow
X%mv=mv
if(allocated(X%Hsample)) deallocate(X%Hsample)
allocate(X%Hsample(nrow))

do i=1,Nrow
    read(unt,*) dummyV
    X%Hsample(i)%Indx=dummyV(IndxCol)
    X%Hsample(i)%DataType=int(dummyV(DataTypeCol),mik)
    X%Hsample(i)%IsEqualTo=dummyV(IsEqualToCol)
    X%Hsample(i)%IsSmallerThan=dummyV(IsSmallerThanCol)
    X%Hsample(i)%IsLargerThan=dummyV(IsLargerThanCol)
    X%Hsample(i)%IsBetween=dummyV(IsBetweenCol)
    if(present(SystErrCol)) then
        X%Hsample(i)%SystErr=int(dummyV(SystErrCol),mik)
    else
        X%Hsample(i)%SystErr=0
    endif
    if(present(QDCol)) then
        X%Hsample(i)%QD=dummyV(QDCol)
    else
        X%Hsample(i)%QD=mv
    endif
enddo
close(unt)
X%n0=count(X%Hsample(:)%DataType==0)
X%n1=count(X%Hsample(:)%DataType==1)
X%n2=count(X%Hsample(:)%DataType==2)
X%n3=count(X%Hsample(:)%DataType==3)
foo=maxval(X%Hsample(:)%SystErr)
X%nSystErr=foo(1)
X%nQD=count(X%Hsample(:)%QD/=mv)
if(present(nickname)) X%nickname=nickname

end subroutine ReadHistoricData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine GetLogLkh_Histo(Hbay,DistPar,DtoPeakPar,SystErrPar,LL, feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute the log-likelihood of historical data
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:01/08/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. Hbay, object of HistoricBayType containg sample, model, priors, etc.
!^*    2. DistPar, parameters of the distribution of peak discharges
!^*    3. <DtoPeakPar>, parameters of the daily-to-peak ratio, ie mean and
!^*    stdev of the ratio. only requested if at least one QD is provided
!^*    4. <SystErrPar>, systematic errors, only requested if at least one
!^*    historical data is affected by a systematic (eg RC) error
!^* OUT
!^*    1.LL; log-likelihood
!^*    2.feas; are parameters feasible?
!^*    3.isnull; Is likelihood==0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use Distribution_tools, only:GetParNumber, GetPdf, GetCdf, GAUSS

type(HistoricBayType), intent(in)::Hbay
real(mrk), intent(in)::DistPar(:)
real(mrk), intent(in), optional::DtoPeakPar(:),SystErrPar(:)
real(mrk), intent(out)::LL
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::npar,i
real(mrk)::par(size(DistPar)),lkh,cdf1,cdf2,Qpeak

err=0;mess='';feas=.true.;isnull=.false.;LL=undefRN

! Size & presence checks
if( Hbay%X%nSystErr>0 .and. .not.present(SystErrPar) ) then
    err=1;mess='GetLogLkh_Histo:fatal: <SystErrPar> not present but nSystErr/=0'
    feas=.false.;return
endif
if( Hbay%X%nQD>0 .and. .not.present(DtoPeakPar) ) then
    err=2;mess='GetLogLkh_Histo:fatal: <DtoPeakPar> not present but nQD/=0'
    feas=.false.;return
endif
if( present(DtoPeakPar) ) then
    if(size(DtoPeakPar)/=2) then
        err=3;mess='GetLogLkh_Histo:fatal: <DtoPeakPar> should be of size 2 (mean/std of daily-to-peak ratio)'
        feas=.false.;return
    endif
endif
CALL GetParNumber(Hbay%distID, npar, err, mess)
If(err>0) then
    mess='GetLogLkh_Histo:'//trim(mess)
    feas=.false.;return
endif
if(size(DistPar)/=npar) then
    feas=.false.;err=4;mess="GetLogLkh_Histo:FATAL:size mismatch [DistPar]";
    return
endif

LL=0._mrk
! First compute likelihood of peak discharges
do i=1,Hbay%X%nTot
    ! Handle systematic error
    par=DistPar
    if (Hbay%X%Hsample(i)%SystErr>0) then
        if(Hbay%X%Hsample(i)%SystErr>size(SystErrPar)) then
            err=5;mess="GetLogLkh_Histo:FATAL:SystErr largen than size(SystErrPar)"
            feas=.false.;return
        else
            select case (trim(Hbay%distID))
            case('GEV','Gumbel')
                ! the following should also makes sense for anyother distribution where par(1:2) are location/scale parameters.
                par(1)=DistPar(1)/SystErrPar(Hbay%X%Hsample(i)%SystErr)
                par(2)=DistPar(2)/SystErrPar(Hbay%X%Hsample(i)%SystErr)
            case('LogNormal')
                ! uses the property: if X~LN(m,s) then aX~LN(m+log(a),s)
                if(SystErrPar(Hbay%X%Hsample(i)%SystErr)<=0._mrk) then;feas=.false.;return;endif
                par(1)=DistPar(1)-log(SystErrPar(Hbay%X%Hsample(i)%SystErr))
            case default
                err=100;mess="GetLogLkh_Histo:FATAL:unsupported distribution"
                feas=.false.;return
            end select
        endif
    elseif(Hbay%X%Hsample(i)%SystErr<0) then
        err=6;mess="GetLogLkh_Histo:FATAL:negative SystErr"
        feas=.false.;return
    endif
    ! Compute contribution to lkh
    select case(Hbay%X%Hsample(i)%DataType)
    case(0)
        call GetPdf(DistId=Hbay%DistID,x=Hbay%X%Hsample(i)%IsEqualTo,par=par,loga=.true.,&
                    pdf=lkh,feas=feas,isnull=isnull,err=err,mess=mess)
        If(err>0) then
            mess='GetLogLkh_Histo:'//trim(mess)
            feas=.false.;return
        endif
        if( (.not.feas) .or. isnull) return
    case(1)
        call GetCdf(DistId=Hbay%DistID,x=Hbay%X%Hsample(i)%IsSmallerThan,par=par,&
                    cdf=cdf1,feas=feas,err=err,mess=mess)
        If(err>0) then
            mess='GetLogLkh_Histo:'//trim(mess)
            feas=.false.;return
        endif
        if(.not.feas) return
        if(cdf1<=0._mrk) then
            isnull=.true.;return
        else
            lkh=log(cdf1)
        endif
    case(2)
        call GetCdf(DistId=Hbay%DistID,x=Hbay%X%Hsample(i)%IsLargerThan,par=par,&
                    cdf=cdf1,feas=feas,err=err,mess=mess)
        If(err>0) then
            mess='GetLogLkh_Histo:'//trim(mess)
            feas=.false.;return
        endif
        if(.not.feas) return
        if(cdf1>=1._mrk) then
            isnull=.true.;return
        else
            lkh=log(1._mrk-cdf1)
        endif
    case(3)
        ! check IsBetween(2)>IsBetween(1)
        if(Hbay%X%Hsample(i)%IsBetween(2)<=Hbay%X%Hsample(i)%IsBetween(1)) then
            err=7;mess="GetLogLkh_Histo:FATAL:Higher bound smaller than lower bound for a type-3 data";
            feas=.false.;return
        endif
        call GetCdf(DistId=Hbay%DistID,x=Hbay%X%Hsample(i)%IsBetween(1),par=par,&
                    cdf=cdf1,feas=feas,err=err,mess=mess)
        If(err>0) then
            mess='GetLogLkh_Histo:'//trim(mess)
            feas=.false.;return
        endif
        if(.not.feas) return
        call GetCdf(DistId=Hbay%DistID,x=Hbay%X%Hsample(i)%IsBetween(2),par=par,&
                    cdf=cdf2,feas=feas,err=err,mess=mess)
        If(err>0) then
            mess='GetLogLkh_Histo:'//trim(mess)
            feas=.false.;return
        endif
        if(.not.feas) return
        if( (cdf2-cdf1)<=0._mrk) then
            isnull=.true.;return
        else
            lkh=log(cdf2-cdf1)
        endif
    case default
        err=8;mess="GetLogLkh_Histo:FATAL:unknown DataType [should be between 0 and 3]";
        feas=.false.;return
    end select

    LL=LL+lkh

    ! Add contribution of QD to likelihood (if any)
    if(Hbay%X%nQD>0) then
        if(Hbay%X%Hsample(i)%QD/=Hbay%X%mv) then
            if(Hbay%X%Hsample(i)%DataType/=0) then
                err=-9;mess="GetLogLkh_Histo:WARNING:QD values ignored since corresponding Qpeak is of type 0";
                !feas=.false.;return
            elseif(Hbay%X%Hsample(i)%IsEqualTo==Hbay%X%mv) then
                err=-10;mess="GetLogLkh_Histo:WARNING:QD values ignored since corresponding Qpeak is missing";
                !feas=.false.;return
            else
                Qpeak=Hbay%X%Hsample(i)%IsEqualTo
                call GetPdf(DistId=GAUSS,x=Hbay%X%Hsample(i)%QD,&
                            par=(/Qpeak*DtoPeakPar(1),Qpeak*DtoPeakPar(2)/),&
                            loga=.true.,pdf=lkh,feas=feas,isnull=isnull,err=err,mess=mess)
                If(err>0) then
                    mess='GetLogLkh_Histo:'//trim(mess)
                    feas=.false.;return
                endif
                if( (.not.feas) .or. isnull) return
                LL=LL+lkh
            endif
        endif
    endif
enddo

end subroutine GetLogLkh_Histo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine GetLogPrior_Histo(DistPar,DtoPeakPar,SystErrPar,&
                                  DistPrior,DtoPeakPrior,SystErrPrior,&
                                  prior,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Get prior
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 02/08/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. DistPar, parameters of the distribution of peak discharges
!^*    2. <DtoPeakPar>, parameters of the daily-to-peak ratio, ie mean and
!^*    stdev of the ratio. only requested if at least one QD is provided
!^*    3. <SystErrPar>, systematic errors, only requested if at least one
!^*    historical data is affected by a systematic (eg RC) error
!^*    4. DistPrior, priors for DistPar
!^*    5. <DtoPeakPrior>, priors for DtoPeakPar
!^*    6. <SystErrPrior>, priors for SystErrPar
!^* OUT
!^*    1.prior; log-pdf of the prior
!^*    2.feas; are parameters feasible?
!^*    3.isnull; Is prior==0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use BayesianEstimation_tools, only:priorlistType
use distribution_tools, only:GetParNumber,GetPdf

real(mrk), intent(in)::DistPar(:)
real(mrk), intent(in), optional::DtoPeakPar(:),SystErrPar(:)
type(PriorListType), intent(in)::DistPrior(:)
type(PriorListType), intent(in), optional::DtoPeakPrior(:),SystErrPrior(:)
real(mrk), intent(out)::prior
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,npar
real(mrk)::pr

err=0;mess='';feas=.true.;isnull=.false.;prior=undefRN

! Size & presence checks
if( present(DtoPeakPar) ) then
    if(.not.present(DtoPeakPrior)) then
        err=1;mess='GetLogPrior_Histo:fatal: <DtoPeakPrior> not present but DtoPeakPar present'
        feas=.false.;return
    else
        if(size(DtoPeakPar)/=size(DtoPeakPrior)) then
            err=2;mess='GetLogPrior_Histo:fatal: size mismatch <DtoPeakPrior> / <DtoPeakPar>'
            feas=.false.;return
        endif
    endif
endif
if( present(SystErrPar)) then
    if(.not.present(SystErrPrior)) then
        err=3;mess='GetLogPrior_Histo:fatal: <SystErrPrior> not present but SystErrPar present'
        feas=.false.;return
    else
        if(size(SystErrPar)/=size(SystErrPrior)) then
            err=4;mess='GetLogPrior_Histo:fatal: size mismatch <SystErrPrior> / <SystErrPar>'
            feas=.false.;return
        endif
    endif
endif
if(size(DistPar)/=size(DistPrior)) then
    err=5;mess='GetLogPrior_Histo:fatal: size mismatch <DistPrior> / <DistPar>'
    feas=.false.;return
endif

! Start computing here
prior=0._mrk
! contribution of DistPar
do i=1,size(DistPar)
    ! check size
    if(.not.allocated(DistPrior(i)%par)) then
        err=6;mess='GetLogPrior_Histo:fatal: unallocated <Prior%Par>'
        feas=.false.;return
    endif
    call GetParNumber(DistID=DistPrior(i)%dist, npar=npar, err=err, mess=mess)
    if(err>0) then
        mess="GetLogPrior_Histo: "//trim(mess);return
    endif
    if(size(DistPrior(i)%par)/=npar) then
        err=7;mess='GetLogPrior_Histo:fatal: size mismatch <Prior%dist> / <Prior%Par>'
        feas=.false.;return
    endif
    call GetPdf(DistId=DistPrior(i)%dist,x=DistPar(i),par=DistPrior(i)%par,&
                loga=.true.,pdf=pr,feas=feas,isnull=isnull,err=err,mess=mess)
    If(err>0) then
        mess='GetLogPrior_Histo:'//trim(mess)
        feas=.false.;return
    endif
    if( (.not.feas) .or. isnull) return
    prior=prior+pr
enddo

! contribution of DtoPeakPar
if(present(DtoPeakPar)) then
    do i=1,size(DtoPeakPar)
        ! check size
        if(.not.allocated(DtoPeakPrior(i)%par)) then
            err=6;mess='GetLogPrior_Histo:fatal: unallocated <Prior%Par>'
            feas=.false.;return
        endif
        call GetParNumber(DistID=DtoPeakPrior(i)%dist, npar=npar, err=err, mess=mess)
        if(err>0) then
            mess="GetLogPrior_Histo: "//trim(mess);return
        endif
        if(size(DtoPeakPrior(i)%par)/=npar) then
            err=7;mess='GetLogPrior_Histo:fatal: size mismatch <Prior%dist> / <Prior%Par>'
            feas=.false.;return
        endif
        call GetPdf(DistId=DtoPeakPrior(i)%dist,x=DtoPeakPar(i),par=DtoPeakPrior(i)%par,&
                    loga=.true.,pdf=pr,feas=feas,isnull=isnull,err=err,mess=mess)
        If(err>0) then
            mess='GetLogPrior_Histo:'//trim(mess)
            feas=.false.;return
        endif
        if( (.not.feas) .or. isnull) return
        prior=prior+pr
    enddo
endif

! contribution of SystErrPar
if(present(SystErrPar)) then
    do i=1,size(SystErrPar)
        ! check size
        if(.not.allocated(SystErrPrior(i)%par)) then
            err=6;mess='GetLogPrior_Histo:fatal: unallocated <Prior%Par>'
            feas=.false.;return
        endif
        call GetParNumber(DistID=SystErrPrior(i)%dist, npar=npar, err=err, mess=mess)
        if(err>0) then
            mess="GetLogPrior_Histo: "//trim(mess);return
        endif
        if(size(SystErrPrior(i)%par)/=npar) then
            err=7;mess='GetLogPrior_Histo:fatal: size mismatch <Prior%dist> / <Prior%Par>'
            feas=.false.;return
        endif
        call GetPdf(DistId=SystErrPrior(i)%dist,x=SystErrPar(i),par=SystErrPrior(i)%par,&
                    loga=.true.,pdf=pr,feas=feas,isnull=isnull,err=err,mess=mess)
        If(err>0) then
            mess='GetLogPrior_Histo:'//trim(mess)
            feas=.false.;return
        endif
        if( (.not.feas) .or. isnull) return
        prior=prior+pr
    enddo
endif

end subroutine GetLogPrior_Histo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine GetLogPost_Histo(DistPar,DtoPeakPar,SystErrPar,& !parameters
                                 Hbay,& ! data, model & priors in object of type HistoricBayType
                                 post,feas,isnull,err,mess) !outputs
!^**********************************************************************
!^* Purpose: Compute log-posterior
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 02/08/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^*    1. DistPar, parameters of the distribution of peak discharges
!^*    2. <DtoPeakPar>, parameters of the daily-to-peak ratio, ie mean and
!^*    stdev of the ratio. only requested if at least one QD is provided
!^*    3. <SystErrPar>, systematic errors, only requested if at least one
!^*    historical data is affected by a systematic (eg RC) error
!^*    4. Hbay, object of HistoricBayType containing sample, model, priors, etc.
!^* OUT
!^*    1.post; log-pdf of the posterior
!^*    2.feas; are parameters feasible?
!^*    3.isnull; Is likelihood==0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
type(HistoricBayType), intent(in)::Hbay
real(mrk), intent(in)::DistPar(:)
real(mrk), intent(in), optional::DtoPeakPar(:),SystErrPar(:)
real(mrk), intent(out)::post
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk)::prior,lkh

err=0;mess='';feas=.true.;isnull=.false.;post=undefRN

! Get prior
call GetLogPrior_Histo(DistPar=DistPar,DtoPeakPar=DtoPeakPar,SystErrPar=SystErrPar,&
                       DistPrior=Hbay%DistPrior,&
                       DtoPeakPrior=Hbay%DtoPeakPrior,&
                       SystErrPrior=Hbay%SystErrPrior,&
                       prior=prior,feas=feas,isnull=isnull,err=err,mess=mess)
If(err>0) then
    mess='GetLogPost_Histo:'//trim(mess)
    feas=.false.;return
endif
if( (.not.feas) .or. isnull) return

! Get likelihood
call GetLogLkh_Histo(Hbay=Hbay,&
                     DistPar=DistPar,DtoPeakPar=DtoPeakPar,SystErrPar=SystErrPar,&
                     LL=lkh,feas=feas,isnull=isnull,err=err,mess=mess)
If(err>0) then
    mess='GetLogPost_Histo:'//trim(mess)
    feas=.false.;return
endif
if( (.not.feas) .or. isnull) return

! Get posterior
post=prior+lkh

end subroutine GetLogPost_Histo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Posterior_wrapper(x,feas,isnull,fx,fAux,err,mess)
!^**********************************************************************
!^* Purpose: wrapper to GetLogPost_Histo to comply with MCMC interface
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 02/08/2011
!^**********************************************************************
!^* Comments: mv handled out of here
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, teta
!^* OUT
!^*    1.feas; is teta feasible?
!^*    2.isnull; Is posterior==0?
!^*    3.fx; posterior value
!^*    4.fAux; Auxiliaries (by-product of posterior computation if any)
!^*    5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    6.mess, error message
!^**********************************************************************

real(mrk),intent(in)::x(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
! locals
real(mrk), allocatable:: DistPar(:),DtoPeakPar(:),SystErrPar(:)

err=0;feas=.true.;isnull=.false.;fx=undefRN;mess=''

if( Jbay%X%nSystErr==0 .and. Jbay%X%nQD==0) then
    if(size(x)==Jbay%npar) then
       call GetLogPost_Histo(DistPar=x(1:Jbay%npar),& !parameters
                          Hbay=Jbay,& ! data, model & priors in object of type HistoricBayType
                          post=fx,feas=feas,isnull=isnull,err=err,mess=mess) !outputs
    else
        call GetLogPost_Histo(DistPar=x(1:Jbay%npar),DtoPeakPar=x( (Jbay%npar+1):(Jbay%npar+2) ),& !parameters
                          Hbay=Jbay,& ! data, model & priors in object of type HistoricBayType
                          post=fx,feas=feas,isnull=isnull,err=err,mess=mess) !outputs
    endif
elseif(Jbay%X%nSystErr==0) then
    call GetLogPost_Histo(DistPar=x(1:Jbay%npar),DtoPeakPar=x( (Jbay%npar+1):(Jbay%npar+2) ),& !parameters
                          Hbay=Jbay,& ! data, model & priors in object of type HistoricBayType
                          post=fx,feas=feas,isnull=isnull,err=err,mess=mess) !outputs
elseif(Jbay%X%nQD==0) then
    if(size(x)==(Jbay%npar+Jbay%X%nSystErr)) then
        call GetLogPost_Histo(DistPar=x(1:Jbay%npar),&
                          SystErrPar=x( (Jbay%npar+1):(Jbay%npar+Jbay%X%nSystErr) ),& !parameters
                          Hbay=Jbay,& ! data, model & priors in object of type HistoricBayType
                          post=fx,feas=feas,isnull=isnull,err=err,mess=mess) !outputs
    else
        call GetLogPost_Histo(DistPar=x(1:Jbay%npar),&
                          DtoPeakPar=x( (Jbay%npar+1):(Jbay%npar+2) ),&
                          SystErrPar=x( (Jbay%npar+3):(Jbay%npar+2+Jbay%X%nSystErr) ),& !parameters
                          Hbay=Jbay,& ! data, model & priors in object of type HistoricBayType
                          post=fx,feas=feas,isnull=isnull,err=err,mess=mess) !outputs
    endif
else
    call GetLogPost_Histo(DistPar=x(1:Jbay%npar),&
                          DtoPeakPar=x( (Jbay%npar+1):(Jbay%npar+2) ),&
                          SystErrPar=x( (Jbay%npar+2+1):(Jbay%npar+2+Jbay%X%nSystErr) ),& !parameters
                          Hbay=Jbay,& ! data, model & priors in object of type HistoricBayType
                          post=fx,feas=feas,isnull=isnull,err=err,mess=mess) !outputs
endif

if(present(fAux)) fAux=UndefRN

end subroutine Posterior_wrapper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HistoricBay_EstimPar(Hbay,& ! Object "HistoricBayType", containing data, model, priors, etc...
        ! Tuning of the MCMC sampler
        start,startStd,&
        nAdapt,nCycles,&
        MinMoveRate,MaxMoveRate,&
        DownMult,UpMult,&
        OutFile, &
        ! error handling
        err,mess)
!^**********************************************************************
!^* Purpose: Parameter estimation for Historical data with Bayes-MCMC
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 02/08/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.Hbay, Object "HistoricBayType", containing data, model, priors, etc...
!^*    2.start, starting point
!^*    3.startStd, starting std for each marginal jump distribution
!^*    4.nAdapt, number of iterations before adapting the jump std
!^*    5.nCycles, number of adaptations (total number of iteration is hence nAdapt*nCycles
!^*        If nCycles<1, user will be asked stop/continue
!^*    6.MinMoveRate, objective move rate (lower bound)
!^*    7.MaxMoveRate, objective move rate (upper bound)
!^*    8.DownMult, value multiplying the jump std whan move rate is too low
!^*    9.UpMult, value multiplying the jump std whan move rate is too high
!^*    10.OutFile, address of output file for MCMC samples
!^* OUT
!^*    1.
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************

use MCMCStrategy_tools, only:Adaptive_Metro_OAAT
use distribution_tools, only:GetParNumber

type(HistoricBayType), intent(in)::Hbay
real(mrk), intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
real(mrk), intent(inout)::start(:), startStd(:)
integer(mik), intent(in)::nAdapt,nCycles
character(*), intent(in)::OutFile
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::n0,n,npar,i, foo(1)
logical, allocatable::mvlist(:)
real(mrk)::fx

err=0;mess=''

! ----------------------------------------
! Part I: handle mv
! ----------------------------------------
n0=Hbay%X%nTot
if(allocated(mvlist)) deallocate(mvlist)
allocate(mvlist(n0))
mvlist= .not.(Hbay%X%Hsample(:)%DataType==0 .and. Hbay%X%Hsample(:)%IsEqualTo==Hbay%X%mv)
n=count(mvlist)
if(n==0) then
    err=1;mess='HistoricBay_EstimPar:FATAL: no non-missing values!';return
endif

! ----------------------------------------
! part II: populate global variable Jbay
! ----------------------------------------
! Model general setup
Jbay%nickname=Hbay%nickname
Jbay%distID=Hbay%distID
CALL GetParNumber(Jbay%distID, npar, err, mess)
if(err>0) then
    mess="HistoricBay_EstimPar: "//trim(mess);return
endif
Jbay%npar=npar
! Data
Jbay%X%nTot=n
Jbay%X%mv=Hbay%X%mv
if(allocated(Jbay%X%Hsample)) deallocate(Jbay%X%Hsample)
allocate(Jbay%X%Hsample(n))
Jbay%X%Hsample=pack(Hbay%X%Hsample,mvlist)
Jbay%X%n0=count(Jbay%X%Hsample(:)%DataType==0)
Jbay%X%n1=count(Jbay%X%Hsample(:)%DataType==1)
Jbay%X%n2=count(Jbay%X%Hsample(:)%DataType==2)
Jbay%X%n3=count(Jbay%X%Hsample(:)%DataType==3)
foo=maxval(Jbay%X%Hsample(:)%SystErr)
Jbay%X%nSystErr=foo(1)
Jbay%X%nQD=count(Jbay%X%Hsample(:)%QD/=Jbay%X%mv)
Jbay%X%nickname=Hbay%X%nickname
! priors for DistPar
if (Jbay%npar/=Hbay%npar) then
    err=1;mess='HistoricBay_EstimPar:FATAL: size mismatch <Jbay%npar> / <Hbay%npar>'
    return
endif
if(associated(Jbay%DistPrior)) deallocate(Jbay%DistPrior)
allocate(Jbay%DistPrior(Jbay%npar))
do i=1,Jbay%npar
    Jbay%DistPrior(i)%dist=Hbay%DistPrior(i)%dist
    npar=size(Hbay%DistPrior(i)%par)
    if(allocated(Jbay%DistPrior(i)%par)) deallocate(Jbay%DistPrior(i)%par)
    allocate(Jbay%DistPrior(i)%par(npar))
    Jbay%DistPrior(i)%par=Hbay%DistPrior(i)%par
enddo

! priors for DtoPeak
!if (Jbay%X%nQD/=Hbay%X%nQD) then
!    err=2;mess='HistoricBay_EstimPar:FATAL: size mismatch <Jbay%X%nQD> / <Hbay%X%nQD>'
!    return
!endif
do i=1,2
    Jbay%DtoPeakPrior(i)%dist=Hbay%DtoPeakPrior(i)%dist
    npar=size(Hbay%DtoPeakPrior(i)%par)
    if(allocated(Jbay%DtoPeakPrior(i)%par)) deallocate(Jbay%DtoPeakPrior(i)%par)
    allocate(Jbay%DtoPeakPrior(i)%par(npar))
    Jbay%DtoPeakPrior(i)%par=Hbay%DtoPeakPrior(i)%par
enddo

! priors for SystErr
if (Jbay%X%nSystErr/=Hbay%X%nSystErr) then
    err=3;mess='HistoricBay_EstimPar:FATAL: size mismatch <Jbay%X%nSystErr> / <Hbay%X%nSystErr>'
    return
endif
if(associated(Jbay%SystErrPrior)) deallocate(Jbay%SystErrPrior)
if(Jbay%X%nSystErr>0) then
    allocate(Jbay%SystErrPrior(Jbay%X%nSystErr))
    do i=1,Jbay%X%nSystErr
        Jbay%SystErrPrior(i)%dist=Hbay%SystErrPrior(i)%dist
        npar=size(Hbay%SystErrPrior(i)%par)
        if(allocated(Jbay%SystErrPrior(i)%par)) deallocate(Jbay%SystErrPrior(i)%par)
        allocate(Jbay%SystErrPrior(i)%par(npar))
        Jbay%SystErrPrior(i)%par=Hbay%SystErrPrior(i)%par
    enddo
endif


! ----------------------------------------
! part III: MCMC
! ----------------------------------------

call Adaptive_Metro_OAAT(f=Posterior_wrapper,x=start,&
        fx=fx,std=startStd,&
        nAdapt=nAdapt,nCycles=nCycles,&
        MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
        DownMult=DownMult,UpMult=UpMult,&
        OutFile=OutFile,err=err,mess=mess)

end subroutine HistoricBay_EstimPar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HistoricBay_Predictive(File,Nburn,Nread,Nslim,& ! IN: Read properties
                Nrep, & ! IN: number of reps for EACH MCMC sample + data dist
                Hbay,& ! Info about the inference
                modal,pred,predQD,& ! OUT: modal par estimate & replicates from the predictive
                err, mess)
!^**********************************************************************
!^* Purpose: Analyse MCMC samples & get predictive
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 02/08/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. File, containing MCMC samples
!^*    2. Nburn, number of lines to discard
!^*    3. Nread, number of lines to read
!^*    4. Nslim, only one line every Nslim will be used
!^*    5. Nrep, for each line, Nrep deviates will be generated
!^*    6. Hbay, object "HistoricBayType" with all requested info
!^* OUT
!^*    1.modal, modal parameter estimate
!^*    2.pred, deviates from the predictive (of Qpeak)
!^*    3.<predQD>, deviates from the predictive of QD
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetParNumber, Generate, GAUSS
use utilities_dmsl_kit, only:getspareunit

integer(mik), intent(in)::Nburn, Nread, Nslim, Nrep
character(*),intent(in)::File
type(HistoricBayType), intent(in)::Hbay
real(mrk), intent(out)::modal(:)
real(mrk), pointer::pred(:)
real(mrk), pointer, optional::predQD(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::unt,errcode,n,k, ml(1),npred,i,j
real(mrk), allocatable::Temp(:,:), mcmc(:,:)
real(mrk)::ratio
logical::feas

err=0;mess=''

call getspareunit(unt,err=err,message=mess)
if(err>0) then
    mess="HistoricBay_Predictive: "//trim(mess);return
endif

open(unit=unt,file=File,status='old')
read(unt,*) !headers
do i=1,Nburn ! burn
    read(unt,*,iostat=errcode)
    if(errcode/=0) then ! problem reading file
        err=1
        mess='HistoricBay_Predictive:FATAL:problem reading file during burnin - cannot proceed!'
        return
    endif
enddo

! allocate Temp & read
if(allocated(Temp)) deallocate(Temp)
allocate(Temp(Nread,Hbay%npar+2+Hbay%X%nSystErr+1))
do i=1,Nread ! read used lines
    k=i
    read(unt,*,iostat=errcode) Temp(i,:)
    if(errcode/=0) then ! problem reading file
        err=-1;k=k-1
        mess='HistoricBay_Predictive:WARNING:problem reading file before [Nread] lines!'
        EXIT
    endif
enddo
close(unt)

! Slim
n=k/Nslim
if(allocated(mcmc)) deallocate(mcmc)
allocate(mcmc(n,Hbay%npar+2+Hbay%X%nSystErr+1))
mcmc=Temp(1:n:nSlim,:)

! get modal estimate
ml=Maxloc(mcmc(:,Hbay%npar+2+Hbay%X%nSystErr+1))
if(size(modal)/=Hbay%npar+2+Hbay%X%nSystErr) then
    err=1; mess='HistoricBay_Predictive:SEVERE:size mismath [modal]'
    modal=undefRN
else
    modal=mcmc(ml(1),1:Hbay%npar+2+Hbay%X%nSystErr)
endif

! get deviates from predictive
npred=n*nrep
if(associated(pred)) deallocate(pred)
allocate(pred(npred))
if(present(predQD))then
    if(associated(predQD)) deallocate(predQD)
    allocate(predQD(npred))
endif
k=0
do i=1,n
    do j=1,nrep
        k=k+1
        call Generate(DistId=Hbay%distID,par=mcmc(i,1:Hbay%npar),&
                gen=pred(k),feas=feas,err=err,mess=mess)
        If(err>0) then
            mess='HistoricBay_Predictive: '//trim(mess);return
        endif
        if(.not. feas) then ! shouldn't happen - MCMC should NOT generate infeasible par!
            err=1;
            mess='HistoricBay_Predictive: unfeasible parameter!';return
        endif
        if(present(predQD)) then
            call Generate(DistId=GAUSS,par=mcmc(i,(Hbay%npar+1):(Hbay%npar+2)),&
                gen=ratio,feas=feas,err=err,mess=mess)
            If(err>0) then
                mess='HistoricBay_Predictive: '//trim(mess);return
            endif
            if(.not. feas) then ! shouldn't happen - MCMC should NOT generate infeasible par!
                err=1;
                mess='HistoricBay_Predictive: unfeasible parameter!';return
            endif
            predQD(k)=pred(k)*ratio
        endif
    enddo
enddo

end subroutine HistoricBay_Predictive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FatalExit(message)
! exit program after a fatal error
character(*),intent(in)::message
write(*,*) trim(message)
read(*,*)
STOP
end subroutine FatalExit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module HistoricBay_tools

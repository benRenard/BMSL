module DimReduction_tools

!~**********************************************************************
!~* Purpose: PCA & other approaches
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!~**********************************************************************
!~* Last modified: 09/11/2012
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
public :: PCA,KernelPCA,S_MECA,S_MECA_LOOCV,S_MECA_ChooseNdim

type, public::pcaType
    real(mrk), allocatable::eigval(:),eigvec(:,:),proj(:,:),Xcenter(:),Xscale(:)
end type pcaType
type, public::mecaType
    real(mrk), allocatable::eigval(:),lambda(:,:),tau(:,:),psi(:,:),tauHat(:,:)
    real(mrk), allocatable::Ytilde(:,:),Yhat(:,:) ! reconstructed with PCA and ECA respectively
    real(mrk), allocatable::Ycenter(:),Yscale(:),Phicenter(:),Phiscale(:) ! values used to center and/or scale Y and Phi
end type mecaType
type, public::LOOCVtype
    real(mrk), allocatable::Y(:,:) ! original data
    real(mrk), allocatable::Yhat(:,:) ! reconstructed with ECA 
    real(mrk), allocatable::tcor(:),scor(:) ! temporal & spatial correlations between Yhat & Y
end type LOOCVtype

Contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PCA(X,DoCenter,DoStd,option,nDim,folder,pcaRes,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:09/11/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. X, data matrix
!^*		2. DoCenter, perform centering?
!^*		3. DoStd, perform Standardazing?
!^*		4. [option], computational variant: 'standard' or 'SVD', default 'standard'
!^*		   'standard' more efficient for "vertical" X (ie many rows, few columns)
!^*		   'SVD' more efficient for "horizontal" X (ie few rows, many columns)
!^*		5. [ndim], number of dimensions to keep; default = 4
!^*		6. [folder], result written 2 file if folder is provided
!^* OUT
!^*		1.pcaRes, object of type pcaType containing the results
!^*		2.err
!^*		3.mess
!^**********************************************************************
use numerix_dmsl_kit, only:getmeanvar,indexx_qsort
use linalg_dmsl_kit, only:svd_dcmp,getEigValsVecsRealSym
real(mrk), intent(inout)::X(:,:)
character(*), intent(in),optional::option,folder
integer(mik),intent(in),optional::ndim
logical, intent(in), optional::DoCenter,DoStd
type(pcaType), intent(out)::pcaRes
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
real(mrk), parameter::epsilon=0.000001_mrk
logical, parameter::DoC_def=.false.,DoS_def=.false.
character(250), parameter::opt_def='standard'
integer(mik), parameter::ND_def=4
logical::DoC,DoS
character(250)::opt
integer(mik)::N,M,i,indx(size(X,dim=1)),r,keep,ND
real(mrk)::mu(size(X,dim=2)),sig(size(X,dim=2)),A(size(X,dim=2),size(X,dim=1)),&
            W(size(X,dim=1)),V(size(X,dim=1),size(X,dim=1)),&
            C(size(X,dim=2),size(X,dim=2)),E(size(X,dim=2))
real(mrk),allocatable::U1(:,:),delta(:,:),V1(:,:),Y(:,:)

! init
err=0;mess=''
! Get Dimensions
N=size(X,dim=1);M=size(X,dim=2)
! Handle optional variables
if(present(DoCenter)) then;DoC=DoCenter;else;DoC=DoC_def;endif
if(present(DoStd)) then;DoS=DoStd;else;DoS=DoS_def;endif
if(present(option)) then;opt=option;else;opt=opt_def;endif
if(present(Ndim)) then;ND=Ndim;else;ND=ND_def;endif
! Center/Standardize if requested
if(allocated(pcaRes%Xcenter)) deallocate(pcaRes%Xcenter);allocate(pcaRes%Xcenter(M))
if(allocated(pcaRes%Xscale)) deallocate(pcaRes%Xscale);allocate(pcaRes%Xscale(M))
if(DoS .or. DoC) then
    call getmeanvar(x=X,mean=mu,var=sig,package=2,method='f',err=err,message=mess)
    if(err/=0) then;mess='PCA:'//trim(mess);return;endif
    if(DoS) then;sig=sqrt(sig);else;sig=1._mrk;endif
    if(.not.DoC) mu=0._mrk
    do i=1,M
        X(:,i)=(X(:,i)-mu(i))/sig(i)
    enddo
    pcaRes%Xcenter=mu;pcaRes%Xscale=sig
else
    pcaRes%Xcenter=0._mrk;pcaRes%Xscale=1._mrk
endif

select case(trim(opt))
    case('standard')
        ! Get Covariance matrix
        C=matmul(transpose(X),X)/real(n,mrk)
        ! eigenvalue decomposition
        call getEigValsVecsRealSym(a2evex=C,doEigVex=.true.,emeth=2,sortEigs=-1,eigvals=E,err=err,message=mess)
        if(err/=0) then;mess='PCA:'//trim(mess);return;endif
        r=count(E>=epsilon);keep=min(r,ND)
        if(allocated(pcaRes%eigval)) deallocate(pcaRes%eigval);allocate(pcaRes%eigval(r))
        pcaRes%eigval=E(1:r)
        if(allocated(pcaRes%eigvec)) deallocate(pcaRes%eigvec);allocate(pcaRes%eigvec(M,keep))
        pcaRes%eigvec=C(:,1:keep)
        if(allocated(pcaRes%proj)) deallocate(pcaRes%proj);allocate(pcaRes%proj(N,keep))
        pcaRes%proj=matmul(X,pcaRes%eigvec)
    case('SVD','svd')
        ! perform SVD
        A=transpose(X)/sqrt(n-1._mrk)
        call svd_dcmp(a=A,w=W,v=V,err=err,message=mess)
        if(err/=0) then;mess='PCA:'//trim(mess);return;endif
        ! Sort W & get rank
        call indexx_qsort(arr=W,indx=indx,ascnd=.false.,err=err,message=mess)
        if(err/=0) then;mess='PCA:'//trim(mess);return;endif
        r=count(W>=epsilon);keep=min(r,ND)
        ! re-order and truncate
        if(allocated(U1)) deallocate(U1);allocate(U1(M,keep))
        if(allocated(V1)) deallocate(V1);allocate(V1(N,keep))
        if(allocated(delta)) deallocate(delta);allocate(delta(keep,keep))
        if(allocated(Y)) deallocate(Y);allocate(Y(N,keep))
        U1=A(:,indx(1:keep));V1=V(:,indx(1:keep));
        delta=0._mrk
        do i=1,keep
            delta(i,i)=W(indx(i))
        enddo
        ! projection
        Y=sqrt(n-1._mrk)*MatMul(V1,delta)
        if(allocated(pcaRes%eigval)) deallocate(pcaRes%eigval);allocate(pcaRes%eigval(r))
        pcaRes%eigval=W(indx(1:r))**2
        if(allocated(pcaRes%eigvec)) deallocate(pcaRes%eigvec);allocate(pcaRes%eigvec(M,keep))
        pcaRes%eigvec=U1
        if(allocated(pcaRes%proj)) deallocate(pcaRes%proj);allocate(pcaRes%proj(N,keep))
        pcaRes%proj=Y
    case default
        err=4;mess='PCA:Fatal:unknown option';return
    end select
if(present(folder)) call PCA_WriteRes(folder,pcaRes,err,mess)
end subroutine PCA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine PCA_WriteRes(folder,pcaRes,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:09/11/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. folder
!^*		1.pcaRes, object of type pcaType containing the results
!^* OUT
!^*		1.err
!^*		2.mess
!^**********************************************************************
use utilities_dmsl_kit, only:getspareunit
character(*),intent(in)::folder
type(pcaType), intent(in)::pcaRes
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
character(250),parameter::eigvalFile='Eigval.txt',eigvecFile='Eigvec.txt',projFile='Proj.txt'
integer(mik)::unt,i,n

call getspareunit(unt,err,mess)
if(err/=0) then;mess='PCA_WriteRes:'//trim(mess);return;endif
open(unit=unt,file=trim(folder)//trim(eigvalFile),status='replace')
write(unt,'(3A14)') 'i','Eigval','Cumulative'
do i=1,size(pcaRes%eigval)
    write(unt,'(I14,2e14.6)') i,pcaRes%eigval(i),sum(pcaRes%eigval(1:i))/sum(pcaRes%eigval)
enddo
close(unt)
call getspareunit(unt,err,mess)
if(err/=0) then;mess='PCA_WriteRes:'//trim(mess);return;endif
open(unit=unt,file=trim(folder)//trim(eigvecFile),status='replace')
n=size(pcaRes%eigvec,2)
do i=1,size(pcaRes%eigvec,1)
    write(unt,'(<n>e14.6)') pcaRes%eigvec(i,:)
enddo
close(unt)
call getspareunit(unt,err,mess)
if(err/=0) then;mess='PCA_WriteRes:'//trim(mess);return;endif
open(unit=unt,file=trim(folder)//trim(projFile),status='replace')
n=size(pcaRes%proj,2)
do i=1,size(pcaRes%proj,1)
    write(unt,'(<n>e14.6)') pcaRes%proj(i,:)
enddo
close(unt)

end subroutine PCA_WriteRes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine KernelPCA(X,kernel,kernelpar,DoCenter,DoStd,nDim,folder,pcaRes,err,mess)
!^**********************************************************************
!^* Purpose: Kernel pca
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:09/11/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. X, data matrix
!^*		2. Kernel, e.g. 'Kernel_Gauss' or 'Kernel_Polynomial'
!^*		3. [Kernelpar], kernel parameters - see specific kernel subs for default values
!^*		4. DoCenter, perform centering?
!^*		5. DoStd, perform Standardazing?
!^*		6. [ndim], number of dimensions to keep; default = 4
!^*		7. [folder], result written 2 file if folder is provided
!^* OUT
!^*		1.pcaRes, object of type pcaType containing the results
!^*		2.err
!^*		3.mess
!^**********************************************************************
use numerix_dmsl_kit, only:getmeanvar
use linalg_dmsl_kit, only:getEigValsVecsRealSym

real(mrk), intent(inout)::X(:,:)
character(*), intent(in)::kernel
real(mrk), intent(in),optional::kernelpar(:)
character(*), intent(in),optional::folder
integer(mik),intent(in),optional::ndim
logical, intent(in), optional::DoCenter,DoStd
type(pcaType), intent(out)::pcaRes
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
real(mrk), parameter::epsilon=0.000001_mrk
logical, parameter::DoC_def=.false.,DoS_def=.false.
integer(mik), parameter::ND_def=4
logical::DoC,DoS
integer(mik)::N,M,i,r,keep,ND
real(mrk)::mu(size(X,dim=2)),sig(size(X,dim=2)),E(size(X,dim=1)),&
        G(size(X,dim=1),size(X,dim=1)),ColSum(size(X,dim=1)),TotSum,&
        ones(size(X,dim=1)),J(size(X,dim=1),size(X,dim=1))
!A(size(X,dim=2),size(X,dim=1)),&
!            W(size(X,dim=1)),V(size(X,dim=1),size(X,dim=1)),&
!            C(size(X,dim=2),size(X,dim=2)),E(size(X,dim=2))
real(mrk),allocatable::delta(:,:)!,V1(:,:),Y(:,:)!U1(:,:),

! init
err=0;mess=''
! Get Dimensions
N=size(X,dim=1);M=size(X,dim=2)
! Handle optional variables
if(present(DoCenter)) then;DoC=DoCenter;else;DoC=DoC_def;endif
if(present(DoStd)) then;DoS=DoStd;else;DoS=DoS_def;endif
if(present(Ndim)) then;ND=Ndim;else;ND=ND_def;endif
! Center/Standardize if requested
if(DoS .or. DoC) then
    call getmeanvar(x=X,mean=mu,var=sig,package=2,method='f',err=err,message=mess)
    if(err/=0) then;mess='PCA:'//trim(mess);return;endif
    if(DoS) then;sig=sqrt(sig);else;sig=1._mrk;endif
    if(.not.DoC) mu=0._mrk
    do i=1,M
        X(:,i)=(X(:,i)-mu(i))/sig(i)
    enddo
endif

! Apply Kernel
call ApplyKernel(kernel=kernel,par=kernelpar,X=X,Gramian=G,err=err,mess=mess)
if(err/=0) then;mess='KernelPCA:'//trim(mess);return;endif
! Normalisation of Gramian (ie centering into unknown transformed space)
ColSum=sum(G,dim=1)/real(N,mrk)
TotSum=sum(ColSum)/real(N,mrk)
ones=1._mrk
J=matmul(reshape(ones,(/N,1/)),reshape(ColSum,(/1,N/)))
G=G-J-transpose(J)+TotSum
! Apply eigenvalue decomposition to the Gramian
call getEigValsVecsRealSym(a2evex=G,doEigVex=.true.,emeth=2,sortEigs=-1,eigvals=E,err=err,message=mess)
if(err/=0) then;mess='KernelPCA:'//trim(mess);return;endif
r=count(E>=epsilon);keep=min(r,ND)
if(allocated(pcaRes%eigval)) deallocate(pcaRes%eigval);allocate(pcaRes%eigval(r))
pcaRes%eigval=E(1:r)
if(allocated(pcaRes%eigvec)) deallocate(pcaRes%eigvec);allocate(pcaRes%eigvec(N,keep))
pcaRes%eigvec=G(:,1:keep)
if(allocated(pcaRes%proj)) deallocate(pcaRes%proj);allocate(pcaRes%proj(N,keep))
if(allocated(delta)) deallocate(delta);allocate(delta(keep,keep))
delta=0._mrk
do i=1,keep
    delta(i,i)=sqrt(E(i))
enddo
! projection
pcaRes%proj=MatMul(pcaRes%eigvec,delta)
! Write if requested
if(present(folder)) call PCA_WriteRes(folder,pcaRes,err,mess)
end subroutine KernelPCA
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ApplyKernel(kernel,par,X,Gramian,err,mess)
!^**********************************************************************
!^* Purpose: Apply a Kernel to derive the Gramian matrix of X
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / Started @ Columbia University
!^**********************************************************************
!^* Last modified:12/11/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Kernel, e.g. 'Kernel_Gauss' or 'Kernel_Polynomial'
!^*		2. [par], kernel parameters
!^*		1. X, data matrix
!^* OUT
!^*		1.Gramian, the Gramian matrix
!^*		2.err
!^*		3.mess
!^**********************************************************************
character(*), intent(in)::kernel
real(mrk),intent(in), optional::par(:)
real(mrk),intent(in)::X(:,:)
real(mrk),intent(out)::Gramian(size(X,dim=1),size(X,dim=1))
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
real(mrk),parameter::parGauss_def=1._mrk
integer(mik)::N,i,j
real(mrk)::param

err=0;mess=''
N=size(X,dim=1)
select case(trim(kernel))
case('Kernel_Gauss')
    if(present(par)) then
        if(size(par)==1) then
            param=par(1)
        else
            param=parGauss_def
        endif
    else
        param=parGauss_def
    endif
    do i=1,N
        do j=1,i
            Gramian(i,j)=exp(-(dot_product((X(i,:)-X(j,:)),(X(i,:)-X(j,:))))/(2._mrk*param**2))
            Gramian(j,i)=Gramian(i,j)
        enddo
    enddo
case default
    err=1;mess='ApplyKernel:Fatal:unknown kernel';return
end select

end subroutine ApplyKernel 
 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine S_MECA(Y,Phi,Y_DoCenter,Y_DoStd,Phi_DoCenter,Phi_DoStd,PCAoption,nDim,folder,mecaRes,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:08/10/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Y, target data matrix
!^*		2. Phi, Explanatory data matrix
!^*		3. Y_DoCenter, perform Y centering?
!^*		4. Y_DoStd, perform Y Standardizing?
!^*		5. Phi_DoCenter, perform Phi centering?
!^*		6. Phi_DoStd, perform Phi Standardazing?
!^*		7. [option], computational variant for PCA: 'standard' or 'SVD', default 'standard'
!^*		   'standard' more efficient for "vertical" X (ie many rows, few columns)
!^*		   'SVD' more efficient for "horizontal" X (ie few rows, many columns)
!^*		8. [ndim], number of dimensions to keep; default = 4
!^*		9. [folder], result written 2 file if folder is provided
!^* OUT
!^*		1.mecaRes, object of type pcaType containing the results
!^*		2.err
!^*		3.mess
!^**********************************************************************
use numerix_dmsl_kit, only:getmeanvar
use EmpiricalStats_tools, only:GetPearsonCorr
real(mrk), intent(inout)::Y(:,:),Phi(:,:)
character(*), intent(in),optional::PCAoption,folder
integer(mik),intent(in),optional::ndim
logical, intent(in), optional::Y_DoCenter,Y_DoStd,Phi_DoCenter,Phi_DoStd
type(mecaType), intent(out)::mecaRes
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::Nt,Nx,Ns,i,s,M
type(pcaType)::YPCA
real(mrk)::mu(size(Phi,dim=2)),sig(size(Phi,dim=2)),Phihat(size(Phi,dim=2))
logical, parameter::DoC_def=.false.,DoS_def=.false.
logical::DoC,DoS

err=0;mess=''
Nt=size(Y,dim=1);Nx=size(Y,dim=2);Ns=size(Phi,dim=2)
if(size(Phi,dim=1)/=Nt) then
    err=1;mess='S_MECA:Fatal:Size mismatch [Y,Phi]';return
endif
! Y PCA
call PCA(Y,Y_DoCenter,Y_DoStd,PCAoption,ndim,folder,YPCA,err,mess)
if(err/=0) then;mess='S_MECA:Fatal:'//trim(mess);return;endif
if(allocated(mecaRes%Ycenter)) deallocate(mecaRes%Ycenter);allocate(mecaRes%Ycenter(Nx))
if(allocated(mecaRes%Yscale)) deallocate(mecaRes%Yscale);allocate(mecaRes%Yscale(Nx))
mecares%Ycenter=YPCA%Xcenter;mecares%Yscale=YPCA%Xscale
if(allocated(mecaRes%eigval)) deallocate(mecaRes%eigval);allocate(mecaRes%eigval(size(YPCA%eigval)))
mecaRes%eigval=YPCA%eigval
if(allocated(mecaRes%lambda)) deallocate(mecaRes%lambda);allocate(mecaRes%lambda(Nx,ndim))
mecaRes%lambda=YPCA%eigvec
if(allocated(mecaRes%tau)) deallocate(mecaRes%tau);allocate(mecaRes%tau(Nt,ndim))
mecaRes%tau=YPCA%proj
if(allocated(mecaRes%psi)) deallocate(mecaRes%psi);allocate(mecaRes%psi(Ns,ndim))
! Phi Constrained PCA
if(present(Phi_DoCenter)) then;DoC=Phi_DoCenter;else;DoC=DoC_def;endif
if(present(Phi_DoStd)) then;DoS=Phi_DoStd;else;DoS=DoS_def;endif
if(allocated(mecaRes%Phicenter)) deallocate(mecaRes%Phicenter);allocate(mecaRes%Phicenter(Ns))
if(allocated(mecaRes%Phiscale)) deallocate(mecaRes%Phiscale);allocate(mecaRes%Phiscale(Ns))
if(DoS .or. DoC) then
    call getmeanvar(x=Phi,mean=mu,var=sig,package=2,method='f',err=err,message=mess)
    if(err/=0) then;mess='S_MECA:'//trim(mess);return;endif
    if(DoS) then;sig=sqrt(sig);else;sig=1._mrk;endif
    if(.not.DoC) mu=0._mrk
    do i=1,size(Phi,dim=2)
        Phi(:,i)=(Phi(:,i)-mu(i))/sig(i)
    enddo
    mecares%Phicenter=mu;mecares%Phiscale=sig;
else
    mecares%Phicenter=0._mrk;mecares%Phiscale=1._mrk;
endif
do i=1,ndim
    do s=1,Ns
        mecaRes%psi(s,i)=dot_product(Phi(:,s),mecaRes%tau(:,i))/dot_product(mecaRes%tau(:,i),mecaRes%tau(:,i))
    enddo
enddo
! Reconstruct Tau
if(allocated(mecaRes%tauHat)) deallocate(mecaRes%tauHat);allocate(mecaRes%tauHat(Nt,ndim))
call S_MECA_ReconstructTau(Phi,mecaRes%psi,mecaRes%TauHat,err,mess)
if(err/=0) then;mess='S_MECA:Fatal:'//trim(mess);return;endif
! Ytilde & Yhat
if(allocated(mecaRes%Ytilde)) deallocate(mecaRes%Ytilde);allocate(mecaRes%Ytilde(Nt,Nx))
if(allocated(mecaRes%Yhat)) deallocate(mecaRes%Yhat);allocate(mecaRes%Yhat(Nt,Nx))
call S_MECA_ReconstructY(mecaRes%lambda,mecaRes%tau,mecaRes%Ytilde,err,mess)
if(err/=0) then;mess='S_MECA:Fatal:'//trim(mess);return;endif
call S_MECA_ReconstructY(mecaRes%lambda,mecaRes%tauHat,mecaRes%Yhat,err,mess)
if(err/=0) then;mess='S_MECA:Fatal:'//trim(mess);return;endif

end subroutine S_MECA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine S_MECA_ReconstructTau(Phi,Psi,TauHat,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:08/10/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Phi, Explanatory data matrix
!^*		2. Psi
!^* OUT
!^*		1.TauHat, reconstructed temporal patterns
!^*		2.err
!^*		3.mess
!^**********************************************************************
use linalg_dmsl_kit, only:choles_invrt
real(mrk), intent(in)::Phi(:,:), Psi(:,:)
real(mrk), intent(out)::TauHat(size(Phi,dim=1),size(Psi,dim=2))
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
real(mrk)::A(size(Psi,dim=2),size(Psi,dim=2)),B(size(Psi,dim=2))
integer(mik)::i,j,t,Nk,Nt
logical::PD

err=0;mess='';TauHat=UndefRN
Nk=size(Psi,dim=2);Nt=size(Phi,dim=1)
! Construct & invert A
do i=1,Nk
    do j=i,Nk
        A(i,j)=dot_product(Psi(:,i),Psi(:,j))
        A(j,i)=A(i,j)
    enddo
enddo
call choles_invrt(ainv=A,posDefinite=PD,err=err,message=mess)
if(err/=0) then;mess='S_MECA_ReconstructTau:'//trim(mess);return;endif
if(.not.PD) then;err=1;mess='S_MECA_ReconstructTau:Fatal:[A] non pos-def';endif
do t=1,Nt
    ! construct B
    do j=1,Nk
        B(j)=dot_product(Phi(t,:),Psi(:,j))
    enddo
    TauHat(t,:)=matmul(B,A)
enddo
end subroutine S_MECA_ReconstructTau
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine S_MECA_ReconstructY(lambda,tau,Y,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:08/10/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. lambda, Y-spatial components
!^*		2. tau, Y-temporal components
!^* OUT
!^*		1.Y, reconstructed data
!^*		2.err
!^*		3.mess
!^**********************************************************************
real(mrk), intent(in)::lambda(:,:), tau(:,:)
real(mrk), intent(out)::Y(size(tau,dim=1),size(lambda,dim=1))
integer(mik), intent(out)::err
character(*), intent(out)::mess
! locals
integer(mik)::Nt,Nx,t,x

err=0;mess='';Y=UndefRN
Nt=size(tau,dim=1);Nx=size(lambda,dim=1)

do t=1,Nt
    do x=1,Nx
        Y(t,x)=dot_product(lambda(x,:),tau(t,:))
    enddo
enddo 
end subroutine S_MECA_ReconstructY
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine S_MECA_LOOCV(Y,Phi,Y_DoCenter,Y_DoStd,Phi_DoCenter,Phi_DoStd,PCAoption,nDim,folder,LOOCVres,err,mess)
!^**********************************************************************
!^* Purpose: Leave-One-Out Cross Validation for S-MECA
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:20/12/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Y, target data matrix
!^*		2. Phi, Explanatory data matrix
!^*		3. Y_DoCenter, perform Y centering?
!^*		4. Y_DoStd, perform Y Standardizing?
!^*		5. Phi_DoCenter, perform Phi centering?
!^*		6. Phi_DoStd, perform Phi Standardazing?
!^*		7. [option], computational variant for PCA: 'standard' or 'SVD', default 'standard'
!^*		   'standard' more efficient for "vertical" X (ie many rows, few columns)
!^*		   'SVD' more efficient for "horizontal" X (ie few rows, many columns)
!^*		8. [ndim], number of dimensions to keep; default = 4
!^*		9. [folder], result written 2 file if folder is provided
!^* OUT
!^*		1.LOOCVres, object of type LOOCV type containing the results
!^*		2.err
!^*		3.mess
!^**********************************************************************
use numerix_dmsl_kit, only:getmeanvar
use EmpiricalStats_tools, only:GetPearsonCorr
use Utilities_dmsl_kit, only:terminateRowColMat
real(mrk), intent(inout)::Y(:,:),Phi(:,:)
character(*), intent(in),optional::PCAoption,folder
integer(mik),intent(in),optional::ndim
logical, intent(in), optional::Y_DoCenter,Y_DoStd,Phi_DoCenter,Phi_DoStd
type(LOOCVtype), intent(out)::LOOCVres
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::Nt,Nx,Ns,Nd,t,onetont(size(Y,dim=1)),onetonx(size(Y,dim=2)),onetons(size(Phi,dim=2))
real(mrk)::Y0(size(Y,dim=1)-1,size(Y,dim=2)),Phi0(size(Phi,dim=1)-1,size(Phi,dim=2))
real(mrk)::Ystar(1,size(Y,dim=2)),Phistar(1,size(Phi,dim=2))
real(mrk), allocatable:: TauHat(:,:)
logical::mask(size(Y,dim=1))
type(mecatype)::mecaRes

err=0;mess=''
Nt=size(Y,dim=1);Nx=size(Y,dim=2);Ns=size(Phi,dim=2)
if(size(Phi,dim=1)/=Nt) then
    err=1;mess='S_MECA_LOOCV:Fatal:Size mismatch [Y,Phi]';return
endif

if(allocated(LOOCVres%Y)) deallocate(LOOCVres%Y);allocate(LOOCVres%Y(Nt,Nx))
if(allocated(LOOCVres%Yhat)) deallocate(LOOCVres%Yhat);allocate(LOOCVres%Yhat(Nt,Nx))

onetont=(/(t,t=1,Nt)/);onetonx=(/(t,t=1,Nx)/);onetons=(/(t,t=1,Ns)/)
do t=1,Nt
    mask=(onetont/=t)
    call terminateRowColMat(m0=Y,m1=Y0,activeRow=mask,activeCol=onetonx>=0,err=err,message=mess)
    if(err/=0) then;mess='S_MECA_LOOCV:Fatal:'//trim(mess);return;endif
    call terminateRowColMat(m0=Phi,m1=Phi0,activeRow=mask,activeCol=onetons>=0,err=err,message=mess)
    if(err/=0) then;mess='S_MECA_LOOCV:Fatal:'//trim(mess);return;endif
    call S_MECA(Y0,Phi0,Y_DoCenter,Y_DoStd,Phi_DoCenter,Phi_DoStd,PCAoption,nDim,folder,mecaRes,err,mess)
    if(err/=0) then;mess='S_MECA_LOOCV:Fatal:'//trim(mess);return;endif
    Ystar(1,:)=(Y(t,:)-mecares%Ycenter)/mecares%Yscale
    Phistar(1,:)=(Phi(t,:)-mecares%Phicenter)/mecares%Phiscale
    ! Reconstruct Tau
    Nd=size(mecaRes%Tau,dim=2)
    if(allocated(TauHat)) deallocate(TauHat);allocate(TauHat(1,Nd))
    call S_MECA_ReconstructTau(Phistar,mecaRes%psi,TauHat,err,mess)
    if(err/=0) then;mess='S_MECA_LOOCV:Fatal:'//trim(mess);return;endif
    ! Yhat
    LOOCVres%Y(t,:)=Y(t,:)
    call S_MECA_ReconstructY(mecaRes%lambda,TauHat,LOOCVres%Yhat(t,:),err,mess)
    if(err/=0) then;mess='S_MECA_LOOCV:Fatal:'//trim(mess);return;endif
enddo

if(allocated(LOOCVres%tcor)) deallocate(LOOCVres%tcor);allocate(LOOCVres%tcor(Nx))
if(allocated(LOOCVres%scor)) deallocate(LOOCVres%scor);allocate(LOOCVres%scor(Nt))
do t=1,Nt
    call GetPearsonCorr(LOOCVres%Y(t,:),LOOCVres%Yhat(t,:),LOOCVres%scor(t),err,mess)
    if(err/=0) then;mess='S_MECA_LOOCV:Fatal:'//trim(mess);return;endif
enddo
do t=1,Nx
    call GetPearsonCorr(LOOCVres%Y(:,t),LOOCVres%Yhat(:,t),LOOCVres%tcor(t),err,mess)
    if(err/=0) then;mess='S_MECA_LOOCV:Fatal:'//trim(mess);return;endif
enddo


end subroutine S_MECA_LOOCV
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine S_MECA_ChooseNdim(Y,Phi,Y_DoCenter,Y_DoStd,Phi_DoCenter,Phi_DoStd,&
                           PCAoption,nDimMax,folder,Ndim,cormatrix,err,mess)
!^**********************************************************************
!^* Purpose: choose the number of components by LOOCV (maximise median correlation between Y and Yhat)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:21/12/2013
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Y, target data matrix
!^*		2. Phi, Explanatory data matrix
!^*		3. Y_DoCenter, perform Y centering?
!^*		4. Y_DoStd, perform Y Standardizing?
!^*		5. Phi_DoCenter, perform Phi centering?
!^*		6. Phi_DoStd, perform Phi Standardazing?
!^*		7. [option], computational variant for PCA: 'standard' or 'SVD', default 'standard'
!^*		   'standard' more efficient for "vertical" X (ie many rows, few columns)
!^*		   'SVD' more efficient for "horizontal" X (ie few rows, many columns)
!^*		8. [ndim], number of dimensions to keep; default = 4
!^*		9. [folder], result written 2 file if folder is provided
!^* OUT
!^*		1.LOOCVres, object of type LOOCV type containing the results
!^*		2.err
!^*		3.mess
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalQuantile
real(mrk), intent(inout)::Y(:,:),Phi(:,:)
character(*), intent(in),optional::PCAoption,folder
integer(mik),intent(in)::ndimmax
logical, intent(in), optional::Y_DoCenter,Y_DoStd,Phi_DoCenter,Phi_DoStd
real(mrk), intent(out)::cormatrix(size(Y,dim=2),NdimMax)
integer(mik), intent(out)::ndim,err
character(*), intent(out)::mess
!locals
integer(mik)::n,foo(1)
type(LOOCVtype)::LOOCVres
real(mrk)::cormed(NdimMax)

err=0;mess=''
do n=1,NdimMax
    call S_MECA_LOOCV(Y,Phi,Y_DoCenter,Y_DoStd,Phi_DoCenter,Phi_DoStd,PCAoption,n,folder,LOOCVres,err,mess)
    if(err/=0) then;mess='S_MECA_ChooseNdim:'//trim(mess);return;endif
    cormatrix(:,n)=LOOCVres%tcor
    call GetEmpiricalQuantile(p=0.5_mrk,x=LOOCVres%tcor,q=cormed(n),err=err,mess=mess)
    if(err/=0) then;mess='S_MECA_ChooseNdim:'//trim(mess);return;endif
enddo
foo=maxloc(cormed)
ndim=foo(1)

end subroutine S_MECA_ChooseNdim
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module DimReduction_tools
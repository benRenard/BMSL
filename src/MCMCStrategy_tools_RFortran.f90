module MCMCStrategy_tools

!~**********************************************************************
!~* Purpose: MCMC strategies
!~**********************************************************************
!~* Programmer: Ben Renard,Xun SUN, Irstea
!~**********************************************************************
!~* Last modified: 09/02/2012
!~**********************************************************************
!~* Comments: 
!~*        - For the optional component
!~*        fAux: Auxiliary information
!~*        UseRfortran: show figures in the R or not?
!~*        showPar: vector of index of parameter, that will be show in R
!~*        ranges: ylim associated to the above parameters
!~*        OutMat: pointer of real, same information as in the OutFile
!~*        BurnFactor: In the mixted algorithm, the covariance matrix for
!~*             the Gaussian Jump is computed based on the first run's 
!~*             parameters with burning first "BurnFactor" of these parameters.
!~*        scaleGJ: scale for the Gaussian Jump, default value is 2.4
!~*        GRIndex: pointer of real, GRIndex values
!~*        GRBlocSize: A predefined bloc size to show the Guassian Jump,
!~*             default value is equal to nAdapt2
!~*        GRBurnFactor: percentage of the total length of parameters used
!~*             to compute de Gelman-Rubin criterion
!~*        - In the multi chain case, some inputs and outputs may have
!~*             one more dimension for the number of chain 
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1.Adaptive_Metro: Metropolis–Hastings algorithm with Adaptive
!~*		2.Adaptive_Metro_OAAT: Gibbs sampling algorithm with Adaptive
!~*		3.Metro_Mix: 3 steps algorithm:
!~*		        - run with algrothms 1.
!~*		        - Use the above results (After burning) to find
!~*		          the covariance matrix, then run with algorithm 2.
!~*		        - Use the same covariance matrix and the scale found above
!~*		          to run with Gaussian jumps
!~*		4.AdaptiveMS_Metro: Same algorithm as 1. Multi-start points case.
!~*		5.AdaptiveMS_Metro_OAAT: Same algorithm as 2. Multi-start points case.
!~*		6.MetroMS_Mix: Same algorithm as 3. Multi-start points case.
!~*		7.MetroMS_Mix2: 3 steps algorithm:
!~*		        - Start with a single points using algorithm 1.
!~*		        - Use the above results (After burning) to find
!~*		          the covariance matrix, then run with algorithm 2. 
!~*		        - After the first 2 steps, n re-starting points are randomly
!^*               selected from the previous results. Then run multi-chain
!^*               with these re-starting points with Gaussian jumps.
!~*		8.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use MCMC_tools

implicit none
Private
public :: Adaptive_Metro_OAAT, Adaptive_Metro, Metro_Mix,&
           AdaptiveMS_Metro_OAAT, AdaptiveMS_Metro, MetroMS_Mix,&
            MetroMS_Mix2,Adaptive_Metro_OAAT_v2

logical::UseRF

Contains


subroutine Adaptive_Metro(f,x,fx,fAux,covar,scale,nAdapt,nCycles,&
				MinMoveRate,MaxMoveRate,DownMult,UpMult,&
				UseRFortran,&
				showPar,ranges,&
				OutFile,OutMat,headers,err,mess)

!^**********************************************************************
!^* Purpose: Metropolis with adaptive scale factor
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified:09/11/2008
!^**********************************************************************
!^* Comments: Only the scale factor is adapted - var is choleskized
!^* out of this sub
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
use utilities_dmsl_kit, only:number_string
Use RFortran

real(mrk),intent(inout)::x(:),covar(:,:),scale
integer(mik),intent(in)::nAdapt,nCycles
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in),optional::OutFile,headers(:)
logical,intent(in),optional::UseRFortran
integer(mik),intent(in),optional::showpar(:)
real(mrk),intent(in),optional::ranges(:,:)
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
real(mrk),pointer,optional::outMat(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
logical::feas,isnull,again,ok
character(100), allocatable::head(:), fAuxHead(:)
logical::move
real(mrk)::moverate
integer(mik)::n,i,j,p,k,compt
real(mrk),allocatable::parlist(:,:)
real(mrk)::postlist(nAdapt),postlist2(nCycles)

interface
    subroutine f(x,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface

if (present(UseRFortran)) then
    UseRF=UseRFortran
else
    UseRF=.true.
endif

!Handle Output file headers
n=size(x,dim=1)
if(present(OutFile)) then
	open(unit=1,file=OutFile,status='replace')
	if(present(headers)) then ! trust user - no size check at the moment
		k=size(headers)
		write(1,'(<k>A14)') headers
	else
		if(allocated(head)) deallocate(head)
		allocate(head(n))
		do i=1,n
			head(i)='par'//trim(number_string(i))
		enddo
		if(present(fAux)) then
			p=size(fAux,dim=1)
			if(allocated(fAuxHead)) deallocate(fAuxHead)
			allocate(fAuxHead(p))
			do i=1,p
				fAuxhead(i)=' fAux'//trim(number_string(i))
			enddo
			write(1,'(<n+1+p>A14)') head(:),'Objectivefunction',fAuxHead(:)
			deallocate(fAuxHead)
		else
			write(1,'(<n+1>A14)') head(:),'Objectivefunction'
		endif
		deallocate(head)
	endif
endif

! Check starting value
call f(x=x,feas=feas,isnull=isnull,fx=fx,fAux=fAux,err=err,mess=mess)
If(err>0) then
	mess='Adaptive_Metro'//trim(mess)
	if(present(OutFile)) close(1)
    return
endif
if((.not.feas).or.isnull) then
	err=1;mess='Adaptive_Metro:Fatal:Unfeasible starting point'
	if(present(OutFile)) close(1)
    return
endif
if(present(OutFile)) then
	if(present(fAux)) then
		write (1,'(<n+1+p>e14.6)') x,fx,fAux
	else
		write (1,'(<n+1>e14.6)') x,fx
	endif
endif

! Associate outMat
if (present(outMat)) then
    if (present(fAux)) then
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1+size(fAux,dim=1)))
        outMat(1,1:n)=x(:)
        outMat(1,n+1)=fx
        outMat(1,(n+2):(n+1+size(fAux,dim=1)))=fAux
    else
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1))
        outMat(1,1:n)=x(:)
        outMat(1,n+1)=fx
    endif
endif

		    
!allocate parlist
    if (present(showpar)) then
        !allocate parlist
            if(allocated(parlist)) deallocate(parlist)
            allocate(parlist(nCycles,size(showpar)))
        !check ranges
        if (present(ranges)) then
	        if (size(ranges,dim=2)/=2_mik) then
	            err=1;mess="AdaptiveMS_Metro_OAAT:ranges should include 2 real numbers"
	            if(present(OutFile)) close(1)
                return
	        endif
	        if (size(ranges,dim=1)/=size(showpar)) then
	            err=1;mess="AdaptiveMS_Metro_OAAT:ranges size mismatch"
	            if(present(OutFile)) close(1)
                return
	        endif
	    endif
    endif


if (UseRF) then
    ok = Rinit();ok=Reval('X11()')

    again=.true.;compt=0;ok=Reval("mrl<-c()")
    do while(again)
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
		    call Metro_GaussJump(f,x,fx,fAux,(scale**2)*covar,move,err,mess)
		    If(err>0) then
			    mess='Adaptive_Metro'//trim(mess)
			    if(present(OutFile)) close(1)
                return
		    endif
		    postlist(i)=fx
		    !Write2File
		    if(present(OutFile)) then
			    if(present(fAux)) then
				    write (1,'(<n+1+p>e14.6)') x,fx,fAux
			    else
				    write (1,'(<n+1>e14.6)') x,fx
			    endif
		    endif
		    if(present(outMat)) then
		        outMat(compt*nAdapt+i+1,1:n)=x
		        outMat(compt*nAdapt+i+1,n+1)=fx
		        if (present(fAux)) then
		            outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)))=fAux
		        endif
		    endif
		    !Update Move Rates
		    if(move) moverate=moverate+1._mrk/nAdapt
	    enddo
	    postlist2(compt+1)=fx
    	
	    if (present(showpar)) then
		    do j=1,size(showpar)
		        parlist(compt+1,j)=x(showpar(j))
		    enddo
	    endif
    	
	    if (present(showpar)) then
            ok=Rput('k',size(showpar))
            ok=reval('layout(matrix(c(rep(1,k),rep(2,k),rep(3,k),rep(4:(k+3),rep(3,k))),k*3,2))')
        else
            ok=Reval('layout(seq(1,3))')
        endif
    !    ok=Reval('layout(c(1, 2))')
        ok=Rput('mr', moverate)
	    ok=Rput('compt', compt+1)
	    ok=Reval('mrl[compt]<-mr')
        ok=Reval('matplot(mrl,type="l", xlab="Iteration(nCycles)", ylab="Move Rate", main="Move Rates (From Begin)", ylim=c(0, 1),cex.lab=2,cex.axis=2,cex.main=3)') 
        ok=Rput('post', postlist)
        ok=Reval('matplot(post,type="l", xlab="Iteration(nAdapt)", ylab="Log-posterior", main="Log-posterior (Bloc)",cex.lab=2,cex.axis=2,cex.main=3)') 
        ok=Rput('post2', postlist2) 
        ok=Rput('compt',compt)
        ok=Reval('matplot(post2[1:(1+compt)],type="l", xlab="Iteration(nCycles)", ylab="Log-posterior", main="Log-posterior (From Begin)",cex.lab=2,cex.axis=2,cex.main=3)')   
    	
        if (present(showpar)) then
            ok=Rput('pars', parlist)
            ok=Rput('compt',compt)
            ok=Rput('showpar',showpar)
            do k=1,size(showpar)
                if (present(ranges)) then
                    ok=Rput('rmin',ranges(k,1))
                    ok=Rput('rmax',ranges(k,2))
                    ok=Rput('k',k)
                    ok=Reval('matplot(pars[1:(compt+1),k],type="l", xlab="Iteration(nCycles)", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                else
                    ok=Rput('k',k)
    !                ok=Reval('rmin=quantile(pars[1:(compt+1),k],0.02)')
    !                ok=Reval('rmax=quantile(pars[1:(compt+1),k],0.98)')
                    ok=Reval('matplot(pars[1:(compt+1),k],type="l", xlab="Iteration(nCycles)", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",cex.lab=2,cex.axis=2,cex.main=3)')
                endif
            enddo
        endif
    	
	    ! Adapt Jump standard deviation
	    if(moverate<=MinMoveRate) scale=scale*DownMult
	    if(moverate>=MaxMoveRate) scale=scale*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if (present(OutFile)) close(1)
    ok = Rclose()  ! End R session

    if (present(showPar)) deallocate(parlist)

else

    !ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0
    do while(again)
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
		    call Metro_GaussJump(f,x,fx,fAux,(scale**2)*covar,move,err,mess)
		    If(err>0) then
			    mess='Adaptive_Metro'//trim(mess)
			    if(present(OutFile)) close(1)
                return
		    endif
		    postlist(i)=fx
		    !Write2File
		    if(present(OutFile)) then
			    if(present(fAux)) then
				    write (1,'(<n+1+p>e14.6)') x,fx,fAux
			    else
				    write (1,'(<n+1>e14.6)') x,fx
			    endif
		    endif
		    if(present(outMat)) then
		        outMat(compt*nAdapt+i+1,1:n)=x
		        outMat(compt*nAdapt+i+1,n+1)=fx
		        if (present(fAux)) then
		            outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)))=fAux
		        endif
		    endif
		    !Update Move Rates
		    if(move) moverate=moverate+1._mrk/nAdapt
	    enddo
    !    ok=Reval('layout(c(1, 2))')
    !    ok=Rput('mr', moverate)
    !	 ok=Rput('compt', compt+1)
    !	 ok=Reval('mrl[compt]<-mr')
    !    ok=Reval('matplot(mrl,type="l", xlab="Parameter", ylab="Move Rate", main="Move Rates", ylim=c(0, 1))')
    !    ok=Rput('post', postlist)
    !    ok=Reval('matplot(post,type="l", xlab="Iteration", ylab="Log-posterior", main="Log-posterior")')

	    ! Adapt Jump standard deviation
	    if(moverate<=MinMoveRate) scale=scale*DownMult
	    if(moverate>=MaxMoveRate) scale=scale*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if(present(OutFile)) close(1)
    !ok = Rclose()  ! End R session
endif

end subroutine Adaptive_Metro

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Adaptive_Metro_OAAT(f,x,fx,fAux,std,nAdapt,nCycles,&
				MinMoveRate,MaxMoveRate,DownMult,UpMult,&
				UseRFortran,&
				showPar,ranges,&
				OutFile,outMat,headers,err,mess)

!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard, University of Newcastle
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
use utilities_dmsl_kit, only:number_string
Use RFortran

real(mrk),intent(inout)::x(:),std(:)
integer(mik),intent(in)::nAdapt,nCycles
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in),optional::OutFile,headers(:)
logical,intent(in),optional::UseRFortran
integer(mik),intent(in),optional::showpar(:)
real(mrk),intent(in),optional::ranges(:,:)
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
real(mrk),pointer,optional::outMat(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
logical::feas,isnull,again,ok
character(100), allocatable::head(:), fAuxHead(:)
logical,allocatable::move(:)
real(mrk),allocatable::moverate(:),parlist(:,:)
integer(mik)::n,i,j,p,k,compt
real(mrk)::postlist(nAdapt),postlist2(nCycles)

interface
    subroutine f(x,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface


if (present(UseRFortran)) then
    UseRF=UseRFortran
else
    UseRF=.true.
endif


!Handle Output file headers
n=size(x,dim=1)
if(present(OutFile)) then
	open(unit=1,file=OutFile,status='replace')
	if(present(headers)) then ! trust user - no size check at the moment
		k=size(headers)
		write(1,'(<k>A14)') headers
	else
		if(allocated(head)) deallocate(head)
		allocate(head(n))
		do i=1,n
			head(i)='par'//trim(number_string(i))
		enddo
		if(present(fAux)) then
			p=size(fAux,dim=1)
			if(allocated(fAuxHead)) deallocate(fAuxHead)
			allocate(fAuxHead(p))
			do i=1,p
				fAuxhead(i)=' fAux'//trim(number_string(i))
			enddo
			write(1,'(<n+1+p>A14)') head(:),'Objectivefunction',fAuxHead(:)
			deallocate(fAuxHead)
		else
			write(1,'(<n+1>A14)') head(:),'Objectivefunction'
		endif
		deallocate(head)
	endif
endif

! Check starting value
call f(x=x,feas=feas,isnull=isnull,fx=fx,fAux=fAux,err=err,mess=mess)
If(err>0) then
	mess='Adaptive_Metro_OAAT'//trim(mess)
	if(present(OutFile)) close(1)
    return
endif
if((.not.feas).or.isnull) then
	err=1;mess='Adaptive_Metro_OAAT:Fatal:Unfeasible starting point'
	if(present(OutFile)) close(1)
    return
endif
if(present(OutFile)) then
	if(present(fAux)) then
		write (1,'(<n+1+p>e14.6)') x,fx,fAux
	else
		write (1,'(<n+1>e14.6)') x,fx
	endif
endif

! Associate outMat
if (present(outMat)) then
    if (present(fAux)) then
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1+size(fAux,dim=1)))
        outMat(1,1:n)=x(:)
        outMat(1,n+1)=fx
        outMat(1,(n+2):(n+1+size(fAux,dim=1)))=fAux
    else
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1))
        outMat(1,1:n)=x(:)
        outMat(1,n+1)=fx
    endif
endif

! Allocate MoveRate stuff
if(allocated(move)) deallocate(move)
allocate(move(n))
if(allocated(moverate)) deallocate(moverate)
allocate(moverate(n))

!allocate parlist
if (present(showpar)) then
    !allocate parlist
        if(allocated(parlist)) deallocate(parlist)
        allocate(parlist(nCycles,size(showpar)))
    !check ranges
    if (present(ranges)) then
        if (size(ranges,dim=2)/=2_mik) then
            err=1;mess="AdaptiveMS_Metro_OAAT;ranges should include 2 real numbers"
            if(present(OutFile)) close(1)
            return
        endif
        if (size(ranges,dim=1)/=size(showpar)) then
            err=1;mess="AdaptiveMS_Metro_OAAT;ranges size mismatch"
            if(present(OutFile)) close(1)
            return
        endif
    endif
endif
    
    
if (UseRF) then
    ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0
    do while(again)
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
		    call Metro_OAAT_GaussJump(f,x,fx,fAux,std,move,err,mess)
		    If(err>0) then
			    mess='Adaptive_Metro_OAAT'//trim(mess)
			    if(present(OutFile)) close(1)
                return
		    endif
		    postlist(i)=fx
		    !Write2File
		    if(present(OutFile)) then
			    if(present(fAux)) then
				    write (1,'(<n+1+p>e14.6)') x,fx,fAux
			    else
				    write (1,'(<n+1>e14.6)') x,fx
			    endif
		    endif
		    if(present(outMat)) then
		        outMat(compt*nAdapt+i+1,1:n)=x
		        outMat(compt*nAdapt+i+1,n+1)=fx
		        if (present(fAux)) then
		            outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)))=fAux
		        endif
		    endif
		    !Update Move Rates
		    where(move) moverate=moverate+1._mrk/nAdapt
	    enddo
	    postlist2(compt+1)=fx
    	
	    if (present(showpar)) then
		    do j=1,size(showpar)
		        parlist(compt+1,j)=x(showpar(j))
		    enddo
	    endif
    	
	    if (present(showpar)) then
            ok=Rput('k',size(showpar))
            ok=reval('layout(matrix(c(rep(1,k),rep(2,k),rep(3,k),rep(4:(k+3),rep(3,k))),k*3,2))')
        else
            ok=Reval('layout(seq(1,3))')
        endif
        
    !    ok=Reval('layout(c(1, 2))')
        ok=Rput('mr', moverate)
        ok=Reval('matplot(mr,type="l", xlab="Parameter", ylab="Move Rate", main="Move Rates", ylim=c(0, 1),cex.lab=2,cex.axis=2,cex.main=3)') 
        ok=Rput('post', postlist)
        ok=Reval('matplot(post,type="l", xlab="Iteration(nAdapt)", ylab="Log-posterior", main="Log-posterior (Bloc)",cex.lab=2,cex.axis=2,cex.main=3)') 
        ok=Rput('post2', postlist2) 
        ok=Rput('compt',compt)
        ok=Reval('matplot(post2[1:(1+compt)],type="l", xlab="Iteration(nCycles)", ylab="Log-posterior", main="Log-posterior (From Begin)",cex.lab=2,cex.axis=2,cex.main=3)')   
    	
        if (present(showpar)) then
            ok=Rput('pars', parlist)
            ok=Rput('compt',compt)
            ok=Rput('showpar',showpar)
            do k=1,size(showpar)
                if (present(ranges)) then
                    ok=Rput('rmin',ranges(k,1))
                    ok=Rput('rmax',ranges(k,2))
                    ok=Rput('k',k)
                    ok=Reval('matplot(pars[1:(compt+1),k],type="l", xlab="Iteration(nCycles)", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                else
                    ok=Rput('k',k)
    !                ok=Reval('rmin=quantile(pars[1:(compt+1),k],0.02)')
    !                ok=Reval('rmax=quantile(pars[1:(compt+1),k],0.98)')
                    ok=Reval('matplot(pars[1:(compt+1),k],type="l", xlab="Iteration(nCycles)", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",cex.lab=2,cex.axis=2,cex.main=3)')
                endif
            enddo
        endif

	    ! Adapt Jump standard deviation
	    where(moverate<=MinMoveRate) std=std*DownMult
	    where(moverate>=MaxMoveRate) std=std*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if(present(OutFile)) close(1)
    ok = Rclose()  ! End R session

else
    !ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0;!ok=Reval("mrl<-c()")
    do while(again)
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
		    call Metro_OAAT_GaussJump(f,x,fx,fAux,std,move,err,mess)
		    If(err>0) then
			    mess='Adaptive_Metro_OAAT'//trim(mess)
			    if(present(OutFile)) close(1)
                return
		    endif
		    postlist(i)=fx
		    !Write2File
		    if(present(OutFile)) then
			    if(present(fAux)) then
				    write (1,'(<n+1+p>e14.6)') x,fx,fAux
			    else
				    write (1,'(<n+1>e14.6)') x,fx
			    endif
		    endif
		    if(present(outMat)) then
		        outMat(compt*nAdapt+i+1,1:n)=x
		        outMat(compt*nAdapt+i+1,n+1)=fx
		        if (present(fAux)) then
		            outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)))=fAux
		        endif
		    endif
		    !Update Move Rates
		    where(move) moverate=moverate+1._mrk/nAdapt
	    enddo
    !    ok=Reval('layout(c(1, 2))')
    !    ok=Rput('mr', moverate)
    !    ok=Reval('matplot(mr,type="l", xlab="Parameter", ylab="Move Rate", main="Move Rates", ylim=c(0, 1))')
    !    ok=Rput('post', postlist)
    !    ok=Reval('matplot(post,type="l", xlab="Iteration", ylab="Log-posterior", main="Log-posterior")')

	    ! Adapt Jump standard deviation
	    where(moverate<=MinMoveRate) std=std*DownMult
	    where(moverate>=MaxMoveRate) std=std*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if(present(OutFile)) close(1)
    !ok = Rclose()  ! End R session
endif

deallocate(move)
deallocate(moverate)
if (present(showPar)) deallocate(parlist)

end subroutine Adaptive_Metro_OAAT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine Adaptive_Metro_OAAT_v2(f,x,fx,fAux,std,nAdapt,nCycles,&
				MinMoveRate,MaxMoveRate,DownMult,UpMult,&
				UseRFortran,&
				showPar,ranges,&
				OutFile,outMat,headers,err,mess)

!^**********************************************************************
!^* Purpose: exactly as Adaptive_Metro_OAAT, but with the second version 
!^*          of the interface to f 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon, coded @ Columbia Uni
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
use utilities_dmsl_kit, only:number_string
Use RFortran

real(mrk),intent(inout)::x(:),std(:)
integer(mik),intent(in)::nAdapt,nCycles
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in),optional::OutFile,headers(:)
logical,intent(in),optional::UseRFortran
integer(mik),intent(in),optional::showpar(:)
real(mrk),intent(in),optional::ranges(:,:)
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
real(mrk),pointer,optional::outMat(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
logical::feas,isnull,again,ok
character(100), allocatable::head(:), fAuxHead(:)
logical,allocatable::move(:)
real(mrk),allocatable::moverate(:),parlist(:,:)
integer(mik)::n,i,j,p,k,compt
real(mrk)::postlist(nAdapt),postlist2(nCycles)

interface
    subroutine f(x,k,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        integer(mik),intent(in)::k
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface


if (present(UseRFortran)) then
    UseRF=UseRFortran
else
    UseRF=.true.
endif


!Handle Output file headers
n=size(x,dim=1)
if(present(OutFile)) then
	open(unit=1,file=OutFile,status='replace')
	if(present(headers)) then ! trust user - no size check at the moment
		k=size(headers)
		write(1,'(<k>A14)') headers
	else
		if(allocated(head)) deallocate(head)
		allocate(head(n))
		do i=1,n
			head(i)='par'//trim(number_string(i))
		enddo
		if(present(fAux)) then
			p=size(fAux,dim=1)
			if(allocated(fAuxHead)) deallocate(fAuxHead)
			allocate(fAuxHead(p))
			do i=1,p
				fAuxhead(i)=' fAux'//trim(number_string(i))
			enddo
			write(1,'(<n+1+p>A14)') head(:),'Objectivefunction',fAuxHead(:)
			deallocate(fAuxHead)
		else
			write(1,'(<n+1>A14)') head(:),'Objectivefunction'
		endif
		deallocate(head)
	endif
endif

! Check starting value
call f(x=x,k=0,feas=feas,isnull=isnull,fx=fx,fAux=fAux,err=err,mess=mess)
If(err>0) then
	mess='Adaptive_Metro_OAAT_v2'//trim(mess)
	if(present(OutFile)) close(1)
    return
endif
if((.not.feas).or.isnull) then
	err=1;mess='Adaptive_Metro_OAAT_v2:Fatal:Unfeasible starting point'
	if(present(OutFile)) close(1)
    return
endif
if(present(OutFile)) then
	if(present(fAux)) then
		write (1,'(<n+1+p>e14.6)') x,fx,fAux
	else
		write (1,'(<n+1>e14.6)') x,fx
	endif
endif

! Associate outMat
if (present(outMat)) then
    if (present(fAux)) then
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1+size(fAux,dim=1)))
        outMat(1,1:n)=x(:)
        outMat(1,n+1)=fx
        outMat(1,(n+2):(n+1+size(fAux,dim=1)))=fAux
    else
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1))
        outMat(1,1:n)=x(:)
        outMat(1,n+1)=fx
    endif
endif

! Allocate MoveRate stuff
if(allocated(move)) deallocate(move)
allocate(move(n))
if(allocated(moverate)) deallocate(moverate)
allocate(moverate(n))

!allocate parlist
if (present(showpar)) then
    !allocate parlist
        if(allocated(parlist)) deallocate(parlist)
        allocate(parlist(nCycles,size(showpar)))
    !check ranges
    if (present(ranges)) then
        if (size(ranges,dim=2)/=2_mik) then
            err=1;mess="Adaptive_Metro_OAAT_v2;ranges should include 2 real numbers"
            if(present(OutFile)) close(1)
            return
        endif
        if (size(ranges,dim=1)/=size(showpar)) then
            err=1;mess="Adaptive_Metro_OAAT_v2;ranges size mismatch"
            if(present(OutFile)) close(1)
            return
        endif
    endif
endif
    
    
if (UseRF) then
    ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0
    do while(again)
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
		    call Metro_OAAT_GaussJump_v2(f,x,fx,fAux,std,move,err,mess)
		    If(err>0) then
			    mess='Adaptive_Metro_OAAT_v2'//trim(mess)
			    if(present(OutFile)) close(1)
                return
		    endif
		    postlist(i)=fx
		    !Write2File
		    if(present(OutFile)) then
			    if(present(fAux)) then
				    write (1,'(<n+1+p>e14.6)') x,fx,fAux
			    else
				    write (1,'(<n+1>e14.6)') x,fx
			    endif
		    endif
		    if(present(outMat)) then
		        outMat(compt*nAdapt+i+1,1:n)=x
		        outMat(compt*nAdapt+i+1,n+1)=fx
		        if (present(fAux)) then
		            outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)))=fAux
		        endif
		    endif
		    !Update Move Rates
		    where(move) moverate=moverate+1._mrk/nAdapt
	    enddo
	    postlist2(compt+1)=fx
    	
	    if (present(showpar)) then
		    do j=1,size(showpar)
		        parlist(compt+1,j)=x(showpar(j))
		    enddo
	    endif
    	
	    if (present(showpar)) then
            ok=Rput('k',size(showpar))
            ok=reval('layout(matrix(c(rep(1,k),rep(2,k),rep(3,k),rep(4:(k+3),rep(3,k))),k*3,2))')
        else
            ok=Reval('layout(seq(1,3))')
        endif
        
    !    ok=Reval('layout(c(1, 2))')
        ok=Rput('mr', moverate)
        ok=Reval('matplot(mr,type="l", xlab="Parameter", ylab="Move Rate", main="Move Rates", ylim=c(0, 1),cex.lab=2,cex.axis=2,cex.main=3)') 
        ok=Rput('post', postlist)
        ok=Reval('matplot(post,type="l", xlab="Iteration(nAdapt)", ylab="Log-posterior", main="Log-posterior (Bloc)",cex.lab=2,cex.axis=2,cex.main=3)') 
        ok=Rput('post2', postlist2) 
        ok=Rput('compt',compt)
        ok=Reval('matplot(post2[1:(1+compt)],type="l", xlab="Iteration(nCycles)", ylab="Log-posterior", main="Log-posterior (From Begin)",cex.lab=2,cex.axis=2,cex.main=3)')   
    	
        if (present(showpar)) then
            ok=Rput('pars', parlist)
            ok=Rput('compt',compt)
            ok=Rput('showpar',showpar)
            do k=1,size(showpar)
                if (present(ranges)) then
                    ok=Rput('rmin',ranges(k,1))
                    ok=Rput('rmax',ranges(k,2))
                    ok=Rput('k',k)
                    ok=Reval('matplot(pars[1:(compt+1),k],type="l", xlab="Iteration(nCycles)", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                else
                    ok=Rput('k',k)
    !                ok=Reval('rmin=quantile(pars[1:(compt+1),k],0.02)')
    !                ok=Reval('rmax=quantile(pars[1:(compt+1),k],0.98)')
                    ok=Reval('matplot(pars[1:(compt+1),k],type="l", xlab="Iteration(nCycles)", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",cex.lab=2,cex.axis=2,cex.main=3)')
                endif
            enddo
        endif

	    ! Adapt Jump standard deviation
	    where(moverate<=MinMoveRate) std=std*DownMult
	    where(moverate>=MaxMoveRate) std=std*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if(present(OutFile)) close(1)
    ok = Rclose()  ! End R session

else
    !ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0;!ok=Reval("mrl<-c()")
    do while(again)
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
		    call Metro_OAAT_GaussJump_v2(f,x,fx,fAux,std,move,err,mess)
		    If(err>0) then
			    mess='Adaptive_Metro_OAAT_v2'//trim(mess)
			    if(present(OutFile)) close(1)
                return
		    endif
		    postlist(i)=fx
		    !Write2File
		    if(present(OutFile)) then
			    if(present(fAux)) then
				    write (1,'(<n+1+p>e14.6)') x,fx,fAux
			    else
				    write (1,'(<n+1>e14.6)') x,fx
			    endif
		    endif
		    if(present(outMat)) then
		        outMat(compt*nAdapt+i+1,1:n)=x
		        outMat(compt*nAdapt+i+1,n+1)=fx
		        if (present(fAux)) then
		            outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)))=fAux
		        endif
		    endif
		    !Update Move Rates
		    where(move) moverate=moverate+1._mrk/nAdapt
	    enddo
    !    ok=Reval('layout(c(1, 2))')
    !    ok=Rput('mr', moverate)
    !    ok=Reval('matplot(mr,type="l", xlab="Parameter", ylab="Move Rate", main="Move Rates", ylim=c(0, 1))')
    !    ok=Rput('post', postlist)
    !    ok=Reval('matplot(post,type="l", xlab="Iteration", ylab="Log-posterior", main="Log-posterior")')

	    ! Adapt Jump standard deviation
	    where(moverate<=MinMoveRate) std=std*DownMult
	    where(moverate>=MaxMoveRate) std=std*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if(present(OutFile)) close(1)
    !ok = Rclose()  ! End R session
endif

deallocate(move)
deallocate(moverate)
if (present(showPar)) deallocate(parlist)

end subroutine Adaptive_Metro_OAAT_v2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Metro_Mix(f,x,fx,fAux,std,&
                nAdapt1,nCycles1,BurnFactor,&
                nAdapt2,nCycles2,&
                nSim,scaleGJ,GRBlocSize,&
				MinMoveRate,MaxMoveRate,DownMult,UpMult,&
				UseRFortran,&
				showPar,ranges,&
				OutFile1,OutFile2,OutFile3,&
				outMat1, outMat2, outMat3, &
				headers,err,mess)
!^**********************************************************************
!^* Purpose: 
!^**********************************************************************
!^* Programmer: Ben Renard & Xun SUN, Irstea
!^**********************************************************************
!^* Last modified:26/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* OUT
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* INOUT
!^*     1.
!^*     2.
!^**********************************************************************
use numerix_dmsl_kit, only:normaldev,getvar
use linalg_dmsl_kit, only:choles_dcmp
use utilities_dmsl_kit, only:number_string
use RFortran

real(mrk),intent(inout)::x(:),std(:)
real(mrk),intent(inout),optional::scaleGJ
integer(mik),intent(in),optional::GRBlocSize
integer(mik),intent(in)::nAdapt1,nCycles1,nAdapt2,nCycles2,nSim
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor
character(*), intent(in),optional::Outfile1,Outfile2,Outfile3,headers(:)
real(mrk),pointer, optional::outMat1(:,:),outMat2(:,:),outMat3(:,:)
logical,intent(in),optional::UseRFortran
integer(mik),intent(in),optional::showpar(:)
real(mrk),intent(in),optional::ranges(:,:)
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
real(mrk),allocatable::parlist(:,:),postlist2(:)
real(mrk),pointer::MCMCsamples(:,:)=>NULL()
integer(mik)::n,nFaux,nAdapt3,nCycles3
character(100), allocatable::head(:), fAuxHead(:)
real(mrk)::covar(size(x),size(x)),scale
integer(mik)::i,k,p,compt,comptSim,ok
logical::move,pos,again

interface
    subroutine f(x,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface

if (present(UseRFortran)) then
    UseRF=UseRFortran
else
    UseRF=.true.
endif

if (present(GRBlocSize)) then
    nAdapt3=GRBlocSize
else
    nAdapt3=nAdapt2
endif

nCycles3=Int(nSim/nAdapt3,mik)
    
err=0;mess=''
n=size(x)

if(allocated(postlist2)) deallocate(postlist2)
allocate(postlist2(nCycles1+nCycles2+nCycles3))
!allocate parlist
if (present(showpar)) then
    !allocate parlist
        if(allocated(parlist)) deallocate(parlist)
        allocate(parlist(nCycles1+nCycles2+nCycles3,size(showpar)))
    !check ranges
    if (present(ranges)) then
        if (size(ranges,dim=2)/=2_mik) then
            err=1;mess="MetroMS_Mix:ranges should include 2 real numbers"
            return
        endif
        if (size(ranges,dim=1)/=size(showpar)) then
            err=1;mess="MetroMS_Mix:ranges size mismatch"
            return
        endif
    endif
endif

!MCMC 1
call Adaptive_Metro_OAAT(f=f,x=x,fx=fx,fAux=fAux,std=std,nAdapt=nAdapt1,nCycles=nCycles1,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				UseRfortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				OutFile=trim(OutFile1),outMat=MCMCsamples,headers=headers,err=err,mess=mess)
				
If(err>0) then
	mess='Metro_Mix'//trim(mess);return
endif

if (present(outMat1)) then
    if (associated(outMat1)) deallocate(outMat1)
    allocate(outMat1(size(MCMCsamples,dim=1),size(MCMCsamples,dim=2)))
    outMat1=MCMCsamples
endif

! MCMC 2
covar=getvar(x=MCMCsamples(NINT(BurnFactor*nAdapt1*nCycles1+1):(nAdapt1*nCycles1+1),1:n),&
			package=2, method="c")
			
! Choleskize covariance
call choles_dcmp(A=covar, posDefinite=pos, err=err, message=mess)

postlist2(1:nCycles1)=MCMCsamples(2:(nAdapt1*nCycles1+1):nAdapt1,n+1)

if (present(showpar)) then
    do k=1,size(showpar)
	    parlist(1:nCycles1,k)=MCMCsamples(2:(nAdapt1*nCycles1+1):nAdapt1,showpar(k))
    enddo
endif

if (present(scaleGJ)) then
    scale=scaleGJ
else
    scale=2.4_mrk/sqrt(1._mrk*(n))			
endif

call Adaptive_Metro(f=f,x=x,fx=fx,fAux=fAux,covar=covar,scale=scale,&
		nAdapt=nAdapt2,nCycles=nCycles2,&
		MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
		DownMult=DownMult,UpMult=UpMult,&
		UseRfortran=UseRF,&
		showPar=showPar,ranges=ranges,&
		OutFile=trim(OutFile2),outMat=MCMCsamples,headers=headers,err=err,mess=mess)	
If(err>0) then
	mess='Metro_Mix'//trim(mess);return
endif

if (present(outMat2)) then
    if (associated(outMat2)) deallocate(outMat2)
    allocate(outMat2(size(MCMCsamples,dim=1),size(MCMCsamples,dim=2)))
    outMat2=MCMCsamples
endif

postlist2((nCycles1+1):(nCycles1+nCycles2))=MCMCsamples(2:(nAdapt2*nCycles2+1):nAdapt2,n+1)

if (present(showpar)) then
    do k=1,size(showpar)
	    parlist((nCycles1+1):(nCycles1+nCycles2),k)=MCMCsamples(2:(nAdapt2*nCycles2+1):nAdapt2,showpar(k))
    enddo
endif

! MCMC 3
! associate outMat
if (present(outMat3)) then
    if (present(fAux)) then
        if(associated(outMat3)) deallocate(outMat3)
        allocate(outMat3(nSim,n+1+size(fAux,dim=1)))
    else
        if(associated(outMat3)) deallocate(outMat3)
        allocate(outMat3(nSim,n+1))
    endif
endif

if (present(OutFile3)) then
    open(unit=1,file=trim(OutFile3))
    if(present(headers)) then ! trust user - no size check at the moment
		k=size(headers)
    	write((1),'(<k>A14)') headers
 	else
		if(allocated(head)) deallocate(head)
		allocate(head(n))
		do i=1,n
			head(i)='par'//trim(number_string(i))
		enddo
		if(present(fAux)) then
			p=size(fAux,dim=1)
			if(allocated(fAuxHead)) deallocate(fAuxHead)
			allocate(fAuxHead(p))
			do i=1,p
				fAuxhead(i)=' fAux'//trim(number_string(i))
			enddo
			write((1),'(<n+1+p>A14)') head(:),'Objectivefunction',fAuxHead(:)
	    	deallocate(fAuxHead)
		else
   			write((1),'(<n+1>A14)') head(:),'Objectivefunction'
		endif
		deallocate(head)
	endif
endif

if (UseRF) then
    ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0;comptsim=0
    ok=Rput("nCycles",nCycles3)
    do while(again)
        do i=1,nAdapt3
            if (present(fAux)) then
                nfaux=size(faux)
    	        call Metro_GaussJump(f=f,&
	                x=x,fx=fx,fAux=fAux,&
	                covar=(scale**2)*covar,move=move,err=err,mess=mess)
	            If(err>0) then
			        mess='Metro_Mix:'//trim(mess)
			        if(present(OutFile3)) close(1)
                    return
		        endif
	            write(1,'(<n+1+nfaux>e14.6)') x,fx,fAux
	        else
	            call Metro_GaussJump(f=f,&
	                x=x,fx=fx,fAux=fAux,&
	                covar=(scale**2)*covar,move=move,err=err,mess=mess)
	            If(err>0) then
			        mess='Metro_Mix:'//trim(mess)
			        if(present(OutFile3)) close(1)
                    return
		        endif
	            write(1,'(<n+1>e14.6)') x,fx
	        endif
        	            
        	if(present(outMat3) .and. again) then
	            outMat3(comptSim+1,1:n)=x
	            outMat3(comptSim+1,n+1)=fx
	            if (present(fAux)) then
	                outMat3(comptSim+1,(n+1+1):(n+1+nfaux))=fAux(:)
	            endif
	        endif
	        comptSim=comptSim+1
	        if (comptSim>=nSim) again=.false.
        enddo
        postlist2(nCycles1+nCycles2+compt+1)=fx
        
        if (present(showpar)) then
		    do i=1,size(showpar)
		        parlist(nCycles1+nCycles2+compt+1,i)=x(showpar(i))
		    enddo
	    endif
    	
        if (present(showpar)) then
            ok=Rput('k',size(showpar))
            ok=reval('layout(matrix(c(rep(1,k),rep(2:(k+1),rep(1,k))),k*1,2))')
        else
    !        ok=Reval('layout(c(1,2))')
        endif
        
        ok=Rput('post2', postlist2)
        ok=Rput('compt',(compt+nCycles1+nCycles2))
        ok=Reval('yminFB=quantile(post2[1:(1+compt)],0.02)')
        ok=Reval('ymaxFB=quantile(post2[1:(1+compt)],0.98)')
        ok=Reval('matplot(post2[1:(1+compt)],type="l", ylim=c(yminFB,ymaxFB), xlab="Iteration(nCycles)", ylab="Log-posterior", main="Log-posterior (From Begin)",cex.lab=2,cex.axis=2,cex.main=3)')

        if (present(showpar)) then
            ok=Rput('pars', parlist)
            ok=Rput('compt',(compt+nCycles1+nCycles2))
            ok=Rput('showpar',showpar)
            do k=1,size(showpar)
                if (present(ranges)) then
                    ok=Rput('rmin',ranges(k,1))
                    ok=Rput('rmax',ranges(k,2))
                    ok=Rput('k',k)
                    ok=Reval('matplot(pars[1:(compt+1),k],type="l", xlab="Iteration(nCycles)", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                else
                    ok=Rput('k',k)
                    ok=Reval('rmin=quantile(pars[1:(compt+1),k],0.02)')
                    ok=Reval('rmax=quantile(pars[1:(compt+1),k],0.98)')
                    ok=Reval('matplot(pars[1:(compt+1),k],type="l", xlab="Iteration(nCycles)", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                endif
            enddo
        endif
        compt=compt+1
    enddo

    if (present(OutFile3)) close(1)
    ok = Rclose()  ! End R session

else

    !ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0;comptsim=0
    !ok=Rput("nCycles",nCycles3)
    do while(again)
        do i=1,nAdapt3
            if (present(fAux)) then
                nfaux=size(faux,dim=1)
    	        call Metro_GaussJump(f=f,&
	                x=x,fx=fx,fAux=fAux,&
	                covar=(scale**2)*covar,move=move,err=err,mess=mess)
	            If(err>0) then
			        mess='Metro_Mix:'//trim(mess)
			        if(present(OutFile3)) close(1)
                    return
		        endif
	            write(1,'(<n+1+nfaux>e14.6)') x,fx,fAux
	        else
	            call Metro_GaussJump(f=f,&
	                x=x,fx=fx,fAux=fAux,&
	                covar=(scale**2)*covar,move=move,err=err,mess=mess)
	            If(err>0) then
			        mess='Metro_Mix:'//trim(mess)
			        if(present(OutFile3)) close(1)
                    return
		        endif
	            write(1,'(<n+1>e14.6)') x,fx
	        endif
	        
	        if(present(outMat3) .and. again) then
	            outMat3(comptSim+1,1:n)=x
	            outMat3(comptSim+1,n+1)=fx
	            if (present(fAux)) then
	                outMat3(comptSim+1,(n+1+1):(n+1+size(faux,dim=1)))=fAux(:)
	            endif
	        endif
	        comptSim=comptSim+1
	        if (comptSim>=nSim) again=.false.
        enddo
        
        compt=compt+1
    enddo

    if (present(OutFile3)) close(1)
    !ok = Rclose()  ! End R session

endif

deallocate(postlist2)
if (present(showPar)) deallocate(parlist)
deallocate(MCMCsamples)

end subroutine Metro_Mix


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine AdaptiveMS_Metro_OAAT(f,x,fx,fAux,std,nAdapt,nCycles,&
				MinMoveRate,MaxMoveRate,DownMult,UpMult,&
				UseRFortran,&
				showPar,ranges,&
				GRIndex,GRBurnFactor,&
				Outdoc,Filename,outMat,headers,err,mess)

!^**********************************************************************
!^* Purpose: Multi start points chain
!^**********************************************************************
!^* Programmer: Ben Renard, Xun SUN, Irstea
!^**********************************************************************
!^* Last modified:30/01/2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.
!^*		2.showpar, indice of parameters, of which mcmc chain will be shown
!^*		3.ranges, range associated to these parameters
!^* OUT
!^*		1.
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^* INOUT
!^*		1.
!^*		2.
!^*		3.
!^**********************************************************************
use utilities_dmsl_kit, only:number_string
Use RFortran

real(mrk),intent(inout)::x(:,:),std(:,:)
integer(mik),intent(in)::nAdapt,nCycles
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in),optional::OutDoc,filename,headers(:)
logical,intent(in),optional::UseRFortran
integer(mik),intent(in),optional::showpar(:)
real(mrk),intent(in),optional::ranges(:,:)
real(mrk),pointer::fx(:)
real(mrk),pointer,optional::GRIndex(:,:),outMat(:,:,:)
real(mrk),intent(in),optional::GRBurnFactor
real(mrk),intent(out),optional::fAux(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
logical::feas,isnull,again,ok
character(100), allocatable::head(:), fAuxHead(:)
logical,allocatable::move(:,:)
!real(mrk)::temp
real(mrk),allocatable::moverate(:,:),postlist(:,:),postlist2(:,:),parlist(:,:,:),phi(:,:,:),GRIndice(:,:)
real(mrk)::GRBF
character(250)::OutFile
integer(mik)::n,i,j,p,k,compt,n_chain
character(1)::nc


interface
    subroutine f(x,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface


if (present(UseRFortran)) then
    UseRF=UseRFortran
else
    UseRF=.true.
endif

if (present(GRBurnfactor)) then
    GRBF=GRBurnfactor
else
    GRBF=0.5_mrk
endif

!number of chains
n_chain=size(x,dim=2)
if (size(std,dim=2)/=n_chain) then
    err=1;mess="AdaptiveMS_Metro_OAAT:std size mismatch"
    return
endif

if(allocated(postlist)) deallocate(postlist)
allocate(postlist(nAdapt,n_chain))

if(allocated(postlist2)) deallocate(postlist2)
allocate(postlist2(n_chain,nCycles))

if(associated(fx)) deallocate(fx)
allocate(fx(n_chain))



!Handle Output file headers
n=size(x,dim=1)
if(present(outDoc) .and. present(Filename)) then
    do j=1,n_chain
        write(nc,'(I1)') j
    	open(unit=j,file=trim(OutDoc)//"C"//nc//"_"//trim(Filename),status='replace')
    enddo
	if(present(headers)) then ! trust user - no size check at the moment
		k=size(headers)
		do j=1,n_chain
    		write(j,'(<k>A14)') headers
	    enddo
	else
		if(allocated(head)) deallocate(head)
		allocate(head(n))
		do i=1,n
			head(i)='par'//trim(number_string(i))
		enddo
		if(present(fAux)) then
			p=size(fAux,dim=1)
			if(allocated(fAuxHead)) deallocate(fAuxHead)
			allocate(fAuxHead(p))
			do i=1,p
				fAuxhead(i)=' fAux'//trim(number_string(i))
			enddo
			do j=1,n_chain
    			write(j,'(<n+1+p>A14)') head(:),'Objectivefunction',fAuxHead(:)
	        enddo
			deallocate(fAuxHead)
		else
		    do j=1,n_chain
    			write(j,'(<n+1>A14)') head(:),'Objectivefunction'
		    enddo
		endif
		deallocate(head)
	endif
endif
! Check starting value
do j=1,n_chain
    if (present(fAux)) then
        call f(x=x(:,j),feas=feas,isnull=isnull,fx=fx(j),fAux=fAux(:,j),err=err,mess=mess)
    else
        call f(x=x(:,j),feas=feas,isnull=isnull,fx=fx(j),err=err,mess=mess)
    endif
    If(err>0) then
    	mess='AdaptiveMS_Metro_OAAT'//trim(mess);return
    endif
    if((.not.feas).or.isnull) then
        write(nc,'(I1)') j
	    err=1;mess='AdaptiveMS_Metro_OAAT:Fatal:Unfeasible starting point '//nc;return
    endif
    if(present(outDoc) .and. present(Filename)) then
	    if(present(fAux)) then
		    write (j,'(<n+1+p>e14.6)') x(:,j),fx(j),fAux(:,j)
    	else
	    	write (j,'(<n+1>e14.6)') x(:,j),fx(j)
	    endif
    endif
enddo

! Associate outMat
if (present(outMat)) then
    if (present(fAux)) then
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1+size(fAux,dim=1),n_chain))
        outMat(1,1:n,:)=x(:,:)
        outMat(1,n+1,:)=fx(:)
        outMat(1,(n+2):(n+1+size(fAux,dim=1)),:)=fAux(:,:)
    else
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1,n_chain))
        outMat(1,1:n,:)=x(:,:)
        outMat(1,n+1,:)=fx(:)
    endif
endif


! Allocate MoveRate stuff
if(allocated(move)) deallocate(move)
allocate(move(n,n_chain))
if(allocated(moverate)) deallocate(moverate)
allocate(moverate(n,n_chain))
if(allocated(phi)) deallocate(phi)
allocate(phi(n,nAdapt*nCycles,n_chain))
if(allocated(GRIndice)) deallocate(GRIndice)
allocate(GRIndice(n,nCycles))


!allocate parlist
if (present(showpar)) then
    !allocate parlist
        if(allocated(parlist)) deallocate(parlist)
        allocate(parlist(nCycles,size(showpar),n_chain))
    !check ranges
    if (present(ranges)) then
        if (size(ranges,dim=2)/=2_mik) then
            err=1;mess="AdaptiveMS_Metro_OAAT;ranges should include 2 real numbers"
            return
        endif
        if (size(ranges,dim=1)/=size(showpar)) then
            err=1;mess="AdaptiveMS_Metro_OAAT;ranges size mismatch"
            return
        endif
    endif
endif
    

if (UseRF) then
    ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0
    do while(again)
     
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
	        do j=1,n_chain
	            if (present(fAux)) then
    		        call Metro_OAAT_GaussJump(f=f,x=x(:,j),fx=fx(j),fAux=fAux(:,j),stdev=std(:,j),move=move(:,j),err=err,mess=mess)
    		    else
    		        call Metro_OAAT_GaussJump(f=f,x=x(:,j),fx=fx(j),stdev=std(:,j),move=move(:,j),err=err,mess=mess)
    		    endif
		        If(err>0) then
			        mess='AdaptiveMS_Metro_OAAT'//trim(mess);return
		        endif
		        postlist(i,j)=fx(j)
		        phi(:,(compt*nAdapt+i),j)=x(:,j)
		    !Write2File
		        if(present(outDoc) .and. present(Filename)) then
			        if(present(fAux)) then
				        write (j,'(<n+1+p>e14.6)') x(:,j),fx(j),fAux(:,j)
			        else
				        write (j,'(<n+1>e14.6)') x(:,j),fx(j)
			        endif
		        endif
		        if(present(outMat)) then
		            outMat(compt*nAdapt+i+1,1:n,j)=x(:,j)
		            outMat(compt*nAdapt+i+1,n+1,j)=fx(j)
		            if (present(fAux)) then
		                outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)),j)=fAux(:,j)
		            endif
		        endif
		        !Update Move Rates
		        where(move(:,j)) moverate(:,j)=moverate(:,j)+1._mrk/nAdapt
		    enddo
	    enddo
	    postlist2(:,compt+1)=fx
    		
        do k=1,n
    	    call get_GRIndice(estimand=phi(k,nint(GRBF*(compt+1)*nAdapt):((compt+1)*nAdapt),:),indice=GRIndice(k,(compt+1)),err=err,mess=mess)
	        If(err>0) then
		        mess='AdaptiveMS_Metro_OAAT'//trim(mess);return
	        endif
	    enddo
    	
	    if (present(showpar)) then
		    do j=1,size(showpar)
		        parlist(compt+1,j,:)=x(showpar(j),:)
		    enddo
	    endif
    	
        if (present(showpar)) then
            ok=Rput('k',size(showpar))
            ok=reval('layout(matrix(c(rep(1,k),rep(2,k),rep(3,k),rep(4,k),rep(5:(k+4),rep(2,k))),k*2,3))')
        else
            ok=Reval('layout(matrix(1:4,2,2))')
        endif
        ok=Rput('mr', moverate)
        ok=Reval('matplot(mr[,1],type="l", xlab="Parameter", ylab="Move Rate", main="Move Rates (Bloc)", ylim=c(0, 1),cex.lab=2,cex.axis=2,cex.main=3)')
        do j=2,n_chain
        ok=Rput('j',j)
        ok=Reval('lines(mr[,j],col=j)')
        enddo
        
        ok=Rput('GR', GRIndice)
        ok=Rput('compt',compt)
        ok=Rput('ymax',max(1.6_mrk,0.5_mrk+maxval(GRIndice(:,compt+1))))
        ok=Reval('matplot(GR[1,1:(compt+1)],type="l",ylim=c(0.8,ymax), xlab="Iteration(nCycles) for each parameter", ylab="Gelman-Rubin criterion", main="Gelman-Rubin criterion", cex.lab=2,cex.axis=2,cex.main=3,col=3)')
        do k=2,n
        ok=Rput('k',k)
        ok=Reval('lines(GR[k,1:(compt+1)],col=k+2)')
        enddo
        ok=Reval('abline(h=1,lwd=2)')
        ok=Reval('abline(h=1.2,lwd=2,col=2)')
        !ok=Reval('legend("topright",c("Each line for one parameter estimation"),lty=c(1),pch=c(NA),col=c(1),cex=1.25,lwd=2)')
        
        ok=Rput('post', postlist)
        ok=Rput('yminbloc',minval(postlist))
        ok=Rput('ymaxbloc',maxval(postlist))
        ok=Reval('matplot(post[,1],type="l", ylim=c(yminbloc,ymaxbloc), xlab="Iteration(nAdapt) for each chain", ylab="Log-posterior", main="Log-posterior (Bloc)",cex.lab=2,cex.axis=2,cex.main=3)')
        do j=2,n_chain
        ok=Rput('j',j)
        ok=Reval('lines(post[,j],col=j)')
        enddo
        
        ok=Rput('post2', postlist2)
        ok=Rput('compt',compt)
        ok=Reval('yminFB=quantile(post2[,1:(1+compt)],0.02)')
        ok=Reval('ymaxFB=quantile(post2[,1:(1+compt)],0.98)')
        ok=Reval('matplot(post2[1,1:(1+compt)],type="l", ylim=c(yminFB,ymaxFB), xlab="Iteration(nCycles) for each chain", ylab="Log-posterior", main="Log-posterior (From Begin)",cex.lab=2,cex.axis=2,cex.main=3)')
        do j=2,n_chain
        ok=Rput('j',j)
        ok=Reval('lines(post2[j,1:(1+compt)],col=j)')
        enddo
        
        
        if (present(showpar)) then
            ok=Rput('pars', parlist)
            ok=Rput('compt',compt)
            ok=Rput('showpar',showpar)
            do k=1,size(showpar)
                if (present(ranges)) then
                    ok=Rput('rmin',ranges(k,1))
                    ok=Rput('rmax',ranges(k,2))
                    ok=Rput('k',k)
                    ok=Reval('matplot(pars[1:(compt+1),k,1],type="l", xlab="Iteration(nCycles) for each chain", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                else
                    ok=Rput('k',k)
                    ok=Reval('rmin=quantile(pars[1:(compt+1),k,],0.02)')
                    ok=Reval('rmax=quantile(pars[1:(compt+1),k,],0.98)')
                    ok=Reval('matplot(pars[1:(compt+1),k,1],type="l", xlab="Iteration(nCycles) for each chain", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                endif
                do j=2,n_chain
                    ok=Rput('j',j)
                    ok=Reval('lines(pars[1:(compt+1),k,j],col=j)')
                enddo
            enddo
        endif
        

	    ! Adapt Jump standard deviation
	    where(moverate<=MinMoveRate) std=std*DownMult
	    where(moverate>=MaxMoveRate) std=std*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if(present(GRIndex)) then
        if(associated(GRIndex)) deallocate(GRIndex)
        allocate(GRIndex(n,nCycles))
        GRIndex=GRIndice
    endif
    if(present(outDoc) .and. present(Filename)) then
    do j=1,n_chain
        close(j)
    enddo
    endif
    ok = Rclose()  ! End R session

else

    !ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0;!ok=Reval("mrl<-c()")
    do while(again)
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
	        do j=1,n_chain
	            if (present(fAux)) then
    		        call Metro_OAAT_GaussJump(f=f,x=x(:,j),fx=fx(j),fAux=fAux(:,j),stdev=std(:,j),move=move(:,j),err=err,mess=mess)
		        else
		            call Metro_OAAT_GaussJump(f=f,x=x(:,j),fx=fx(j),stdev=std(:,j),move=move(:,j),err=err,mess=mess)
		        endif
		        If(err>0) then
			        mess='AdaptiveMS_Metro_OAAT'//trim(mess);return
		        endif
		        postlist(i,j)=fx(j)
		        phi(:,(compt*nAdapt+i),j)=x(:,j)
		    !Write2File
		        if(present(outDoc) .and. present(Filename)) then
			        if(present(fAux)) then
				        write (j,'(<n+1+p>e14.6)') x(:,j),fx(j),fAux(:,j)
			        else
				        write (j,'(<n+1>e14.6)') x(:,j),fx(j)
			        endif
		        endif
		        if(present(outMat)) then
		            outMat(compt*nAdapt+i+1,1:n,j)=x(:,j)
		            outMat(compt*nAdapt+i+1,n+1,j)=fx(j)
		            if (present(fAux)) then
		                outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)),j)=fAux(:,j)
		            endif
		        endif
		        !Update Move Rates
		        where(move(:,j)) moverate(:,j)=moverate(:,j)+1._mrk/nAdapt
		    enddo
	    enddo
    	
	    do k=1,n
    	    call get_GRIndice(estimand=phi(k,nint(GRBF*(compt+1)*nAdapt):((compt+1)*nAdapt),:),indice=GRIndice(k,(compt+1)),err=err,mess=mess)
	        If(err>0) then
		        mess='AdaptiveMS_Metro_OAAT'//trim(mess);return
	        endif
	    enddo
        !ok=Reval('layout(c(1, 2))')
        !ok=Rput('mr', moverate)
        !ok=Reval('matplot(mr(:,1),type="l", xlab="Parameter", ylab="Move Rate", main="Move Rates", ylim=c(0, 1))')
        !do j=2,n_chain
        !ok=Rput('j',j)
        !ok=Reval('lines(mr(:,j),col=j)')
        !enddo
        !ok=Rput('post', postlist)
        !ok=Reval('matplot(post(:,1),type="l", xlab="Iteration", ylab="Log-posterior", main="Log-posterior")')
        !do j=2,n_chain
        !ok=Rput('j',j)
        !ok=Reval('lines(post(:,j),col=j)')
        !enddo

	    ! Adapt Jump standard deviation
	    where(moverate<=MinMoveRate) std=std*DownMult
	    where(moverate>=MaxMoveRate) std=std*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if(present(GRIndex)) then
        if(associated(GRIndex)) deallocate(GRIndex)
        allocate(GRIndex(n,nCycles))
        GRIndex=GRIndice
    endif

    if(present(outDoc) .and. present(Filename)) then
    do j=1,n_chain
        close(j)
    enddo
    endif
    !ok = Rclose()  ! End R session

endif

deallocate(move)
deallocate(moverate)
deallocate(postlist)
deallocate(postlist2)
if (present(showPar)) deallocate(parlist)
deallocate(phi)
deallocate(GRIndice)

end subroutine AdaptiveMS_Metro_OAAT


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine AdaptiveMS_Metro(f,x,fx,fAux,covar,scale,nAdapt,nCycles,&
				MinMoveRate,MaxMoveRate,DownMult,UpMult,&
				UseRFortran,&
				showPar,ranges,&
				GRIndex,GRBurnfactor,&
				outDoc,Filename,outMat,headers,err,mess)

!^**********************************************************************
!^* Purpose: Metropolis with adaptive scale factor (Multi chains)
!^**********************************************************************
!^* Programmer: Ben Renard, Xun SUN, Irstea
!^**********************************************************************
!^* Last modified:31/01/2012
!^**********************************************************************
!^* Comments: Only the scale factor is adapted - var is choleskized
!^* out of this sub
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
use utilities_dmsl_kit, only:number_string
Use RFortran

real(mrk),intent(inout)::x(:,:),covar(:,:,:),scale(:)
integer(mik),intent(in)::nAdapt,nCycles
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in),optional::OutDoc, Filename, headers(:)
logical,intent(in),optional::UseRFortran
integer(mik),intent(in),optional::showpar(:)
real(mrk),intent(in),optional::ranges(:,:),GRBurnfactor
real(mrk),pointer::fx(:)
real(mrk),pointer,optional::GRIndex(:,:),outMat(:,:,:)
real(mrk),intent(out),optional::fAux(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
logical::feas,isnull,again,ok
character(100), allocatable::head(:), fAuxHead(:)
logical,allocatable::move(:)
real(mrk),allocatable::moverate(:),postlist(:,:),postlist2(:,:),parlist(:,:,:),phi(:,:,:),GRIndice(:,:)
integer(mik)::n,i,j,p,k,compt,n_chain
real(mrk)::GRBF
character(1)::nc

interface
    subroutine f(x,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface

if (present(UseRFortran)) then
    UseRF=UseRFortran
else
    UseRF=.true.
endif

if (present(GRBurnfactor)) then
    GRBF=GRBurnfactor
else
    GRBF=0.5_mrk
endif

!number of chains
n_chain=size(x,dim=2)
if (size(scale)/=n_chain) then
    err=1;mess="AdaptiveMS_Metro:scale size mismatch"
    return
endif

if(allocated(postlist)) deallocate(postlist)
allocate(postlist(nAdapt,n_chain))

if(allocated(postlist2)) deallocate(postlist2)
allocate(postlist2(n_chain,nCycles))

if(associated(fx)) deallocate(fx)
allocate(fx(n_chain))


!Handle Output file headers
n=size(x,dim=1)
if(present(outDoc) .and. present(Filename)) then
    do j=1,n_chain
        write(nc,'(I1)') j
    	open(unit=j,file=trim(OutDoc)//"C"//nc//"_"//trim(Filename),status='replace')
    enddo
	if(present(headers)) then ! trust user - no size check at the moment
		k=size(headers)
		do j=1,n_chain
    		write(j,'(<k>A14)') headers
        enddo
	else
		if(allocated(head)) deallocate(head)
		allocate(head(n))
		do i=1,n
			head(i)='par'//trim(number_string(i))
		enddo
		if(present(fAux)) then
			p=size(fAux,dim=1)
			if(allocated(fAuxHead)) deallocate(fAuxHead)
			allocate(fAuxHead(p))
			do i=1,p
				fAuxhead(i)=' fAux'//trim(number_string(i))
			enddo
			do j=1,n_chain
    			write(j,'(<n+1+p>A14)') head(:),'Objectivefunction',fAuxHead(:)
    		enddo
	    	deallocate(fAuxHead)
		else
		    do j=1,n_chain
    			write(j,'(<n+1>A14)') head(:),'Objectivefunction'
	        enddo
		endif
		deallocate(head)
	endif
endif

! Check starting value
do j=1,n_chain
    if (present(fAux)) then
        call f(x=x(:,j),feas=feas,isnull=isnull,fx=fx(j),fAux=fAux(:,j),err=err,mess=mess)
    else
        call f(x=x(:,j),feas=feas,isnull=isnull,fx=fx(j),err=err,mess=mess)
    endif
    If(err>0) then
	    mess='AdaptiveMS_Metro'//trim(mess);return
    endif
    if((.not.feas).or.isnull) then
	    err=1;mess='AdaptiveMS_Metro:Fatal:Unfeasible starting point';return
    endif
    if(present(outDoc) .and. present(Filename)) then
	    if(present(fAux)) then
		    write (j,'(<n+1+p>e14.6)') x(:,j),fx(j),fAux(:,j)
	    else
		    write (j,'(<n+1>e14.6)') x(:,j),fx(j)
	    endif
    endif
enddo

! associate outMat
if (present(outMat)) then
    if (present(fAux)) then
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1+size(fAux,dim=1),n_chain))
        outMat(1,1:n,:)=x(:,:)
        outMat(1,n+1,:)=fx(:)
        outMat(1,(n+2):(n+1+size(fAux,dim=1)),:)=fAux(:,:)
    else
        if(associated(outMat)) deallocate(outMat)
        allocate(outMat(nAdapt*nCycles+1,n+1,n_chain))
        outMat(1,1:n,:)=x(:,:)
        outMat(1,n+1,:)=fx(:)
    endif
endif

! Allocate MoveRate stuff
if(allocated(move)) deallocate(move)
allocate(move(n_chain))
if(allocated(moverate)) deallocate(moverate)
allocate(moverate(n_chain))
if(allocated(phi)) deallocate(phi)
allocate(phi(n,nAdapt*nCycles,n_chain))
if(allocated(GRIndice)) deallocate(GRIndice)
allocate(GRIndice(n,nCycles))


!allocate parlist
if (present(showpar)) then
    !allocate parlist
        if(allocated(parlist)) deallocate(parlist)
        allocate(parlist(nCycles,size(showpar),n_chain))
    !check ranges
    if (present(ranges)) then
        if (size(ranges,dim=2)/=2_mik) then
            err=1;mess="AdaptiveMS_Metro_OAAT:ranges should include 2 real numbers"
            return
        endif
        if (size(ranges,dim=1)/=size(showpar)) then
            err=1;mess="AdaptiveMS_Metro_OAAT:ranges size mismatch"
            return
        endif
    endif
endif


if (UseRF) then
    ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0;
    ok=Rput("n_chain",n_chain)
    ok=Rput("nCycles",nCycles)
    ok=Reval("mrl<-matrix(0,nCycles,n_chain)")
    do while(again)
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
	        do j=1,n_chain
	            if (present(fAux)) then
		            call Metro_GaussJump(f=f,x=x(:,j),fx=fx(j),fAux=fAux(:,j),&
		                        covar=(scale(j)**2)*covar(:,:,j),move=move(j),err=err,mess=mess)
		        else
		            call Metro_GaussJump(f=f,x=x(:,j),fx=fx(j),&
		                        covar=(scale(j)**2)*covar(:,:,j),move=move(j),err=err,mess=mess)
		        endif
		        If(err>0) then
			        mess='AdaptiveMS_Metro'//trim(mess);return
		        endif
		        postlist(i,j)=fx(j)
		        phi(:,(compt*nAdapt+i),j)=x(:,j)
		        !Write2File
		        if(present(outDoc) .and. present(Filename)) then
			        if(present(fAux)) then
				        write (j,'(<n+1+p>e14.6)') x(:,j),fx(j),fAux(:,j)
			        else
				        write (j,'(<n+1>e14.6)') x(:,j),fx(j)
			        endif
		        endif
		        if(present(outMat)) then
		            outMat(compt*nAdapt+i+1,1:n,j)=x(:,j)
		            outMat(compt*nAdapt+i+1,n+1,j)=fx(j)
		            if (present(fAux)) then
		                outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)),j)=fAux(:,j)
		            endif
		        endif
		        !Update Move Rates
		        if(move(j)) moverate(j)=moverate(j)+1._mrk/nAdapt
	        enddo
	    enddo
        postlist2(:,compt+1)=fx

	    do k=1,n
    	    call get_GRIndice(estimand=phi(k,nint(GRBF*(compt+1)*nAdapt):((compt+1)*nAdapt),:),indice=GRIndice(k,(compt+1)),err=err,mess=mess)
	        If(err>0) then
		        mess='AdaptiveMS_Metro:'//trim(mess);return
	        endif
	    enddo
    	
	    if (present(showpar)) then
		    do j=1,size(showpar)
		        parlist(compt+1,j,:)=x(showpar(j),:)
		    enddo
	    endif
    	
        if (present(showpar)) then
            ok=Rput('k',size(showpar))
            ok=reval('layout(matrix(c(rep(1,k),rep(2,k),rep(3,k),rep(4,k),rep(5:(k+4),rep(2,k))),k*2,3))')
        else
            ok=Reval('layout(matrix(1:4,2,2))')
        endif
        ok=Rput('mr', moverate)
        ok=Rput('compt', compt+1)
        ok=Reval('mrl[compt,]<-mr')
        ok=Reval('matplot(mrl[1:compt,1],type="l", xlab="Iteration(nCycles) for each chain", ylab="Move Rate", main="Move Rates (From Begin)", ylim=c(0, 1),cex.lab=2,cex.axis=2,cex.main=3)')
        do j=2,n_chain
        ok=Rput('j',j)
        ok=Reval('lines(mrl[1:compt,j],col=j)')
        enddo
        
        ok=Rput('GR', GRIndice)
        ok=Rput('compt',compt)
        ok=Rput('ymax',max(1.6_mrk,0.5_mrk+maxval(GRIndice(:,compt+1))))
        ok=Reval('matplot(GR[1,1:(compt+1)],type="l",ylim=c(0.8,ymax), xlab="Iteration(nCycles) for each parameter", ylab="Gelman-Rubin criterion", main="Gelman-Rubin criterion", cex.lab=2,cex.axis=2,cex.main=3,col=3)')
        do k=2,n
        ok=Rput('k',k)
        ok=Reval('lines(GR[k,1:(compt+1)],col=k+2)')
        enddo
        ok=Reval('abline(h=1,lwd=2)')
        ok=Reval('abline(h=1.2,lwd=2,col=2)')
        !ok=Reval('legend("topright",c("for each parameter"),lty=c(1),pch=c(NA),col=c(1),cex=1.25,lwd=2)')
        
        ok=Rput('post', postlist)
        ok=Rput('yminbloc',minval(postlist))
        ok=Rput('ymaxbloc',maxval(postlist))
        ok=Reval('matplot(post[,1],type="l", ylim=c(yminbloc,ymaxbloc), xlab="Iteration for each chain", ylab="Log-posterior", main="Log-posterior (Bloc)",cex.lab=2,cex.axis=2,cex.main=3)')
        do j=2,n_chain
        ok=Rput('j',j)
        ok=Reval('lines(post[,j],col=j)')
        enddo
        
        ok=Rput('post2', postlist2)
        ok=Rput('compt',compt)
        ok=Reval('yminFB=quantile(post2[,1:(1+compt)],0.02)')
        ok=Reval('ymaxFB=quantile(post2[,1:(1+compt)],0.98)')
        ok=Reval('matplot(post2[1,1:(1+compt)],type="l", ylim=c(yminFB,ymaxFB), xlab="Iteration(nCycles) for each chain", ylab="Log-posterior", main="Log-posterior (From Begin)",cex.lab=2,cex.axis=2,cex.main=3)')
        do j=2,n_chain
        ok=Rput('j',j)
        ok=Reval('lines(post2[j,1:(1+compt)],col=j)')
        enddo

        
        if (present(showpar)) then
            ok=Rput('pars', parlist)
            ok=Rput('compt',compt)
            ok=Rput('showpar',showpar)
            do k=1,size(showpar)
                if (present(ranges)) then
                    ok=Rput('rmin',ranges(k,1))
                    ok=Rput('rmax',ranges(k,2))
                    ok=Rput('k',k)
                    ok=Reval('matplot(pars[1:(compt+1),k,1],type="l", xlab="Iteration(nCycles) for each chain", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                else
                    ok=Rput('k',k)
                    ok=Reval('rmin=quantile(pars[1:(compt+1),k,],0.02)')
                    ok=Reval('rmax=quantile(pars[1:(compt+1),k,],0.98)')
                    ok=Reval('matplot(pars[1:(compt+1),k,1],type="l", xlab="Iteration(nCycles) for each chain", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                endif
                do j=2,n_chain
                    ok=Rput('j',j)
                    ok=Reval('lines(pars[1:(compt+1),k,j],col=j)')
                enddo
            enddo
        endif
        

	    ! Adapt Jump standard deviation
	    where(moverate<=MinMoveRate) scale=scale*DownMult
	    where(moverate>=MaxMoveRate) scale=scale*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if(present(GRIndex)) then
        if(associated(GRIndex)) deallocate(GRIndex)
        allocate(GRIndex(n,nCycles))
        GRIndex=GRIndice
    endif

    if(present(outDoc) .and. present(Filename)) then
        do j=1,n_chain
            close(j)
        enddo
    endif
    ok = Rclose()  ! End R session

else
    !ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0
    do while(again)
	    !Start sampling
	    moverate=0._mrk
	    do i=1,nAdapt
	        do j=1,n_chain
		        if (present(fAux)) then
    		        call Metro_GaussJump(f=f,x=x(:,j),fx=fx(j),fAux=fAux(:,j),&
    		                    covar=(scale(j)**2)*covar(:,:,j),move=move(j),err=err,mess=mess)
	            else
	                call Metro_GaussJump(f=f,x=x(:,j),fx=fx(j),&
    		                    covar=(scale(j)**2)*covar(:,:,j),move=move(j),err=err,mess=mess)
	            endif
		        If(err>0) then
			        mess='AdaptiveMS_Metro'//trim(mess);return
		        endif
		        postlist(i,j)=fx(j)
		        phi(:,(compt*nAdapt+i),j)=x(:,j)
		        !Write2File
		        if(present(outDoc) .and. present(Filename)) then
			        if(present(fAux)) then
				        write (j,'(<n+1+p>e14.6)') x(:,j),fx(j),fAux(:,j)
			        else
				        write (j,'(<n+1>e14.6)') x(:,j),fx(j)
			        endif
		        endif
		        if(present(outMat)) then
		            outMat(compt*nAdapt+i+1,1:n,j)=x(:,j)
		            outMat(compt*nAdapt+i+1,n+1,j)=fx(j)
		            if (present(fAux)) then
		                outMat(compt*nAdapt+i+1,(n+1+1):(n+1+size(fAux,dim=1)),j)=fAux(:,j)
		            endif
		        endif
		        !Update Move Rates
		        if(move(j)) moverate(j)=moverate(j)+1._mrk/nAdapt
	        enddo
	    enddo
    	
	    do k=1,n
    	    call get_GRIndice(estimand=phi(k,nint(GRBF*(compt+1)*nAdapt):((compt+1)*nAdapt),:),indice=GRIndice(k,(compt+1)),err=err,mess=mess)
	        If(err>0) then
		        mess='AdaptiveMS_Metro:'//trim(mess);return
	        endif
	    enddo

	    ! Adapt Jump standard deviation
	    where(moverate<=MinMoveRate) scale=scale*DownMult
	    where(moverate>=MaxMoveRate) scale=scale*UpMult
	    compt=compt+1
	    if(nCycles>=1) then
		    if(compt>=nCycles) again=.false.
	    else
	        write(*,*) 'Again?'
		    read(*, *) again
	    endif
    enddo

    if(present(GRIndex)) then
        if(associated(GRIndex)) deallocate(GRIndex)
        allocate(GRIndex(n,nCycles))
        GRIndex=GRIndice
    endif

    if(present(outDoc) .and. present(Filename)) then
        do j=1,n_chain
            close(j)
        enddo
    endif
    !ok = Rclose()  ! End R session

endif


deallocate(move)
deallocate(moverate)
deallocate(postlist)
deallocate(postlist2)
if (present(showPar)) deallocate(parlist)
deallocate(phi)
deallocate(GRIndice)


end subroutine AdaptiveMS_Metro

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine MetroMS_Mix(f,x,fx,fAux,std,&
                nAdapt1,nCycles1,BurnFactor,&
                nAdapt2,nCycles2,&
                nSim,scaleGJ,GRBlocSize,&
				MinMoveRate,MaxMoveRate,DownMult,UpMult,&
				UseRFortran,&
				showPar,ranges,&
				GRIndex,GRBurnfactor,&
				OutDoc,Filename1,Filename2,Filename3,headers,err,mess)
!^**********************************************************************
!^* Purpose: Multi chains mixed algorithme with different starting points
!^**********************************************************************
!^* Programmer: Ben Renard & Xun SUN, Irstea
!^**********************************************************************
!^* Last modified: 31/01/2012
!^**********************************************************************
!^* Comments: There are three different algorithms inside, and the GR Index
!^*         is computed separately with the values of each algorithm.
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* OUT
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* INOUT
!^*     1.
!^*     2.
!^**********************************************************************
use numerix_dmsl_kit, only:normaldev,getvar
use linalg_dmsl_kit, only:choles_dcmp
use utilities_dmsl_kit, only:number_string
use RFortran

real(mrk),intent(inout)::x(:,:),std(:,:)
real(mrk),intent(inout),optional::scaleGJ(:)
integer(mik),intent(in),optional::GRBlocSize
integer(mik),intent(in)::nAdapt1,nCycles1,nAdapt2,nCycles2,nSim
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor
character(*), intent(in),optional::OutDoc,filename1,filename2,filename3,headers(:)
logical,intent(in),optional::UseRFortran
integer(mik),intent(in),optional::showpar(:)
real(mrk),intent(in),optional::ranges(:,:),GRBurnfactor
real(mrk),pointer::fx(:)
real(mrk),pointer,optional::GRIndex(:,:)
real(mrk),intent(out),optional::fAux(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
real(mrk),allocatable::scale(:),covar(:,:,:),postlist2(:,:),parlist(:,:,:),phi(:,:,:),GR3(:,:),GR(:,:)
real(mrk),pointer::MCMCsamples(:,:,:)=>NULL(),GR1(:,:)=>NULL(),GR2(:,:)=>NULL()
integer(mik)::n,nFaux,n_chain,nAdapt3,nCycles3
character(100), allocatable::head(:), fAuxHead(:)
integer(mik)::i,j,k,p,compt,comptSim,ok
logical,allocatable::move(:)
real(mrk)::GRBF
logical::pos,again
character(1)::nc

interface
    subroutine f(x,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface

if (present(UseRFortran)) then
    UseRF=UseRFortran
else
    UseRF=.true.
endif

if (present(GRBurnfactor)) then
    GRBF=GRBurnfactor
else
    GRBF=0.5_mrk
endif

if (present(GRBlocSize)) then
    nAdapt3=GRBlocSize
else
    nAdapt3=nAdapt2
endif

nCycles3=Int(nSim/nAdapt3,mik)

err=0;mess=''
n=size(x,dim=1)
n_chain=size(x,dim=2)


if(allocated(postlist2)) deallocate(postlist2)
allocate(postlist2(n_chain,nCycles1+nCycles2+nCycles3))
if(allocated(phi)) deallocate(phi)
allocate(phi(n,nCycles1*nAdapt1+nCycles2*nAdapt2+nCycles3*nAdapt3,n_chain))

!allocate parlist
if (present(showpar)) then
    !allocate parlist
        if(allocated(parlist)) deallocate(parlist)
        allocate(parlist(nCycles1+nCycles2+nCycles3,size(showpar),n_chain))
    !check ranges
    if (present(ranges)) then
        if (size(ranges,dim=2)/=2_mik) then
            err=1;mess="MetroMS_Mix:ranges should include 2 real numbers"
            return
        endif
        if (size(ranges,dim=1)/=size(showpar)) then
            err=1;mess="MetroMS_Mix:ranges size mismatch"
            return
        endif
    endif
endif

    
!MCMC 1
call AdaptiveMS_Metro_OAAT(f=f,x=x,fx=fx,fAux=fAux,std=std,nAdapt=nAdapt1,nCycles=nCycles1,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				UseRfortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				GRIndex=GR1,GRBurnfactor=GRBurnfactor,&
				Outdoc=outdoc,Filename=Filename1,outMat=MCMCsamples,headers=headers,err=err,mess=mess)
				
If(err>0) then
	mess='MetroMS_Mix'//trim(mess);return
endif

! MCMC 2    
if(allocated(covar)) deallocate(covar)
allocate(covar(n,n,n_chain))

do j=1,n_chain
    do i=1,nAdapt1*nCycles1
	    phi(:,i,j)=MCMCsamples(i+1,1:n,j)
    enddo

    covar(:,:,j)=getvar(x=MCMCsamples(NINT(BurnFactor*nAdapt1*nCycles1+1):(nAdapt1*nCycles1+1),1:n,j),&
			package=2, method="c")
    ! Choleskize covariance
    call choles_dcmp(A=covar(:,:,j), posDefinite=pos, err=err, message=mess)
    postlist2(j,1:nCycles1)=MCMCsamples(2:(nAdapt1*nCycles1+1):nAdapt1,n+1,j)
enddo

if (present(showpar)) then
    do k=1,size(showpar)
	    parlist(1:nCycles1,k,:)=MCMCsamples(2:(nAdapt1*nCycles1+1):nAdapt1,showpar(k),:)
    enddo
endif

if(allocated(scale)) deallocate(scale)
allocate(scale(N_chain))

if (present(scaleGJ)) then
    scale=scaleGJ
else
    scale=2.4_mrk/sqrt(1._mrk*(n))			
endif

call AdaptiveMS_Metro(f=f,x=x,fx=fx,fAux=fAux,covar=covar,scale=scale,&
		nAdapt=nAdapt2,nCycles=nCycles2,&
		MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
		DownMult=DownMult,UpMult=UpMult,&
		UseRfortran=UseRF,&
		showPar=showPar,ranges=ranges,&
		GRIndex=GR2,GRBurnfactor=GRBurnfactor,&
		Outdoc=outdoc,Filename=Filename2,outMat=MCMCsamples,headers=headers,err=err,mess=mess)	
If(err>0) then
	mess='MetroMS_Mix'//trim(mess);return
endif


! MCMC 3
do j=1,n_chain
    do i=1,nAdapt2*nCycles2
	    phi(:,(i+nAdapt1*nCycles1),j)=MCMCsamples(i+1,1:n,j)
    enddo
    postlist2(j,(nCycles1+1):(nCycles1+nCycles2))=MCMCsamples(2:(nAdapt2*nCycles2+1):nAdapt2,n+1,j)
enddo

if (present(showpar)) then
    do k=1,size(showpar)
	    parlist((nCycles1+1):(nCycles1+nCycles2),k,:)=MCMCsamples(2:(nAdapt2*nCycles2+1):nAdapt2,showpar(k),:)
    enddo
endif

if(allocated(move)) deallocate(move)
allocate(move(n_chain))
if(allocated(GR3)) deallocate(GR3)
allocate(GR3(n,nCycles3))
if(allocated(GR)) deallocate(GR)
allocate(GR(n,(nCycles1+nCycles2+nCycles3)))

if(present(outDoc) .and. present(Filename3)) then
    do j=1,n_chain
        write(nc,'(I1)') j
        open(unit=(j),file=trim(OutDoc)//"C"//nc//"_"//trim(Filename3))
    enddo

	if(present(headers)) then ! trust user - no size check at the moment
		k=size(headers)
		do j=1,n_chain
    		write((j),'(<k>A14)') headers
        enddo
	else
		if(allocated(head)) deallocate(head)
		allocate(head(n))
		do i=1,n
			head(i)='par'//trim(number_string(i))
		enddo
		if(present(fAux)) then
			p=size(fAux,dim=1)
			if(allocated(fAuxHead)) deallocate(fAuxHead)
			allocate(fAuxHead(p))
			do i=1,p
				fAuxhead(i)=' fAux'//trim(number_string(i))
			enddo
			do j=1,n_chain
    			write((j),'(<n+1+p>A14)') head(:),'Objectivefunction',fAuxHead(:)
    		enddo
	    	deallocate(fAuxHead)
		else
		    do j=1,n_chain
    			write((j),'(<n+1>A14)') head(:),'Objectivefunction'
	        enddo
		endif
		deallocate(head)
	endif
endif

if (UseRF) then
    ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0;comptsim=0
    ok=Rput("n_chain",n_chain)
    ok=Rput("nCycles",nCycles3)
    do while(again)
        do i=1,nAdapt3
            do j=1,n_chain
                if (present(fAux)) then
                    call Metro_GaussJump(f=f,&
	                    x=x(:,j),fx=fx(j),fAux=fAux(:,j),&
	                    covar=(scale(j)**2)*covar(:,:,j),move=move(j),err=err,mess=mess)
	            else
	                call Metro_GaussJump(f=f,&
	                    x=x(:,j),fx=fx(j),&
	                    covar=(scale(j)**2)*covar(:,:,j),move=move(j),err=err,mess=mess)
	            endif
	            If(err>0) then
			        mess='MetroMS_Mix:'//trim(mess);return
		        endif
                phi(:,(nAdapt1*nCycles1+nAdapt2*nCycles2+compt*nAdapt3+i),j)=x(:,j)   
    		            
	            if(present(outDoc) .and. present(Filename3)) then        
                    if (present(fAux)) then
                        nfaux=size(faux,dim=1)
	                    write((j),'(<n+1+nfaux>e14.6)') x(:,j),fx(j),fAux(:,j)
	                else
	                    write((j),'(<n+1>e14.6)') x(:,j),fx(j)
	                endif
	            endif
            enddo
            comptSim=comptSim+1
            if (comptSim>=nSim) again=.false.
        enddo
        postlist2(:,nCycles1+nCycles2+compt+1)=fx
        
        do k=1,n
    	    call get_GRIndice(estimand=phi(k,nint(GRBF*(nAdapt1*nCycles1+nAdapt2*nCycles2+(compt+1)*nAdapt3)):(nAdapt1*nCycles1+nAdapt2*nCycles2+(compt+1)*nAdapt3),:),&
    	                indice=GR3(k,(compt+1)),err=err,mess=mess)
	        If(err>0) then
		        mess='MetroMS_Mix:'//trim(mess);return
	        endif
	    enddo
        
        if (present(showpar)) then
		    do j=1,size(showpar)
		        parlist(nCycles1+nCycles2+compt+1,j,:)=x(showpar(j),:)
		    enddo
	    endif
    	
	    if (present(showpar)) then
            ok=Rput('k',size(showpar))
            ok=reval('layout(matrix(c(rep(1,k),rep(2,k),rep(3:(k+2),rep(2,k))),k*2,2))')
        else
            ok=Reval('layout(c(1,2))')
        endif
        
        GR(:,1:nCycles1)=GR1
        GR(:,(nCycles1+1):(nCycles1+nCycles2))=GR2
        GR(:,(nCycles1+nCycles2+1):(nCycles1+nCycles2+nCycles3))=GR3
        ok=Rput('Nc1', nCycles1)
        ok=Rput('Nc2', nCycles1+nCycles2)
	    ok=Rput('GR', GR)
        ok=Rput('compt',(compt+nCycles1+nCycles2))
        ok=Rput('ymax',max(1.6_mrk,0.5_mrk+maxval(GR(:,compt+nCycles1+nCycles2+1))))
        ok=Reval('matplot(GR[1,1:(compt+1)],type="l",ylim=c(0.8,ymax), xlab="Iteration(nCycles) for each parameter", ylab="Gelman-Rubin criterion", main="Gelman-Rubin criterion", cex.lab=2,cex.axis=2,cex.main=3,col=3)')
        ok=Reval('abline(v=Nc1)')
        ok=Reval('abline(v=Nc2)')
        do k=2,n
        ok=Rput('k',k)
        ok=Reval('lines(GR[k,1:(compt+1)],col=k+2)')
        enddo
        ok=Reval('abline(h=1,lwd=2)')
        ok=Reval('abline(h=1.2,lwd=2,col=2)')
        !ok=Reval('legend("topright",c("for each parameter"),lty=c(1),pch=c(NA),col=c(1),cex=1.25,lwd=2)')

        ok=Rput('post2', postlist2)
        ok=Rput('compt',(compt+nCycles1+nCycles2))
        ok=Reval('yminFB=quantile(post2[,1:(1+compt)],0.02)')
        ok=Reval('ymaxFB=quantile(post2[,1:(1+compt)],0.98)')
        ok=Reval('matplot(post2[1,1:(1+compt)],type="l", ylim=c(yminFB,ymaxFB), xlab="Iteration(nCycles) for each chain", ylab="Log-posterior", main="Log-posterior (From Begin)",cex.lab=2,cex.axis=2,cex.main=3)')
        do j=2,n_chain
        ok=Rput('j',j)
        ok=Reval('lines(post2[j,1:(1+compt)],col=j)')
        enddo

        if (present(showpar)) then
            ok=Rput('pars', parlist)
            ok=Rput('compt',(compt+nCycles1+nCycles2))
            ok=Rput('showpar',showpar)
            do k=1,size(showpar)
                if (present(ranges)) then
                    ok=Rput('rmin',ranges(k,1))
                    ok=Rput('rmax',ranges(k,2))
                    ok=Rput('k',k)
                    ok=Reval('matplot(pars[1:(compt+1),k,1],type="l", xlab="Iteration(nCycles) for each chain", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                else
                    ok=Rput('k',k)
                    ok=Reval('rmin=quantile(pars[1:(compt+1),k,],0.02)')
                    ok=Reval('rmax=quantile(pars[1:(compt+1),k,],0.98)')
                    ok=Reval('matplot(pars[1:(compt+1),k,1],type="l", xlab="Iteration(nCycles) for each chain", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                endif
                do j=2,n_chain
                    ok=Rput('j',j)
                    ok=Reval('lines(pars[1:(compt+1),k,j],col=j)')
                enddo
            enddo
        endif
        
        compt=compt+1
     
    enddo

else
    again=.true.;compt=0;comptsim=0
    do while(again)
        do i=1,nAdapt3
            do j=1,n_chain
                if (present(fAux)) then
                    call Metro_GaussJump(f=f,&
	                    x=x(:,j),fx=fx(j),fAux=fAux(:,j),&
	                    covar=(scale(j)**2)*covar(:,:,j),move=move(j),err=err,mess=mess)
	            else
	                call Metro_GaussJump(f=f,&
	                    x=x(:,j),fx=fx(j),&
	                    covar=(scale(j)**2)*covar(:,:,j),move=move(j),err=err,mess=mess)
	            endif
	            If(err>0) then
			        mess='MetroMS_Mix:'//trim(mess);return
		        endif
                phi(:,(nAdapt1*nCycles1+nAdapt2*nCycles2+compt*nAdapt3+i),j)=x(:,j)   
    		            
	            if(present(outDoc) .and. present(Filename3)) then        
                    if (present(fAux)) then
                        nfaux=size(faux,dim=1)
	                    write((j),'(<n+1+nfaux>e14.6)') x(:,j),fx(j),fAux(:,j)
	                else
	                    write((j),'(<n+1>e14.6)') x(:,j),fx(j)
	                endif
	            endif
            enddo
            comptSim=comptSim+1
            if (comptSim>=nSim) again=.false.
        enddo
        
        do k=1,n
    	    call get_GRIndice(estimand=phi(k,nint(GRBF*(nAdapt1*nCycles1+nAdapt2*nCycles2+(compt+1)*nAdapt3)):(nAdapt1*nCycles1+nAdapt2*nCycles2+(compt+1)*nAdapt3),:),&
    	                indice=GR3(k,(compt+1)),err=err,mess=mess)
	        If(err>0) then
		        mess='MetroMS_Mix:'//trim(mess);return
	        endif
	    enddo
        compt=compt+1
     
    enddo

endif

if(present(GRIndex)) then
    if(associated(GRIndex)) deallocate(GRIndex)
    allocate(GRIndex(n,nCycles1+nCycles2+nCycles3))
    GRIndex=GR
endif

if(present(outDoc) .and. present(Filename3)) then
    do j=1,n_chain
        close(j)
    enddo
endif

deallocate(move)
deallocate(scale)
deallocate(covar)
deallocate(postlist2)
if (present(showPar)) deallocate(parlist)
deallocate(phi)
deallocate(GR)
deallocate(GR3)
deallocate(MCMCsamples)
deallocate(GR1)
deallocate(GR2)

end subroutine MetroMS_Mix


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MetroMS_Mix2(f,x,fx,fAux,std,&
                nAdapt1,nCycles1,BurnFactor,&
                nAdapt2,nCycles2,&
                nSim,scaleGJ,GRBlocSize,&
				MinMoveRate,MaxMoveRate,DownMult,UpMult,&
				n_chain,&
				UseRFortran,&
				showPar,ranges,&
				GRIndex,GRBurnfactor,&
				OutDoc,Filename1,Filename2,Filename3,&
				OutMat1, OutMat2, OutMat3,&
				headers,err,mess)
!^**********************************************************************
!^* Purpose: Multi chains mixed algorithme with ONE starting point
!^**********************************************************************
!^* Programmer: Xun SUN, Irstea
!^**********************************************************************
!^* Last modified:09/02/2012
!^**********************************************************************
!^* Comments: Starts with a single points. After OAAT and Adaptive_Metro, 
!^*      n re-starting points are randomly selected from the previous results.
!^*      Then go with Gaussian jumps.
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     XXX. OutDoc: Output Directory where MCMC files will be written
!^*     XXX. Filename1: MCMC file for Adaptive-OAAT step
!^*     XXX. Filename2: MCMC file for Adaptive-Metro step
!^*     XXX. Filename3: MCMC file for production step; a prefix 'Ck_' will be added for the kth chain
!^* OUT
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* INOUT
!^*     1.
!^*     2.
!^**********************************************************************
use numerix_dmsl_kit, only:normaldev,getvar,uniran
use linalg_dmsl_kit, only:choles_dcmp
use utilities_dmsl_kit, only:number_string
use RFortran

real(mrk),intent(inout)::x(:),std(:)
real(mrk),intent(inout),optional::scaleGJ
integer(mik),intent(in),optional::GRBlocSize
integer(mik),intent(in)::nAdapt1,nCycles1,nAdapt2,nCycles2,nSim,n_chain
real(mrk),intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult,BurnFactor
character(*), intent(in),optional::OutDoc,filename1,filename2,filename3,headers(:)
real(mrk),pointer,optional::outMat1(:,:),outMat2(:,:),outMat3(:,:,:) !outMat contains the 3 dimentional matrix, the additional dimension is n_chain
!character(*), intent(in)::
logical,intent(in),optional::UseRFortran
integer(mik),intent(in),optional::showpar(:)
real(mrk),intent(in),optional::ranges(:,:),GRBurnfactor
real(mrk),pointer::fx(:)
real(mrk),pointer,optional::GRIndex(:,:)
real(mrk),intent(out),optional::fAux(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
real(mrk),allocatable::covar(:,:),MCMCsamples2(:,:),postlist2(:,:),parlist(:,:,:),phi(:,:,:),GR3(:,:),startTemp(:,:),unifR(:)
real(mrk),pointer::MCMCsamples(:,:)=>NULL(),MCMCs2temp(:,:)=>NULL()
integer(mik)::n,nFaux,nAdapt3,nCycles3
character(100), allocatable::head(:), fAuxHead(:)
integer(mik)::i,j,k,p,compt,comptSim,ok
logical,allocatable::move(:)
real(mrk)::GRBF,fxtemp,scale
logical::pos,again
character(1)::nc
character(500)::OutDocTemp

interface
    subroutine f(x,feas,isnull,fx,fAux,err,mess)
        use kinds_dmsl_kit
        real(mrk),intent(in)::x(:)
        logical,intent(out)::feas,isnull
        real(mrk),intent(out)::fx
        real(mrk),intent(out),optional::fAux(:)
        integer(mik),intent(out)::err
        character(*),intent(out)::mess
    end subroutine f
end interface

if (present(UseRFortran)) then
    UseRF=UseRFortran
else
    UseRF=.true.
endif

if (present(GRBurnfactor)) then
    GRBF=GRBurnfactor
else
    GRBF=0.5_mrk
endif

if (present(GRBlocSize)) then
    nAdapt3=GRBlocSize
else
    nAdapt3=nAdapt2
endif

nCycles3=Int(nSim/nAdapt3,mrk)

err=0;mess=''
n=size(x,dim=1)

OutDocTemp=""
if (present(OutDoc)) OutDocTemp=trim(OutDoc)//trim(OutDocTemp)

if(allocated(postlist2)) deallocate(postlist2)
allocate(postlist2(n_chain,nCycles3))
if(allocated(phi)) deallocate(phi)
allocate(phi(n,nCycles3*nAdapt3,n_chain))
if(associated(fx)) deallocate(fx)
allocate(fx(n_chain))

!allocate parlist
if (present(showpar)) then
    !allocate parlist
        if(allocated(parlist)) deallocate(parlist)
        allocate(parlist(nCycles1+nCycles2+nCycles3,size(showpar),n_chain))
    !check ranges
    if (present(ranges)) then
        if (size(ranges,dim=2)/=2_mik) then
            err=1;mess="MetroMS_Mix:ranges should include 2 real numbers"
            return
        endif
        if (size(ranges,dim=1)/=size(showpar)) then
            err=1;mess="MetroMS_Mix:ranges size mismatch"
            return
        endif
    endif
endif
    
    
!MCMC 1
if (present(Filename1)) then
    if (present(fAux)) then
        call Adaptive_Metro_OAAT(f=f,x=x,fx=fxtemp,fAux=fAux(:,1),std=std,nAdapt=nAdapt1,nCycles=nCycles1,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				UseRfortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				OutFile=trim(OutDocTemp)//trim(Filename1),outMat=MCMCsamples ,headers=headers,err=err,mess=mess)        
    else
        call Adaptive_Metro_OAAT(f=f,x=x,fx=fxtemp,std=std,nAdapt=nAdapt1,nCycles=nCycles1,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				UseRfortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				OutFile=trim(OutDocTemp)//trim(Filename1),outMat=MCMCsamples ,headers=headers,err=err,mess=mess)    
    endif

else
    if (present(fAux)) then
        call Adaptive_Metro_OAAT(f=f,x=x,fx=fxtemp,fAux=fAux(:,1),std=std,nAdapt=nAdapt1,nCycles=nCycles1,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				UseRfortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				outMat=MCMCsamples ,headers=headers,err=err,mess=mess)       
    else
        call Adaptive_Metro_OAAT(f=f,x=x,fx=fxtemp,std=std,nAdapt=nAdapt1,nCycles=nCycles1,&
				MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
				DownMult=DownMult,UpMult=UpMult,&
				UseRfortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				outMat=MCMCsamples ,headers=headers,err=err,mess=mess)    
    endif

endif				
If(err>0) then
	mess='MetroMS_Mix'//trim(mess);return
endif
If (present(OutMat1)) then
    if (associated(outMat1)) deallocate(outMat1)
    allocate(outMat1(size(MCMCsamples,dim=1),size(MCMCsamples,dim=2)))
    outMat1=MCMCsamples
endif

! MCMC 2
!if(allocated(MCMCsamples)) deallocate(MCMCsamples)
!allocate(MCMCsamples(nAdapt1*nCycles1,n+1))
    
if(allocated(covar)) deallocate(covar)
allocate(covar(n,n))

!j=1
!open(unit=j,file=trim(OutDoc)//trim(Filename1))
!read(j,*)
!do i=1,nAdapt1*nCycles1
!    read(j,*) MCMCsamples(i,:) !to be checked
!enddo
!close(j)

covar(:,:)=getvar(x=MCMCsamples(NINT(BurnFactor*nAdapt1*nCycles1):nAdapt1*nCycles1,1:n),&
		package=2, method="c")
! Choleskize covariance
call choles_dcmp(A=covar(:,:), posDefinite=pos, err=err, message=mess)


if (present(scaleGJ)) then
    scale=scaleGJ
else
    scale=2.4_mrk/sqrt(1._mrk*(n))			
endif

if (present(Filename2)) then
    if (present(fAux)) then
        call Adaptive_Metro(f=f,x=x,fx=fxtemp,fAux=fAux(:,1),covar=covar,scale=scale,&
		        nAdapt=nAdapt2,nCycles=nCycles2,&
		        MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
		        DownMult=DownMult,UpMult=UpMult,&
				UseRFortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				OutFile=trim(OutDocTemp)//trim(Filename2),outMat=MCMCs2temp,headers=headers,err=err,mess=mess)    
    else
        call Adaptive_Metro(f=f,x=x,fx=fxtemp,covar=covar,scale=scale,&
		        nAdapt=nAdapt2,nCycles=nCycles2,&
		        MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
		        DownMult=DownMult,UpMult=UpMult,&
				UseRFortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				OutFile=trim(OutDocTemp)//trim(Filename2),outMat=MCMCs2temp,headers=headers,err=err,mess=mess)    
    endif

else
    if (present(fAux)) then
        call Adaptive_Metro(f=f,x=x,fx=fxtemp,fAux=fAux(:,1),covar=covar,scale=scale,&
		        nAdapt=nAdapt2,nCycles=nCycles2,&
		        MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
		        DownMult=DownMult,UpMult=UpMult,&
				UseRFortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				outMat=MCMCs2temp,headers=headers,err=err,mess=mess)       
    else
        call Adaptive_Metro(f=f,x=x,fx=fxtemp,covar=covar,scale=scale,&
		        nAdapt=nAdapt2,nCycles=nCycles2,&
		        MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
		        DownMult=DownMult,UpMult=UpMult,&
				UseRFortran=UseRF,&
				showPar=showPar,ranges=ranges,&
				outMat=MCMCs2temp,headers=headers,err=err,mess=mess)   
    endif
endif				
If(err>0) then
	mess='MetroMS_Mix'//trim(mess);return
endif
If (present(OutMat2)) then
    if (associated(outMat2)) deallocate(outMat2)
    allocate(outMat2(size(MCMCs2temp,dim=1),size(MCMCs2temp,dim=2)))
    outMat2=MCMCs2temp
endif


! MCMC 3
if (present(fAux)) then
    if(allocated(MCMCsamples2)) deallocate(MCMCsamples2)
    allocate(MCMCsamples2(nAdapt1*nCycles1+nAdapt2*nCycles2,n+1+size(fAux,dim=1)))
else
    if(allocated(MCMCsamples2)) deallocate(MCMCsamples2)
    allocate(MCMCsamples2(nAdapt1*nCycles1+nAdapt2*nCycles2,n+1))
endif
MCMCsamples2(1:(nAdapt1*nCycles1),:)=MCMCsamples(1:(nAdapt1*nCycles1),:)
MCMCsamples2((nAdapt1*nCycles1+1):(nAdapt1*nCycles1+nAdapt2*nCycles2),:)=MCMCs2temp(1:(nAdapt2*nCycles2),:)    

!j=1
!open(unit=j,file=trim(OutDoc)//trim(Filename1))
!read(j,*)
!do i=1,nAdapt2*nCycles2
!    read(j,*) MCMCsamples2(nAdapt1*nCycles1+i,:)
!enddo
!close(j)

if(allocated(unifR)) deallocate(unifR)
allocate(unifR(n_chain-1))
do j=1,(n_chain-1)
    call uniran(unifR(j),real(nAdapt1*nCycles1+1,mrk),real(nAdapt1*nCycles1+nAdapt2*nCycles2,mrk))
enddo

if(allocated(startTemp)) deallocate(startTemp)
allocate(startTemp(n,n_chain))

if(allocated(move)) deallocate(move)
allocate(move(n_chain))
if(allocated(GR3)) deallocate(GR3)
allocate(GR3(n,nCycles3))

! associate outMat
if (present(outMat3)) then
    if (present(fAux)) then
        if(associated(outMat3)) deallocate(outMat3)
        allocate(outMat3(nSim,n+1+size(fAux,dim=1),n_chain))
!        outMat3(1,1:n,:)=x(:,:)
!        outMat3(1,n+1,:)=fx(:)
!        outMat3(1,(n+2):(n+1+size(fAux,dim=1)),:)=fAux(:,:)
    else
        if(associated(outMat3)) deallocate(outMat3)
        allocate(outMat3(nSim,n+1,n_chain))
!        outMat3(1,1:n,:)=x(:,:)
!        outMat3(1,n+1,:)=fx(:)
    endif
endif

! define startTemp and fx
startTemp(:,1)=x
fx(1)=fxtemp
if (present(fAux)) fAux(:,1)=fAux(:,1)
do j=1,(n_chain-1)
    startTemp(:,j+1)=MCMCsamples2(Nint(unifR(j)),1:n)
    fx(j+1)=MCMCsamples2(Nint(unifR(j)),n+1)
    if (present(fAux)) fAux(:,j+1)=MCMCsamples2(Nint(unifR(j)),(n+1+1):(n+1+size(fAux,dim=1)))
enddo

if(present(Filename3)) then
    do j=1,n_chain
        write(nc,'(I1)') j
        open(unit=(j),file=trim(OutDocTemp)//"C"//nc//"_"//trim(Filename3))
    enddo

	if(present(headers)) then ! trust user - no size check at the moment
		k=size(headers)
		do j=1,n_chain
    		write((j),'(<k>A14)') headers
        enddo
	else
		if(allocated(head)) deallocate(head)
		allocate(head(n))
		do i=1,n
			head(i)='par'//trim(number_string(i))
		enddo
		if(present(fAux)) then
			p=size(fAux,dim=1)
			if(allocated(fAuxHead)) deallocate(fAuxHead)
			allocate(fAuxHead(p))
			do i=1,p
				fAuxhead(i)=' fAux'//trim(number_string(i))
			enddo
			do j=1,n_chain
    			write((j),'(<n+1+p>A14)') head(:),'Objectivefunction',fAuxHead(:)
    		enddo
	    	deallocate(fAuxHead)
		else
		    do j=1,n_chain
    			write((j),'(<n+1>A14)') head(:),'Objectivefunction'
	        enddo
		endif
		deallocate(head)
	endif
endif

if (UseRF) then
    ok = Rinit();ok=Reval('X11()')
    again=.true.;compt=0;comptsim=0
    ok=Rput("n_chain",n_chain)
    ok=Rput("nCycles",nCycles3)
    do while(again)
        do i=1,nAdapt3
            do j=1,n_chain
                if (present(fAux)) then
                    call Metro_GaussJump(f=f,&
	                    x=startTemp(:,j),fx=fx(j),fAux=fAux(:,j),&
	                    covar=(scale**2)*covar,move=move(j),err=err,mess=mess)
	            else
	                call Metro_GaussJump(f=f,&
	                    x=startTemp(:,j),fx=fx(j),&
	                    covar=(scale**2)*covar,move=move(j),err=err,mess=mess)
	            endif
	            If(err>0) then
			        mess='MetroMS_Mix:'//trim(mess);return
		        endif
                phi(:,(compt*nAdapt3+i),j)=startTemp(:,j)   
    		            
	            if(present(Filename3)) then        
                    if (present(fAux)) then
                        nfaux=size(faux,dim=1)
	                    write((j),'(<n+1+nfaux>e14.6)') startTemp(:,j),fx(j),fAux(:,j)
	                else
	                    write((j),'(<n+1>e14.6)') startTemp(:,j),fx(j)
	                endif
	            endif
	            
	            if(present(outMat3) .and. again) then
		            outMat3(comptSim+1,1:n,j)=startTemp(:,j)
		            outMat3(comptSim+1,n+1,j)=fx(j)
		            if (present(fAux)) then
		                nfaux=size(faux,dim=1)
		                outMat3(comptSim+1,(n+1+1):(n+1+nfaux),j)=fAux(:,j)
		            endif
		        endif
            enddo
            comptSim=comptSim+1
            if (comptSim>=nSim) again=.false.
        enddo
        postlist2(:,compt+1)=fx
        
        do k=1,n
    	    call get_GRIndice(estimand=phi(k,nint(GRBF*((compt+1)*nAdapt3)):((compt+1)*nAdapt3),:),&
    	                indice=GR3(k,(compt+1)),err=err,mess=mess)
	        If(err>0) then
		        mess='MetroMS_Mix:'//trim(mess);return
	        endif
	    enddo
        
        if (present(showpar)) then
		    do j=1,size(showpar)
		        parlist(compt+1,j,:)=startTemp(showpar(j),:)
		    enddo
	    endif
    	
	    if (present(showpar)) then
            ok=Rput('k',size(showpar))
            ok=reval('layout(matrix(c(rep(1,k),rep(2,k),rep(3:(k+2),rep(2,k))),k*2,2))')
        else
            ok=Reval('layout(c(1,2))')
        endif
        
	    ok=Rput('GR', GR3)
        ok=Rput('compt',compt)
        ok=Rput('ymax',max(1.6_mrk,0.5_mrk+maxval(GR3(:,compt+1))))
        ok=Reval('matplot(GR[1,1:(compt+1)],type="l",ylim=c(0.8,ymax), xlab="Iteration(nCycles) for each parameter", ylab="Gelman-Rubin criterion", main="Gelman-Rubin criterion", cex.lab=2,cex.axis=2,cex.main=3,col=3)')

        do k=2,n
        ok=Rput('k',k)
        ok=Reval('lines(GR[k,1:(compt+1)],col=k+2)')
        enddo
        ok=Reval('abline(h=1,lwd=2)')
        ok=Reval('abline(h=1.2,lwd=2,col=2)')
        !ok=Reval('legend("topright",c("for each parameter"),lty=c(1),pch=c(NA),col=c(1),cex=1.25,lwd=2)')

        ok=Rput('post2', postlist2)
        ok=Rput('compt',compt)
        ok=Reval('yminFB=quantile(post2[,1:(1+compt)],0.02)')
        ok=Reval('ymaxFB=quantile(post2[,1:(1+compt)],0.98)')
        ok=Reval('matplot(post2[1,1:(1+compt)],type="l", ylim=c(yminFB,ymaxFB), xlab="Iteration(nCycles) for each chain", ylab="Log-posterior", main="Log-posterior (From Begin)",cex.lab=2,cex.axis=2,cex.main=3)')
        do j=2,n_chain
        ok=Rput('j',j)
        ok=Reval('lines(post2[j,1:(1+compt)],col=j)')
        enddo

        if (present(showpar)) then
            ok=Rput('pars', parlist)
            ok=Rput('compt',compt)
            ok=Rput('showpar',showpar)
            do k=1,size(showpar)
                if (present(ranges)) then
                    ok=Rput('rmin',ranges(k,1))
                    ok=Rput('rmax',ranges(k,2))
                    ok=Rput('k',k)
                    ok=Reval('matplot(pars[1:(compt+1),k,1],type="l", xlab="Iteration(nCycles) for each chain", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                else
                    ok=Rput('k',k)
                    ok=Reval('rmin=quantile(pars[1:(compt+1),k,],0.02)')
                    ok=Reval('rmax=quantile(pars[1:(compt+1),k,],0.98)')
                    ok=Reval('matplot(pars[1:(compt+1),k,1],type="l", xlab="Iteration(nCycles) for each chain", ylab=paste("par ",showpar[k],sep=""), main="MCMC chain",ylim=c(rmin,rmax),cex.lab=2,cex.axis=2,cex.main=3)')
                endif
                do j=2,n_chain
                    ok=Rput('j',j)
                    ok=Reval('lines(pars[1:(compt+1),k,j],col=j)')
                enddo
            enddo
        endif
        compt=compt+1     
    enddo
    
else
    again=.true.;compt=0;comptsim=0
    do while(again)
        do i=1,nAdapt3
            do j=1,n_chain
                if (present(fAux)) then
                    call Metro_GaussJump(f=f,&
	                    x=startTemp(:,j),fx=fx(j),fAux=fAux(:,j),&
	                    covar=(scale**2)*covar,move=move(j),err=err,mess=mess)
	            else
	                call Metro_GaussJump(f=f,&
	                    x=startTemp(:,j),fx=fx(j),&
	                    covar=(scale**2)*covar,move=move(j),err=err,mess=mess)
	            endif
	            If(err>0) then
			        mess='MetroMS_Mix:'//trim(mess);return
		        endif
                phi(:,(compt*nAdapt3+i),j)=startTemp(:,j)    
    		            
	            if(present(Filename3)) then        
                    if (present(fAux)) then
                        nfaux=size(faux,dim=1)
	                    write((j),'(<n+1+nfaux>e14.6)') startTemp(:,j),fx(j),fAux(:,j)
	                else
	                    write((j),'(<n+1>e14.6)') startTemp(:,j),fx(j)
	                endif
	            endif
	            
	            if(present(outMat3) .and. again) then
		            outMat3(comptSim+1,1:n,j)=startTemp(:,j)
		            outMat3(comptSim+1,n+1,j)=fx(j)
		            if (present(fAux)) then
		                nfaux=size(faux,dim=1)
		                outMat3(comptSim+1,(n+1+1):(n+1+nfaux),j)=fAux(:,j)
		            endif
		        endif
            enddo
            comptSim=comptSim+1
            if (comptSim>=nSim) again=.false.
        enddo
        
        do k=1,n
    	    call get_GRIndice(estimand=phi(k,nint(GRBF*((compt+1)*nAdapt3)):((compt+1)*nAdapt3),:),&
    	                indice=GR3(k,(compt+1)),err=err,mess=mess)
	        If(err>0) then
		        mess='MetroMS_Mix:'//trim(mess);return
	        endif
	    enddo
        compt=compt+1
     
    enddo

endif

if(present(GRIndex)) then
    if(associated(GRIndex)) deallocate(GRIndex)
    allocate(GRIndex(n,nCycles3))
    GRIndex=GR3
endif

if(present(Filename3)) then
    do j=1,n_chain
        close(j)
    enddo
endif

deallocate(move)
deallocate(covar)
deallocate(postlist2)
if (present(showPar)) deallocate(parlist)
deallocate(phi)
deallocate(GR3)
deallocate(startTemp)
deallocate(unifR)
deallocate(MCMCsamples)
deallocate(MCMCsamples2)
deallocate(MCMCs2temp)
end subroutine MetroMS_Mix2



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!#########
!#Private#
!#########

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_GRIndice(estimand,indice,err,mess)
!^**********************************************************************
!^* Purpose: Computation of Gelman-Rubin criterion
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:  /  /2012
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.estimand(length of chain, number of parallel sequences)
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* OUT
!^*     1.
!^*     2.
!^*     3.
!^*     4.
!^*     5.
!^* INOUT
!^*     1.
!^*     2.
!^**********************************************************************
real(mrk),intent(in)::estimand(:,:)
real(mrk),intent(out)::indice
integer(mik), intent(out)::err
character(100),intent(out)::mess
!local
real(mrk)::B,W,meantot
real(mrk),allocatable::meanInt(:),s2(:),temp(:)
integer(mik)::i,j,n,m


err=0;mess=''
n=size(estimand,dim=1)
m=size(estimand,dim=2)

if (m==1_mik) then
    err=1;mess="get_GRIndice:Only exite one chain"
    return
endif

if(allocated(meanInt)) deallocate(meanInt)
allocate(meanInt(m))
if(allocated(s2)) deallocate(s2)
allocate(s2(m))



do j=1,m
    meanInt(j)=sum(estimand(:,j))/n
enddo  
  
if(allocated(temp)) deallocate(temp)
allocate(temp(n))
do j=1,m
    temp=meanInt(j)-estimand(:,j)
    s2(j)=DOT_PRODUCT(temp,temp)/(n-1._mrk)
enddo

meantot=sum(meanInt)/m

if(allocated(temp)) deallocate(temp)
allocate(temp(m))
temp=meanInt-meantot
B=n/(m-1._mrk)*DOT_PRODUCT(temp,temp)
W=sum(s2)/m

indice=sqrt(((n-1._mrk)/n*W+1._mrk/n*B)/W)

deallocate(meanInt)
deallocate(s2)
deallocate(temp)

end subroutine get_GRIndice

!subroutine ComputeGR(nrow, values, criterion)
!! Computation of Gelman-Rubin criterion
!
!integer, intent(in)::nrow
!REAL(RK), parameter:: GRSIZE=0.5_RK
!real(rk), intent(in):: values(nrow, gibbs%nParallel, gibbs%dim)
!real(rk), intent(out):: criterion(gibbs%dim)
!real(rk):: X(int(nrow*GRSIZE), gibbs%nParallel), 
!moy(gibbs%nParallel),s2(gibbs%nParallel), bigmoy, B, W
!integer i, j, k, m, n
!
!n=int(nrow*GRSIZE)
!m=gibbs%nParallel
!
!do k=1, gibbs%dim
!    X=values(int((1+nrow)*GRSIZE):nrow, :, k)
!    do j=1, m
!        moy(j)=(1.0_rk/(1.0_rk*n))*sum(X(:, j))
!        s2(j)=(1.0_rk/(1.0_rk*(n-1.0_rk)))*sum((X(:, j)-moy(j))**2)
!    enddo
!    bigmoy=(1.0_rk/(1.0_rk*m))*sum(moy)
!    B=((1.0_rk*n)/(1.0_rk*(m-1.0_rk)))*sum((moy-bigmoy)**2)
!    W=(1.0_rk/(1.0_rk*m))*sum(s2)
!    criterion(k)=sqrt((((n-1.0_rk)/(n*1.0_rk))*W+((1.0_rk/(1.0_rk*n)))*B)/W)
!enddo
!
!end subroutine ComputeGR
end module MCMCStrategy_tools
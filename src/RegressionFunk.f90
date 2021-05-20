module RegressionFunk

!~**********************************************************************
!~* Purpose: Libaray of Regression functions
!~**********************************************************************
!~* Programmer: Xun SUN, Ben Renard - Cemagref
!~**********************************************************************
!~* Last modified: 16/09/2011
!~**********************************************************************
!~* Comments: Naming is given by the LINK function (not the inverse link!) 
!~* eg in model log(y)=reg(x|teta)+eps, one should use: 
!~* call InverseLinkFunk(IvsLinkFunk='log',etc. 
!~* call LinkFunk(LinkFunk='log',etc. 
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: Fix Link naming
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. ApplyRegressionFunk(RegressionFunk,covariate,regressionPar,out,err,mess)
!~*		2. ApplyRegionalRegFunk(RegressionFunk,t_cov,s_cov,st_cov,regressionPar,out,err,mess)
!~*		3. InverseLinkFunk(IvsLinkFunk,x,fx,feas,err,mess)
!~*		4. GetRegParNumber(RegressionFunk,covariate,npar,err,mess)
!~*		5. LinkFunk(LinkFunk,x,fx,feas,err,mess)
!~*		6. GetRegCovNumber(RegressionFunk,ncov,err,mess)
!~*		7. GetRegioanlRegCovNumber(RegressionFunk,ncov,err,mess)
!~**********************************************************************
use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: ApplyRegressionFunk, ApplyRegionalRegFunk ,InverseLinkFunk, LinkFunk, GetRegParNumber,&
             GetRegCovNumber, GetRegioanlRegCovNumber !,GetSpatParNb,GetSpatCovNb,ApplySpatRegFunc

Character(100), parameter, PUBLIC:: &
            Linear='Linear',Linaire='Linear',&
            Cosine='Cosinus',Cosinus='Cosinus', &
            id='Identity',identity='Identity',&
            Power_reg="Power",&
            XtraFlo_reg="XtraFlo"

contains
!!!!!!!!!!!!!!!!!!!!!!
subroutine ApplyRegressionFunk(RegressionFunk,covariate,Xtra,regressionPar,out,err,mess)
!^**********************************************************************
!^* Purpose: LOCAL Regression Function Libaray
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified: 30/01/2012
!^**********************************************************************
!^* Comments: Modified by BR to include optional "Xtra", used to pass 
!^* optional arguments that might be needed to run the regression. It 
!^* was required for EXTRAFLO, which uses region-specific regressions 
!^* and therefore needs to know the number of regions. Should be 
!^* back-compatible since argument is optional. In XtraFlo regression, 
!^* also used to pass regression type (1=additive, 2=power).
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.RegressionFunk: Names of Regression functions: 
!^*		    Linear: f(Xs|mus)=mu0+x1*mu1+x2*mu2+..... (nb X= nb mu)
!^*		    Cosinus: f(t|mu0,mu1)=mu0*cos(mu1*t)     (nb X=1, nb mu=2)
!^*		2.Covariate
!^*		3.RegressionPar: mu
!^*		4.[Xtra]: optional arguments
!^*	OUT
!^*	    4.out: Function results
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
character(*),intent(in)::RegressionFunk
real(mrk),intent(in),optional::covariate(:)
real(mrk),intent(in)::regressionPar(:)
real(mrk),intent(in), optional::Xtra(:)
real(mrk),intent(out)::out
integer(mik), intent(out)::err
character(*),intent(out)::mess
err=0;mess='';out=undefRN

select case (trim(RegressionFunk))
case('Linear')
    if(size(covariate)/=size(regressionPar)) then
        err=1;mess="ApplyRegressionFunk: [linear] Covariate and regressionPar size mismatch";return
    endif
    if(size(covariate)/=1 .and. .not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: Covariate not present";return
    endif
	out=DOT_PRODUCT (covariate,regressionPar)
case('LinearRelative')
    if(size(covariate)/=(size(regressionPar)-1)) then
        err=1;mess="ApplyRegressionFunk: [LinearRelative] Covariate and regressionPar size mismatch";return
    endif
    if(size(covariate)/=1 .and. .not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: Covariate not present";return
    endif
    if(size(regressionPar)<2) then
        err=1;mess="ApplyRegressionFunk: [LinearRelative] needs at least 2 regression parameters";return
    endif
	out=regressionpar(1)*(1._mrk+DOT_PRODUCT (covariate,regressionPar(2:)))
case('LinearWithIntercept')
    if(size(covariate)/=(size(regressionPar)-1)) then
        err=1;mess="ApplyRegressionFunk: [LinearWithIntercept] Covariate and regressionPar size mismatch";return
    endif
    if(size(covariate)/=1 .and. .not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: Covariate not present";return
    endif
    if(size(regressionPar)<2) then
        err=1;mess="ApplyRegressionFunk: [LinearWithIntercept] needs at least 2 regression parameters";return
    endif
	out=regressionpar(1)+DOT_PRODUCT (covariate,regressionPar(2:))
case(Power_reg)
    if(.not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: [power]: Covariate not present";return
    endif
    if(size(covariate)/=size(regressionPar)) then
        err=1;mess="ApplyRegressionFunk: [power]: Covariate and regressionPar size mismatch";return
    endif
    if(any(covariate<=0._mrk)) then
        err=1;mess="ApplyRegressionFunk: [power]: negative covariates not allowed with power regression";return
    endif
	out=exp(DOT_PRODUCT (log(covariate),regressionPar))
case(XtraFlo_reg)
    if(.not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: [XtraFlo]: Covariate not present";return
    endif
    if (.not.present(Xtra)) then
        err=1;mess="ApplyRegressionFunk: [XtraFlo]: Xtra not optional for this regression";return
    endif
    if (size(Xtra)/=2) then
        err=1;mess="ApplyRegressionFunk: [XtraFlo]:  optional Xtra should have size 2";return
    endif
    call ApplyExtraFloReg(covariate=covariate,regressionPar=regressionPar,&
              nregion=nint(Xtra(1)),RegType=nint(Xtra(2)),out=out,err=err,mess=mess)
    if(err>0) then
        mess="ApplyRegressionFunk:"//trim(mess);return
    endif
case('Linear_L1')
    if(.not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: Covariate not present";return
    endif
    if( (size(covariate)/=1) .or. (size(regressionPar)/=1) ) then
        err=1;mess="ApplyRegressionFunk:[Cosinus] Covariate or regressionPar size mismatch";return
    endif
    out=regressionPar(1)*covariate(1)
case('Linear_L2')
    if(.not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: Covariate not present";return
    endif
    if( (size(covariate)/=2) .or. (size(regressionPar)/=2) ) then
        err=1;mess="ApplyRegressionFunk:[Cosinus] Covariate or regressionPar size mismatch";return
    endif
    out=regressionPar(1)*covariate(1)+regressionPar(2)*covariate(2)
case('Linear_L3')
    if(.not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: Covariate not present";return
    endif
    if( (size(covariate)/=3) .or. (size(regressionPar)/=3) ) then
        err=1;mess="ApplyRegressionFunk:[Cosinus] Covariate or regressionPar size mismatch";return
    endif
    out=regressionPar(1)*covariate(1)+regressionPar(2)*covariate(2)+regressionPar(3)*covariate(3)
case("Cosinus")
    if(.not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: Covariate not present";return
    endif
    if( (size(covariate)/=1) .or. (size(regressionPar)/=2) ) then
        err=1;mess="ApplyRegressionFunk:[Cosinus] Covariate or regressionPar size mismatch";return
    endif
	out=regressionPar(1)*cos(covariate(1)*regressionPar(2))
case("Exp")
    if(.not.present(covariate))then
        err=1;mess="ApplyRegressionFunk: Covariate not present";return
    endif
    if( (size(covariate)/=1) .or. (size(regressionPar)/=2) ) then
        err=1;mess="ApplyRegressionFunk:[Cosinus] Covariate or regressionPar size mismatch";return
    endif
    out=regressionPar(1)*exp(covariate(1)*regressionPar(2))
case("Identity")
    if(size(regressionPar)/=1) then
        err=1;mess="ApplyRegressionFunk:pb";return
    endif
    out= regressionPar(1)

case default
	mess="ApplyRegressionFunk:UnknownRegressionFunk";err=1;return
end select

end subroutine ApplyRegressionFunk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ApplyRegionalRegFunk(RegressionFunk,t_cov,s_cov,st_cov,regressionPar,out,err,mess)
!^**********************************************************************
!^* Purpose: Regional Regression Function Libaray
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:16/09/2011
!^**********************************************************************
!^* Comments: For a fixed time and a fixed site
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.RegressionFunk: Names of Regression functions: 
!^*     2.t_cov: temporal covariates. (e.g. time, NAO)
!^*     3.s_cov: spatial covariates. (e.g. latitude)
!^*     4.st_cov: spatial-temporal covariates. (e.g. type de temps)
!^*     5.RegressionPar: mu
!^* OUT
!^*	    4.out: Function results
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************

character(*),intent(in)::RegressionFunk
real(mrk),intent(in)::t_cov(:)
!real(mrk),pointer::t_cov(:)
real(mrk),intent(in)::s_cov(:)
!real(mrk),pointer::s_cov(:)
real(mrk),intent(in)::st_cov(:)
!real(mrk),pointer::st_cov(:)
!integer(mik),intent(in)::site,time
real(mrk),intent(in)::regressionPar(:)
real(mrk),intent(out)::out
integer(mik), intent(out)::err
character(*),intent(out)::mess
err=0;mess='';out=undefRN

select case (trim(RegressionFunk))
case("Identity")
    if(size(regressionPar)/=1) then
        err=1;mess="ApplyRegionalRegFunk:parameter number incorrect";return
    endif
    out=regressionPar(1)
case("Regional_1")
   ! if( .not. associated(t_cov)) then
    !    err=1; mess="ApplyRegionalRegFunk:Does not exist t_cov";return
    !else
        if(  (size(regressionPar)/=2) .or. (size(t_cov)/=1)) then
            err=1;mess="ApplyRegionalRegFunk:[Regional_1] Covariate or regressionPar size mismatch";return
        endif
    !endif
    out= regressionPar(1)+ regressionPar(2)*t_cov(1)
case("Regional_2")
    !if( .not. associated(s_cov)) then
    !    err=1; mess="ApplyRegionalRegFunk:Does not exist s_cov";return
    !else
        if(  (size(regressionPar)/=2) .or. (size(s_cov)/=1)) then
            err=1;mess="ApplyRegionalRegFunk:[Regional_2] Covariate or regressionPar size mismatch";return
        endif
    !endif
    out= regressionPar(1)+ regressionPar(2)*s_cov(1) 
case("Regional_3")
    !if( (.not. associated(t_cov)) .or. (.not. associated(s_cov))) then
    !    err=1; mess="ApplyRegionalRegFunk:Does not exist t_cov or s_cov";return
    !else
        if(  (size(regressionPar)/=3) .or. (size(s_cov)/=1) .or. (size(t_cov)/=1)) then
            err=1;mess="ApplyRegionalRegFunk:[Regional_3] Covariate or regressionPar size mismatch";return
        endif
    !endif
    out= regressionPar(1)+ regressionPar(2)*s_cov(1)+regressionPar(3)*t_cov(1)   
case("Regional_4")
    !if( (.not. associated(t_cov)) .or. (.not. associated(s_cov))) then
    !    err=1; mess="ApplyRegionalRegFunk:Does not exist t_cov or s_cov";return
    !else
        if(  (size(regressionPar)/=3) .or. (size(s_cov)/=1) .or. (size(t_cov)/=1)) then
            err=1;mess="ApplyRegionalRegFunk:[Regional_4] Covariate or regressionPar size mismatch";return
        endif
    !endif
    out= regressionPar(1)+ (regressionPar(2)+regressionPar(3)*s_cov(1))*t_cov(1)
case("Regional_5")
    !if( .not. associated(st_cov)) then
    !    err=1; mess="ApplyRegionalRegFunk:Does not exist st_cov";return
    !else
        if(  (size(regressionPar)/=2) .or. (size(st_cov)/=1)) then
            err=1;mess="ApplyRegionalRegFunk:[Regional_5] Covariate or regressionPar size mismatch";return
        endif
    !endif
    out= regressionPar(1)+ regressionPar(2)*st_cov(1)
case("Regional_6")
        if(  (size(regressionPar)/=2) .or. (size(t_cov)/=1)) then
            err=1;mess="ApplyRegionalRegFunk:[Regional_1] Covariate or regressionPar size mismatch";return
        endif
        out=regressionPar(1)*(1._mrk+ regressionPar(2)*t_cov(1))
case("Regional_7")
        if(  (size(regressionPar)/=3) .or. (size(t_cov)/=1)) then
            err=1;mess="ApplyRegionalRegFunk:[Regional_1] Covariate or regressionPar size mismatch";return
        endif
        out=regressionPar(1)*(1._mrk+ regressionPar(2)+regressionPar(3)*t_cov(1))
case("Cosine_1")
        if(  (size(regressionPar)/=2) .or. (size(t_cov)/=1)) then
            err=1;mess="ApplyRegionalRegFunk:[Regional_1] Covariate or regressionPar size mismatch";return
        endif
        out=regressionPar(1)*cos(regressionPar(2)*t_cov(1))
case("ChangePoint_1")
        if(  (size(regressionPar)/=3) .or. (size(t_cov)/=2)) then
            err=1;mess="ApplyRegionalRegFunk:[ChangePoint_1] t_cov or regressionPar size mismatch";return
        endif
        if (t_cov(2)<0) then
            out= regressionPar(1)*t_cov(1)+ regressionPar(2)*t_cov(2)
        else
            out= regressionPar(1)*t_cov(1)+ regressionPar(3)*t_cov(2)
        endif
case("ChangePoint_1b")
        if(  (size(regressionPar)/=3) .or. (size(t_cov)/=2)) then
            err=1;mess="ApplyRegionalRegFunk:[ChangePoint_1b] t_cov or regressionPar size mismatch";return
        endif
        if (t_cov(2)<0) then
            out= regressionPar(1)*(t_cov(1)+ regressionPar(2)*t_cov(2))
        else
            out= regressionPar(1)*(t_cov(1)+ regressionPar(3)*t_cov(2))
        endif
case("ChangePoint_2")
        if(  (size(regressionPar)/=4) .or. (size(t_cov)/=2)) then
            err=1;mess="ApplyRegionalRegFunk:[ChangePoint_1] t_cov or regressionPar size mismatch";return
        endif
        if (t_cov(2)<regressionPar(4)) then
            out= regressionPar(1)*t_cov(1)+ regressionPar(2)*t_cov(2)
        else
            out= regressionPar(1)*t_cov(1)+ regressionPar(3)*t_cov(2)
        endif
case("ChangePoint_3")
        if(  (size(regressionPar)/=5) .or. (size(t_cov)/=3)) then
            err=1;mess="ApplyRegionalRegFunk:[ChangePoint_1] t_cov or regressionPar size mismatch";return
        endif
        if (t_cov(2)<0) then
            if (t_cov(3)<0) then
                out= regressionPar(1)*t_cov(1)+ regressionPar(2)*t_cov(2)+regressionPar(4)*t_cov(3)
            else
                out= regressionPar(1)*t_cov(1)+ regressionPar(2)*t_cov(2)+regressionPar(5)*t_cov(3)
            endif
        else
            if (t_cov(3)<0) then
                out= regressionPar(1)*t_cov(1)+ regressionPar(3)*t_cov(2)+regressionPar(4)*t_cov(3)
            else
                out= regressionPar(1)*t_cov(1)+ regressionPar(3)*t_cov(2)+regressionPar(5)*t_cov(3)
            endif
        endif
case("ChangePoint_4")
        if(  (size(regressionPar)/=5) .or. (size(t_cov)/=3)) then
            err=1;mess="ApplyRegionalRegFunk:[ChangePoint_1] t_cov or regressionPar size mismatch";return
        endif
        if (t_cov(2)<0) then
            if (t_cov(3)<0) then
                out= regressionPar(1)*t_cov(1)+ regressionPar(2)*t_cov(2)
            else
                out= regressionPar(1)*t_cov(1)+ regressionPar(3)*t_cov(2)
            endif
        else
            if (t_cov(3)<0) then
                out= regressionPar(1)*t_cov(1)+ regressionPar(4)*t_cov(2)
            else
                out= regressionPar(1)*t_cov(1)+ regressionPar(5)*t_cov(2)
            endif
        endif
case("ChangePoint_5")
        if(  (size(regressionPar)/=2) .or. (size(t_cov)/=2)) then
            err=1;mess="ApplyRegionalRegFunk:[ChangePoint_1] t_cov or regressionPar size mismatch";return
        endif
        if (t_cov(2)<0) then
            out= regressionPar(1)*t_cov(1)
        else
            out= regressionPar(1)*t_cov(1)+ regressionPar(2)*t_cov(2)
        endif
case default
	mess="ApplyRegionalRegFunk:UnknownRegressionFunk";err=1;return
end select

end subroutine ApplyRegionalRegFunk


pure subroutine GetRegParNumber(RegressionFunk,covariate,Xtra,npar,err,mess)

!^**********************************************************************
!^* Purpose: Get the number of parameters of a regression function
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 30/01/2012
!^**********************************************************************
!^* Comments: new optional Xtra, see sub ApplyRegressionFunk
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.RegressionFunk
!^*		2.covariate vector
!^*		3.[Xtra]
!^* OUT
!^*		1.npar
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*),intent(in)::RegressionFunk
real(mrk),optional, intent(in)::covariate(:), Xtra(:)
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess

err=0;mess='';npar=undefIN

select case (trim(RegressionFunk))
case('Linear',Power_reg)
    if (.not. present(covariate)) then
        err=1;mess="GetRegParNumber: Covariate not present";return
    endif
    npar=size(covariate)
case('LinearRelative','LinearWithIntercept')
    if (.not. present(covariate)) then
        err=1;mess="GetRegParNumber: Covariate not present";return
    endif
    npar=size(covariate)+1
case(XtraFlo_reg)
    if (.not. present(covariate)) then
        err=1;mess="GetRegParNumber: Covariate not present";return
    endif
    if (.not.present(Xtra)) then
        err=1;mess="GetRegParNumber: [XtraFlo]: Xtra not optional for this regression";return
    endif
    if (size(Xtra)/=2) then
        err=1;mess="GetRegParNumber: [XtraFlo]:  optional Xtra should have size 2";return
    endif
    npar=nint(Xtra(1))*(size(covariate))
case('Linear_L1')
    npar=1
case('Linear_L2')
    npar=2
case('Linear_L3')
    npar=3
case("Cosinus")
    npar=2
case("Exp")
    npar=2
case("Identity")
    npar=1
case("Regional_1")
    npar=2
case("Regional_2")
    npar=2
case("Regional_3")
    npar=2
case("Regional_4")
    npar=3
case("Regional_5")
    npar=2
case("Regional_6")
    npar=2
case("Regional_7")
    npar=3
case("Cosine_1")
    npar=2
case("ChangePoint_1")
    npar=3
case("ChangePoint_2")
    npar=4
case("ChangePoint_3")
    npar=5
case("ChangePoint_4")
    npar=5
case("ChangePoint_5")
    npar=2
case default
	mess="GetRegParNumber:UnknownRegressionFunk";err=1;return
end select

end subroutine GetRegParNumber


pure subroutine GetRegCovNumber(RegressionFunk,covariate,ncov,err,mess)

!^**********************************************************************
!^* Purpose: Get the number of covariates of a LOCAL regression function 
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref Lyon
!^**********************************************************************
!^* Last modified:12/09/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.RegressionFunk
!^*		2.
!^* OUT
!^*		1.npar
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*),intent(in)::RegressionFunk
real(mrk),optional,intent(in)::covariate(:)
integer(mik), intent(out)::ncov,err
character(*),intent(out)::mess

err=0;mess='';ncov=undefIN

select case (trim(RegressionFunk))
case('Linear',Power_reg,XtraFlo_reg, 'LinearRelative')
    if (.not. present(covariate)) then
        err=1;mess="GetRegCovNumber: Covariate not present";return
    endif
    ncov=size(covariate)
case('Linear_L1')
    ncov=1
case('Linear_L2')
    ncov=2
case('Linear_L3')
    ncov=3
case("Cosinus")
    ncov=1
case("Exp")
    ncov=1
case("Identity")
    ncov=0
case default
	mess="GetRegCovNumber:UnknownRegressionFunk";err=1;return
end select

end subroutine GetRegCovNumber



pure subroutine GetRegioanlRegCovNumber(RegressionFunk,ncov,err,mess)

!^**********************************************************************
!^* Purpose: Get the number of covariates of a regional regression function
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref Lyon
!^**********************************************************************
!^* Last modified:16/09/2011
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.RegressionFunk
!^*		2.
!^* OUT
!^*		1.npar
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*),intent(in)::RegressionFunk
integer(mik), intent(out)::ncov(3),err
character(*),intent(out)::mess

err=0;mess='';ncov=undefIN

select case (trim(RegressionFunk))
case("Identity")
    ncov=(/0,0,0/)
case("Regional_1")
    ncov=(/1,0,0/)
case("Regional_2")
    ncov=(/0,1,0/)
case("Regional_3")
    ncov=(/1,1,0/)
case("Regional_4")
    ncov=(/1,1,0/)
case("Regional_5")
    ncov=(/0,0,1/)
case("Regional_6")
    ncov=(/1,0,0/)
case("Regional_7")
    ncov=(/1,1,0/)
case("Cosine_1")
    ncov=(/1,1,0/)
case("ChangePoint_1")
    ncov=(/2,0,0/)
case("ChangePoint_2")
    ncov=(/2,0,0/)
case("ChangePoint_3")
    ncov=(/3,0,0/)
case("ChangePoint_4")
    ncov=(/3,0,0/)
case("ChangePoint_5")
    ncov=(/2,0,0/)
case default
	mess="GetRegionalRegCovNumber:UnknownRegressionFunk";err=1;return
end select

end subroutine GetRegioanlRegCovNumber



!!!!!!!!!!!!!!!!!!!!!!
subroutine InverseLinkFunk(IvsLinkFunk,x,fx,feas,err,mess)
!^**********************************************************************
!^* Purpose: Inverse Link Function Library
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified: 06/11/2011
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.InverseLinkFunk: Names of Regression functions: 

!^* OUT
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
character(*),intent(in)::IvsLinkFunk
real(mrk),intent(in)::x
real(mrk),intent(out)::fx
logical,intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess

err=0;mess='';fx=undefRN;feas=.true.

select case (trim(IvsLinkFunk))
case('Identity')
    fx=x
case('Log')
    fx=exp(x)
case('Inverse')
    if (x==0.) then
        feas=.False. ;mess='IvsLinkFunk: denominator=0';return;
    endif
    fx=1./x
case('logit')
    fx=exp(x)/(exp(x)+1.)
case default
	mess="ApplyRegressionFunk:UnknownLinkFunk";err=1;return
end select

end subroutine InverseLinkFunk

!!!!!!!!!!!!!!!!!!!!!!
subroutine LinkFunk(LFunk,x,fx,feas,err,mess)
!^**********************************************************************
!^* Purpose:  Link Function Library
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref
!^**********************************************************************
!^* Last modified: 17/05/2011
!^**********************************************************************
!^* Comments: 
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.LFunk: Names of Link function
!^*		2.x 
!^*		3.fx
!^* OUT
!^*		1.feas
!^*		1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		2.mess, error message
!^**********************************************************************
character(*),intent(in)::LFunk
real(mrk),intent(in)::x
real(mrk),intent(out)::fx
logical,intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess

err=0;mess='';fx=undefRN;feas=.true.

select case (trim(LFunk))
case('Identity')
    fx=x
case('Log')
    fx=log(x)
case('Inverse')
    if (x==0.) then
        feas=.False. ;mess='LinkFunk: denominator=0';return;
    endif
    fx=1./x
case('logit')
    if (x<=0._mrk .or. x>=1._mrk) then
        feas=.False. ;mess='LinkFunk: x out of [0;1]';return;
    endif
    fx=log(x)/(1._mrk-x)
case default
	mess="LinkFunk:UnknownLinkFunk";err=1;return
end select

end subroutine LinkFunk

!!!!!!!!!!!
! Private !
!!!!!!!!!!!

pure subroutine ApplyExtraFloReg(covariate,regressionPar,nregion,RegType,out,err,mess)

!^**********************************************************************
!^* Purpose: Extraflo regression, ie region-specific power regression
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 30/01/2012
!^**********************************************************************
!^* Comments: Regtype = 1 for additive, 2 for power
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1.covariate
!^*		2.regressionPar
!^*		3.nregion
!^* OUT
!^*		1.out
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
real(mrk), intent(in)::covariate(:), regressionPar(:)
integer(mik), intent(in)::nregion,regType
real(mrk), intent(out)::out
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::region,npar,ncov
real(mrk), allocatable::par(:)

err=0;mess='';out=undefRN
!Size check
npar=size(regressionPar)
ncov=size(covariate)
if(npar/=(nregion*(ncov))) then
    err=1;mess="ApplyExtraFloReg:FATAL:Size mismatch[npar,ncov,nregion]";return
endif
! Assume region is the first covariate
region=nint(covariate(1),mik)
! get parameters for this region
if(allocated(par)) deallocate(par)
allocate(par(ncov))
par=regressionPar(region::nregion)
! Apply regression
select case(regtype)
case(2)
    if(any(covariate<=0._mrk)) then
        err=1;mess="ApplyExtraFloReg:FATAL:Negative covariates not allowed for power regression";return
    endif
    out=par(1)*exp(dot_product(log(covariate(2:)),par(2:)))
case(1)
    out=par(1)+dot_product(covariate(2:),par(2:))
case(0)    
    ! WARNING::BOTCH:: quick-and-[very]dirty solution to easily implement a constant-per-region model within EXTRAFLO
    out=par(1)
case default
    err=1;mess="ApplyExtraFloReg:FATAL:Unknown Regression type";return
end select

end subroutine ApplyExtraFloReg


end module RegressionFunk
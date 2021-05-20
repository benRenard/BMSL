program TestExpressionParser

use interpreter ! Brazilian parser
use fparser,only: initf, parsef, evalf, EvalErrType, EvalErrMsg ! German parser

implicit none 

integer,parameter:: Nvar=2,Nsim=1
integer::i
character(250):: formula, varname(Nvar),flag
real(realkind)::res,varval(Nvar),t0,t1

formula='log(x/y)'
varname=(/'x','y'/)
varval=(/1.,0.2/)

write(*,*) '***************************************'
write(*,*) '*****************BRAZIL****************'
write(*,*) '***************************************'

! init / evaluate / destroy
CALL CPU_TIME (t0)
call init (formula, varname, flag)
do i=1,Nsim
    res=evaluate(varval)
enddo
call destroyfunc()
CALL CPU_TIME (t1)
write(*,*) trim(formula)//' evaluated at ', varval
write(*,*) 'Result: ', res
write(*,*) 'CPU time for evaluation: ', t1-t0

write(*,*) 'DONE!'

write(*,*) '***************************************'
write(*,*) '*****************GERMANY***************'
write(*,*) '***************************************'
formula='laog(x/y)'
CALL initf (1)
write(*,*) 'INIT Error: ', EvalErrType, EvalErrMsg()
read(*,*)
CALL CPU_TIME (t0)
call parsef(1,formula,varname)
write(*,*) 'PARSE Error: ', EvalErrType, EvalErrMsg()
read(*,*)
do i=1,Nsim
    res=evalf(1,varval)
enddo
CALL CPU_TIME (t1)
write(*,*) trim(formula)//' evaluated at ', varval
write(*,*) 'EVAL Error: ', EvalErrType, EvalErrMsg()
write(*,*) 'Result: ', res
write(*,*) 'CPU time for evaluation: ', t1-t0
read(*,*)


end program TestExpressionParser
module InOutTools
!~**********************************************************************
!~* Purpose: Some tools to read and write a vector or a matrix in Fortran
!~**********************************************************************
!~* Programmer: Xun SUN, Cemagref
!~**********************************************************************
!~* Last modified:10/12/2010
!~**********************************************************************
!~* Comments: Sometimes there are some pbs with the non-default separator 
!~*        while outputing files.
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. mat_input, read a matrix
!~*		2. mat_output, write a vector or a matrix into a file
!~*		3. 
!~*		4. 
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public ::mat_input, mat_output

contains

!*********
subroutine mat_input(Filename, row, column, header, sep, tab, err, mess)
!optional: row=row number, header= .false., sep=';'

!~**********************************************************************
!~* Purpose: Read a matrix data
!~**********************************************************************
!~* Programmer: Xun SUN, Cemagref
!~**********************************************************************
!~* Last modified: 10/12/2010
!~**********************************************************************
!~* Comments:Each line can not exceed 1000 characters (ifnot modify longeur)
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* In
!~*		1.Filename
!~*		2.row: number of row contained in the file (optional)
!~*		3.column: number of column in the file
!~*		4.header: Whether do the data contain titles (optional, default= .false.)
!~*		5.sep: separator symbol (optional, default= ' ')
!~*	Out
!~*		1.tab: stock the file values (real) in a matrix 
!^*		2.err, error code
!^*		3.mess, error message
!~**********************************************************************
character(*), intent(in)::Filename
integer(mik), intent(in), optional::row
integer(mik)::k ! k plays the role of row,
character(len=1), intent(in), optional:: sep
character(len=1):: separator
integer(mik), intent(in)::  column 
logical, intent(in), optional:: header 
logical:: head
integer(mik), intent(out)::err
real(mrk), pointer::tab(:,:) !output matrice of real numbers
character(*), intent(out)::mess
!locals
integer(mik):: numfile, i, l, n, c, errcode  !i and l are parameters for the loop
integer(mik),parameter:: longeur=2000    ! maximum longeur/2 character of each line
character(longeur)::d
logical::nextValue


err=0;mess=''

if (.not. present(sep)) then
    separator=' '
else
    separator=sep
end if

if (.not. present(header)) then
    head=.false.
else
    head=header
end if

numfile=1
k=0


if (.not. present(row)) then !compute the number of row if that is not given
    if (head) k=k-1
    open(unit=numfile,file=trim(Filename),form='formatted',status='old',iostat=err)
    if(err>0) then
        mess="InOutTools:Can't find the file:" // trim(Filename)
        return
    end if
    do
        
        read(numfile,'(a<longeur>)',iostat=err) d
        if (err<0) exit
        k=k+1
        
    end do
   !print *,k
    
    close(numfile)
else 
    k=row
end if

if(associated(tab)) deallocate(tab)
allocate(tab(k,column))

open(unit=numfile,file=trim(Filename),form='formatted',status='old',iostat=err)
    if(err>0) then
        mess="InOutTools:Can't find the file" // trim(Filename)
        return
    end if

if(head) then
    read(numfile, '(a<longeur>)',iostat=errcode) d !header
    if(errcode/=0) then
        err=1
        mess='InOutTools:Read File error while reading header'
        return
    end if
end if




do i= 1, k
    read(numfile, '(a<longeur>)') d
   
    do l= 1 ,column
    
        !to find the next value after the separator
        c=0
        n=0
        nextValue=.True.
        do while(nextValue)
            if (scan(d, ' ')/=n+1) then 
                n= scan(d, separator)
                nextValue=.False.
            else 
                d(1:(longeur-1))=d(2:longeur)
                c=c+1
                
            end if
            if (c==longeur) then
                nextValue=.False.
                err=1
                mess='InOutTools:Scan line error'
                return
            endif
        end do        
        
        
        if(n/=0) then
        
            read(d(1:(n-1)),*,iostat=errcode) tab(i,l)
            
            if(errcode/=0) then
                err=2
                mess='InOutTools:Read File error while reading the data'
                return
            end if
            d(1:(longeur-n))=d((n+1):longeur)
            
        else
            !print *,d(1:5),'11'
            read(d(1:longeur),*,iostat=errcode) tab(i,l)
            if(errcode/=0) then
                err=1
                mess='InOutTools:Read File error while reading the data'
                return
            end if
            !print tab(i,l)
        end if
        
    end do
end do
close(numfile)
end subroutine mat_input


!*********
subroutine mat_output(Filename, header, sep, mat, vect, err, mess)
!~**********************************************************************
!~* Purpose: Write a matrix or a vector into a file
!~**********************************************************************
!~* Programmer: Xun SUN, Cemagref
!~**********************************************************************
!~* Last modified: 10/12/2010
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* In
!~*		1.Filename
!~*		2.header: header name (optional), whose length have to be equal to colume number
!~*		3.sep: separator symbol (optional, default= ' ')
!~*		4.mat: the matrix output
!~*		5.vect: the vect to output
!       !!! You have to decide which one to be output, matrix or vect??
!~*	Out
!~*		1.
!^*		2.err, error code
!^*		3.mess, error message
!~**********************************************************************
character(*), intent(in)::Filename
character(*), intent(in), optional::header(:)
character(len=1),intent(in),optional::sep
character(len=1)::separator
real(mrk), intent(in), optional::mat(:,:)
real(mrk), intent(in), optional::vect(:)
integer(mik), intent(out)::err
character(*), intent(out)::mess
!locals
integer(mik)::n,m,errcode,i,j
real(mrk),allocatable:: table(:,:)
character(15)::temp1
character(4000)::temp2 ! max length of each line
character(15)::title

err=0;mess=''

!check the presence of mat or vect
if(present(mat)) then
    if (present(vect)) then
        err=1
        mess='InOutTools:You may output only one file, a matrix or a vector??'
        return
    else
        m=size(mat(:,1))
        n=size(mat(1,:))
        if (allocated(table)) deallocate(table)
        allocate(table(m,n))
        table=mat
    end if
else
    if (.not. present(vect)) then
        err=1
        mess='InOutTools:No output matrix or vector'
        return
    else
        m=size(vect)
        n=1
        if (allocated(table)) deallocate(table)
        allocate(table(m,n))
        table(:,1)=vect
    end if
end if

!check the sepertor
if (.not. present(sep)) then
    separator=' '
else
    separator=sep
end if

!Writing file
open(2,file=Filename,status='replace',iostat=errcode)
if(errcode/=0) then
    err=2
    mess='InOutTools:Write File error'
    return
end if

if (present(header)) then
    
    if(size(header)/=n) then
        err=1
        mess='InOutTools:size of header and number of colunm mismatch'
        return
    end if
    
    temp2=''
    title=''
    do j=1,(n-1)
        read(header(j),'(a15)') title
        if (separator==' ') then
            write(temp2, '(a<15*(j-1)>,a15,a1)') temp2,title,separator
        else
            temp2=trim(temp2)//title//separator
        endif
    end do
    read(header(n),'(a15)') title
    if (separator==' ') then
        write(temp2, '(a<15*(n-1)>,a15)') temp2,title
    else
        temp2=trim(temp2)//title
    end if
    write(2,'(A<15*n>)') trim(temp2)

end if
    
do i=1,m
    temp2=''
    do j=1,(n-1)
        
        write(temp1,'(e15.6)') table(i,j)
        temp2=trim(temp2)//temp1//separator
        
    end do    
    write(temp1,'(e15.6)') table(i,n)
    temp2=trim(temp2)//temp1
    write(2,'(A<15*n>)') trim(temp2) 
end do

close(2)
end subroutine mat_output


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine addNewMatrix(Filename, newFile, row, column, header, sep, newMat, NewHeader, err, mess)
!^**********************************************************************
!^* Purpose: not finished
!^**********************************************************************
!^* Programmer: Xun SUN, Cemagref
!^**********************************************************************
!^* Last modified:  /  /2011
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
use utilities_dmsl_kit, only: number_string
character(*), intent(in)::Filename
character(*), intent(in), optional::newFile
integer(mik), intent(in), optional::row
character(len=1), intent(in), optional:: sep
integer(mik), intent(in)::  column 
logical, intent(in), optional:: header 
character(15), intent(in), optional:: NewHeader(:)
integer(mik), intent(out)::err
real(mrk), allocatable, intent(in)::newMat(:,:)
character(*), intent(out)::mess
! locals
character(15*column)::oldheader
character(15),allocatable::NHD(:)
character(150)::outputfile
real(mrk),pointer::tab(:,:)
integer(mik)::newMC,i

if (present(header)) then
    if (header) then
        open(unit=1,file=trim(Filename),form='formatted',status='old')
        read(1,'(a<15*column>)',iostat=err) oldheader
        close(1)
    endif
endif
                
        
call mat_input(Filename=Filename, row=row, column=column, header=header, sep=sep, tab=tab, err=err, mess=mess)
if (err>0) then
    mess='addNewMatrix Error:'//trim(mess);return
endif

if (size(tab,dim=1)/=size(newMat,dim=1)) then
    err=1;mess='addNewMatrix Error: new Matrix does not have the same row number as in the document';return
endif


if(allocated(NHD)) deallocate(NHD)
allocate(NHD(size(newMat,dim=2)))

newMC=size(newMat,dim=2)

NHD=''
if (present(newHeader)) then
    if (size(newHeader)/=size(newMat,dim=2)) then
        err=1; mess='addNewMatrix:newHeader size mismatch with newMat column'
        return
    else
        NHD=newHeader
    endif
else
    do i=1,newMC
        NHD= "C_"//trim(number_string(i))
    enddo
endif

if (present(newFile)) then
    outputfile=trim(newFile)
else
    outputfile=trim(Filename)
endif

open(unit=1,file=trim(outputFile),status='replace')
if (.not. present(header)) then
    do i=1,size(tab,dim=1)
        write(1,'(<column+newMC>e15.6)') tab(i,:), newmat(i,:)
    enddo
else
    if (header==.true.) then
        write(1,'(<column+newMC>A15)') oldheader,' ',NHD
        do i=1,size(tab,dim=1)
            write(1,'(<column+newMC>e15.6)') tab(i,:), newmat(i,:)
        enddo
    else
        do i=1,size(tab,dim=1)
            write(1,'(<column+newMC>e15.6)') tab(i,:), newmat(i,:)
        enddo
    endif
endif
close(1)
end subroutine addNewMatrix

end module InOutTools
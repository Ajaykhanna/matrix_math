program matrixsolveprogram
use matrixtools
IMPLICIT NONE

real(8),dimension(:,:),allocatable :: coeff, ans
real(8),dimension(:),  allocatable :: sol
integer :: ios, row_coeff, col_coeff, row_ans, col_ans, i,j
character(80) :: fname_coeff, fname_ans

101 format(a) ! plain text descriptor
102 format(e12.6) 

! take input of matrix file name
write(*,101) 'matrix file must have (row,col) in first line'
write(*,101) 'enter filename of coefficient matrix'
read(*,*) fname_coeff
write(*,101) 'enter filename of answer matrix'
read(*,*) fname_ans
write(*,*)

! open matrix file
open(unit = 11, file = fname_coeff, status = 'old', action = 'read', iostat = ios)
if (ios .ne. 0) stop('error opening unit 11')
open(unit = 12, file = fname_ans, status = 'old', action = 'read', iostat = ios)
if (ios .ne. 0) stop('error opening unit 12')

! read matrix file
read(11,*) row_coeff, col_coeff
allocate(coeff(row_coeff,col_coeff))
read(12,*) row_ans, col_ans
allocate(ans(row_ans,col_ans))

do i = 1,row_coeff
	read(11,*) coeff(i,:)
enddo
do i = 1,row_ans
	read(12,*) ans(i,:)
enddo

write(*,101) 'echo coefficient input'
call matrixwrite(row_coeff,col_coeff,coeff)
write(*,101) 'echo answer input'
call matrixwrite(row_ans,col_ans,ans)

call matrixsolve(row_coeff,col_coeff,coeff,row_ans,col_ans,ans,sol)

do i = 1,row_ans
	write(*,102) sol(i)
enddo
write(*,*)

endprogram matrixsolveprogram
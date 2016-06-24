program matrixsolve
use matrixtools
IMPLICIT NONE

real(8),dimension(:,:),allocatable :: coeff, ans, lower, upper
real(8),dimension(:),  allocatable :: sol, y
real(8) :: lower_sum, upper_sum
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
allocate(sol(row_ans),y(row_ans))

if ((col_coeff .ne. row_ans) .or. (col_ans .ne. 1)) stop('not sure that is a linear system...')

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

do i = 1,row_coeff
	lower_sum = 0.0d0
	do j = 1,(i - 1)
		lower_sum = lower_sum + lower(i,j) * y(j)
	enddo
	y(i) = (ans(i,1) - lower_sum) !/ lower(i,i)
enddo

do i = row_coeff,1,-1
	upper_sum = 0.0d0
	do j = i, (row_coeff - 1)
		upper_sum = upper_sum + upper(i,j) * sol(j)
	enddo
	sol(i) = (y(i) - upper_sum) / upper(i,i)
enddo

do i = 1,row_ans
	write(*,102) y(i)
enddo
write(*,*)
do i = 1,row_ans
	write(*,102) sol(i)
enddo
write(*,*)

endprogram matrixsolve
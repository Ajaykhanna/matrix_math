program matrixsolveprogram
use matrixtools
IMPLICIT NONE

real(8),dimension(:,:),allocatable :: coeff, sol_matrix, ansmat
real(8),dimension(:),  allocatable :: sol, ans
real(8) :: difference
integer :: ios, row_coeff, col_coeff, row_ans, i,j
integer :: row_ansmat, col_ansmat
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
read(12,*) row_ans
allocate(ans(row_ans))
do i = 1,row_coeff
	read(11,*) coeff(i,:)
enddo
do i = 1,row_ans
	read(12,*) ans(i)
enddo

write(*,101) 'echo coefficient input'
call matrixwrite(row_coeff,col_coeff,coeff)
write(*,*)
write(*,101) 'echo answer input'
call vectorwrite(row_ans,ans)

call matrixsolve(row_coeff,col_coeff,coeff,row_ans,ans,sol)

call vectorwrite(row_coeff,sol)

allocate(sol_matrix(row_coeff,1))
do i = 1,row_coeff
	sol_matrix(i,1) = sol(i)
enddo
call matrixmul(row_coeff,col_coeff,coeff,row_coeff,1,sol_matrix,row_ansmat,col_ansmat,ansmat)
! call matrixwrite(row_ansmat,col_ansmat,ansmat)
do i = 1,row_coeff
	difference = ans(i) - ansmat(i,1)
	if ((difference .gt. 1e-5) .or. (difference .lt. -1e-5)) then
		write(*,'(i3)') i
		stop('solution didnt work')
	endif
enddo
write(*,101) 'solution works'


endprogram matrixsolveprogram
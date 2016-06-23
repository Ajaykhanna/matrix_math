program matrixlu
use matrixtools
IMPLICIT NONE

real(8),allocatable,dimension(:,:) :: mat, lower, upper
real(8) :: det_mat
integer :: ios, row, col, i
character(80) :: fname
! logical :: decomp_check

101 format(a) ! plain text descriptor
102 format(e12.6) ! exponential notation, length = 12, decimal = 6

! take input of matrix file name
write(*,101) 'matrix file must have (row,col) in first line'
write(*,101) 'enter filename of matrix'
read(*,*) fname

! open matrix file
open(unit = 11, file = fname, status = 'old', action = 'read', iostat = ios)
if (ios .ne. 0) stop('error opening unit 11')

! read matrix file
read(11,*) row, col
allocate(mat(row,col))
allocate(lower(row,col),upper(row,col))
do i = 1,row
	read(11,*) mat(i,:)
enddo
write(*,101) 'echo matrix input'
call matrixwrite(row,col,mat)

call matrixdecompose(row,col,mat,lower,upper)
call matrixwrite(row,col,lower)
call matrixwrite(row,col,upper)

call matrixdet(row,col,mat,det_mat)
write(*,102) det_mat

endprogram matrixlu
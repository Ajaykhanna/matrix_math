program matrixdetprogram
use matrixtools
IMPLICIT NONE

real(8),dimension(:,:),allocatable :: mat
real(8) :: det_mat
integer :: ios, row, col, i
character(80) :: fname

101 format(a) ! plain text descriptor
102 format(e12.6) 

! take input of matrix file name
write(*,101) 'matrix file must have (row,col) in first line'
write(*,101) 'enter filename of maticient matrix'
read(*,*) fname
write(*,*)

! open matrix file
open(unit = 11, file = fname, status = 'old', action = 'read', iostat = ios)
if (ios .ne. 0) stop('error opening unit 11')


! read matrix file
read(11,*) row, col
allocate(mat(row,col))
do i = 1,row
	read(11,*) mat(i,:)
enddo

write(*,101) 'echo input'
call matrixwrite(row,col,mat)
do i = 1,10
	call matrixdet(row,col,mat,det_mat)
	write(*,102) det_mat
enddo
write(*,*)




endprogram matrixdetprogram
program matrixdetprogram
use matrixtools
IMPLICIT NONE

real(8),dimension(:,:),allocatable :: mat
real(8) :: det_mat
integer :: ios, row, col
character(80) :: fname

101 format(a) ! plain text descriptor
102 format(e12.6) 

! take input of matrix file name
write(*,101) 'enter filename of maticient matrix'
read(*,*) fname
write(*,*)

call matrixread(11,fname,row,col,mat)

write(*,101) 'echo input'
call matrixwrite(row,col,mat)
call matrixdet(row,col,mat,det_mat)
write(*,102) det_mat





endprogram matrixdetprogram
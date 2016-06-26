program matrixinvprogram
use matrixtools
IMPLICIT NONE

integer :: row, col
real(8),dimension(:,:),allocatable :: mat, invmat
character(80) :: fname


101 format(a) ! plain text descriptor

! take input of matrix file name
write(*,101) 'enter filename of matrix'
read(*,*) fname
write(*,*)

call matrixread(11,fname,row,col,mat)

call matrixinv(row,col,mat,invmat)

call matrixwrite(row,col,invmat)



endprogram matrixinvprogram
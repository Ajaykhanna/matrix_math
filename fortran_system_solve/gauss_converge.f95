program gaussconverge
use gausstools
IMPLICIT NONE

character(80) :: fname_mat, fname_ans
integer :: row, col, row_ans
real(8),dimension(:,:),allocatable :: mat
real(8),dimension(:)  ,allocatable :: ans, sol

101 format(a)

write(*,101) 'input filename of matrix'
read(*,*) fname_mat
write(*,101) 'input filename of answer (vector)'
read(*,*) fname_ans


call matrixread(11,fname_mat,row,col,mat)
if (row .ne. col) then
	stop('this program for square matricies')
endif
call vectorread(12,fname_ans,row_ans,ans)
if (row_ans .ne. row) then
	stop('answer vector is not the correct dimension')
endif

call gausssolve(row,mat,ans,sol)


endprogram gaussconverge
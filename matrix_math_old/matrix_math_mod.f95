program matmath
use matrixtools
IMPLICIT NONE

real(8),allocatable, dimension(:,:) :: mat1, mat2, ansmat
character(80) :: fname1, fname2
integer :: row1, col1, row2, col2, ios, r, c, rowans, colans
integer :: i,j,k

101 format(a) ! plain text descriptor

! accept user input
write(*,101) 'matrix file must have (row,col) in 1st line'
write(*,101) 'enter filename of matrix # 1'
read(*,*) fname1
write(*,101) 'enter filename of matrix # 2'
read(*,*) fname2

! open necessary files
open( unit = 11, file = fname1, status = 'old', action = 'read', iostat = ios)
	if (ios .ne. 0) stop ('error opening unit 11')
open( unit = 12, file = fname2, status = 'old', action = 'read', iostat = ios)
	if (ios .ne. 0) stop ('error opening unit 12')

! read matrix size from input files
read(11,*) row1, col1
read(12,*) row2, col2

! allocate appropriate matrix sizes
allocate(mat1(row1,col1))
allocate(mat2(row2,col2))

! read input files
do r = 1, row1
	read(11,*) mat1(r,:)
enddo
do r = 1, row2
	read(12,*) mat2(r,:)
enddo

! echo matrix # 1
write(*,101)
write(*,101) 'echo matrix # 1'
call matrixwrite(row1,col1,mat1)

! echo matrix # 2
write(*,101)
write(*,101) 'echo matrix # 2'
call matrixwrite(row2,col2,mat2)

! MATRIX MULTIPLICATION
! mat1 * mat2
call matrixmul(row1,col1,mat1,row2,col2,mat2,rowans,colans,ansmat)

write(*,101)
write(*,101) 'echo answer matrix'
call matrixwrite(rowans,colans,ansmat)


endprogram matmath
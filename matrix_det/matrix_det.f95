program matrixdet
use matrixtools
IMPLICIT NONE

real(8),dimension(:,:), allocatable :: mat, minor_mat
real(8) :: minor_det, singular_tol, minor_multiplier, det
integer :: ios, r, c, row, col, minor_size, i, j, minor_row, minor_col
character(80) :: fname

101 format(a) ! plain text descriptor
102 format(e12.6)

! take input of matrix file name
write(*,101) 'matrix file must have (row,col) in 1st line'
write(*,101) 'enter filename of matrix'
read(*,*) fname

! open matrix file
open(unit = 11, file = fname, status = 'old', action = 'read', iostat = ios)
	if (ios .ne. 0) stop('error opening unit 11')

! read matrix file
read(11,*) row, col
allocate(mat(row,col))
allocate(minor_mat(row,col))
do r = 1,row
	read(11,*) mat(r,:)
enddo
write(*,101) 'echo matrix input'
call matrixwrite(row,col,mat)
write(*,101)

! make sure matrix is square
if (row .ne. col) stop('determinant is only defined for square matricies')

minor_size = col
det = 0.0d0
! right now this works w/ 3x3 and nothing else
do while (minor_size .ge. 2)
	do c = 1,col
		minor_mat(:,:) = 0.0d0
		minor_row = 0
		! emloy Laplacian determinate technique
		! calculate minor matricies
		minor_multiplier = mat(1,c) * (-1) ** (1 + c)
		minor_size = minor_size - 1
		! write(*,'(i3)') minor_size
		! minor_size operation is wrong
		do i = 1,row
			if (i .eq. 1) cycle
			minor_row = minor_row + 1
			minor_col = 0
			do j = 1,col
				if (j .eq. c) cycle
				minor_col = minor_col + 1
				minor_mat(minor_row,minor_col) = mat(i,j)
				! write(*,'(4i3)') minor_row, minor_col, i, j
			enddo
		enddo	
		call matrixwrite(minor_row,minor_col,minor_mat)
		minor_det = minor_multiplier * ((minor_mat(1,1) * minor_mat(2,2)) - (minor_mat(1,2) * minor_mat(2,1)))
		det = det + minor_det
		write(*,102) minor_det
		write(*,101)
	enddo
enddo!while

write(*,102) det

singular_tol = 1.0d-10
if ((minor_det .lt. singular_tol) .and. (minor_det .gt. (-1 * singular_tol))) then
	write(*,101) 'matrix is singular'
endif




! call matrixwrite(row,col,mat)

endprogram matrixdet
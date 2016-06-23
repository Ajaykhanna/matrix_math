program matrixsolve
IMPLICIT NONE

real(8),allocatable, dimension(:,:) :: identity, mat, invmat, ansmat, errormat
real(8)                             :: increment, error, slope, initial
integer                             :: ios, row, col, r, c, rowans, colans, i, j
character(80)                       :: fname

interface
	subroutine matrixwrite(row,col,mat)
		integer, intent(in)  :: row, col
		integer              :: r, c
		real(8), intent(in)  :: mat(:,:)
	endsubroutine matrixwrite
	
	subroutine matrixmul(row1,col1,mat1,row2,col2,mat2,ansmat,rowans,colans)
		integer, intent(in)  :: row1, col1, row2, col2
		integer              :: dis1, dis2, sim, i, j, k
		integer, intent(out) :: rowans, colans
		real(8), intent(in),  dimension(:,:)              :: mat1, mat2
		real(8), intent(out), dimension(:,:), allocatable :: ansmat
	endsubroutine matrixmul
endinterface

101 format(a)     ! plain text descriptor
102 format(e15.6) ! exponential, length = 15, decimal = 6


write(*,101) 'matrix file must have (row,col) in 1st line'
write(*,101) 'enter filename of matrix'
read(*,*) fname

open(unit = 11, file = fname, status = 'old', action = 'read', iostat = ios);
	if (ios .ne. 0) stop('error opening unit 11')

read(11,*) row, col
if (row .ne. col) stop('right now this only works for square matricies')
allocate(mat(row,col),identity(row,col),invmat(row,col))
allocate(ansmat(row,col),errormat(row,col))
do r = 1,row
	read(11,*) mat(r,:)
enddo

! build identity matrix
do r = 1,row
	do c = 1,col
		if (r .eq. c) then
			identity(r,c) = 1.0d0
		else
			identity(r,c) = 0.0d0
		endif
	enddo
enddo

! set solution to all 1's
initial = 1.0d0
invmat(:,:) = initial

! calculate and record the total error
call matrixmul(row,col,mat,row,col,invmat,ansmat,rowans,colans)
error = 0.0d0
do r = 1,row
	do c = 1,col
		error = error + (ansmat(r,c) - identity(r,c))
	enddo
enddo
! debug write
! write(*,102) error
! call matrixwrite(rowans,colans,ansmat)

! set "bump" value
increment = 0.5d0
! bump each and record the total error
! assume the error is linear
do r = 1,row
	do c= 1,col
		invmat(r,c) = invmat(r,c) + increment
		call matrixmul(row,col,mat,row,col,invmat,ansmat,rowans,colans)
		! calculate the new error based on this "bump"
		errormat(r,c) = 0.0d0
		do i = 1,row
			do j = 1,col
				errormat(r,c) = errormat(r,c) + (ansmat(r,c) - identity(r,c))
			enddo
		enddo
		! set the error to 0 (linearly)
		slope = (errormat(r,c) - error) / (increment)
		invmat(r,c) = initial - (error / slope)		
	enddo
enddo


		


call matrixwrite(row,col,invmat)
write(*,*)
call matrixmul(row,col,mat,row,col,invmat,ansmat,rowans,colans)
call matrixwrite(rowans,colans,ansmat)

endprogram matrixsolve

subroutine matrixwrite(row,col,mat)
	IMPLICIT NONE
	integer, intent(in)  :: row, col
	integer              :: r, c
	real(8), intent(in)  :: mat(:,:)
	
	101 format(a)     ! plain text descriptor
	102 format(e15.6) ! exponential, length = 15, decimal = 6
	
	do r = 1, row
		write(*,101,advance= 'no') '[ '
		do c = 1, col
			write(*,102,advance = 'no') mat(r,c)
		enddo
		write(*,101) ' ]'
	enddo
endsubroutine matrixwrite

subroutine matrixmul(row1,col1,mat1,row2,col2,mat2,ansmat,rowans,colans)
	IMPLICIT NONE
	integer, intent(in)  :: row1, col1, row2, col2
	integer              :: dis1, dis2, sim, i, j, k
	integer, intent(out) :: rowans, colans
	real(8), intent(in),  dimension(:,:)              :: mat1, mat2
	real(8), intent(out), dimension(:,:), allocatable :: ansmat
	
	101 format(a) ! plain text descriptor
	
	! check if matrix multiplication is allowed
	if (col1 .eq. row2) then
		sim = col1
		dis1 = row1
		dis2 = col2
	elseif (row1 .eq. col2) then
		write(*,101) 'I think you want mat2 * mat1'
		write(*,101) 'I cant do that right now'
		write(*,101) 'Switch mat1 & mat2 in input.'
		stop 'END OF PROGRAM'
	else
		stop('matrix multiplication not possible based on dimensions')
	endif
	
	allocate(ansmat(dis1,dis2))
	rowans = dis1
	colans = dis2
	
	do i = 1, dis1
		do j = 1, dis2
			ansmat(i,j) = 0.0d0
			do k = 1, sim
				! debug write
				! write(*,'(3e15.6,3i3)') mat1(i,k), mat2(k,i), (mat1(i,k) * mat2(k,i)), i, j, k
				ansmat(i,j) = ansmat(i,j) + (mat1(i,k) * mat2(k,j))
			enddo
		enddo
	enddo
endsubroutine
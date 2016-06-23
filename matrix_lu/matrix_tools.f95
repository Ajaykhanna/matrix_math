module matrixtools
IMPLICIT NONE

contains
subroutine matrixwrite(row,col,mat)
	IMPLICIT NONE
	integer, intent(in)  :: row, col
	integer              :: r, c
	real(8), intent(in)  :: mat(:,:)
	
	101 format(a)     ! plain text descriptor
	102 format(e12.6,x) ! exponential, length = 12, decimal = 6, follow w/ space
	
	do r = 1, row
		write(*,101,advance= 'no') '[ '
		do c = 1, col
			write(*,102,advance = 'no') mat(r,c)
		enddo
		write(*,101) ']'
	enddo
	write(*,101) ! blank line
	
endsubroutine matrixwrite

subroutine matrixmul(row1,col1,mat1,row2,col2,mat2,rowans,colans,ansmat)
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
endsubroutine matrixmul

subroutine matrixdecompose(row,col,mat,lower,upper)
	IMPLICIT NONE
	integer,intent(in) :: row, col
	integer :: ansrow, anscol
	integer :: i, j, k, w
	real(8),dimension(:,:),intent(in) :: mat
	real(8),dimension(:,:),intent(out),allocatable :: upper, lower
	real(8),dimension(:,:),allocatable :: ansmat
	real(8) :: asum, bsum
	
	101 format(a) ! plain text descriptor
	
	! Doolittle Algorithm
	! See Notre Dame PDF
	
	if (row .ne. col) stop('LU decomposition is only defined for squqre matricies')
	allocate(lower(row,col),upper(row,col))
	lower(:,:) = 0.0d0
	upper(:,:) = 0.0d0
	do i = 1,row
		lower(i,i) = 1.0d0
	enddo
	
	do i = 1,row
		do j = i,row
			asum = 0.0d0
			do w = 1,(i - 1)
				asum = asum + lower(i,w) * upper(w,j)
			enddo
			upper(i,j) = mat(i,j) - asum
		enddo
		do k = 1,col
			bsum = 0.0d0
			do w = 1,(i - 1)
				bsum = bsum + lower(k,w) * upper(w,i)
			enddo
			lower(k,i) = (mat(k,i) - bsum) / upper(i,i)
		enddo
	enddo
	
	! Multiply L * U to make sure it worked
	! TO-DO: Remove this later
	! call matrixmul(row,col,lower,row,col,upper,ansrow,anscol,ansmat)
	! call matrixwrite(ansrow,anscol,ansmat)
	! do i = 1,col
		! do j = 1,row
			! if (ansmat(i,j) .ne. mat(i,j)) then
				! write(*,'(2i3)') i,j
				! stop('didnt work')
			! endif
		! enddo
	! enddo
	! write(*,101) 'decomposition successful'
	
endsubroutine matrixdecompose

subroutine matrixdet(row,col,mat,det_mat)
	IMPLICIT NONE
	integer,intent(in) :: row, col
	integer :: i
	real(8),dimension(:,:),intent(in) :: mat
	real(8),intent(out) :: det_mat
	real(8),dimension(:,:),allocatable :: lower, upper
	
	if (row .ne. col) stop('determinant only defined for square matricies')
	
	! First, decompose the matrix
	call matrixdecompose(row,col,mat,lower,upper)
	! determinant of a triangular matrix is the multiplication of the "angle"
	! det_lower = 1.0d0 ! definition from Doolittle algorithm
	det_mat = 1.0d0
	do i = 1,row
		det_mat = det_mat * upper(i,i)
	enddo

endsubroutine matrixdet
	
endmodule matrixtools
module matrixtools
IMPLICIT NONE

contains
subroutine matrixwrite(row,col,mat)
	IMPLICIT NONE
	integer, intent(in)  :: row, col
	integer              :: r, c
	real(8), intent(in)  :: mat(:,:)
	
	101 format(a)       ! plain text descriptor
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

subroutine vectorwrite(row,vect)
	IMPLICIT NONE
	integer,intent(in) :: row
	integer :: i
	real(8),dimension(:),intent(in) :: vect
	
	103 format(a,e12.6,a)
	
	do i = 1,row
		write(*,103) '[ ', vect(i), ' ]'
	enddo
	
	write(*,*)
endsubroutine vectorwrite

subroutine vectorread(unit_num,fname,row,vector)
	IMPLICIT NONE
	integer,intent(in) :: unit_num
	integer,intent(out) :: row
	integer :: ios, i
	real(8),dimension(:),intent(out),allocatable :: vector
	character(80),intent(in) :: fname

	! vector file must have (row) in the first line

	open(unit = unit_num, file = fname, status = 'old', action = 'read', iostat = ios)
	if (ios .ne. 0) then
		write(*,'(a,i3,a,a)') 'error opening unit', unit_num, ' -- ', fname
		stop('END PROGRAM')
	endif

	read(unit_num,*) row
	allocate(vector(row))
	do i = 1,row
		read(unit_num,*) vector(i)
	enddo
endsubroutine vectorread
	
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

subroutine nonzerodiag(row,col,mat,ans)
	IMPLICIT NONE
	integer,intent(in) :: row, col
	real(8),dimension(:,:),intent(inout) :: mat
	real(8),dimension(:)  ,intent(inout) :: ans
	integer :: i, j, z, n, x, min_non_zero, max_non_zero
	real(8) :: pivot_tol, store_ans
	real(8),dimension(:),allocatable :: store_row
	logical :: diagonal_check
	logical,dimension(:),allocatable :: done
	integer,dimension(:),allocatable :: non_zero

	if (row .ne. col) then
		stop('nonzerodiag only for square matricies')
	endif

	z = 0
	n = row

	allocate(non_zero(n))
	pivot_tol = 1.0d-5
	diagonal_check = .true.

	do i = 1,n
		non_zero(i) = 0
		do j = 1,n
			if ((mat(i,j) .gt. pivot_tol) .or. (mat(i,j) .lt. (-1)*pivot_tol)) then
				non_zero(i) = non_zero(i) + 1
			endif
		enddo
	enddo

	min_non_zero = n
	max_non_zero = 0
	do i = 1,n
		if (non_zero(i) .lt. min_non_zero) then
			min_non_zero = non_zero(i)
		endif
		if (non_zero(i) .gt. max_non_zero) then
			max_non_zero = non_zero(i)
		endif
	enddo

	z = 0
	allocate(done(n))
	do while (diagonal_check)
		done(:) = .false.
		diagonal_check = .false.
		z = z + 1
		do x = min_non_zero,max_non_zero
			do i = 1,n
				if ((non_zero(i) .eq. x) .and. ((mat(i,i) .lt. pivot_tol) .and. (mat(i,i) .gt. (-1)*pivot_tol))) then
					diagonal_check = .true.
					do j = 1,n
						if ((mat(i,j) .gt. pivot_tol) .or. (mat(i,j) .lt. (-1)*pivot_tol)) then
							if (.not. done(j)) then
								! swap
								if(.not. allocated(store_row)) allocate(store_row(n))
								store_row = mat(i,:)
								store_ans = ans(i)
								mat(i,:)  = mat(j,:)
								ans(i)    = ans(j)
								mat(j,:)  = store_row
								ans(j)    = store_ans
								done(j)   = .true.
							endif
						endif
					enddo
				endif
			enddo
		enddo
	enddo

	! Check if it worked (not really necessary)
	! do i = 1,n
	! 	if ((mat(i,i) .lt. pivot_tol) .and. (mat(i,i) .gt. (-1)*pivot_tol)) then
	! 		stop('error')
	! 	endif
	! enddo

endsubroutine nonzerodiag

subroutine matrixdecompose(row,col,mat,lower,upper)
	IMPLICIT NONE
	integer,intent(in) :: row, col
	integer :: ansrow, anscol
	integer :: i, j, k, w
	real(8),dimension(:,:),intent(in) :: mat
	real(8),dimension(:,:),intent(out),allocatable :: lower, upper
	real(8),dimension(:,:),allocatable :: ansmat
	real(8) :: asum, bsum
	
	101 format(a) ! plain text descriptor
	
	! Doolittle Algorithm
	! See wiki PDF
	
	if (row .ne. col) stop('LU decomposition is only defined for square matricies')
	allocate(lower(row,col),upper(row,col))
	lower(:,:) = 0.0d0
	upper(:,:) = 0.0d0
	do i = 1,row
		lower(i,i) = 1.0d0
	enddo
	
	do i = 1,row
		! write(*,'(i6)') i
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
	! do i = 1,col
	! 	do j = 1,row
	! 		if (((ansmat(i,j) - mat(i,j)) .gt. 1.0d-5) .or. ((ansmat(i,j) - mat(i,j)) .lt. -1.0d-5)) then
	! 			write(*,'(2i3)') i,j
	! 			write(*,'(2e12.6)') ansmat(i,j), mat(i,j)
	! 			stop('decomposition didnt work')
	! 		endif
	! 	enddo
	! enddo
	! write(*,101) 'decomposition successful'
	! write(*,*)
	
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

subroutine matrixsolve(row_coeff,col_coeff,mat,row_ans,ans,sol)
	IMPLICIT NONE
	integer,intent(in) :: row_coeff, col_coeff, row_ans
	integer :: i, j
	real(8),dimension(:,:),intent(inout) :: mat
	real(8),dimension(:),intent(inout) :: ans
	real(8),dimension(:),intent(out),allocatable :: sol
	real(8),dimension(:),allocatable :: y
	real(8) :: lower_sum, upper_sum
	real(8),dimension(:,:),allocatable :: lower, upper

	101 format(a)
	
	if (col_coeff .ne. row_ans) stop('not sure that is a linear system...')
	allocate(sol(row_coeff),y(row_coeff))

	call nonzerodiag(row_coeff,col_coeff,mat,ans)

	call matrixdecompose(row_coeff,col_coeff,mat,lower,upper)
	
	! forward subtitution (lower matrix)
	do i = 1,row_coeff
		lower_sum = 0.0d0
		do j = 1,(i - 1)
			lower_sum = lower_sum + (lower(i,j) * y(j))
		enddo
		y(i) = (ans(i) - lower_sum) / lower(i,i)
	enddo
	
	! backward subtitution (upper matrix)
	do i = row_coeff,1,-1
		upper_sum = 0.0d0
		do j = (i + 1),row_coeff
			upper_sum = upper_sum + upper(i,j) * sol(j)
		enddo
		sol(i) = (y(i) - upper_sum) / upper(i,i)
	enddo

endsubroutine matrixsolve

subroutine loweruppersolve(row_coeff,col_coeff,lower,upper,row_ans,ans,sol)
	IMPLICIT NONE
	integer,intent(in) :: row_coeff, col_coeff, row_ans
	integer :: i, j
	real(8),dimension(:,:),intent(in) :: lower, upper
	real(8),dimension(:),intent(in) :: ans
	real(8),dimension(:),intent(out),allocatable :: sol
	real(8),dimension(:),allocatable :: y
	real(8) :: lower_sum, upper_sum

	if (col_coeff .ne. row_ans) stop('not sure that is a linear system...')
	allocate(sol(row_coeff),y(row_coeff))

	do i = 1,row_coeff
		lower_sum = 0.0d0
		do j = 1,(i - 1)
			lower_sum = lower_sum + (lower(i,j) * y(j))
		enddo
		y(i) = (ans(i) - lower_sum) / lower(i,i)
	enddo
	
	! backward subtitution (upper matrix)
	do i = row_coeff,1,-1
		upper_sum = 0.0d0
		do j = (i + 1),row_coeff
			upper_sum = upper_sum + upper(i,j) * sol(j)
		enddo
		sol(i) = (y(i) - upper_sum) / upper(i,i)
	enddo

endsubroutine loweruppersolve


subroutine matrixread(unit_num,fname,row,col,mat)
	IMPLICIT NONE
	integer,intent(in) :: unit_num
	integer,intent(out) :: row,col
	integer :: ios,i
	real(8),intent(out),dimension(:,:),allocatable :: mat
	character(80),intent(in) :: fname
	
	! matrix file must have (row,col) in first line
	
	open(unit = unit_num, file = fname, status = 'old', action = 'read', iostat = ios)
	if (ios .ne. 0) then
		write(*,'(a,i3,a,a)') 'error opening unit', unit_num, ' -- ', fname
		stop('END PROGRAM')
	endif
	
	read(unit_num,*) row,col
	allocate(mat(row,col))
	do i = 1,row
		read(unit_num,*) mat(i,:)
	enddo
endsubroutine matrixread

subroutine matrixinv(row,col,mat,invmat)
	IMPLICIT NONE
	integer,intent(in) :: row,col
	integer :: i
	real(8),dimension(:,:),intent(in) :: mat
	real(8),dimension(:,:),intent(out),allocatable :: invmat
	real(8),dimension(:),allocatable :: invcol, identitycol
	real(8),dimension(:,:),allocatable :: identity
	real(8) :: det_mat, tol
	
	! NEVER EVER USE AN INVERT TO PERFORM A SOLVE!
	! THE INVERSION CALLS THE SOLVE SUBROUTINE N TIMES!
	
	101 format(a) ! plain text descriptor
	
	tol = 1.0d-2
	
	if (row .ne. col) stop ('inversion only for square matricies')
	
	! IEEE-FPE floating point error
	call matrixdet(row,col,mat,det_mat)
	allocate(invmat(row,col),identity(row,col))
	allocate(invcol(row),identitycol(row))
	invmat(:,:) = 0.0d0
	identity(:,:) = 0.0d0
	do i = 1,row
		identity(i,i) = 1.0d0
	enddo
		
	if ((det_mat .lt. tol) .and. (det_mat .gt. ((-1) * tol))) then
		write(*,101) 'matrix is singular (det = 0)'
		write(*,101) 'singular matricies are not invertable'
		stop('END PROGRAM')
	endif
	
	do i = 1,row
		invcol(:) = 0.0d0
		identitycol(:) = identity(:,i)
		
		! call matrixsolve(row,col,mat,row,identitycol,invcol)
		stop('commented out')
		
		invmat(:,i) = invcol(:)
	enddo
	
endsubroutine matrixinv

subroutine gausssolve(n,mat,ans,sol)
	IMPLICIT NONE
	real(8),dimension(:,:),allocatable,intent(inout)  :: mat
	real(8),dimension(:)  ,allocatable,intent(inout)  :: ans
	real(8),dimension(:)  ,allocatable,intent(out) :: sol
	real(8),dimension(:)  ,allocatable             :: old_sol, difference
	integer,intent(in) :: n
	integer :: i, j, z, max_iteration, x
	real(8) :: multiply_sum, solution_tol, convergence

	101 format(a)

	max_iteration = 500
	solution_tol = 1.0d-10
	allocate(sol(n),old_sol(n))
	allocate(difference(n))
	sol = ans
	convergence = 1.0d0
	z = 0

	call nonzerodiag(n,n,mat,ans)

	z = 0
	do while (convergence .gt. solution_tol)
	! do x = 1,1
		convergence = 0.0d0
		z = z + 1
		if (z .eq. (max_iteration + 1)) then
			write(*,'(a,i4,a)') 'failed to converge in ', max_iteration, ' iterations'
		endif
		do i = 1,n
			old_sol(i) = sol(i)
			multiply_sum = 0.0d0
			do j = 1,n
				if (j .ne. i) then
					multiply_sum = multiply_sum + mat(i,j) * sol(j)
				endif
			enddo
			sol(i) = (ans(i) - multiply_sum) / mat(i,i)
			if (isnan(sol(i))) then
				write(*,'(i4)') i
				stop('sol(i) is nan')
			endif
			difference(i) = (sol(i) - old_sol(i)) / sol(i)
			difference(i) = difference(i) ** 2
			convergence = convergence + difference(i)
		enddo
	enddo

	write(*,'(a,i4)') 'iterations ', z

endsubroutine gausssolve
	
endmodule matrixtools
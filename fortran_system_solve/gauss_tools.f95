module gausstools
IMPLICIT NONE
contains

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

subroutine gausssolve(n,mat,ans,sol)
	IMPLICIT NONE
	real(8),dimension(:,:),allocatable,intent(in)  :: mat
	real(8),dimension(:)  ,allocatable,intent(in)  :: ans
	real(8),dimension(:)  ,allocatable,intent(out) :: sol
	real(8),dimension(:)  ,allocatable             :: old_sol, difference
	integer,intent(in) :: n
	integer :: i, j, z, max_iteration
	real(8) :: multiply_sum, solution_tol, convergence

	101 format(a)

	max_iteration = 500
	solution_tol = 1.0d-10
	allocate(sol(n),old_sol(n))
	allocate(difference(n))
	sol = 1.0d0
	convergence = 1.0d0
	z = 0

	do while (convergence .gt. solution_tol)
		convergence = 0.0d0
		z = z + 1
		do i = 1,n
			old_sol(i) = sol(i)
			multiply_sum = 0.0d0
			do j = 1,n
				if (j .ne. i) then
					multiply_sum = multiply_sum + mat(i,j) * sol(j)
				endif
			enddo
			sol(i) = (ans(i) - multiply_sum) / mat(i,i)
			difference(i) = (sol(i) - old_sol(i)) / sol(i)
			difference(i) = difference(i) ** 2
			convergence = convergence + difference(i)
		enddo
	enddo
endsubroutine gausssolve



endmodule gausstools
program matmath
IMPLICIT NONE

real(8),allocatable, dimension(:,:) :: mat1, mat2, ansmat
character(80) :: fname1, fname2
integer :: row1, col1, row2, col2, ios, r, c, rowans, colans
integer :: i,j,k

101 format(a) ! plain text descriptor
102 format(e15.6) ! exponential, length = 15, decimal = 6

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
call matrixmul(row1,col1,mat1,row2,col2,mat1,ansmat,rowans,colans)

write(*,101)
write(*,101) 'echo answer matrix'
call matrixwrite(rowans,colans,ansmat)


endprogram matmath

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

	


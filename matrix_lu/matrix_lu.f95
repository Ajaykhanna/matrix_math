program matrixlu
use matrixtools
IMPLICIT NONE

real(8),allocatable,dimension(:,:) :: mat, lower, upper, ansmat
real(8) :: first_sum, second_sum, third_sum, fourtH_sum
integer :: ios, row, col, i, j, t, ansrow, anscol
character(80) :: fname
! logical :: decomp_check

101 format(a) ! plain text descriptor

! take input of matrix file name
write(*,101) 'matrix file must have (row,col) in first line'
write(*,101) 'enter filename of matrix'
read(*,*) fname

! open matrix file
open(unit = 11, file = fname, status = 'old', action = 'read', iostat = ios)
if (ios .ne. 0) stop('error opening unit 11')

! read matrix file
read(11,*) row, col
allocate(mat(row,col))
allocate(lower(row,col),upper(row,col))
do i = 1,row
	read(11,*) mat(i,:)
enddo
write(*,101) 'echo matrix input'
call matrixwrite(row,col,mat)
write(*,101)

! make sure matrix is square
if (row .ne. col) stop('LU decomposition is only defined for square matricies')

lower(:,:) = 0.0d0
upper(:,:) = 0.0d0
do i = 1,row
	lower(i,i) = 1.0d0
enddo

! call matrixwrite(row,col,lower)
! write(*,*)
! call matrixwrite(row,col,upper)

upper(1,:) = mat(1,:)
lower(:,1) = mat(:,1) / upper(1,1)
do i = 2,(row - 1)
	first_sum = 0.0d0
	second_sum = 0.0d0
	third_sum = 0.0d0
	do t = 1,(i - 1)
		first_sum = first_sum + (lower(i,t) * upper(t,i))
	enddo
	upper(i,i) = mat(i,i) - first_sum
	do j = (i + 1),row
		do t = 1,(i - 1)
			second_sum = second_sum + (lower(i,t) * upper(t,j))
			third_sum = third_sum + (lower(j,t) * upper(t,i))
		enddo
		upper(i,j) = mat(i,j) - second_sum
		lower(j,i) = (mat(j,i) - third_sum) / upper(i,i)
	enddo
enddo
fourth_sum = 0.0d0
do t = 1,(row - 1)
	fourth_sum = fourth_sum + (lower(row,t) * upper(t,col))
enddo
upper(row,col) = mat(row,col) - fourth_sum

call matrixwrite(row,col,upper)
write(*,*)
call matrixwrite(row,col,lower)

call matrixmul(row,col,lower,row,col,upper,ansrow,anscol,ansmat)
write(*,*)
call matrixwrite(ansrow,anscol,ansmat)
! decomp_check = .false.
do i = 1,col
	do j = 1,row
		if (ansmat(i,j) .ne. mat(i,j)) then
			! decomp_check = .true.
			write(*,'(2i3)') i,j
			stop('didnt work')
		endif
	enddo
enddo
		
endprogram matrixlu
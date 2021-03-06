program gausssolver
use gausstools
IMPLICIT NONE

character(80) :: fname_mat, fname_ans
integer :: row, col, row_ans
real(8),dimension(:,:),allocatable :: mat
real(8),dimension(:)  ,allocatable :: ans, sol

integer :: i, j, k
integer :: s, t, n
real(8) :: unique_tol

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

n = row
allocate(sol(n))
ext(1:n,1:n) = mat
ext(:,(n + 1)) = ans
call matrixwrite(n, n + 1, ext)

do i = 1,n
	do j = (i + 1),(n + 1)
		ext(i,j) = ext(i,j) / ext(i,i)
	enddo
	ext(i,i) = 1.0d0
	do j = 1,n
		if (j .ne. i) then
			do k = (i + 1),(n + 1)
				ext(j,k) = ext(j,k) - (ext(i,k) * ext(j,i))
			enddo
		endif
	enddo
enddo

do i = 1,n
	sol(i) = ext(i,(n + 1))
enddo

! call vectorwrite(n,sol)
call matrixwrite(n,n + 1,ext)



endprogram gausssolver
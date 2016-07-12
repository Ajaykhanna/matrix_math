program transportsolve
use transporttools
use matrixtools
IMPLICIT NONE

character(80)  :: fname, fname_input
integer :: group, material
character(3),dimension(:),allocatable :: names
real(8),dimension(:,:),allocatable :: t_in, tr_in, d_in, a_in, c_in, f_in, nu_in, x_in, r_in
real(8),dimension(:,:,:),allocatable :: s_in

integer :: i,j,ios

integer :: cells
real(8) :: dx

integer,dimension(:),allocatable :: material_index_cells
real(8) :: total_length

real(8),dimension(:,:),allocatable :: t, tr, d, a, c, f, nu, x, r
real(8),dimension(:,:,:),allocatable :: s

real(8),dimension(:,:,:),allocatable :: A_mat
real(8),dimension(:,:,:),allocatable :: phi, phibar, current
real(8),dimension(:,:,:),allocatable :: Q_f, Q_up, Q_down, Q
real(8),dimension(:,:),allocatable :: f_mat, y_mat, A_mat_g
real(8),dimension(:),allocatable :: f_mat_g, y_mat_g
integer :: g, g_prime, z, max_iteration
real(8),dimension(:),allocatable :: k
real(8) :: kerror, phibarerror, epsilon_k, epsilon_phi
real(8) :: fission_term, upscatter_term, downscatter_term
real(8) :: numerator, denominator, sum_phibar, sum_phi
real(8),dimension(:,:),allocatable :: ansmat

character(8)  :: date
character(10) :: time
character(20) :: date_out, time_out
real :: start_cpu, finish_cpu


101 format(a) ! plain text descriptor
102 format(i1)

open(unit = 11, file = 'runtime.out', status = 'replace', action = 'write', iostat = ios)
if (ios .ne. 0) then
	write(*,'(a,i3,a,a)') 'error opening unit', 11, ' -- ', fname
	stop('END PROGRAM')
endif

call date_and_time(date,time)
call cpu_time(start_cpu)
date_out = date(5:6) // '/' // date(7:8) // '/' // date(1:4)
time_out = time(1:2) // ':' // time(3:4) // ':' // time(5:10)

write(*,101) 'execution time'
write(*,'(a,x,a)') date_out, time_out

write(*,101) 'input filename'
read(*,*) fname

write(*,101) 'enter test input filename'
read(*,*) fname_input

call xsread(fname,group,material,names,t_in,tr_in,d_in,a_in,c_in,f_in,nu_in,x_in,s_in,r_in)
call transportinput (fname_input,dx,cells,names,material_index_cells,material)


call xsbuild(material_index_cells,cells,group,t_in,tr_in,d_in,a_in,c_in,f_in,nu_in,x_in,r_in,s_in,t,tr,d,a,c,f,nu,x,r,s)

total_length = real(cells,8) * dx

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

! build A_mat matrix of geometry/material properties
allocate(A_mat((3 * cells + 2),(3 * cells + 2),group))
allocate(A_mat_g((3 * cells + 2),(3 * cells + 2)))
A_mat(:,:,:) = 0.0d0
do g = 1,group
	do i = 1,cells
		A_mat((3 * i - 2),(3 * i    ),g) = r(i,g) * dx
		A_mat((3 * i - 2),(3 * i - 1),g) = -1.0d0
		A_mat((3 * i - 2),(3 * i + 2),g) = 1.0d0

		A_mat((3 * i - 1),(3 * i    ),g) = (-1.0d0) * d(i,g)
		A_mat((3 * i - 1),(3 * i + 1),g) = d(i,g)
		A_mat((3 * i - 1),(3 * i + 2),g) = 0.5d0 * dx

		A_mat((3 * i),(3 * i    ),g) = d(i,g)
		A_mat((3 * i),(3 * i - 1),g) = 0.5d0 * dx
		A_mat((3 * i),(3 * i - 2),g) = (-1.0d0) * d(i,g)
	enddo
enddo
! reflective BC in A_mat matrix
A_mat((3 * cells + 1),1,:) = -1.0d0
A_mat((3 * cells + 1),3,:) = 1.0d0
A_mat((3 * cells + 2),(3 * cells + 1),:) = -1.0d0
A_mat((3 * cells + 2),(3 * cells    ),:) = 1.0d0


! initial guesses
z = 1
phibarerror = 1.0d0
kerror = 1.0d0
allocate(f_mat((3 * cells + 2),group))
allocate(f_mat_g(3 * cells + 2))
allocate(y_mat((3 * cells + 2),group))
allocate(y_mat_g(3 * cells + 2))
f_mat(:,:) = 0.0d0

max_iteration = 100
allocate(k(max_iteration))
allocate(phibar(cells,group,max_iteration))
allocate(phi(cells,group,max_iteration))
allocate(current(cells,group,max_iteration))
allocate(Q_f(cells,group,max_iteration))
allocate(Q_up(cells,group,max_iteration))
allocate(Q_down(cells,group,max_iteration))
allocate(Q(cells,group,max_iteration))
k(1) = 1.0d0
do i = 1,cells
	do j = 1,group
		phibar(i,j,z) = 2.1d1 ! initial guess = 21
	enddo
enddo

! convergence tolerance
epsilon_k   = 1.0d-5
epsilon_phi = 1.0d-5

! stop('solver loop')
! solver loop
do while ((phibarerror .gt. epsilon_phi) .or. (kerror .gt. epsilon_k))
	z = z + 1
	write(*,'(i3)') z
	! check for number of iterations
	if (z .eq. max_iteration + 1) then
		do i = 1,(z - 1)
			write(*,'(e12.6)') k(i)
		enddo
		write(*,'(a,i5,a)') 'failed to converge after', max_iteration, ' iterations'
		stop
	endif
	! calculate Q values
	do g = 1,group
		do i = 1,cells
			fission_term = 0.0d0
			do g_prime = 1,group
				fission_term = fission_term + nu(i,g_prime) * f(i,g_prime) * phibar(i,g_prime,(z - 1))
			enddo
			Q_f(i,g,(z - 1)) = x(i,g) / k(z - 1) * fission_term

			upscatter_term = 0.0d0
			do g_prime = (g + 1),group
				upscatter_term = upscatter_term + s(i,g_prime,g) * phibar(i,g_prime,(z - 1))
			enddo
			Q_up(i,g,(z - 1)) = upscatter_term

			downscatter_term = 0.0d0
			do g_prime = 1,(g - 1)
				downscatter_term = downscatter_term + s(i,g_prime,g) * phibar(i,g_prime,z)
			enddo
			Q_down(i,g,z) = downscatter_term

			Q(i,g,z) = Q_down(i,g,z) + Q_up(i,g,(z - 1)) + Q_f(i,g,(z - 1))
			! turn Qs into f vector
			f_mat((3 * i - 2),g) = Q(i,g,z) * dx
		enddo
		
		A_mat_g = A_mat(:,:,g)
		f_mat_g = f_mat(:,g)
		y_mat_g = y_mat(:,g)
		
		call matrixsolve((3 * cells + 2),(3 * cells + 2),A_mat_g,(3 * cells + 2),f_mat_g,y_mat_g)

		write(*,'(i3,a,i3)') z, ' , ', g
		if (isnan(y_mat_g(1))) stop('y is nan')
		y_mat(:,g) = y_mat_g

		! parse y vector for required functions
		do i = 1,(cells + 1)
			if (i .le. cells) then
				phi(i,g,z) = y_mat((3 * i - 2),g)
				current(i,g,z) = y_mat((3 * i - 1),g)
				phibar(i,g,z) = y_mat((3 * i),g)
			else
				phi(i,g,z) = y_mat((3 * i - 2),g)
				current(i,g,z) = y_mat((3 * i - 1),g)
			endif
		enddo
	enddo

	! update the guess for k
	numerator = 0.0d0
	denominator = 0.0d0
	do g_prime = 1,group
		do i = 1,cells
			numerator   = numerator   + nu(i,g_prime) * f(i,g_prime) * phibar(i,g_prime,z) * dx
			denominator = denominator + nu(i,g_prime) * f(i,g_prime) * phibar(i,g_prime,(z - 1)) * dx
		enddo
	enddo
	write(*,*) numerator, denominator
	k(z) = k(z - 1) * (numerator / denominator)

	! make sure k is positive
	if (k(z) .lt. 0) then
		stop('negative k value')
	endif

	! normalized calculated functions
	sum_phibar = 0.0d0
	sum_phi = 0.0d0
	do g = 1,group
		do i = 1,cells
			sum_phibar = sum_phibar + phibar(i,g,z) * dx
			sum_phi    = sum_phi    + phi(i,g,z) * dx
		enddo
	enddo
	phibar(:,:,z) = (phibar(:,:,z) * total_length) / sum_phibar
	phi(:,:,z)    = (phi(:,:,z)    * total_length) / sum_phi

	! calculate iteration error for the next iteration
	phibarerror = 0.0d0 ! TO-DO: FIX THIS SHIT! (involves taking the norm)
	kerror = abs(1 - (k(z - 1) / k(z)))

enddo

write(*,101) 'I FOUND K!!!'
! write(*,101) '(maybe)'
write(*,'(e12.6)') k(z)

call cpu_time(finish_cpu)
write(*,101)
write(*,'(e12.6)') finish_cpu


endprogram transportsolve
program transportsolve
use transporttools
use matrixtools
IMPLICIT NONE

character(80)  :: fname
integer :: group, material
character(3),dimension(:),allocatable :: names
real(8),dimension(:,:),allocatable :: t_in, tr_in, d_in, a_in, c_in, f_in, nu_in, x_in, r_in
real(8),dimension(:,:,:),allocatable :: s_in

integer :: i,j

integer :: ios, cells, pins, npin, pin_num, pos, line_num, total_parts, beginning
integer,dimension(:),allocatable :: pinmap, parts, parts_sum
real(8),dimension(:),allocatable :: dimensions, temp_dimensions
character(3),dimension(:),allocatable :: material_list, temp_material_list
real(8) :: dx
character(100) :: line, label
character(80)  :: fname_test

integer :: regions, region_counter
integer,dimension(:),allocatable :: material_index_cells
real(8) :: capital_H, mesh_tol
real(8),dimension(:),allocatable :: x_coords, mesh
character(3),dimension(:),allocatable :: material_region, material_cells

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

101 format(a) ! plain text descriptor
102 format(i1)

write(*,101) 'input filename'
read(*,*) fname
call xsread(fname,group,material,names,t_in,tr_in,d_in,a_in,c_in,f_in,nu_in,x_in,s_in,r_in)

! call matrixwrite(material,group,r)


write(*,101) 'enter test input filename'
read(*,*) fname_test

open(unit = 11, file = fname_test, status = 'old', action = 'read', iostat = ios)
if (ios .ne. 0) then
	write(*,'(a,i3,a,a)') 'error opening unit', 11, ' -- ', fname_test
	stop('END PROGRAM')
endif
open(unit = 12, file = 'A_mat_g.csv', status = 'replace', action = 'write', iostat = ios)
if (ios .ne. 0) then
	write(*,'(a,i3,a,a)') 'error opening unit', 12, ' -- ', 'A_mat_g.csv'
	stop('END PROGRAM')
endif
open(unit = 13, file = 'f_mat_g.csv', status = 'replace', action = 'write', iostat = ios)
if (ios .ne. 0) then
	write(*,'(a,i3,a,a)') 'error opening unit', 13, ' -- ', 'A_mat_g.csv'
	stop('END PROGRAM')
endif

line_num = 0
total_parts = 0
allocate(dimensions(1))
allocate(temp_dimensions(1))
allocate(material_list(1))
allocate(temp_material_list(1))
do
	line_num = line_num + 1
	read(11,101) line
	pos = scan(line, ' ')
	label = line(1:pos)
	line = line(pos + 1:)
	pos = scan(line,'!')
	if (pos .ne. 0) then
		line = line(1:pos)
	endif

	select case (label)
	case ('!')
		! comment
		cycle
	case ('')
		! blank line
		cycle
	case ('cells')
		read(line,*) cells
	case ('dx')
		read(line,*) dx
	case ('pins')
		read(line,*) pins
		allocate(parts(pins))
		allocate(parts_sum(pins))
	case ('npin')
		read(line,*) npin
		allocate(pinmap(npin))
		pinmap(:) = 0
	case ('pinmap')
		beginning = 1
		do while (pinmap(npin) .eq. 0)
			read(11,*) pinmap(beginning:)
			do i = 1,npin
				if (pinmap(i) .eq. 0) then
					beginning = i
				endif
			enddo
		enddo
	case ('pin')
		read(line,*) pin_num
	case ('parts')
		read(line,*) parts(pin_num)
		total_parts = total_parts + parts(pin_num)
	case ('dimension')
		deallocate(temp_dimensions)
		allocate(temp_dimensions(total_parts - parts(pin_num)))
		temp_dimensions = dimensions
		deallocate(dimensions)
		allocate(dimensions(total_parts))
		dimensions(1:(total_parts - parts(pin_num))) = temp_dimensions
		read(11,*) dimensions(total_parts - parts(pin_num) + 1:total_parts)
	case ('material')
		deallocate(temp_material_list)
		allocate(temp_material_list(total_parts - parts(pin_num)))
		temp_material_list = material_list
		deallocate(material_list)
		allocate(material_list(total_parts))
		material_list(1:(total_parts - parts(pin_num))) = temp_material_list
		read(11,*) material_list(total_parts - parts(pin_num) + 1:total_parts)
	case ('eof')
		write(*,101) 'END OF FILE'
		write(*,101) 
		exit
	case default
		write(*,'(a,i3,2a)') 'unknown input - line #',line_num,' - ',label
	endselect
enddo

parts_sum = 0
do i = 1,pins
	do j = 1,i
		parts_sum(i) = parts_sum(i) + parts(j)
	enddo
enddo
! call vectorwrite(pins,real(parts_sum,8))

regions = 0
do i = 1,npin
	regions = regions + parts(pinmap(i))
enddo
! write(*,'(i4)') regions
allocate(x_coords(regions))
allocate(material_region(regions))

region_counter = 1
do i = 1,npin
	x_coords(region_counter:(region_counter + parts(pinmap(i)) - 1)) =  &
		& x_coords(region_counter - 1) + dimensions((parts_sum(pinmap(i)) - parts(pinmap(i)) + 1):parts_sum(pinmap(i)))
	material_region(region_counter:(region_counter + parts(pinmap(i)) - 1)) = &
		& material_list((parts_sum(pinmap(i)) - parts(pinmap(i)) + 1):parts_sum(pinmap(i)))
	region_counter = region_counter + parts(pinmap(i))
	! write(*,'(2i3)') (parts_sum(pinmap(i)) - parts(pinmap(i)) + 1),parts_sum(pinmap(i))
	! write(*,101) material_list(parts_sum(pinmap(i)) - parts(pinmap(i)) + 1:parts_sum(pinmap(i)))
enddo
capital_H = x_coords(regions)

allocate(mesh(cells))
allocate(material_cells(cells))
allocate(material_index_cells(cells))
mesh(1) = dx
do i = 2,cells
	mesh(i) = mesh(i - 1) + dx
enddo

mesh_tol = 1.0d-3
do i = 1,cells
	do j = 1,regions
		if ((mesh(i) - x_coords(j)) .lt. mesh_tol) then
			material_cells(i) = material_region(j)
			! write(*,101) material_cells(i)
			exit
		endif
	enddo
enddo

material_index_cells = 0
do i = 1,cells
	do j = 1,material
		if (material_cells(i) .eq. names(j)) then
			material_index_cells(i) = j
			! write(*,'(i3)') material_index_cells(i)
			exit
		endif
	enddo
enddo

allocate(t(cells,group))
allocate(tr(cells,group))
allocate(d(cells,group))
allocate(a(cells,group))
allocate(c(cells,group))
allocate(f(cells,group))
allocate(nu(cells,group))
allocate(x(cells,group))
allocate(r(cells,group))
allocate(s(cells,group,group))
do i = 1,cells
	t(i,:) = t_in(material_index_cells(i),:)
	tr(i,:) = tr_in(material_index_cells(i),:)
	d(i,:) = d_in(material_index_cells(i),:)
	a(i,:) = a_in(material_index_cells(i),:)
	c(i,:) = c_in(material_index_cells(i),:)
	f(i,:) = f_in(material_index_cells(i),:)
	nu(i,:) = nu_in(material_index_cells(i),:)
	x(i,:) = x_in(material_index_cells(i),:)
	r(i,:) = r_in(material_index_cells(i),:)
	s(i,:,:) = s_in(material_index_cells(i),:,:)
enddo
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

do while ((phibarerror .gt. epsilon_phi) .or. (kerror .gt. epsilon_k))
	z = z + 1
	write(*,'(i3)') z
	! check for number of iterations
	if (z .eq. 101) then
		stop('failure to converge after 100 iterations')
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
			! write(*,'(2(e12.6,x))') phibar(i,g,(z - 1)), phibar(i,g,z)
			! turn Qs into f vector
			f_mat((3 * i - 2),g) = Q(i,g,z) * dx
		enddo
		
		A_mat_g = A_mat(:,:,g)
		f_mat_g = f_mat(:,g)
		y_mat_g = y_mat(:,g)
		









		do i = 1, (3 * cells + 2)
			do j = 1, (3 * cells + 2)
				write(12,'(e12.6,a)',advance = 'no') A_mat_g(i,j),','
			enddo
			write(12,101) 
		enddo

		do i = 1, (3 * cells + 2)
			write(13,'(e12.6)') f_mat_g(i)
		enddo

		! stop('here')






























		! the slowdown occurs in matrix decomposition
		! call matrixwrite((3 * cells + 2),(3 * cells + 2),A_mat_g)
		call matrixsolve((3 * cells + 2),(3 * cells + 2),A_mat_g,(3 * cells + 2),f_mat_g,y_mat_g)
		write(*,'(i3,a,i3)') z, ' , ', g
		if (isnan(y_mat_g(1))) stop('y is nan')
		y_mat(:,g) = y_mat_g

		do i = 1,(3 * cells + 2)
			! write(*,'(e12.6)') y_mat(i,g)
		enddo

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
			! write(*,'(3(e12.6,x))') phibar(i,g_prime,z), phibar(i,g_prime,(z - 1)), dx
		enddo
	enddo
	write(*,*) numerator
	write(*,*) denominator
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
	phibar(:,:,z) = (phibar(:,:,z) * capital_H) / sum_phibar
	phi(:,:,z)    = (phi(:,:,z)    * capital_H) / sum_phi

	! calculate iteration error for the next iteration
	phibarerror = 0.0d0 ! TO-DO: FIX THIS SHIT! (involves taking the norm)
	kerror = abs(1 - (k(z - 1) / k(z)))

enddo

write(*,101) 'I FOUND K!!!'
write(*,101) '(maybe)'
write(*,'(e12.6)') k(z)


endprogram transportsolve
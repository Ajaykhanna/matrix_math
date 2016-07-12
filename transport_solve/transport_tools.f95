module transporttools
IMPLICIT NONE

contains
subroutine xsread(fname, group, material, names, t, tr, d, a, c, f, nu, x, s, r)
IMPLICIT NONE
	integer,intent(out) :: group, material
	integer :: mat_num, i, j, ios, pos, line_num
	character(3),dimension(:),allocatable,intent(out) :: names
	real(8),dimension(:,:),allocatable,intent(out) :: t, tr, d, a, c, f, nu, x, r
	real(8),dimension(:,:,:),allocatable,intent(out) :: s
	character(80),intent(in) :: fname
	character(100) :: line, label

	101 format(a) ! plain text

	open(unit = 11, file = fname, status = 'old', action = 'read', iostat = ios)
	if (ios .ne. 0) then
		write(*,'(a,i3,a,a)') 'error opening unit', 11, ' -- ', fname
		stop('END PROGRAM')
	endif
	group = 0
	material = 0
	! READING
	line_num = 0
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
		case ('group')
			read(line,*) group
			if (material .ne. 0) then
				call allocatexs(t,tr,d,a,c,f,nu,x,r,s,material,group)
			endif
		case ('material')
			read(line,*) material
			allocate(names(material))
			if (group .ne. 0) then
				call allocatexs(t,tr,d,a,c,f,nu,x,r,s,material,group)
			endif
		case ('mat')
			read(line,*) mat_num
		case ('name')
			read(line,*) names(mat_num)
		case ('t')
			do i = 1,group
				read(11,*) t(mat_num,i)
			enddo
		case ('tr')
			do i = 1,group
				read(11,*) tr(mat_num,i)
			enddo
		case ('d')
			do i = 1,group
				read(11,*) d(mat_num,i)
			enddo
		case ('a')
			do i = 1,group
				read(11,*) a(mat_num,i)
			enddo
		case ('c')
			do i = 1,group
				read(11,*) c(mat_num,i)
			enddo
		case ('f')
			do i = 1,group
				read(11,*) f(mat_num,i)
			enddo
		case ('nu')
			do i = 1,group
				read(11,*) nu(mat_num,i)
			enddo
		case ('x')
			do i = 1,group
				read(11,*) x(mat_num,i)
			enddo
		case ('s')
			do i = 1,group
				read(11,*) s(mat_num,i,:)
			enddo
		case ('eof')
			write(*,101) 'END OF FILE - xsread'
			exit
		case default
			write(*,'(a,i3,3a)') 'unknown input xsread - line #',line_num,' - ',label,line
		endselect
	enddo

	tr(:,:) = t(:,:)
	d(:,:) = 1 / (3 * tr(:,:))

	do i = 1,material
		do j = 1,group
			r(i,j) = t(i,j) - s(i,j,j)
		enddo
	enddo

	close(unit = 11)

endsubroutine xsread

subroutine transportinput(fname,dx,cells,names,material_index_cells,material)
	IMPLICIT NONE
	integer :: ios, line_num, total_parts, pin_num, pos, beginning
	integer :: i
	integer,intent(out) :: cells
	integer :: pins, npin
	integer,dimension(:),allocatable :: pinmap, parts
	! integer,dimension(:),allocatable :: parts!, parts_sum
	real(8),dimension(:),allocatable :: dimensions
	real(8),dimension(:),allocatable :: temp_dimensions
	character(3),dimension(:),allocatable :: temp_material_list
	character(3),dimension(:),allocatable :: material_list
	character(3),dimension(:),allocatable,intent(in) :: names
	integer,intent(in) :: material
	integer,dimension(:),allocatable,intent(out) :: material_index_cells
	real(8) :: dx
	character(100) :: line, label
	character(80)  :: fname

	101 format(a)

	open(unit = 11, file = fname, status = 'old', action = 'read', iostat = ios)
	if (ios .ne. 0) then
		write(*,'(a,i3,a,a)') 'error opening unit', 11, ' -- ', fname
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
			! allocate(parts_sum(pins))
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
			write(*,101) 'END OF FILE - transportinput'
			exit
		case default
			write(*,'(a,i3,2a)') 'unknown input transportinput - line #',line_num,' - ',label
		endselect
	enddo

	close(unit = 11)

	call meshbuild(cells,dx,pins,npin,pinmap,dimensions,material_list,parts,names,material,material_index_cells)

endsubroutine transportinput

subroutine meshbuild(cells,dx,pins,npin,pinmap,dimensions,material_list,parts,names,material,material_index_cells)
	IMPLICIT NONE
	
	integer,intent(in) :: cells, pins, npin, material
	integer,dimension(:),allocatable,intent(in) :: pinmap,parts
	real(8),intent(in) :: dx
	real(8),dimension(:),allocatable,intent(in) :: dimensions
	character(3),dimension(:),allocatable,intent(in) :: material_list,names

	integer,dimension(:),allocatable,intent(out) :: material_index_cells

	integer,dimension(:),allocatable :: parts_sum
	integer :: regions, region_counter, i, j
	real(8),dimension(:),allocatable :: x_coords, mesh
	character(3),dimension(:),allocatable :: material_region, material_cells
	real(8) :: mesh_tol

	allocate(parts_sum(pins))
	parts_sum = 0
	do i = 1,pins
		do j = 1,i
			parts_sum(i) = parts_sum(i) + parts(j)
		enddo
	enddo

	regions = 0
	do i = 1,npin
		regions = regions + parts(pinmap(i))
	enddo
	allocate(x_coords(regions))
	allocate(material_region(regions))

	region_counter = 1
	do i = 1,npin
		x_coords(region_counter:(region_counter + parts(pinmap(i)) - 1)) =  &
			& x_coords(region_counter - 1) + dimensions((parts_sum(pinmap(i)) - parts(pinmap(i)) + 1):parts_sum(pinmap(i)))
		material_region(region_counter:(region_counter + parts(pinmap(i)) - 1)) = &
			& material_list((parts_sum(pinmap(i)) - parts(pinmap(i)) + 1):parts_sum(pinmap(i)))
		region_counter = region_counter + parts(pinmap(i))
	enddo
	! total_length = x_coords(regions)

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
				exit
			endif
		enddo
	enddo

	material_index_cells = 0
	do i = 1,cells
		do j = 1,material
			if (material_cells(i) .eq. names(j)) then
				material_index_cells(i) = j
				exit
			endif
		enddo
	enddo

endsubroutine meshbuild

subroutine xsbuild(material_index_cells,cells,group,t_in,tr_in,d_in,a_in,c_in,f_in,nu_in,x_in,r_in,s_in,t,tr,d,a,c,f,nu,x,r,s)
	IMPLICIT NONE

	integer,dimension(:),allocatable,intent(in) :: material_index_cells
	integer,intent(in) :: cells, group
	real(8),dimension(:,:),allocatable,intent(in) :: t_in,tr_in,d_in,a_in,c_in,f_in,nu_in,x_in,r_in
	real(8),dimension(:,:,:),allocatable,intent(in) :: s_in
	real(8),dimension(:,:),allocatable,intent(out) :: t,tr,d,a,c,f,nu,x,r
	real(8),dimension(:,:,:),allocatable,intent(out) :: s

	integer :: i

	call allocatexs(t,tr,d,a,c,f,nu,x,r,s,cells,group)
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

endsubroutine xsbuild

subroutine allocatexs(t,tr,d,a,c,f,nu,x,r,s,dim1,dim2)
	IMPLICIT NONE
	real(8),dimension(:,:),allocatable,intent(inout) :: t, tr, d, a, c, f, nu, x, r
	real(8),dimension(:,:,:),allocatable,intent(inout) :: s
	integer,intent(in) :: dim1,dim2

	allocate(t(dim1,dim2)) ! total
	allocate(tr(dim1,dim2)) ! transport
	allocate(d(dim1,dim2)) ! diffusion
	allocate(a(dim1,dim2)) ! abssorption
	allocate(c(dim1,dim2)) ! capture
	allocate(f(dim1,dim2)) ! fission
	allocate(nu(dim1,dim2)) ! nu_f
	allocate(x(dim1,dim2))  ! chi
	allocate(r(dim1,dim2)) ! removal
	allocate(s(dim1,dim2,dim2)) ! scattering (2d)

endsubroutine allocatexs


endmodule transporttools
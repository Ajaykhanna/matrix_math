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
			write(*,101) 'END OF FILE'
			write(*,101)
			exit
		case default
			write(*,'(a,i3,3a)') 'unknown input - line #',line_num,' - ',label,line
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

subroutine transportinput(fname,cells,dx,pins,npin,pinmap,dimensions,material_list,parts)
	IMPLICIT NONE
	integer :: ios, line_num, total_parts, pin_num, pos, beginning
	integer :: i
	integer,intent(out) :: cells, pins, npin
	integer,dimension(:),allocatable,intent(out) :: pinmap, parts
	! integer,dimension(:),allocatable :: parts!, parts_sum
	real(8),dimension(:),allocatable,intent(out) :: dimensions
	real(8),dimension(:),allocatable :: temp_dimensions
	character(3),dimension(:),allocatable :: temp_material_list
	character(3),dimension(:),allocatable,intent(out) :: material_list
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
			write(*,101) 'END OF FILE'
			write(*,101) 
			exit
		case default
			write(*,'(a,i3,2a)') 'unknown input - line #',line_num,' - ',label
		endselect
	enddo

	close(unit = 11)

endsubroutine transportinput

subroutine allocatexs(t,tr,d,a,c,f,nu,x,r,s,material,group)
	IMPLICIT NONE
	real(8),dimension(:,:),allocatable,intent(inout) :: t, tr, d, a, c, f, nu, x, r
	real(8),dimension(:,:,:),allocatable,intent(inout) :: s
	integer,intent(in) :: material,group

	allocate(t(material,group)) ! total
	allocate(tr(material,group)) ! transport
	allocate(d(material,group)) ! diffusion
	allocate(a(material,group)) ! abssorption
	allocate(c(material,group)) ! capture
	allocate(f(material,group)) ! fission
	allocate(nu(material,group)) ! nu_f
	allocate(x(material,group))  ! chi
	allocate(r(material,group)) ! removal
	allocate(s(material,group,group)) ! scattering (2d)
endsubroutine allocatexs


endmodule transporttools
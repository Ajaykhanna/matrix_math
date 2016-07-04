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
	102 format(i1) ! single integer

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
			endif
		case ('material')
			read(line,*) material
			allocate(names(material))
			if (group .ne. 0) then
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


endmodule transporttools
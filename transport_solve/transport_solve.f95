program transportsolve
use transporttools
use matrixtools
IMPLICIT NONE

character(80)  :: fname
integer :: group, material
character(3),dimension(:),allocatable :: names
real(8),dimension(:,:),allocatable :: t, tr, d, a, c, f, nu, x, r
real(8),dimension(:,:,:),allocatable :: s

integer :: ios, cells, pins, npin, pin_num, pos, line_num, total_parts, i, beginning
integer,dimension(:),allocatable :: pinmap, parts
real(8),dimension(:),allocatable :: dimension, temp_dimension
character(3),dimension(:),allocatable :: material_list, temp_material_list
real(8) :: dx
character(100) :: line, label
character(80)  :: fname_test

101 format(a) ! plain text descriptor
102 format(i1)

write(*,101) 'input filename'
read(*,*) fname
call xsread(fname,group,material,names,t,tr,d,a,c,f,nu,x,s,r)

! call matrixwrite(material,group,r)


write(*,101) 'enter test input filename'
read(*,*) fname_test

open(unit = 11, file = fname_test, status = 'old', action = 'read', iostat = ios)
if (ios .ne. 0) then
	write(*,'(a,i3,a,a)') 'error opening unit', 11, ' -- ', fname_test
	stop('END PROGRAM')
endif

line_num = 0
total_parts = 0
allocate(dimension(1))
allocate(temp_dimension(1))
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
		deallocate(temp_dimension)
		allocate(temp_dimension(total_parts - parts(pin_num)))
		temp_dimension = dimension
		deallocate(dimension)
		allocate(dimension(total_parts))
		dimension(1:(total_parts - parts(pin_num))) = temp_dimension
		read(11,*) dimension(total_parts - parts(pin_num) + 1:total_parts)
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


endprogram transportsolve
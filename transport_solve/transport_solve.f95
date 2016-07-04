program transportsolve
use transporttools
use matrixtools
IMPLICIT NONE

character(80)  :: fname
integer :: group, material
character(3),dimension(:),allocatable :: names
real(8),dimension(:,:),allocatable :: t, tr, d, a, c, f, nu, x, r
real(8),dimension(:,:,:),allocatable :: s

integer :: ios, cells, pins, npin, pin_num, pos, line_num, total_parts, i
integer,dimension(:),allocatable :: pinmap, parts
real(8),dimension(:),allocatable :: dimension, temp_dimension
character(3),dimension(:),allocatable :: material_list
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
		do while (pinmap(npin) .eq. 0)
			read(11,*) pinmap
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
		! write(*,*) total_parts - parts(pin_num) + 1 , total_parts
		read(11,101) line
		read(line,*) dimension(total_parts - parts(pin_num) + 1:total_parts)
		! do i = 1,total_parts
		! 	write(*,'(e12.6)') dimension(i)
		! enddo
		! write(*,*)
	case ('material')
		read(11,*) material_list(1:parts(pin_num))
	case ('eof')
		write(*,101) 'END OF FILE'
		write(*,101) 
		exit
	case default
		write(*,'(a,i3,2a)') 'unknown input - line #',line_num,' - ',label
	endselect
enddo



endprogram transportsolve
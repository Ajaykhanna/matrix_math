module fileio
IMPLICIT NONE
contains

subroutine search(str,int)
IMPLICIT NONE
	character(:),intent(in) :: str
	integer,intent(in) :: int
	integer :: i

	do i = 1,
		

real(8) :: chi_tol, chi_sum
logical :: chi_check


chi_tol = 1.0d-5
chi_check = .false.
chi_sum
do i = 1,material
	do j = 1,group
		chi_sum = chi_sum + x(i,j)
	enddo
	if (((chi_sum - 1.0d0) .gt. chi_tol) .or. ((chi_sum - 1.0d0) .lt. (-1) * chi_tol)) then
		chi_check = .true.
	endif
enddo
if (chi_check) then
	stop('chi check failed')
endif

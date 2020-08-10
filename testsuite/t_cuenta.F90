program test1
use grid_temp_mod
implicit none
	integer :: iunit
	open(newunit=iunit,file="f_test.csv",status="old",action="read")
	print *,cuenta(iunit)
end program test1

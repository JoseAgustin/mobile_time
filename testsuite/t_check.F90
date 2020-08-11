!  Test for check subroutine
program test_check
use grid_temp_mobile, only :check
integer ::istatus =0

        print *,"Status= ",istatus
	call check(istatus)

end program test_check

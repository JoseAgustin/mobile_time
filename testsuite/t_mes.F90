!  test mes
program test_mes
use grid_temp_mobile, only: mes
implicit none
    integer:: imes
    character(len=2):: cmes
    do imes=1,12
	write(cmes,'(I2.2)') imes
	write(6,'(A3,x)') mes(cmes)
    end do
end program test_mes

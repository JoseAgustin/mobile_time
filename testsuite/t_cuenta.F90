program test1
implicit none
integer :: iunit
integer :: cuenta
open(newunit=iunit,file="f_test.csv",status="old",action="read")

 print *,cuenta(iunit)

end program test1

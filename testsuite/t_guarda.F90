!  program test subroutine guarda_malla
program test_guarda_malla
use grid_temp_mobile
implicit none
  emision=100./24.
  eday= 4.
  nic=28*34
  call system("mkdir -v data")
  call guarda_malla
  call system("rm -rv data")
end program test_guarda_malla

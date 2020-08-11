!  program test subroutine guarda_malla_nc
program test_guarda_malla_nc
use grid_temp_mobile
implicit none
  emision=100./24.
  eday= 4.
  nic=28*34
  call system("ln -s ../data .")
  call lee_atributos
  call guarda_malla_nc
  call system("rm data emission_1990*.nc")
end program test_guarda_malla_nc

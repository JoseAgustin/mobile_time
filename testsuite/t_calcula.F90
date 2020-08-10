!  test calcula_emision
program test_calcula_emision
use grid_temp_mobile
implicit none
    call system("ln -s ../data .")
    call lee_atributos
    call lee_actividades
    call lee_factor_emision
    call calcula_emision
call system("rm data")
end program test_calcula_emision

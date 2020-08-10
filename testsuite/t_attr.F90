!  test lee_atributos
program test_lee_atributos
use grid_temp_mobile, only: lee_atributos,logs
implicit none
    call system("ln -s ../data .")
    call lee_atributos
    call system("rm data")
end program test_lee_atributos

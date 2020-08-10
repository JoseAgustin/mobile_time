!  test lee_actividades
program test_lee_actividades
use grid_temp_mobile, only: lee_actividades,logs
implicit none
    call system("ln -s ../data .")
    call lee_actividades
    call system("rm data")
end program test_lee_actividades

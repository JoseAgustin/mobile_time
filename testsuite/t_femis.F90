!  test lee_factor_emision
program test_lee_factor_emision
use grid_temp_mobile, only: lee_factor_emision,logs
implicit none
    call system("ln -s ../data .")
    call lee_factor_emision
    call system("rm data")
end program test_lee_factor_emision

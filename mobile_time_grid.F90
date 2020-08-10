!   __ _ _ __(_) __| |
!  / _` | '__| |/ _` |
! | (_| | |  | | (_| |
!  \__, |_|  |_|\__,_|
!  |___/           _     _ _    _
!  _ __ ___   ___ | |__ (_) |  | |_ ___ _ __ ___  _ __
! | '_ ` _ \ / _ \| '_ \| | |  | __/ _ \ '_ ` _ \| '_ \
! | | | | | | (_) | |_) | | |  | ||  __/ | | | | | |_) |
! |_| |_| |_|\___/|_.__/|_|_|___\__\___|_| |_| |_| .__/
!                          |_____|               |_|
!>  @brief Program to obtain the temporal distribution over CDMX
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
program grid_mobil_temp
use grid_temp_mobile

  call logs("STARTS RUNNING     ")

  call lee_atributos

  call lee_actividades

  call lee_factor_emision

  call calcula_emision

  call guarda_malla

  call guarda_malla_nc

end program grid_mobil_temp

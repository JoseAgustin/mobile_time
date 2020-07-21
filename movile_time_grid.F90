!     program movile_time_grid.F90
!>  @brief Variables for the identification of temporal profiles in a grid
!>
!>   |  ID_time_period | Time period |
!>   |:---:|:---   |
!>   |  1 | Working Day |
!>   |  2 | Saturday |
!>   |  3 | Sunday   |
!>   | 11 | Week     |
!>   | 12 | Year     |
!>
!>   |  itype  | Description    | itype | Description |
!>   |:---:    |:---            |:---: |:---          |
!>   | 11     | Automoviles     |  15  | Otros buses  |
!>   | 12     | Ligeros         |  16  | Medianos    |
!>   | 13     | Microbuses      |  17  | Pesados     |
!>   | 14     | Ruta 100        |  18  | Camiones Municipales |
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
!>
module grid_temp_mod
!> Number of grids in the mesh
integer,parameter :: n_grid=952
!> Number of rows in results data
integer,parameter :: n_rows=75889
!> Number of viality segments
integer,parameter :: nint=11848
integer idcell(nint),igeo(nint),idsrc(nint),isrc(nint)
!> GRID ID from the used mesh
integer :: grid_id(n_grid)
!> ID of period time
integer :: ID_time_period(n_rows)
!> Source type
integer :: itype(n_rows)
!> Viality  ID from src_td
integer :: isource(n_rows)
!>     autn     Number of cars in a specifil viality and hour         c
real,dimension(24,n_rows) ::autn
!> longitud coordinate for the mesh
real :: long(n_grid)
!> latitude coordinate for the mesh
real :: lat(n_grid)
real cutla(nint),cw(nint)
common /intersec/ cutla,cw,idcell,idsrc,igeo,isrc
common /cellattr/ grid_id,long,lat
common /srctd/ autn, isource, ID_time_period,itype

contains
!> @brief Program to convert UTM coordinates to lat lon.
!> @author SCRAM EPA
!> @date 1990
!>https://www.epa.gov/scram/air-quality-dispersion-modeling-related-model-support-programs#concor
!
  subroutine  utm_2_ll(utmx,utmy,utmz,lat,long)
  implicit none
  integer,intent (IN):: utmZ
  real,intent (IN):: utmx,utmy
  real,intent (OUT) ::lat,long
  real:: utmym,dlong,dlongp,dlongx
  real::  DEGRAD,latx
!
!     THIS SUBROUTINE CONVERTS from UTM coordinates
!
  DATA DEGRAD/0.017453/
      
  latx = utmy / 111.
  dlongx = (utmx - 500.) / ((111.226 + 0.0053 * latx) * &
   (COS (DEGRAD * latx)))
      utmym = utmy - (3187. * SIN (2. * DEGRAD * latx) * &
      (1. -COS (DEGRAD * dlongx)))
      lat = (utmym - 2.41 - 0.00903 * latx * latx) / 110.270
      dlongp = (utmx - 500.) / ((111.226 + 0.0053 * lat) * &
    (COS (DEGRAD * lat)))
      long = -(180 - (6 * utmZ) + 3 - dlongp )
      
  end subroutine utm_2_ll
!>  @brief Reads file grid cell attributes IDgrid, utmx, utmy
!>
!> Converts the UMTx and UTMy to longitude, latitude coordinates.
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine lee_atributos
  implicit none
  integer:: iunit      !  Unit ID for the cell_attr.txt file
  integer :: i,idum
  real :: utmx,utmy
  open(newunit=iunit,file="cell_attr.txt",ACTION="READ")
  do i=1,n_grid
   read (iunit,*)idum,grid_id(i),utmx,utmy
    call utm_2_ll(utmx,utmy,14,lat(i),long(i))
  end do
  close(iunit)
  i=1
  !write(6,*) grid_id(i),lat(i),long(i)
  i=n_grid
  !write(6,*) grid_id(i),lat(i),long(i)
  write(6,*) "Finish reading cell_attr.txt"
end subroutine lee_atributos
!>  @brief Reads for each grid the numbero fo vehicles per category by hour
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine lee_actividades
  implicit none
  integer:: iunit      !  Unit ID for thesrc_td.txt file
  integer :: idum,i,j
  open(newunit=iunit,file="src_td.txt",ACTION="READ")
  do j=1,n_rows
  read (iunit,*)idum,idum,idum,isource(j),&
        ID_time_period(j),itype(j), &
        idum,idum,(autn(i,j),i=1,24)!,(vel(i,j),i=1,24)
  end do
  close(iunit)
  write(6,*) "Finish reading src_td.txt"
  open(newunit=iunit,file="intersection.txt",ACTION="READ")
  do  j=1,nint
   read(iunit,*)idum,idcell(j),isrc(j),igeo(j),idsrc(j), &
             cutla(j),cw(j)
  end do
  write(6,*) "Finish reading intersection.txt"
  close(iunit)

end subroutine lee_actividades
!>  @brief Make a single grid code with vehicular activity
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine genera_malla
  implicit none

end subroutine genera_malla
!>  @brief Stores the mesh in a file
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine guarda_malla

end subroutine guarda_malla

end module grid_temp_mod
!>  @brief Program to obtain the temporal distribution over CDMX
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
program grid_movil_temp
use grid_temp_mod

  call lee_atributos

  call lee_actividades

  call genera_malla

  call guarda_malla
end program

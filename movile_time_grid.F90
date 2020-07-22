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
integer,parameter :: nef = 5, nef2 =5,ntypd= 3, nhr=24
integer,parameter :: nfe = 7  ! EPA
integer,parameter :: nic=952 ;!> Number of rows in results data
integer,parameter :: ntd=75889; !> Number of viality segments
integer,parameter :: nint=11848; !> Number of viality lengths
integer,parameter :: natt=5554!> GRID ID from the used mesh from intersection
integer :: idcell(nint) ;!> GEOmetry ID for viality from intersection
integer :: igeo(nint)   ;!> GEOmetry ID for viality from intersection
integer :: idsrc(nint)  ;!> GEOmetry ID for viality from intersection
integer :: isrc(nint)   ;!> GRID ID from the used mesh
integer :: grid_id(nic) ;!> ID of period time
integer :: ID_time_period(ntd);!> Source type
integer :: itype(ntd)   ;!> Viality  ID from src_td
integer :: isource(ntd)
integer ::idsrc2(natt),igeo2(natt),kstype(natt),vel(nhr,ntd)
!>  Number of cars in a specifil viality and hour
real,dimension(nhr,ntd) :: autn;!> total number of vehicles per hour in each grid
real :: long(nic)    ;!> latitude coordinate for the mesh
real :: lat(nic); !> fraction
real :: cutla(nint) ;!> fraction
real :: cw(nint)   ;!> street length
real :: slen(natt,2) ;!> number of vehicles
real :: scar(nic,8,1:2)
real :: vele(nef,nfe),ehc(nef,nfe),eco(nef,nfe),eno(nef,nfe)
real :: vele2(nef2,7),ehc2(nef2,7),eco2(nef2,7),eno2(nef2,7)
real :: fcor(nhr),start(nhr)

common /intersec/ cutla,cw,idcell,idsrc,igeo,isrc
common /cellattr/ grid_id,long,lat
common /srctd/    autn, isource, ID_time_period,itype
common /facttuv/  vele,ehc,eco,eno!,eso
common /factsec/  vele2,ehc2,eco2,eno2
common /miscell/  scar,fcor,start!,vv,et,fcorr,ffr,emision
common /srcattr/  slen,kstype,igeo2,idsrc2

contains
!> @brief Program to convert UTM coordinates to lat lon.
!>
!>https://www.epa.gov/scram/air-quality-dispersion-modeling-related-model-support-programs#concor
!> @author SCRAM EPA
!> @date 1990
!> @param  utmy Coordinate in axis _y_ in km
!> @param  utmx Coordinate in axis _x_ in km
!> @param  utmz Zone
!> @param  lat Coordinate latitude in decimal degrees
!> @param  long Coordinate longitude in decimal degrees
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
  integer:: iunit      !  Unit ID for the cell_attr.csv file
  integer :: i,idum,j
  real :: utmx,utmy,dum

  open(newunit=iunit,file="data/src_attr.csv",ACTION="READ")
  do  j=1,natt
    read(iunit,*)idum,igeo2(j),idsrc2(j),idum,idum,idum,idum&
    ,kstype(j),slen(j,1),dum,slen(j,2)
  end do
  close(iunit)
  write(6,140)
  open(newunit=iunit,file="data/cell_attr.csv",ACTION="READ")
  do i=1,nic
   read (iunit,*)idum,grid_id(i),utmx,utmy
    call utm_2_ll(utmx,utmy,14,lat(i),long(i))
  end do
  close(iunit)
  i=1
  !write(6,*) grid_id(i),lat(i),long(i)
  i=nic
  !write(6,*) grid_id(i),lat(i),long(i)
  write(6,150)
140 format(9X,'******  END READING src_attr.txt',9X,'******')
150 format(9X,'******  END READING cell_attr.txt',8X,'******')

end subroutine lee_atributos
!>  @brief Reads for each grid the numbero fo vehicles per category by hour
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine lee_actividades
  implicit none
  integer:: iunit      !  Unit ID for thesrc_td.csv file
  integer :: idum,i,j
  open(newunit=iunit,file="data/src_td.csv",ACTION="READ")
  do j=1,ntd
  read (iunit,*)idum,idum,idum,isource(j),&
        ID_time_period(j),itype(j), &
        idum,idum,(autn(i,j),i=1,nhr),(vel(i,j),i=1,nhr)
  end do
  close(iunit)
  write(6,160)
  open(newunit=iunit,file="data/intersection.csv",ACTION="READ")
  do  j=1,nint
   read(iunit,*)idum,idcell(j),isrc(j),igeo(j),idsrc(j), &
             cutla(j),cw(j)
  end do
  write(6,150)
  close(iunit)
150 format(9X,'******  END READING intersection.txt',5X,'******')
160 format(9X,'******  END READING src_td.csv',11X,'******')
end subroutine lee_actividades
!>  @brief Reads emissions factor from EPA, and for cool start
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine lee_factor_emision
implicit none
integer:: iunit      !  Unit ID for file reading
integer :: i,j
character*25  title
open(newunit=iunit,file="data/factepa.txt",ACTION="READ")

!..
!    -----------  Reading    factepa.txt      unit  15    ----------
!..
  do  j=1,nef
    read(iunit,'(a25)')title
    read(iunit,'(a2)')title
    read(iunit,'(a2)')title
    read(iunit,'(a2)')title
    do  i =1,nfe
      read(iunit,*)vele(j,i),ehc(j,i),eco(j,i),eno(j,i) !,eso(j,i)
    end do
  end do
  close (iunit)
  write(6,130)
!..
!    -----------  Reading    factsec.txt     cold start    ----------
!..
open(newunit=iunit,file="data/factsec.txt",ACTION="READ")

  do  j=1,5
    read(iunit,'(a25)')title
    read(iunit,'(a2)')title
    read(iunit,'(a2)')title
    read(iunit,'(a2)')title
    do  i =1,7
      read(iunit,*)vele2(j,i),ehc2(j,i),eco2(j,i),eno2(j,i)
    end do
  end do
  close (iunit)
  write(6,140)
!    -----------  Reading    factvar.dat        unit  19  ----------
!..
open(newunit=iunit,file="data/factvar.dat",ACTION="READ")

  read(iunit,'(a25)')title
  do  i=1,nhr
    read(iunit,*)fcor(i)
  end do
  close (iunit)
  write(6,150)
!    -----------  Reading    fraarran.dat      unit  18   ----------
!..
open(newunit=iunit,file="data/fraarran.dat",ACTION="READ")
  read(iunit,'(a25)') title
  do  i=1,nhr
    read(iunit,*) start(i)
  end do
  close (iunit)
  write(6,160)
!..
130 format(9X,'******  END READING factepa.txt',10X,'******')
140 format(9X,'******  END READING factsec.txt',10X,'******')
150 format(9X,'******  END READING fraarran.dat',9X,'******')
160 format(9X,'******  END READING factvar.dat',10X,'******')

end subroutine lee_factor_emision
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

  call lee_factor_emision

  call genera_malla

  call guarda_malla
end program

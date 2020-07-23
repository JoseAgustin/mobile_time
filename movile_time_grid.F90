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
!>   |  nef    | Description |
!>   |:---:    |:---             |
!>  | 1   | Automovil    |
!>  | 2   | Ligeros   |
!>  | 3  |  Microbuses   |
!>  | 4  |  Ruta 100, Pesados & Municipales TUV   |
!>  | 5   | Otros buses   |
!>  | 6  |  Medianos   |
!>
!>   |  itype  | Description    | itype | Description |
!>   |:---:    |:---            |:---: |:---          |
!>   | 11     | Automoviles     |  15  | Otros buses  |
!>   | 12     | Ligeros         |  16  | Medianos    |
!>   | 13     | Microbuses      |  17  | Pesados     |
!>   | 14     | Ruta 100        |  18  | Camiones Municipales |
!>
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
!>
module grid_temp_mod
!> number species in emission factors from EPA
integer,parameter :: nef = 5 ;!> number species in emission factors cold start from EPA
integer,parameter :: nef2 =5 ;!> number of period 1= wee, 2=saturday 30sunday
integer,parameter :: ntypd= 3 ;!> number of hours per day
integer,parameter :: nhr=24  ; !> number species in emission file
integer,parameter :: nspc = 5 ; !> number of speeds per EF specie
integer,parameter :: nfe = 7  ; !> Number of grids in the mesh
integer,parameter :: nic=28*34;!> Number of rows in results data 952 
integer,parameter :: ntd=75889; !> Number of viality segments
integer,parameter :: nint=11848; !> Number of viality lengths
integer,parameter :: natt=5554 ; !> cell ID from the used mesh from intersection
integer :: idcell(nint) ;!> Source Type (1 line 2 Area) from intersection
integer :: igeo(nint)   ;!> ID for viality segment from intersection
integer :: idsrc(nint)  ;!> GEOmetry ID for viality from intersection
integer :: isrc(nint)   ;!> GRID ID from the used mesh
integer :: grid_id(nic) ;!> ID of period time
integer :: ID_time_period(ntd);!> Source type
integer :: itype(ntd)   ;!> Viality  ID from src_td
integer :: isource(ntd) ;!> Source identification from file src_attr.txt
integer :: idsrc2(natt) ;!> Type of source (1 Line  2 Area) from file src_attr.txt
integer :: igeo2(natt)  ;!> Source classification
integer :: kstype(natt)
!>  Number of cars in a specifil viality and hour
real,dimension(nhr,ntd) :: autn;!> total number of vehicles per hour in each grid
real :: long(nic)    ;!> latitude coordinate for the mesh
real :: lat(nic); !> Cut-leng or cut-area of source segment
real :: cutla(nint) ;!> Relative weight of the viality in the grid
real :: cw(nint)   ;!> Linear source lenght km
real :: slen(natt,2) ;!> number of vehicles
real :: scar(nic,8,1:2) ;!> velocity in emissions factos
real :: vele(nef,nfe)  ;!> VOC emission factor
real :: ehc(nef,nfe)  ;!> CO emission factor
real :: eco(nef,nfe)  ;!> NO emission factor
real :: eno(nef,nfe)  ;!> velocity in emissions factos cold start
real :: vele2(nef2,7) ;!> VOC emission factor cold start
real :: ehc2(nef2,7)  ;!> CO emission factor cold start
real :: eco2(nef2,7)  ;!> NO emission factor cold start
real :: eno2(nef2,7)  ;!> Correction factor for start mode
real :: fcor(nhr)     ;!> Fraction of cars that starts cold
real :: start(nhr)    ;!>Emission factor at speed vel
real :: efact(nspc)   ;!>Emission factor cold start at speed vel
real :: efac2(nspc)   ;!> Street velocity from src_td
real :: vel(nhr,ntd)  ;!> Total emission per day
real :: eday(nic,nspc,ntypd,8,1:2); !> Movil emision per cell,hour,specie,day type
real :: emision(nic,nhr,nspc,ntypd)

common /intersec/ cutla,cw,idcell,idsrc,igeo,isrc
common /cellattr/ grid_id,long,lat
common /srctd/    autn, isource, ID_time_period,itype
common /facttuv/  vele,ehc,eco,eno!,eso
common /factsec/  vele2,ehc2,eco2,eno2
common /miscell/  scar,fcor,start,emision!,vv,et,fcorr,ffr
common /srcattr/  slen,kstype,igeo2,idsrc2
common /computs/ efact,efact2,eday
 
contains
!        _                 ____      _ _
!  _   _| |_ _ __ ___     |___ \    | | |
! | | | | __| '_ ` _ \      __) |   | | |
! | |_| | |_| | | | | |    / __/    | | |
!  \__,_|\__|_| |_| |_|___|_____|___|_|_|
!                    |_____|   |_____|
!> @brief Program to convert UTM coordinates to lat lon.
!>
!>https://www.epa.gov/scram/air-quality-dispersion-modeling-related-model-support-programs#concor
!> @author SCRAM EPA
!> @date 1990
!> @param  utmy Coordinate in axis _y_ in km
!> @param  utmx Coordinate in axis _x_ in km
!> @param  utmz Zone for the UTMx and UTMy coordinates
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
!   _
!  | | ___  ___
!  | |/ _ \/ _ \
!  | |  __/  __/
!  |_|\___|\___|
!       _        _ _           _
!  __ _| |_ _ __(_) |__  _   _| |_ ___  ___
! / _` | __| '__| | '_ \| | | | __/ _ \/ __|
!| (_| | |_| |  | | |_) | |_| | || (_) \__ \
! \__,_|\__|_|  |_|_.__/ \__,_|\__\___/|___/
!
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
!  _
! | | ___  ___
! | |/ _ \/ _ \
! | |  __/  __/
! |_|\___|\___|
!             _   _       _     _           _
!   __ _  ___| |_(_)_   _(_) __| | __ _  __| | ___  ___
!  / _` |/ __| __| \ \ / / |/ _` |/ _` |/ _` |/ _ \/ __|
! | (_| | (__| |_| |\ V /| | (_| | (_| | (_| |  __/\__ \
!  \__,_|\___|\__|_| \_/ |_|\__,_|\__,_|\__,_|\___||___/
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
!  _
! | | ___  ___
! | |/ _ \/ _ \
! | |  __/  __/
! |_|\___|\___|
!  __            _                            _     _
! / _| __ _  ___| |_ ___  _ __  ___ _ __ ___ (_)___(_) ___  _ __
!| |_ / _` |/ __| __/ _ \| '__|/ _ \ '_ ` _ \| / __| |/ _ \| '_ \
!|  _| (_| | (__| || (_) | |  |  __/ | | | | | \__ \ | (_) | | | |
!|_|  \__,_|\___|\__\___/|_|___\___|_| |_| |_|_|___/_|\___/|_| |_|
!                         |_____|
!>  @brief Reads emissions factor from EPA, and for cool start
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine lee_factor_emision
implicit none
integer:: iunit      !  Unit ID for file reading
integer :: i,j
character(len=25)::  header
open(newunit=iunit,file="data/factepa.txt",ACTION="READ")

!..
!    -----------  Reading    factepa.txt      unit  15    ----------
!..
  do  j=1,nef
    read(iunit,'(a25)')header
    read(iunit,'(a2)')header
    read(iunit,'(a2)')header
    read(iunit,'(a2)')header
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
    read(iunit,'(a25)')header
    read(iunit,'(a2)')header
    read(iunit,'(a2)')header
    read(iunit,'(a2)')header
    do  i =1,7
      read(iunit,*)vele2(j,i),ehc2(j,i),eco2(j,i),eno2(j,i)
    end do
  end do
  close (iunit)
  write(6,140)
!    -----------  Reading    factvar.dat        unit  19  ----------
!..
open(newunit=iunit,file="data/factvar.dat",ACTION="READ")

  read(iunit,'(a25)')header
  do  i=1,nhr
    read(iunit,*)fcor(i)
  end do
  close (iunit)
  write(6,150)
!    -----------  Reading    fraarran.dat      unit  18   ----------
!..
open(newunit=iunit,file="data/fraarran.dat",ACTION="READ")
  read(iunit,'(a25)') header
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
!                                                     _ _
!   __ _  ___ _ __   ___ _ __ __ _     _ __ ___   __ _| | | __ _
!  / _` |/ _ \ '_ \ / _ \ '__/ _` |   | '_ ` _ \ / _` | | |/ _` |
! | (_| |  __/ | | |  __/ | | (_| |   | | | | | | (_| | | | (_| |
!  \__, |\___|_| |_|\___|_|  \__,_|___|_| |_| |_|\__,_|_|_|\__,_|
!  |___/                         |_____|
!>  @brief Make a single array with emissions
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine genera_malla
  implicit none
  integer:: l,i,n,m
  integer:: indx
  real :: sl,fcorr,ffr
  real :: temp
  do n = ntd ,1,-1                   !Main LOOP initialization
    do m = 1,nint
      if(isource(n).eq.idsrc(m).and.isrc(m).lt.2 .and.igeo(m).lt.2) then
      indx = idcell(m)
      do l=1,nhr
      ! for EPA
      efact(1)= emisfac2(itype(n),vel(l,n),vele ,ehc )
      efact(2)= emisfac2(itype(n),vel(l,n),vele ,eco )
      efact(3)= emisfac2(itype(n),vel(l,n),vele ,eno )
      ! Emissions factor for cold start
      efac2(1)= emisfac2(itype(n),vel(l,n),vele2,ehc2)
      efac2(2)= emisfac2(itype(n),vel(l,n),vele2,eco2)
      efac2(3)= emisfac2(itype(n),vel(l,n),vele2,eno2)
      if (itype(n).ge. 14 ) then
        efact(4) =efact(1)
        efac2(4) =efac2(1)
        efact(1) = 0.0
        efac2(1) = 0.0
      else
        efact(4) = 0.0
        efac2(4) = 0.0
      end if
!..
!     ----------   Localization of the viality lenght     ----------
!..
       call viality(igeo(m),idsrc(m),igeo2,idsrc2,kstype,slen &
                ,fcor,start,l,fcorr,ffr,sl)
       sl = sl*cw(m)
!     ----------   Computation of the emissions for each specie
        do i = 1,nspc ! 1 HC(gasoline), 2 CO, 3 NOx, 4 HC(Diesel),5 SO2
!..
          temp =(autn(l,n)*sl*efact(i)*(1.0-ffr)+ &
                     autn(l,n)*sl*efac2(i)*ffr)*fcorr
!..        Xing  is 7.5% of the emission of segment viality
         if (isrc(m).eq.2) temp=temp*0.075
!..
            emision(indx,l,i,ID_time_period(n))= temp/3600.0 &
                              + emision(indx,l,i,ID_time_period(n))
!..
!    ----------  Computation of dayly emissions  EDAY     ----------
!
            eday(indx,i,ID_time_period(n),itype(n)-10,2)= temp/1000 &
                      + eday(indx,i,ID_time_period(n),itype(n)-10,2)
        end do    ! i nspc
      end do      ! m nint
    end if
    end do        ! ntd
  end do
end subroutine genera_malla
!                            _                           _ _
!   __ _ _   _  __ _ _ __ __| | __ _     _ __ ___   __ _| | | __ _
!  / _` | | | |/ _` | '__/ _` |/ _` |   | '_ ` _ \ / _` | | |/ _` |
! | (_| | |_| | (_| | | | (_| | (_| |   | | | | | | (_| | | | (_| |
!  \__, |\__,_|\__,_|_|  \__,_|\__,_|___|_| |_| |_|\__,_|_|_|\__,_|
!  |___/                           |_____|
!>  @brief Stores the mesh in a file
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine guarda_malla
implicit none
integer :: i,j,l,iday,irec
integer :: iunit
  write(6,180)
    open (newunit=iunit,file='data/movil.dat', &
    status='unknown',access='direct',form='unformatted' &
    ,recl=nic*4)
    irec = 0
    do iday=1,ntypd
      do l = 1,nhr
        do i = 1,nspc
          irec = irec +1
          write(iunit,rec=irec)(emision(j,l,i,iday),j=1,nic)
        end do   !  i specie
      end do      !  l
    end do
180 format(7X,'      Wrinting output file for GrADS')
end subroutine guarda_malla
!                            _
!   __ _ _   _  __ _ _ __ __| | __ _
!  / _` | | | |/ _` | '__/ _` |/ _` |
! | (_| | |_| | (_| | | | (_| | (_| |
!  \__, |\__,_|\__,_|_|  \__,_|\__,_|
!  |___/          _ _
! _ __ ___   __ _| | | __ _     _ __   ___
!| '_ ` _ \ / _` | | |/ _` |   | '_ \ / __|
!| | | | | | (_| | | | (_| |   | | | | (__
!|_| |_| |_|\__,_|_|_|\__,_|___|_| |_|\___|
!                         |_____|
!>  @brief Stores the mesh in a netcdf file
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine guarda_malla_nc
use netcdf
implicit none
integer, parameter :: NDIMS=6,nx=28,ny=34, zlev=1
integer :: i,j,k,l,iday,ispc,ncid,it
integer :: iit
integer :: dimids2(2),dimids3(3),dimids4(4)
integer :: id_unlimit ;!>latitude ID in netcdf file
integer ::id_varlat ;!>longitude ID in netcdf file
integer ::id_varlong ;!>pollutant emission ID in netcdf file
integer :: id_var(nspc)
integer,dimension(NDIMS):: dim,id_dim ;!> longitude coordinates
real,dimension(nx,ny)::xlong ;!> latitude coordinates
real,dimension(nx,ny)::xlat  ;!> emissions values
real,dimension(nx,ny,1,1)::emis
character (len=19),dimension(NDIMS) ::sdim
character(len=19) :: current_date
character(len= 8) :: date
character(len=10) :: time
character(len=24) :: hoy
character(len=20) :: FILE_NAME
character(len=19),dimension(1,1)::Times
character(len=11),dimension(5)::ename ;!> Emissions long name
character(len=26),dimension(5):: cname
data sdim /"Time               ","DateStrLen         ","west_east          ",&
&          "south_north        ","bottom_top         ","emissions_zdim_stag"/
ename=(/'E_VOC       ','E_CO        ','E_NOx       ',&
        'E_VOC_diesel','E_SO2       '/)
cname=(/'VOC from gasoline vehicles','Carbon Monoxide emissions ', &
        'Nitrogen oxides emissions ','VOC from diesel vehicles  ', &
        'Sulfur dioxide emissions  '/)
  write(6,180)
  FILE_NAME="emission_CDMX1990.nc"
  call check( nf90_create(path =FILE_NAME,cmode = NF90_CLASSIC_MODEL,ncid = ncid) )
  !     Define dimensiones
  xlong=reshape(long,(/nx,ny/))
  xlat=reshape(lat,(/nx,ny/))
  call date_and_time(date,time)
  hoy=date(7:8)//'-'//mes(date(5:6))//'-'//date(1:4)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
  print *,hoy

  dim=(/1,19,nx,ny,1,zlev/)
  current_date="1990-01-12_00:00:00"
  call check( nf90_def_dim(ncid,sdim(1), NF90_UNLIMITED, id_dim(1)) )
  do i=2,NDIMS
      call check( nf90_def_dim(ncid, sdim(i), dim(i), id_dim(i)) )
  end do
  dimids2 = (/id_dim(2),id_dim(1)/)
  dimids3 = (/id_dim(3),id_dim(2),id_dim(1) /)
  dimids4 = (/id_dim(3),id_dim(4),id_dim(6),id_dim(1)/)

  write(6,181)
  call check( nf90_put_att(ncid, NF90_GLOBAL, "TITLE","Emissions from TUV study 1990"))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "START_DATE",current_date))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION",nx))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION",ny))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "BOTTOM-TOP_GRID_DIMENSION",1))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "DX",2*1000))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "DY",2*1000))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "CEN_LAT",xlat(nx/2,ny/2)))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "CEN_LON",xlong(nx/2,ny/2)))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "TRUELAT1",17.5))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "TRUELAT2",29.5))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "MOAD_CEN_LAT",xlat(nx/2,ny/2)))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "STAND_LON",xlong(nx/2,ny/2)))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "POLE_LAT",90.))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "POLE_LON",0.))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "GRIDTYPE","C"))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "GMT",12.))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "JULYR",1990))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "JULDAY",1))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "MAP_PROJ",1))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "MMINLU","USGS"))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "MECHANISM","None"))
  call check( nf90_put_att(ncid, NF90_GLOBAL, "CREATION_DATE",hoy))
!  Define las variables
  call check( nf90_def_var(ncid, "Times", NF90_CHAR, dimids2,id_unlimit ) )
  call check( nf90_def_var(ncid, "XLONG", NF90_REAL,(/id_dim(3),id_dim(4),id_dim(1)/),id_varlong) )
  call check( nf90_def_var(ncid, "XLAT", NF90_REAL,(/id_dim(3),id_dim(4),id_dim(1)/),id_varlat ) )

! Assign  attributes
  call check( nf90_put_att(ncid, id_varlong, "FieldType", 104 ) )
  call check( nf90_put_att(ncid, id_varlong, "MemoryOrder", "XYZ") )
  call check( nf90_put_att(ncid, id_varlong, "description", "LONGITUDE, WEST IS NEGATIVE") )
  call check( nf90_put_att(ncid, id_varlong, "units", "degree_east"))
  call check( nf90_put_att(ncid, id_varlong, "axis", "X") )
  call check( nf90_put_att(ncid, id_varlat, "FieldType", 104 ) )
  call check( nf90_put_att(ncid, id_varlat, "MemoryOrder", "XYZ") )
  call check( nf90_put_att(ncid, id_varlat, "description", "LATITUDE, SOUTH IS NEGATIVE") )
  call check( nf90_put_att(ncid, id_varlat, "units", "degree_north"))
  call check( nf90_put_att(ncid, id_varlat, "axis", "X") )
!  Attributos para cada variable
  do i=1,nspc
   call crea_attr(ncid,4,dimids4,ename(i),cname(i),"g km^-2 s^-1",id_var(i))
  end do
!   Terminan definiciones
  call check( nf90_enddef(ncid) )
  do iday=1,ntypd
  write(6,182) iday
  write(current_date(09:10),'(I2.2)') iday+4
  tiempo: do it=1,nhr
    iit=it+24*(iday-1)
    write(current_date(12:13),'(I2.2)') it-1
    Times(1,1)=current_date(1:19)
    write(6,'(A,x,I2.2)')'TIEMPO: ', iit
    call check( nf90_put_var(ncid, id_unlimit,Times,start=(/1,iit/)) )
    call check( nf90_put_var(ncid, id_varlong,xlong,start=(/1,1,iit/)) )
    call check( nf90_put_var(ncid, id_varlat,xlat,start=(/1,1,iit/)) )
   
    do ispc=1,nspc
      do i=1,nx
        do j=1,ny
          k=i+28*(j-1)
          emis(i,j,1,1)=emision(k,it,ispc,iday)
        end do
      end do
      call check( nf90_put_var(ncid, id_var(ispc),emis,start=(/1,1,1,iit/)))
    end do
   end do TIEMPO
  end do
  call check( nf90_close(ncid) )
180 format(7X,'      Wrinting in output file for netcdf')
181 format(7X,'      Atributos Globales NF90_GLOBAL')
182 format(7X,'      Guarda variables dia: ',I2.2)
end subroutine guarda_malla_nc
!        _       _ _ _
! __   _(_) __ _| (_) |_ _   _
! \ \ / / |/ _` | | | __| | | |
!  \ V /| | (_| | | | |_| |_| |
!   \_/ |_|\__,_|_|_|\__|\__, |
!                        |___/
!>  @brief Localization of the viality and depending of his type it
!>  assings  the value of start or fcorr, and longitud
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/21/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
!>  @param ig Source Type (1 line 2 Area)
!>  @param id viality source identification
!>  @param ig2 Source Type (1 line 2 Area) from intersections
!>  @param isrc viality source identification from intersections
!>  @param kstype Source classification
!>  @param slen  Viality length in km
!>  @param fcor correction factor for start mode
!>  @param start car fraction staring cool
!>  @param nh  time hour for set the EF
!>  @param fcorr  correction factor (out)
!>  @param ffr  correction factor (out)
!>  @param sl Viality length in km (out)
  Subroutine viality(ig,id,ig2,isrc,kstype,slen &
                    ,fcor,start,nh,fcorr,ffr,sl)
  integer,intent(in):: ig,id,isrc(:),kstype(:)
  integer,intent(in):: ig2(:),nh
  real,intent(in) :: fcor(:),start(:)
  real,intent(in) :: slen(:,:)
  real,intent(out):: fcorr,ffr, sl
  integer:: i,flag
!..
  flag=0
  do i = size(isrc),1,-1
    if(id .eq.isrc(i) .and. ig.eq.ig2(i)) then
      sl = slen(i,ig)
      if(kstype(i).gt.10) then
        flag = 1
        fcorr = fcor(nh)
        ffr   = start(nh)
        if(kstype(i).gt.20) then
          ffr= start(nh)
          fcorr =1.0
          exit
        end if
        exit
      else
        flag=1
        fcorr = 1.0
        ffr   = 0.0
        exit
      end if
    end if
  end do
!..
  if (flag .eq. 0 )  then
    write(6,200) id
    stop
  end if
  if(sl.eq.0) then
    write(6,201) sl,ig,slen(i,1),slen(i,2),i
    stop
  end if
!      for area sources is 10% of their cutleng.
  if(ig.eq.2) sl =sl*0.10
 200      format('Invalid Cell',I6)
 201      format('Invalid Longitud',f4.0,'At igeo ',I3,&
                'len(1)=',f6.4,'Len(2)=',f6.4,I4)
  return
!*********************************************************************
!*********             END OF SUBROUTINE VIALITY             *********
!*********************************************************************
  end
!                 _      __            ____
!   ___ _ __ ___ (_)___ / _| __ _  ___|___ \
!  / _ \ '_ ` _ \| / __| |_ / _` |/ __| __) |
! |  __/ | | | | | \__ \  _| (_| | (__ / __/
!  \___|_| |_| |_|_|___/_|  \__,_|\___|_____|
!
!>  @brief Emission factor computation for velociity and EF array
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/22/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
!>  @param ncartype Type of vehicle
!>  @param velocity speed in viality
!>  @param vem velocity array from emission factors file
!>  @param comp Emission factor for vel and specie
  real function emisfac2(ncartype,velocity,vem,comp)
  integer :: ncartype
  integer :: icar, i
  real:: velocity,vem(:,:),comp(:,:)
!    ncartype  Type of vehicle   icar    Type of vehicle
!      11       Automoviles       1      Vehiculos ligeros a Gas
!      12       Ligeros           2      Camionetas ligeras a Gaso
!      13       Microbuses        2      Camionetas ligeras a Gaso
!      13       Microbuses        5      Camiones Pesados a Gasoli
!      14       Ruta 100          3      Camiones ligeros a diesel
!      15       Otros camiones    3      Camiones ligeros a diesel
!      16       Medianos          3      Camiones ligeros a diesel
!      17       Pesados           4      Vehiculos pesados a diesel
!      18       Camiones mpales.  3      Camiones ligeros a diesel
!
  icar = 0
  if( ncartype .eq. 11) icar =1
  if( ncartype .eq. 12) icar =2
  if( ncartype .eq. 13) icar =5
  if( ncartype .eq. 14) icar =3
  if( ncartype .eq. 15) icar =3
  if( ncartype .eq. 16) icar =3
  if( ncartype .eq. 17) icar =4
  if( ncartype .eq. 18) icar =3
  if (icar .eq. 0) then
    write(6,*)'Invalid car type',ncartype
    stop
  end if
  !..
  i=0
  if(velocity.le.vem(icar,2)) i=1
  if(velocity.gt.vem(icar,2).and.velocity.le.vem(icar,3)) i=2
  if(velocity.gt.vem(icar,3).and.velocity.le.vem(icar,4)) i=3
  if(velocity.gt.vem(icar,4).and.velocity.le.vem(icar,5)) i=4
  if(velocity.gt.vem(icar,5).and.velocity.le.vem(icar,6)) i=5
  if(velocity.gt.vem(icar,6)) i=5
  if (i.eq. 0) then
    write(6,200)velocity
    stop
  end if
!..
    emisfac2= comp(icar,i) +(velocity-vem(icar,i)) * &
        (comp(icar,i+1)-comp(icar,i))/(vem(icar,i+1)-vem(icar,i))
!..
200  format('Invalid velocity',E10.1)
  return
!*********************************************************************
!*********             END OF FUNCTION EMISFAC2              *********
!*********************************************************************
  end
!       _               _
!   ___| |__   ___  ___| | __
!  / __| '_ \ / _ \/ __| |/ /
! | (__| | | |  __/ (__|   <
!  \___|_| |_|\___|\___|_|\_\
!>  @brief Verifies no error in netcdf function call
!>  @param status NetCDF functions return a non-zero status codes on error.
subroutine check(status)
use netcdf
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 2
    end if
end subroutine check
!  _ __ ___   ___  ___
! | '_ ` _ \ / _ \/ __|
! | | | | | |  __/\__ \
! |_| |_| |_|\___||___/
!
!>  @brief Returns the month in characters from month number
!>   @author  Jose Agustin Garcia Reynoso
!>   @date  07/13/2020
!>   @version  2.2
!>   @copyright Universidad Nacional Autonoma de Mexico 2020
!>   @param  num number of the month
character(len=3)function mes(num)
    character*2 num
    select case (num)
    case('01');mes='Jan'
    case('02');mes='Feb'
    case('03');mes='Mar'
    case('04');mes='Apr'
    case('05');mes='May'
    case('06');mes='Jun'
    case('07');mes='Jul'
    case('08');mes='Aug'
    case('09');mes='Sep'
    case('10');mes='Oct'
    case('11');mes='Nov'
    case('12');mes='Dec'
    case default
        print *,"   **************************"
        print *,"Month:",num," does not exists!!"
        print *,"   **************************"
        stop  "End program, review namelist_emiss.nml"
    end select
    return
end function
!   ___ _ __ ___  __ _     __ _| |_| |_ _ __
!  / __| '__/ _ \/ _` |   / _` | __| __| '__|
! | (__| | |  __/ (_| |  | (_| | |_| |_| |
!  \___|_|  \___|\__,_|___\__,_|\__|\__|_|
!                    |_____|
!>  @brief Creates attributes for each variable in the netcdf file
!>   @author  Jose Agustin Garcia Reynoso
!>   @date  07/13/2020
!>   @version  2.2
!>   @copyright Universidad Nacional Autonoma de Mexico 2020
!>   @param ncid netcdf file ID
!>   @param idm number of items in dimids array
!>   @param dimids ID dimensons array
!>   @param svar variable name
!>   @param cname description variable name
!>   @param cunits units of the variable
!>   @param id_var variable ID
subroutine crea_attr(ncid,idm,dimids,svar,cname,cunits,id_var)
use netcdf
    implicit none
    integer , INTENT(IN) ::ncid,idm
    integer, INTENT(out) :: id_var
    integer, INTENT(IN),dimension(idm):: dimids
    character(len=*), INTENT(IN)::svar,cname,cunits
    character(len=50) :: cvar
    cvar="Emissions rate of "//trim(cname)

    call check( nf90_def_var(ncid, svar, NF90_REAL, dimids,id_var ) )
    ! Assign  attributes
    call check( nf90_put_att(ncid, id_var, "FieldType", 104 ) )
    call check( nf90_put_att(ncid, id_var, "MemoryOrder", "XYZ") )
    call check( nf90_put_att(ncid, id_var, "description", Cvar) )
    call check( nf90_put_att(ncid, id_var, "units", cunits))
    call check( nf90_put_att(ncid, id_var, "stagger", "Z") )
    call check( nf90_put_att(ncid, id_var, "coordinates", "XLONG XLAT") )
    ! print *,"Entro a Attributos de variable",dimids,id,jd
    return
end subroutine crea_attr
end module grid_temp_mod
!             _     _
!   __ _ _ __(_) __| |
!  / _` | '__| |/ _` |
! | (_| | |  | | (_| |
!  \__, |_|  |_|\__,_|
!  |___/                _ _    _
!  _ __ ___   _____   _(_) |  | |_ ___ _ __ ___  _ __
! | '_ ` _ \ / _ \ \ / / | |  | __/ _ \ '_ ` _ \| '_ \
! | | | | | | (_) \ V /| | |  | ||  __/ | | | | | |_) |
! |_| |_| |_|\___/ \_/ |_|_|___\__\___|_| |_| |_| .__/
!                         |_____|               |_|
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

  call guarda_malla_nc

end program

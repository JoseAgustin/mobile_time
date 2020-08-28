!     module mod_mobile_time_grid.F90
!>  @brief Variables for the identification of temporal profiles in a grid
!>
!>   |  ID_time_period | Time period |
!>   |:---:            |:---   |
!>   |  1 | Working Day|
!>   |  2 | Saturday |
!>   |  3 | Sunday   |
!>   | 11 | Week     |
!>   | 12 | Year     |
!>
!>   |  nef    | Description |
!>   |:---:    |:---         |
!>   | 1   | Automovil       |
!>   | 2   | Ligeros         |
!>   | 3   | Microbuses      |
!>   | 4   | Ruta 100, Pesados & Municipales TUV   |
!>   | 5   | Otros buses     |
!>   | 6   | Medianos        |
!>
!>   |  veh_type |SCC  |Description | veh_type|SCC| Description |
!>   |:---:|:---       |:---        |:---: |:---  | :---        |
!>   | 11  |2201001330 |Automoviles |  15  |2230075330| Otros buses |
!>   | 12  |2201040330 |Ligeros     |  16  |2201040330| Medianos    |
!>   | 13  |2230070270 |Microbuses  |  17  |2230074330| Pesados     |
!>   | 14  |2230001000 |Ruta 100    |  18  |2230001000| Camiones Municipales |
!>
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
!>
module grid_temp_mobile
!> number species in emission factors from EPA
integer,parameter :: nef   =5
!> number species in emission factors cold start from EPA
integer,parameter :: nef2  =5 ;!> Type of day 1 week day, 2 Saturday, 3 Sunday
integer,parameter :: ntypd =3 ;!> number of hours per day
integer,parameter :: nhr   =24;!> number species in emission file
integer,parameter :: nspc  = 5;!> number of speeds per EF specie
integer,parameter :: nfe   = 7;!> number of vehicle types 11 to 18
integer,parameter :: nveht = 8  ;!> Number of grids in the mesh 952
integer :: nic !=28*34
!> Number of rows in results data
integer,parameter :: ntd=75889;!> Number of viality segments
integer,parameter :: nint=11848;!> Number of viality lengths
integer :: natt=5554 ;!> Grid ID from the used mesh from intersection
integer :: id_grid_INT(nint)
!> Source type geometry (1 line 2 Area) from intersection
integer :: geometry_type(nint) ;!> Viality segment ID from intersection
integer :: id_source_INT(nint)
!> Geometry source type from intersection 2nd viality or area
integer :: geo_type_INT(nint)  ;!> GRID ID for each cell in the emissions mesh
integer,allocatable  :: id_grid(:) ;!> Period time ID 1-weekday 2-saturday 3-Sunday
integer :: ID_time_period(ntd);!> Vehicular source type ID= 11 to 18
integer :: veh_type(ntd)   ;!> Viality segment ID from src_td
integer :: id_source_TD(ntd) ;!> Source identification from file src_attr.txt
integer,allocatable  :: id_source_ATT(:)
!> Geometry source type (1 Line  2 Area) from file src_attr.txt
integer,allocatable :: geometry_type2(:);!> Source classification ID viality 1,5,6
integer,allocatable  :: source_type(:)
!>  Number of cars in a specific viality and hour
real,dimension(nhr,ntd) :: veh_number;!> total number of vehicles per hour in each grid
real,allocatable  :: long(:)     ;!> latitude coordinate for the mesh
real,allocatable  :: lat(:)      ;!> Cut-lenght or cut-area of source segment
real :: cutla(nint)   ;!> Relative weight of the viality in the grid
real :: r_weight(nint) ;!> Linear source 1 lenght km or 2 surface area km2 size
real,allocatable  :: source_size(:,:)  ;!> speed in emissions factos
real :: ef_speed(nef,nfe) ;!> VOC emission factor
real :: ef_hc(nef,nfe)  ;!> CO emission factor
real :: ef_co(nef,nfe)  ;!> NO emission factor
real :: ef_no(nef,nfe)  ;!> SO2 emission factor
real :: ef_so2(nef,nfe) ;!> speed in emissions factos cold start engine
real :: ef_speed_cold(nef2,nfe) ;!> VOC emission factor cold start engine
real :: ef_hc_cold(nef2,nfe)  ;!> CO emission factor cold start engine
real :: ef_co_cold(nef2,nfe)  ;!> NO emission factor cold start engine
real :: ef_no_cold(nef2,nfe)  ;!> Correction factor for start engine mode
real :: fcor(nhr)             ;!> Fraction of cold engine cars
real :: f_cold_engine_car(nhr) ;!> Emission factor for specific specie
real :: emiss_factor(nspc)     ;!> Emission factor cold start for specific specie
real :: emis_fact_cold(nspc)   ;!> cars speed in each viality from src_td
real :: veh_speed(nhr,ntd)     ;!> Total emission per cell, specie and day
real :: eday(28*34,nspc,ntypd)
!> Total emission per cell, specie, day and vehicle type
real :: edve(28*34,nspc,ntypd,nveht)
!> Mobile source emision per cell,hour,specie,day type
real :: emision(28*34,nhr,nspc,ntypd)
 !> Emision per cell,hour,specie,day and vehicle type
real :: emi_veh(28*34,nhr,nspc,ntypd,nveht)
!> Number of vehicles per grid, hour and type
real :: no_veh(28*34,nhr,nveht+1,ntypd)
common /intersec/ cutla,r_weight,id_grid_INT,id_source_INT,geometry_type,geo_type_INT
common /srcattr/ natt
common /cellattr/ nic
common /srctd/ veh_number, id_source_TD, ID_time_period,veh_type,veh_speed
common /factemision/ ef_speed,ef_hc,ef_co,ef_no
common /factsec/ ef_speed_cold,ef_hc_cold,ef_co_cold,ef_no_cold
common /miscell/ fcor,f_cold_engine_car
common /computs/ emiss_factor,emis_fact_cold,eday,edve,emision,emi_veh,no_veh

public :: logs,lee_atributos,lee_actividades,lee_factor_emision,calcula_emision, &
         guarda_malla,guarda_malla_nc
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
!> @date 90239
!> @param  utmy Coordinate in axis _y_ in km
!> @param  utmx Coordinate in axis _x_ in km
!> @param  utmz Zone for the UTMx and UTMy coordinates
!> @param  lat Coordinate latitude in decimal degrees
!> @param  long Coordinate longitude in decimal degrees
  subroutine  utm_2_ll(utmx,utmy,utmz,lat,long)
  implicit none
  integer,intent(IN):: utmZ
  real,intent (IN)  :: utmx,utmy
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
!>  @brief Reads file grid cell attributes ID_grid, utmx, utmy
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
    natt=cuenta(iunit)
    allocate(geometry_type2(natt),id_source_ATT(natt))
    allocate(source_type(natt),source_size(natt,2))
    do  j=1,natt
        read(iunit,*)idum,geometry_type2(j),id_source_ATT(j),idum,idum,idum,idum&
        ,source_type(j),source_size(j,1),dum,source_size(j,2)
    end do
    close(iunit)
    call logs("END READING src_attr.txt")
    open(newunit=iunit,file="data/cell_attr.csv",ACTION="READ")
    nic=cuenta(iunit)
    allocate(id_grid(nic),lat(nic),long(nic))
    do i=1,nic
        read (iunit,*)idum,id_grid(i),utmx,utmy
        call utm_2_ll(utmx,utmy,14,lat(i),long(i))
    end do
    close(iunit)
    i=1
    !write(6,*) id_grid(i),lat(i),long(i)
    i=nic
    !write(6,*) id_grid(i),lat(i),long(i)
    call logs("END READING cell_attr.txt")
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
  read (iunit,*)idum,idum,idum,id_source_TD(j),&
        ID_time_period(j),veh_type(j), &
        idum,idum,(veh_number(i,j),i=1,nhr),(veh_speed(i,j),i=1,nhr)
  end do
  close(iunit)
  call logs("END READING intersection.txt")
  open(newunit=iunit,file="data/intersection.csv",ACTION="READ")
  do  j=1,nint
   read(iunit,*)idum,id_grid_INT(j),geo_type_INT(j),geometry_type(j),id_source_INT(j), &
             cutla(j),r_weight(j)
  end do
    call logs("END READING src_td.csv")
  close(iunit)
end subroutine lee_actividades
!  _
! | | ___  ___
! | |/ _ \/ _ \
! | |  __/  __/
! |_|\___|\___|
!  __            _                            _     _
! / _| __ _  ___| |_ ___  _ __   ___ _ __ ___ (_)___(_) ___  _ __
!| |_ / _` |/ __| __/ _ \| '__| / _ \ '_ ` _ \| / __| |/ _ \| '_ \
!|  _| (_| | (__| || (_) | |   |  __/ | | | | | \__ \ | (_) | | | |
!|_|  \__,_|\___|\__\___/|_|____\___|_| |_| |_|_|___/_|\___/|_| |_|
!                          |_____|
!>  @brief Reads emissions factor from EPA, and for cold engine car start
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
      read(iunit,*)ef_speed(j,i),ef_hc(j,i),ef_co(j,i),ef_no(j,i),ef_so2(j,i)
    end do
  end do
  close (iunit)
  call logs('END READING factepa.txt')
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
      read(iunit,*)ef_speed_cold(j,i),ef_hc_cold(j,i),ef_co_cold(j,i),ef_no_cold(j,i)
    end do
  end do
  close (iunit)
  call logs('END READING factsec.txt')
!    -----------  Reading    factvar.dat        unit  19  ----------
!..
open(newunit=iunit,file="data/factvar.dat",ACTION="READ")

  read(iunit,'(a25)')header
  do  i=1,nhr
    read(iunit,*)fcor(i)
  end do
  close (iunit)
  call logs('END READING fraarran.dat')
!    -----------  Reading    fraarran.dat      unit  18   ----------
!..
open(newunit=iunit,file="data/fraarran.dat",ACTION="READ")
  read(iunit,'(a25)') header
  do  i=1,nhr
    read(iunit,*) f_cold_engine_car(i)
  end do
  close (iunit)
  call logs('END READING factvar.dat')
!..
end subroutine lee_factor_emision
!            _            _                         _     _
!   ___ __ _| | ___ _   _| | __ _     ___ _ __ ___ (_)___(_) ___  _ __
!  / __/ _` | |/ __| | | | |/ _` |   / _ \ '_ ` _ \| / __| |/ _ \| '_ \
! | (_| (_| | | (__| |_| | | (_| |  |  __/ | | | | | \__ \ | (_) | | | |
!  \___\__,_|_|\___|\__,_|_|\__,_|___\___|_| |_| |_|_|___/_|\___/|_| |_|
!                               |_____|
!>  @brief computes a single array with emissions using information from grid and EF
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/23/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine calcula_emision
  implicit none
  integer:: ntime,isp,n,m
  integer:: indx
  real :: sl,fcorr,ffr
  real :: temp
  call logs("STARTS EMISSIONS COMPUTATIONS")
  no_veh=0. ! number of vehicles per grid hour
 do n = ntd ,1,-1                   !Main LOOP initialization
    do m = 1,nint
        if(id_source_TD(n).eq.id_source_INT(m).and.geo_type_INT(m).lt.2&
        .and. geometry_type(m) .lt.2 ) then
            indx = id_grid_INT(m)
!..
!     ----------   Localization of the viality lenght     ----------
!..
                call viality(geometry_type(m),id_source_INT(m),geometry_type2,&
                     id_source_ATT,source_type,source_size,fcor,f_cold_engine_car,&
                     ntime,fcorr,ffr,sl)
                sl = sl*r_weight(m)
            do ntime=1,nhr
            ! for EPA
            emiss_factor(1)= emisfac2(veh_type(n),veh_speed(ntime,n),ef_speed ,ef_hc )
            emiss_factor(2)= emisfac2(veh_type(n),veh_speed(ntime,n),ef_speed ,ef_co )
            emiss_factor(3)= emisfac2(veh_type(n),veh_speed(ntime,n),ef_speed ,ef_no )
            emiss_factor(5)= emisfac2(veh_type(n),veh_speed(ntime,n),ef_speed ,ef_so2)
            ! Emissions factor for cold start
            emis_fact_cold(1)= emisfac2(veh_type(n),veh_speed(ntime,n),ef_speed_cold,ef_hc_cold)
            emis_fact_cold(2)= emisfac2(veh_type(n),veh_speed(ntime,n),ef_speed_cold,ef_co_cold)
            emis_fact_cold(3)= emisfac2(veh_type(n),veh_speed(ntime,n),ef_speed_cold,ef_no_cold)
            emis_fact_cold(5)= emiss_factor(5)
            if (veh_type(n).ge. 14 ) then
                emiss_factor(4)   = emiss_factor(1)
                emis_fact_cold(4) = emis_fact_cold(1)
                emiss_factor(1)   = 0.0
                emis_fact_cold(1) = 0.0
            else
                emiss_factor(4)   = 0.0
                emis_fact_cold(4) = 0.0
            end if
!    Number of vehicles all for index 1 and from 2 to 9 for the others
            no_veh(indx,ntime,1,ID_time_period(n))=veh_number(ntime,n)&
                   + no_veh(indx,ntime,veh_type(n)-10,ID_time_period(n))
            no_veh(indx,ntime,veh_type(n)-9,ID_time_period(n))=veh_number(ntime,n) &
                        + no_veh(indx,ntime,veh_type(n)-10,ID_time_period(n))
!     ----------   Computation of the emissions for each specie
                do isp = 1,nspc ! 1 HC(gasoline), 2 CO, 3 NOx, 4 HC(Diesel), 5 SO2
!..
                  temp =(veh_number(ntime,n)*sl*emiss_factor(isp)*(1.0-ffr)+ &
                         veh_number(ntime,n)*sl*emis_fact_cold(isp)*ffr)*fcorr
!..        Xing  is 7.5% of the emission of segment viality
                 if (geo_type_INT(m).eq.2) temp=temp*0.075
!..
                    emision(indx,ntime,isp,ID_time_period(n))= temp/3600.0 &
                            + emision(indx,ntime,isp,ID_time_period(n))
                    emi_veh(indx,ntime,isp,ID_time_period(n),veh_type(n)-10)= temp/3600.0 &
                            + emi_veh(indx,ntime,isp,ID_time_period(n),veh_type(n)-10)
!..
!    ----------  Computation of dayly emissions  EDAY     ----------
!
                    eday(indx,isp,ID_time_period(n))= temp/3600.0 &
                            + eday(indx,isp,ID_time_period(n))
                    edve(indx,isp,ID_time_period(n),veh_type(n)-10)= temp/3600.0 &
                            + edve(indx,isp,ID_time_period(n),veh_type(n)-10)
                end do    ! i nspc
            end do  ! ntime
        end if
    end do  ! m nint
end do ! ntd
end subroutine calcula_emision
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
real:: emis(nic)
  call logs("Wrinting output file for GrADS")
    open (newunit=iunit,file='data/movil.dat', &
    status='unknown',access='direct',form='unformatted' &
    ,recl=nic*4)
    irec = 0
    do iday=1,ntypd
      do l = 1,nhr
        do i = 1,nspc
          do j=1,nic
            if(eday(j,i,iday).gt.0)then
              emis(j)=emision(j,l,i,iday)/eday(j,i,iday)
            else
              emis(j)=0.0
            end if
          end do
          irec = irec +1
          write(iunit,rec=irec)(emis(j),j=1,nic)
        end do   !  i specie
      end do      !  l
    end do
    open (newunit=iunit,file='data/movil_day.dat', &
    status='unknown',access='direct',form='unformatted' &
    ,recl=nic*4)
    irec = 0
    do iday=1,ntypd
        do i = 1,nspc
          irec = irec +1
          write(iunit,rec=irec)(0.010416*eday(j,i,iday),j=1,nic)
        end do   !  i specie
    end do
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
!>  @brief Stores the emissions mesh in a netcdf file
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/20/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
subroutine guarda_malla_nc
use netcdf
implicit none
integer, parameter :: NDIMS=6,nx=28,ny=34, vtype=9
integer :: i,j,k,l,m,iday,ispc,ncid,it,ivt
integer :: iit
integer :: dimids(3),dimids2(2),dimids3(3),dimids4(4),dimidsc(2),dimiscc(2)
integer :: id_unlimit ;!> id_varlat latitude ID in netcdf file
integer :: id_varlat  ;!> id_varlong longitude ID in netcdf file
integer :: id_varlong ;!> id_scc ID of the SCC for each vehicle type
integer :: id_varscc  ;!> id_vardesc Vehicular types description
integer :: id_vardesc ;!> pollutant TP and Total emission ID in netcdf file
integer :: id_var(nspc*2) ;!> vehicules type id in netcdf file
integer :: id_nveh     ;!> dimensions values
integer,dimension(NDIMS):: dim;!> dimensions ID
integer,dimension(NDIMS):: id_dim ;!>xlong longitude coordinates
real,dimension(nx,ny)::xlong ;!>xlat latitude coordinates
real,dimension(nx,ny)::xlat  ;!>tprof temporal profiles
real,dimension(nx,ny,vtype,1)::tprof ;!>emis_day dayly emissions
real,dimension(nx,ny,vtype):: emis_day;!> num cars per grid, hour and type
real,dimension(nx,ny,vtype,1):: numero_veh
character (len=19),dimension(NDIMS) ::sdim
character(len=19) :: current_date
character(len=16) :: vehiclename
character(len= 8) :: date
character(len=10) :: time
character(len=24) :: hoy, fecha_creado
character(len=26) :: FILE_NAME
character(len=19),dimension(1,1)::Times
character(len=19),dimension(vtype,1)::cveh_type
character(len=10),dimension(vtype,1)::scc
character(len=11),dimension(2*nspc):: ename ;!> Emissions long name
character(len=26),dimension(2*nspc):: cname
character(len=42),dimension(ntypd):: title
character(len=250)::summary
data title /"Movile temporal profiles for weekdays V4.0",&
            "Movile temporal profiles for Saturday V4.0",&
            "Movile temporal profiles for Sunday V4.0  "/
data sdim /"Time               ","DateStrLen         ","west_east          ",&
&          "south_north        ","emissions_zdim_stag","vehicle_type       "/
ename=(/'TP_VOC       ','TP_CO        ','TP_NO        ','TP_VOC_diesel',&
        'TP_SO2       ','E_VOC        ','E_CO         ','E_NO         ',&
        'E_VOC_diesel ','E_SO2        '/)
cname=(/'VOC gasoline vehicle      ','Carbon Monoxide           ', &
        'Nitrogen Oxide            ','VOC Diesel                ', &
        'Sulfur Dioxide            ','VOC gasoline vehicle      ', &
        'Carbon Monoxide           ', &
        'Nitrogen monoxide         ','VOC Diesel vehicles       ', &
        'Sulfur dioxide emissions  '/)
 summary="Using information from 11,848 street segments with activity data "//&
"(type, number, speed, street length, engine cold start fraction) for computing"//&
" a daily and hourly emission per grid. A ratio is used for calculating a temporal"//&
" profile per grid"
cveh_type(1,1)="All_Types_vehicles ";cveh_type(2,1)="Automoviles     ,11"
cveh_type(3,1)="Ligeros         ,12";cveh_type(4,1)="Microbuses      ,13"
cveh_type(5,1)="Ruta_100        ,14";cveh_type(6,1)="Otros_buses     ,15"
cveh_type(7,1)="Medianos        ,16";cveh_type(8,1)="Pesados         ,17"
cveh_type(9,1)="Camion_Municipal,18"
!>   |Type |   SCC    |Description      |Type |  SCC     |Description |
!>   | --- | ---| ---| --- | ---| ---|
!>   | 11 |2201001330 |Automoviles     |  15  |2230075330|  Otros buses |
!>   | 12 |2201040330 |Ligeros         |  16  |2230060330|  Medianos    |
!>   | 13 |2230070270 |Microbuses      |  17  |2230074330|  Pesados     |
!>   | 14 |2230001000 |Ruta 100        |  18  |2230001000| Camiones Municipales |
scc(1,1)="2201001330";scc(2,1)="2201001330";scc(3,1)="2201040330";scc(4,1)="2230070270"
scc(5,1)="2230001000";scc(6,1)="2230075330";scc(7,1)="2230060330";scc(8,1)="2230074330"
scc(9,1)="2230001000"

  call logs("Wrinting output file for netcdf")
  dim=(/1,19,nx,ny,10,vtype/)
  call date_and_time(date,time)
  fecha_creado=date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T'//time(1:2)//':'//time(3:4)//':00Z'
  hoy=date(7:8)//'-'//mes(date(5:6))//'-'//date(1:4)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
  call logs("   HOY: "//hoy)
  xlong=reshape(long,(/nx,ny/))
  xlat=reshape(lat,(/nx,ny/))
  do iday=1,ntypd
    current_date="1990-01-05_00:00:00"
    write(current_date(09:10),'(I2.2)') iday+4
    FILE_NAME="emission_"//current_date(1:13)//".nc"
    write(6,182) iday, FILE_NAME
    call check( nf90_create(path =FILE_NAME,cmode = NF90_CLASSIC_MODEL,ncid = ncid) )
!     Define dimensiones
    call check( nf90_def_dim(ncid,sdim(1), NF90_UNLIMITED, id_dim(1)) )
    do i=2,NDIMS
        call check( nf90_def_dim(ncid, sdim(i), dim(i), id_dim(i)) )
    end do
    dimids  = (/id_dim(3),id_dim(4),id_dim(6)/) ! nx,ny,vtype
    dimidsc = (/id_dim(2),id_dim(6)/)           ! Str(19),vtype
    dimiscc = (/id_dim(5),id_dim(6)/)           ! Str(10),time
    dimids2 = (/id_dim(2),id_dim(1)/)           ! Str(19), time
    dimids3 = (/id_dim(3),id_dim(4),id_dim(1)/) ! nx,ny, time
    dimids4 = (/id_dim(3),id_dim(4),id_dim(6),id_dim(1)/)! nx,ny,vtype,time
    if (iday.eq.1) call logs("Atributos Globales NF90_GLOBAL")
    call check( nf90_put_att(ncid, NF90_GLOBAL, "TITLE",title(iday)))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "START_DATE",current_date))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION",nx))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION",ny))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "BOTTOM-TOP_GRID_DIMENSION",1))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "DAY","FRD"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "DX",2000))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "DY",2000))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "CEN_LAT",xlat(nx/2,ny/2)))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "CEN_LON",xlong(nx/2,ny/2)))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "TRUELAT1",17.5))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "TRUELAT2",29.5))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "MOAD_CEN_LAT",xlat(nx/2,ny/2)))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "STAND_LON",xlong(nx/2,ny/2)))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "POLE_LAT",90.))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "POLE_LON",0.))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "GRIDTYPE","C"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "GMT",-6.))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "JULYR",1990))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "JULDAY",5))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "MAP_PROJ",1))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "standard_parallel","(/17.5,29.5/)"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "grid_mapping_name","lambert_conformal_conic"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "MMINLU","USGS"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "MECHANISM","NONE"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "creator_institution", &
    "Centro de Ciencias de la Atmosfera, UNAM"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "creator_type","institution"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "contributor_name",&
    "Agustin Garcia, agustin@atmosfera.unam.mx"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "contributor_role","Researcher"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "cdm_data_type","Grid"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "acknowledgment","CCA, UNAM"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "publisher_institution","CCA, UNAM"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "publisher_url","www.atmosfera.unam.mx"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "publisher_type","institution"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "product_version","1.0"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "CREATION_DATE",hoy))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "date_issued",fecha_creado))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "date_created",fecha_creado))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "date_modified",fecha_creado))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "date_metadata_modified",fecha_creado))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_start","1990-01-05T00:00:00Z"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_end","1990-01-08T00:00:00Z"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_duration","P3D"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_resolution","PT1H"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_units","degrees_east"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_units","degrees_north"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_max",maxval(xlat)))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_min",minval(xlat)))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_max",maxval(xlong)))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_min",minval(xlong)))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "geospatial_bounds_crs","EPSG:4979"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "id","temporal_frac_using_EPA_EF_1993"))
    call check( nf90_put_att(ncid, NF90_GLOBAL, "summary",summary))
!  Define las variables
    call check( nf90_def_var(ncid, "Times", NF90_CHAR, dimids2,id_unlimit ) )
    call check( nf90_def_var(ncid, "Type",  NF90_CHAR, dimidsc,id_vardesc ) )
    call check( nf90_def_var(ncid, "SCC" ,  NF90_CHAR, dimiscc,id_varscc  ) )
    call check( nf90_def_var(ncid, "XLONG", NF90_REAL, dimids3,id_varlong) )
    call check( nf90_def_var(ncid, "XLAT" , NF90_REAL, dimids3,id_varlat ) )
! Assign  attributes
    call check( nf90_put_att(ncid, id_varlong, "FieldType", 104 ) )
    call check( nf90_put_att(ncid, id_varlong, "MemoryOrder", "XYZ") )
    call check( nf90_put_att(ncid, id_varlong, "description", "LONGITUDE, WEST IS NEGATIVE") )
    call check( nf90_put_att(ncid, id_varlong, "standard_name", "grid_longitude") )
    call check( nf90_put_att(ncid, id_varlong, "units", "degrees_east"))
    call check( nf90_put_att(ncid, id_varlong, "axis", "X") )
    call check( nf90_put_att(ncid, id_varlat, "FieldType", 104 ) )
    call check( nf90_put_att(ncid, id_varlat, "MemoryOrder", "XYZ") )
    call check( nf90_put_att(ncid, id_varlat, "description", "LATITUDE, SOUTH IS NEGATIVE") )
    call check( nf90_put_att(ncid, id_varlat, "standard_name", "grid_latitude") )
    call check( nf90_put_att(ncid, id_varlat, "units", "degrees_north"))
    call check( nf90_put_att(ncid, id_varlat, "axis", "Y") )
    call check( nf90_put_att(ncid, id_vardesc, "description", "Vehicular types from 1 to 9") )
    call check( nf90_put_att(ncid, id_varscc,"description","Source Clasification Code for each vehicle type") )

!  Attributos para cada perfil temporal
    do i=1,nspc
     call crea_attr(ncid,0,dimids4,ename(i),cname(i),"1",id_var(i)) !adimensional=1
    end do
!  Attributos para cada emision diaria
    do i=1+nspc,nspc+nspc
     call crea_attr(ncid,1,dimids,ename(i),cname(i),"g km^-2 s^-1",id_var(i))
    end do
!  Attributos para numero de vehiculos
      vehiclename="traffic"
      call crea_attr(ncid,2,dimids4,vehiclename,vehiclename,"number",id_nveh)
!   Terminan definiciones
    call check( nf90_enddef(ncid) )

  call check( nf90_put_var(ncid, id_vardesc,cveh_type,start=(/1,1/)) )
  call check( nf90_put_var(ncid, id_varscc,scc,start=(/1,1/)) )
  tiempo: do it=1,nhr
    iit=it
    write(current_date(12:13),'(I2.2)') it-1
    Times(1,1)=current_date(1:19)
    !write(6,'(A,x,I2.2)')'TIEMPO: ', iit
    call check( nf90_put_var(ncid, id_unlimit,Times,start=(/1,it/)) )
    call check( nf90_put_var(ncid, id_varlong,xlong,start=(/1,1,it/)) )
    call check( nf90_put_var(ncid, id_varlat,xlat,start=(/1,1,it/)) )
    do ivt=1,vtype
        do i=1,nx
            do j=1,ny
              k=i+28*(j-1)
            numero_veh(i,j,ivt,1)=no_veh(k,iit,ivt,iday)
            end do
        end do
    end do
    call check(nf90_put_var(ncid,id_nveh,numero_veh,start=(/1,1,1,it/)))
    do ispc=1,nspc
      do i=1,nx
        do j=1,ny
          k=i+28*(j-1)
            if(eday(k,ispc,iday).gt.0) then
                tprof(i,j,1,1)=emision(k,iit,ispc,iday)/eday(k,ispc,iday)
            else
                tprof(i,j,1,1)=0.0
            end if
            do m=1,nveht
                if(edve(k,ispc,iday,m).gt.0) then
                 tprof(i,j,m+1,1)=emi_veh(k,iit,ispc,iday,m)/edve(k,ispc,iday,m)
                else
                 tprof(i,j,m+1,1)=0
                end if
            end do
          if (it.eq. 1) then
            emis_day(i,j,1)=eday(k,ispc,iday)/4./real(nhr)! grid cell 4 km^2
            do m=1,nveht
            emis_day(i,j,m+1)=edve(k,ispc,iday,m)/4./real(nhr)! grid cell 4 km^2
            end do
          end if
        end do
      end do
      call check( nf90_put_var(ncid, id_var(ispc),tprof,start=(/1,1,1,it/)))
      if (it.eq. 1)call check(nf90_put_var(ncid,id_var(ispc+nspc),emis_day,start=(/1,1,1/)))
    end do
   end do TIEMPO
   call check( nf90_close(ncid) )
  end do !   day
182 format(5X,'      Guarda variables dia: ',I2.2,x,A26)
end subroutine guarda_malla_nc
!        _       _ _ _
! __   _(_) __ _| (_) |_ _   _
! \ \ / / |/ _` | | | __| | | |
!  \ V /| | (_| | | | |_| |_| |
!   \_/ |_|\__,_|_|_|\__|\__, |
!                        |___/
!>  @brief Localization of the viality and depending of his type it
!>  assings  the value of cold start fraction (ffr), fcorr, and longitud
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/21/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
!>  @param ig Source Type (1 line 2 Area)
!>  @param id viality source identification
!>  @param ig2 Source Type (1 line 2 Area) from intersections
!>  @param isrc viality source identification from intersections
!>  @param kstype Source classification
!>  @param slen  source size in km or km^2, for 1-line, 2-Area
!>  @param fcor correction factor for start mode
!>  @param f_cold_engine_car car fraction staring cool
!>  @param nh  time hour for set the EF
!>  @param fcorr  correction factor of ffr (out)
!>  @param ffr  car fraction staring factor (out)
!>  @param sl Viality length in km (out)
  Subroutine viality(ig,id,ig2,isrc,kstype,slen &
                    ,fcor,f_cold_engine_car,nh,fcorr,ffr,sl)
  integer,intent(in):: ig,id,isrc(:),kstype(:)
  integer,intent(in):: ig2(:),nh
  real,intent(in) :: fcor(:),f_cold_engine_car(:)
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
        ffr   = f_cold_engine_car(nh)
        if(kstype(i).gt.20) then
          ffr= f_cold_engine_car(nh)
          fcorr =1.0
          exit
        end if
        exit
      else
        flag  = 1
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
 200      format('xxx ERROR: Invalid Cell',I6)
 201      format('XXX ERROR: Invalid Longitud',f4.0,'At geometry_type ',I3,&
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
!>  @param velocity vehicle speed in viality
!>  @param vem Velocity array from emission factors file
!>  @param comp Emission factor per speed and specie from emission factors file
  real function emisfac2(ncartype,velocity,vem,comp)
  integer :: ncartype
  integer :: icar, i
  real:: velocity,vem(:,:),comp(:,:)
!> | ncartype |Type of vehicle | icar  |   Type of vehicle |
!> |   ---    |---             | ---   |--- |
!> |    11    |   Automoviles  |  1 |  Vehiculos ligeros a Gasolina  |
!> |    12    |   Ligeros      |  2 |  Camionetas ligeras a Gasolina |
!> |    13    |   Microbuses   |  2 |  Camionetas ligeras a Gasolina |
!> |    13    |   Microbuses   |  5 |  Camiones Pesados a Gasolina |
!> |    14    |   Ruta 100     |  3 |  Camiones ligeros a Diesel   |
!> |    15    |  Otros camiones|  3 |  Camiones ligeros a Diesel   |
!> |    16    |   Medianos     |  3 |  Camiones ligeros a Diesel   |
!> |    17    |   Pesados      |  4 |  Vehiculos pesados a Diesel  |
!> |    18    |  Camiones municipales| 3| Camiones ligeros a Diesel|
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
    write(6,190) ncartype
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
!.. linear interpolation using vehicular speed
    emisfac2= comp(icar,i) +(velocity-vem(icar,i)) * &
        (comp(icar,i+1)-comp(icar,i))/(vem(icar,i+1)-vem(icar,i))
!..
190  format('XXXX ERROR: Invalid car type:',x,I3)
200  format('XXXX ERROR: Invalid speed',E10.1)

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
!> @param status NetCDF functions return a non-zero status codes on error.
!> @copyright 1993-2020 University Corporation for Atmospheric Research/Unidata
! Copyright 1993-2020 University Corporation for Atmospheric Research/Unidata
!
!Portions of this software were developed by the Unidata Program at the
!University Corporation for Atmospheric Research.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors
! may be used to endorse or promote products derived from this software without
! specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
! OF THE POSSIBILITY OF SUCH DAMAGE.
subroutine check(status)
use netcdf
    integer, intent (in) :: status
    if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 2
    end if
end subroutine check
!>  @brief count the number of rowns in a file
!>   @author  Jose Agustin Garcia Reynoso
!>   @date  07/13/2020
!>   @version  2.2
!>   @copyright Universidad Nacional Autonoma de Mexico 2020
!>   @param  iunit file unit where the count has to be made
integer function  cuenta(iunit)
    implicit none
    integer,intent(IN) :: iunit
    integer :: io
    cuenta = 0
    DO
        READ(iunit,*,iostat=io)
        IF (io/=0) EXIT
        cuenta = cuenta + 1
    END DO
    rewind(iunit)
    return
end

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
    character*2,intent(IN):: num
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
!                               _   _
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
!>   @param ifl number of items in dimids array
!>   @param dimids ID dimensons array
!>   @param svar short variable name
!>   @param cname description variable name
!>   @param cunits units of the variable
!>   @param id_var variable ID
subroutine crea_attr(ncid,ifl,dimids,svar,cname,cunits,id_var)
use netcdf
    implicit none
    integer , INTENT(IN) ::ncid,ifl
    integer, INTENT(out) :: id_var
    integer, INTENT(IN),dimension(:):: dimids
    character(len=*), INTENT(IN)::svar,cname,cunits
    character(len=50) :: cvar
    if (ifl.eq.0)cvar="temporal_profile "//trim(cname)
    if (ifl.eq.1) cvar="Flux "//trim(cname)
    if (ifl.eq.2) cvar="Number vehicles type "//trim(cname)
    call check( nf90_def_var(ncid, svar, NF90_REAL, dimids,id_var ) )
    ! Assign  attributes
    call check( nf90_put_att(ncid, id_var, "FieldType", 104 ) )
    call check( nf90_put_att(ncid, id_var, "MemoryOrder", "XYZ") )
    call check( nf90_put_att(ncid, id_var, "standard_name", cvar) )
    call check( nf90_put_att(ncid, id_var, "description", cvar) )
    call check( nf90_put_att(ncid, id_var, "units", cunits))
    call check( nf90_put_att(ncid, id_var, "stagger", "Z") )
    call check( nf90_put_att(ncid, id_var, "coordinates", "XLONG XLAT") )
    call check( nf90_put_att(ncid, id_var, "coverage_content_type","modelResult"))
    ! print *,"Entro a Attributos de variable",dimids,id,jd
    return
end subroutine crea_attr
!>  @brief display log during different program stages
!>   @author  Jose Agustin Garcia Reynoso
!>   @date  08/08/2020
!>   @version  2.2
!>   @copyright Universidad Nacional Autonoma de Mexico 2020
!>   @param texto text to be displayed
subroutine logs(texto)
    implicit none
    character(len=*),intent(in):: texto
    write(6,'(3x,6("*"),x,A35,x,"******")') texto
end subroutine
end module grid_temp_mobile

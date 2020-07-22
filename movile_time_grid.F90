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
integer,parameter :: nef = 5, nef2 =5,ntypd= 3, nhr=24, nspc = 5
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
integer ::idsrc2(natt),igeo2(natt),kstype(natt)
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
real :: fcor(nhr),start(nhr); !Emission factor at speed vel
real :: efact(nspc) ; !Emission factor cold start at speed vel
real :: efac2(nspc)  ;!> Street velocity from src_td
real :: vel(nhr,ntd) ;!> Total emission per day
real :: eday(nic,nspc,ntypd,8,1:2) !> Emissions per cell
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
!>  @brief Localization of the viality and depending of his type it
!>  assings  the value of start or fcorr, and longitud
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/21/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
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
!>  @brief Emission factor computation for velociity and EF array
!>  @author Jose Agustin Garcia Reynoso
!>  @date 07/22/2020
!>  @version  1.0
!>  @copyright Universidad Nacional Autonoma de Mexico
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

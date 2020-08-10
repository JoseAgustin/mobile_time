! Test subroutine utm_2_ll
program utm2ll
use grid_temp_mobile, only:utm_2_ll
implicit none
integer:: mz=14
real :: mx,my
real :: latitud, longitud

mx=465.00
my=2115.0

  call utm_2_ll(mx,my,mz,latitud,longitud)
  write(6,*)mx,my," = ",latitud, longitud
mx=519.0
  call utm_2_ll(mx,my,mz,latitud,longitud)
  write(6,*)mx,my," = ",latitud, longitud
my=2181.0
  call utm_2_ll(mx,my,mz,latitud,longitud)
  write(6,*)mx,my," = ",latitud, longitud
mx=465.00
  call utm_2_ll(mx,my,mz,latitud,longitud)
  write(6,*)mx,my," = ",latitud, longitud

end program utm2ll

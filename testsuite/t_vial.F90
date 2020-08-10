!  test viality
program test_viality
use grid_temp_mobile
implicit none
    integer:: ntime,isp,n,m,indx
    real :: sl,fcorr,ffr
   call logs("Prueba subroutina viality   ")
!
    call system("ln -s ../data .")
    call lee_atributos
    call lee_actividades
    call lee_factor_emision
    do n=2115,2145,3
      do m= 1,nint
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
             write(6,100)indx,sl,r_weight(m),sl/r_weight(m)
          end if
        end do
      end do
call system("rm data")
100 format(I7,x,f6.3,x,f5.3,x,f6.3)
end program test_viality

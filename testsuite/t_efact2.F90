!  test emisfac2
program test_emisfac2
use grid_temp_mobile
implicit none
    integer:: ntime,isp,n,m,indx
!
    call system("ln -s ../data .")
    call lee_atributos
    call lee_actividades
    call lee_factor_emision
    do n=15,22
      do m= 1,nint
        if(id_source_TD(n).eq.id_source_INT(m).and.geo_type_INT(m).lt.2&
        .and. geometry_type(m) .lt.2 ) then
            indx = id_grid_INT(m)
            do ntime=1,2
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
          end do
          write(6,'(I2,5(F7.3,x))') veh_type(n),emiss_factor
          write(6,'(I2,5(F7.3,x))') ntime-1,emis_fact_cold
          exit
          end if
        end do
      end do
call system("rm data")
end program test_emisfac2

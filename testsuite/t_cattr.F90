!  Test for crea_attr subroutine
program test_crea_attr
use netcdf
use grid_temp_mobile, only :check,crea_attr,logs
integer, parameter :: NDIMS=6,nx=28,ny=34, vtype=9
integer :: ncid,id_unlimit
integer :: dimids2(2),dimids4(4),id_var
integer,dimension(NDIMS):: dim,id_dim
character (len=19),dimension(NDIMS) ::sdim
character(len=11):: ename ='TP_VOC     '
character(len=26):: cname ='VOC gasoline vehicle      '
data sdim /"Time               ","DateStrLen         ","west_east          ",&
&          "south_north        ","emissions_zdim_stag","vehicle_type       "/

    dim=(/1,19,nx,ny,10,vtype/)

    call check( nf90_create(path ="test.nc",cmode = NF90_CLASSIC_MODEL,ncid = ncid) )
    call check( nf90_def_dim(ncid,sdim(1), NF90_UNLIMITED, id_dim(1)) )
    do i=2,NDIMS
        call check( nf90_def_dim(ncid, sdim(i), dim(i), id_dim(i)) )
    end do
    dimids2 = (/id_dim(2),id_dim(1)/)
    dimids4 = (/id_dim(3),id_dim(4),id_dim(6),id_dim(1)/)
    call check( nf90_put_att(ncid, NF90_GLOBAL, "TITLE","Test 11 attributes creation"))
    call check( nf90_def_var(ncid, "Times", NF90_CHAR, dimids2,id_unlimit ) )
    call logs("Entra a crea_attr")
    call crea_attr(ncid,1,dimids4,ename,cname,"g km^-2 s^-1",id_var)
    call check( nf90_enddef(ncid) )
    call check( nf90_close(ncid) )
    call system("rm test.nc")
end program test_crea_attr

  integer function  cuenta(iunit)
  implicit none
  integer :: iunit
  integer :: nl,io
  nl = 0 
  DO
    READ(iunit,*,iostat=io)
    IF (io/=0) EXIT
    nl = nl + 1
  END DO
  rewind(iunit)
  cuenta=nl
 return
 end

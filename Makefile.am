#
# This file is processed by automake (produces file Makefile.in)
#
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = mobile_time_grid.exe

# Additional files to be cleaned with 'make clean'
CLEANFILES = *.o *.mod 
# Compiler flags
  FC =$(NC_FC)
  AM_FCFLAGS  = $(NC_INC)
  AM_FCFLAGS += $(NC_LIB)

mobile_time_grid_exe_SOURCES = mobile_time_grid.F90

EXTRA_DIST = README.md assets _config.yml Doc/Makefile Doc/doxygen_mobile_time.conf

.PHONY: $(EXTRA_DIST)
$(EXTRA_DIST):

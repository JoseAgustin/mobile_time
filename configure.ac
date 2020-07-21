# Process this file with autoconf to produce a configure script.
# To produce files: configure, Makefiles and Scripts
#    $ automake
#    $ autoconf
#
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([movil_time_grid], [1.0], [agustin@atmosfera.unam.mx])

# Define our M4 macro directory
AC_CONFIG_AUX_DIR([autoconf])
AM_INIT_AUTOMAKE([1.11 foreign])

# Set default prefix (where directory bin is created)
# This is the default top directory of the installation
AC_PREFIX_DEFAULT(`pwd`)

# Checks for programs.
AC_LANG(Fortran)
AC_PROG_FC
# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_LANG_PUSH(Fortran)
AC_FC_FREEFORM

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# Write configuration on the screen
AC_MSG_NOTICE([---------------------------------------------------------])
AC_MSG_NOTICE([Configuration complete - $PACKAGE_NAME-$PACKAGE_VERSION $MSG_PARALLEL])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([Fortran compiler:           FC=$FC])
AC_MSG_NOTICE([Fortran flags:              FCFLAGS=$FCFLAGS])
AC_MSG_NOTICE([Install prefix:             --prefix=$prefix])
AC_MSG_NOTICE([Executables install prefix: --exec_prefix=$exec_prefix])
AC_MSG_NOTICE([Binary directory:           --bindir=$bindir])
AC_MSG_NOTICE([---------------------------------------------------------])

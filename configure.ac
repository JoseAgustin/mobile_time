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
AC_PROG_MAKE_SET
AC_PROG_CPP

# Set language for configuration checks
AC_LANG(Fortran)
AC_PROG_FC
# Make precious variables for ./configure --help
AC_ARG_VAR(NETCDF_ROOT,[Root of netCDF installation or] (compile-time))
AC_ARG_VAR(NETCDF_INC,[  Location of netCDF headers] (compile-time))
AC_ARG_VAR(NC_LIB,[  Location of netCDF library] (compile-time))

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_LANG_PUSH(Fortran)
AC_FC_FREEFORM
if test -z "${NETCDF_ROOT}"; then
  AC_CHECK_PROG([NETCDF_CONFIG],[nf-config],[nf-config])
        # If not found, check for nc-config
  AS_IF([test -z "${NETCDF_CONFIG}"],
     [AC_CHECK_PROG([NETCDF_CONFIG],[nc-config],[nc-config])])
        NC_LIB="`${NETCDF_CONFIG} --flibs` ";
        NC_INC="-I`${NETCDF_CONFIG} --includedir`";
        NC_FC="`${NETCDF_CONFIG} --fc` "
        AC_PATH_PROG([NETCDF_ROOT], [ncdump])
        NETCDF_ROOT=`dirname $NETCDF_ROOT`   # Strip /ncdump
        NETCDF_ROOT=`dirname $NETCDF_ROOT`;  # Strip /bin
    else
        NETCDF_CONFIG="${NETCDF_ROOT}/bin/nf-config"
        NC_FC="`${NETCDF_CONFIG} --fc` "
        NC_LIB="`${NETCDF_CONFIG} --flibs` "
        NC_INC="-I`${NETCDF_CONFIG} --includedir` ";
    fi

AC_MSG_NOTICE([setting netcdf root directory NETCDF_ROOT=$NETCDF_ROOT])
AC_MSG_CHECKING([for the existence of the netcdf root directory])
AS_IF([test -d $NETCDF_ROOT],[AC_MSG_RESULT([$NETCDF_ROOT])],
    [AC_MSG_ERROR([directory $NETCDF_ROOT not found])])

# Set the netcdf compiler flags
AC_MSG_CHECKING([for netCDF include directory])
AS_IF([test -z "$NETCDF_INC"],[NETCDF_INC=$NETCDF_ROOT/include])
AC_MSG_RESULT($NETCDF_INC)
AC_MSG_CHECKING([for the existence of the netcdf include directory])
AS_IF([test -d $NETCDF_INC],[AC_MSG_RESULT([ok])],
    [AC_MSG_ERROR([directory $NETCDF_INC not found])])
AC_MSG_CHECKING([for file netcdf.mod in the include directory])
AS_IF([test -f $NETCDF_INC/netcdf.mod],
   [AC_MSG_RESULT(ok)],
   [AC_MSG_ERROR([File: netcdf.mod was not found in $NETCDF_INC
   Probably the Fortran version of NetCDF is not correctly installed])])
AS_IF([ test -z "{$NC_INC}"],[NC_INC=-I$NETCDF_INC])
# Set the netcdf linker flags
AC_MSG_CHECKING([for netCDF library directory])
AS_IF([test "X$NC_LIB" = "X"],[NC_LIB=$NETCDF_ROOT/lib])
AC_MSG_RESULT([$NC_LIB])
AS_IF([test ! -d $NETCDF_ROOT/lib ],[AC_MSG_ERROR([No such directory: $NC_LIB])])
AS_IF([test -n $NC_LIB ],[AC_MSG_CHECKING([for -lnetcdff]);
    AS_IF([test -f $NETCDF_ROOT/lib/libnetcdff.a],
      [NC_LIB="$NC_LIB -lnetcdff -lnetcdf";AC_MSG_RESULT([yes])],
      [NC_LIB="$NC_LIB -lnetcdf";AC_MSG_RESULT([no])
   ])
])
#
# NC_INC and NC_LIB contain the flags for compilation and linking and are
# exported to Makefiles

AC_SUBST(NC_INC)
AC_SUBST(NC_LIB)
AC_SUBST(NC_FC)
# Checks for library functions.

AC_LANG_POP(Fortran)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# Write configuration on the screen
AC_MSG_NOTICE([---------------------------------------------------------])
AC_MSG_NOTICE([Configuration complete - $PACKAGE_NAME-$PACKAGE_VERSION $MSG_PARALLEL])
AC_MSG_NOTICE([])
AS_IF([test "x$NC_VERSION" != "x"],AC_MSG_NOTICE([Using $NC_VERSION]))
AC_MSG_NOTICE([Fortran compiler:           FC=$FC])
AC_MSG_NOTICE([Fortran flags:              FCFLAGS=$FCFLAGS])
AC_MSG_NOTICE([Root directory of netcdf:   NETCDF=$NETCDF_ROOT])
AC_MSG_NOTICE([Compiler flags for netcdf:  NC_INC=$NC_INC])
AC_MSG_NOTICE([Linker flags for netcdf:    NC_LIB=$NC_LIB])
AC_MSG_NOTICE([Fortran compiler netcdf:    NC_FC=$NC_FC])
AC_MSG_NOTICE([Install prefix:             --prefix=$prefix])
AC_MSG_NOTICE([Executables install prefix: --exec_prefix=$exec_prefix])
AC_MSG_NOTICE([Binary directory:           --bindir=$bindir])
AC_MSG_NOTICE([---------------------------------------------------------])

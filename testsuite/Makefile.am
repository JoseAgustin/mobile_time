#
# This file is processed by automake (produces file Makefile.in)
#
# Compiler flags -I$(top_srcdir)
FC =$(NC_FC)
AM_FCFLAGS  = $(NC_INC) -I$(top_srcdir)
AM_FCFLAGS += $(NC_LIB)

FCFLAGS = $(AM_FCFLAGS)

TESTS = test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11
TESTS += test12 test13 test14
check_PROGRAMS= test1 test2 test3 test4 test5 test6 test7 test8 test9 test10
check_PROGRAMS += test11 test12 test13 test14
test1_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_cuenta.F90
test2_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_mes.F90
test3_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_utm2ll.F90
test4_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_logs.F90
test5_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_activ.F90
test6_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_attr.F90
test7_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_femis.F90
test8_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_vial.F90
test9_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_efact2.F90
test10_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_check.F90
test11_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_cattr.F90
test12_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_calcula.F90
test13_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_guarda.F90
test14_SOURCES= $(top_builddir)/mod_mobile_time_grid.F90 t_guarda_nc.F90

# Additional files to be cleaned with 'make clean'
CLEANFILES = *.o *.mod package.m4 testsuite

# The `:;' works around a Bash 3.2 bug when the output is not writeable.
$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	:;{ \
	echo '# Signature of the current package.' && \
 	echo 'm4_define([AT_PACKAGE_NAME],' && \
 	echo '  [@PACKAGE_NAME@])' && \
 	echo 'm4_define([AT_PACKAGE_TARNAME],' && \
 	echo '  [@PACKAGE_TARNAME@])' && \
 	echo 'm4_define([AT_PACKAGE_VERSION],' && \
 	echo '  [@PACKAGE_VERSION@])' && \
	echo 'm4_define([AT_PACKAGE_STRING],' && \
	echo '  [@PACKAGE_STRING@])' && \
	echo 'm4_define([AT_PACKAGE_BUGREPORT],' && \
	echo '  [@PACKAGE_BUGREPORT@])'; \
	echo 'm4_define([AT_PACKAGE_URL],' && \
	echo '  [@PACKAGE_URL@])'; \
	}>'$(srcdir)/package.m4'

     
EXTRA_DIST = testsuite.at $(srcdir)/package.m4 $(TESTSUITE)
TESTSUITE = $(srcdir)/testsuite
     
check-local: atconfig  $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' $(TESTSUITEFLAGS)

maintainer-clean-local:
	rm -rf aclocal.m4 configure Makefile.in autotools m4
     
installcheck-local: atconfig  $(TESTSUITE)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' \
	$(TESTSUITEFLAGS)
     
clean-local:
	test ! -f '$(TESTSUITE)' || \
	$(SHELL) '$(TESTSUITE)' --clean
     
AUTOM4TE = $(SHELL) $(srcdir)/build-aux/missing --run autom4te
AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/package.m4
	$(AUTOTEST) -I '$(srcdir)' -o $@.tmp $@.at
	mv $@.tmp $@

.PHONY: $(EXTRA_DIST)
$(EXTRA_DIST):

AC_PREREQ(2.62)
AC_INIT([vorticus],[0.1])
AC_CANONICAL_SYSTEM
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror])
AM_SILENT_RULES([yes])
AM_PROG_AR

AC_LANG([C++])
AC_PROG_CC
AC_PROG_CXX
AC_OPENMP

AX_BLAS([],[AC_MSG_ERROR([vorticus requires BLAS.])])
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AC_SUBST(AM_CPPFLAGS, "$BOOST_CPPFLAGS  $OPENMP_CXXFLAGS -I$srcdir/.. $AM_CPPFLAGS" )
AC_SUBST(AM_CXXFLAGS, "$OPENMP_CXXFLAGS $AM_CXXFLAGS" )
AC_SUBST(AM_LDFLAGS,  "$BOOST_LDFLAGS   $OPENMP_CXXFLAGS $AM_LDFLAGS" )

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 src/makefile
                 src/vrm/makefile
                 src/fmm/makefile
                 src/math/makefile
		         src/misc/makefile
                 src/geometry/makefile])

LT_PREREQ([2.2])
LT_INIT([disable-shared])
AC_OUTPUT


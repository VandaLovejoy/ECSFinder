AC_INIT(SISSIz, 3.0, sissiz.vbc5@univie.ac.at, SISSIz)
AM_INIT_AUTOMAKE([ -Wall foreign])
AM_PROG_AR
AC_PROG_CC 


if test "$GCC" = yes; then
  AC_DEFINE(UNUSED, __attribute__ ((unused)), avoid warning about unused variables)
else
  AC_DEFINE(UNUSED,, only for gcc)
fi

AC_CANONICAL_HOST
is_linux=no
is_mac=no

case "${host_os}" in
    linux*)
        is_linux=yes
        ;;
    darwin*)
        is_mac=yes
        ;;
    *)
        AC_MSG_NOTICE(["Other OS: $host_os"])
        ;;
esac

AM_CONDITIONAL([DISABLE_OPENMP], [test "$is_mac" = "yes"])
AC_MSG_NOTICE($host_os)
AM_COND_IF([DISABLE_OPENMP], [AC_MSG_NOTICE(OPENMP disabled for Mac OS users)],[AC_MSG_NOTICE(OPENMP enabled)])
          
AC_SUBST(VERSION)
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)


AC_HEADER_STDC
AC_CHECK_HEADERS(malloc.h strings.h unistd.h)
AC_C_CONST
AC_TYPE_SIZE_T
AC_INLINE

AC_CHECK_FUNCS(strdup strstr strchr erand48)

AC_PROG_CXX

AC_PROG_RANLIB

AC_PROG_INSTALL

AC_CONFIG_HEADERS(config.h) 

AC_OUTPUT([Makefile librna/Makefile sissiz/Makefile levmar/Makefile src/Makefile phyml/Makefile examples/Makefile])





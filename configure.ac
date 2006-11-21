AC_INIT(cg-dwf, 1.3.3,[avp@lns.mit.edu])
AC_CONFIG_AUX_DIR(autoconfdir)

dnl removed so that powerpc-bgl-blrts-gnu will work
dnl AC_CANONICAL_HOST

AC_PROG_CC_C_O([bgl-gcc gcc])

AC_LANG([C])
AC_PROG_RANLIB
AC_CHECK_TOOL(AR, ar, [ar])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([cg-dwf-config-internal.h])

AC_ARG_ENABLE(cg-dwf,
	AC_HELP_STRING([--enable-cg-dwf=<cpuarch> <cpuarch>=sse,bluelight,altivec]),
	[ac_ARCH=${enableval}],
	[ac_ARCH="sse"]
)

AC_ARG_WITH(qmp,
  AC_HELP_STRING(
    [--with-qmp=DIR],
    [Build QDP++ on top of QMP, where QMP is installed in DIR.]
  ),
  [QMP_HOME="${with_qmp}"]
)

case ${ac_ARCH} in
sse)
	AC_MSG_NOTICE([Configuring for SSE architecture])
	DWF_CFLAGS="-Wall -msse2 -DUSE_QMP2"
	AC_DEFINE([CG_DWF_ARCH_SSE],[],[Set to indicate target arch is SSE])
	;;
bluelight)
	AC_MSG_NOTICE([Configuring for BlueLight architecture])
	DWF_CFLAGS="-Wall -mbluelight -DUSE_QMP2"
	AC_DEFINE([CG_DWF_ARCH_BLUELIGHT],[],[Set to indicate target arch is bluelight])
	;;
altivec)
	AC_MSG_NOTICE([Configuring for AltiVec architecture])
	DWF_CFLAGS="-Wall -maltivec -DUSE_QMP2"
	AC_DEFINE([CG_DWF_ARCH_ALTIVEC],[],[Set to indicate target arch is AltiVec])
	;;
*)
	AC_MSG_ERROR([Unknown cpu architecture. Use --enable-cpu-arch to select a known one])
	;;
esac

AC_SUBST(DWF_CPUARCH, [${ac_ARCH}])
AC_SUBST(DWF_CFLAGS)
	
AC_MSG_NOTICE([ Checking for QMP ])

if test "X${QMP_HOME}X" = "XX" ; then
   AC_PATH_PROG([QMP_CONFIG], [qmp-config], [])
else
   AC_PATH_PROG([QMP_CONFIG], [qmp-config], [], [${QMP_HOME}/bin:${PATH}])
fi

if test "X${QMP_CONFIG}X" != "XX" ; then
   AC_MSG_NOTICE([Found QMP configuration program ${QMP_CONFIG}])
   AC_SUBST(QMP_HOME,    "${QMP_HOME}")
   AC_SUBST(QMP_CFLAGS,  "`${QMP_CONFIG} --cflags`")
   AC_MSG_NOTICE([QMP compile flags: ${QMP_CFLAGS}])
   AC_SUBST(QMP_LDFLAGS, "`${QMP_CONFIG} --ldflags`")
   AC_MSG_NOTICE([QMP linking flags: ${QMP_LDFLAGS}])
   AC_SUBST(QMP_LIBS,    "`${QMP_CONFIG} --libs`")
   AC_MSG_NOTICE([QMP libraries flags: ${QMP_LIBS}])
else 
   AC_MSG_WARN([QMP configuration program qmp-config not found.])
   AC_MSG_WARN([Set environment variables QMP_CFLAGS QMP_LDFAGS QMP_LIBS
        before configure])
fi
dnl verify we can compile and link against QMP, if needed
PAC_QMP_COMPILE_CC_FUNC(
  [${QMP_CFLAGS}],
  [${QMP_LDFLAGS}],
  [${QMP_LIBS}],
  ,
  ,
  [qmp_compile_ok=yes],
  [qmp_compile_ok=no]
)

AC_MSG_CHECKING([if we can compile  of a simple QMP C program])
if test "X${qmp_compile_ok}X" = "XyesX" ; then
   AC_MSG_RESULT(yes)
else
   AC_MSG_RESULT(no)
   AC_MSG_ERROR([Cannot compile a basic QMP C++ program!
        Check QMP_CFLAGS, QMP_LDFLAGS, QMP_LIBS.])
fi

AM_CONDITIONAL(ARCH_SSE, [test "X${ac_ARCH}X" = "XsseX"])
AM_CONDITIONAL(ARCH_BLUELIGHT, [test "X${ac_ARCH}X" = "XbluelightX"])
AM_CONDITIONAL(ARCH_ALTIVEC, [test "X${ac_ARCH}X" = "XaltivecX"])



AC_CONFIG_FILES([cg-dwf-config],[chmod u+x cg-dwf-config])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT

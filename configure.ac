######################################################################################################
###
###  This file is based on that in the rjags package (version 4-10) with minor modifications
###  [configure stuff relating to modules removed as RoBMA installs within the package directory
###  and minimum JAGS version requirement]
###  Matthew Denwood, 24th February 2021
###  Original file is Copyright (C) Martyn Plummer, licensed under GPL-2
###
######################################################################################################


AC_INIT([RoBMA], [2])

if test -z "${R_HOME}"; then
   AC_MSG_ERROR("R_HOME is not defined")
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

# Try using pkg-config to set compiler flags

PKG_CONFIG_OK="FALSE"

AC_PATH_PROG([PKG_CONFIG], [pkg-config], [], [$PATH:/usr/local/bin])
if test [ -n "${PKG_CONFIG}" ] ; then
   if `${PKG_CONFIG} --exists jags`; then
       if `${PKG_CONFIG} --exists --atleast-version=4.0.0 jags`; then
         JAGS_LIBS=`${PKG_CONFIG} --libs jags`
         JAGS_CFLAGS=`${PKG_CONFIG} --cflags jags`
         JAGS_LIBDIR=`${PKG_CONFIG} --variable libdir jags`
         if test [ -z "${JAGS_LIBS}" -a -z "${JAGS_CFLAGS}" ] ; then
            AC_MSG_WARN([Cannot get information about jags ${PACKAGE_VERSION} from pkg-config])
         else
            AC_MSG_NOTICE([Setting compile and link flags according to pkg-config])
            AC_MSG_NOTICE([Compile flags are $JAGS_CFLAGS])
            AC_MSG_NOTICE([Link flags are $JAGS_LIBS])
            PKG_CONFIG_OK="TRUE"
         fi
      else
         jags_installed_version=`${PKG_CONFIG} --modversion jags`
         AC_MSG_ERROR([pkg-config found JAGS ${jags_installed_version} but rjags requires JAGS >= 4.x.y])
      fi
   else
      AC_MSG_WARN([pkg-config file for jags ${PACKAGE_VERSION} unavailable])
      AC_MSG_WARN([Consider adding the directory containing `jags.pc`])
      AC_MSG_WARN([to the PKG_CONFIG_PATH environment variable])
   fi
fi


########################## Begin legacy configuration #####################
if test "${PKG_CONFIG_OK}" != TRUE; then

AC_MSG_NOTICE([Attempting legacy configuration of rjags])

# Set JAGS_PREFIX
#
# For installations in a non-standard location, set the prefix for
# JAGS manually, using either option --with-jags-prefix or environment
# variable JAGS_PREFIX. Otherwise we get the prefix from the jags
# executable.

AC_ARG_WITH([jags-prefix],
            AS_HELP_STRING([--with-jags-prefix=PREFIX_PATH],
                           [Prefix used when JAGS was installed]),
            [jags_prefix_path=$withval], [jags_prefix_path=""])
if test [ -n "$jags_prefix_path" ] ; then
   JAGS_PREFIX=${jags_prefix_path}
else
   if test [ -z "${JAGS_PREFIX}" ] ; then
   AC_PATH_PROG(jags_bin_path, [jags])
   if test -n "$jags_bin_path"; then
      jags_bin_path=`AS_DIRNAME(["$jags_bin_path"])`
      JAGS_PREFIX=`AS_DIRNAME(["$jags_bin_path"])`
   else
       AC_MSG_ERROR("automatic detection of JAGS failed. Please use pkg-config to locate the JAGS library. See the INSTALL file for details.")
       fi
   fi
fi

# Set include path with option --with-jags-includedir or environment
# variable JAGS_INCLUDEDIR. By default, this is determined from
# JAGS_PREFIX.

AC_ARG_WITH([jags-includedir],
            AS_HELP_STRING([--with-jags-includedir=INCLUDE_PATH],
                           [the value of includedir used when JAGS was installed]),
            [jags_include_path=$withval], [jags_include_path=""])
if test [ -n "$jags_include_path" ] ; then
   JAGS_INCLUDEDIR=${jags_include_path}
else
   if test [ -z "${JAGS_INCLUDEDIR}" ] ; then
      JAGS_INCLUDEDIR=${JAGS_PREFIX}/include
   fi
fi
JAGS_PKGINCLUDEDIR=${JAGS_INCLUDEDIR}/JAGS

AC_LANG_PUSH(C++)
OCPPFLAGS=${CPPFLAGS}
CPPFLAGS=-I${JAGS_PKGINCLUDEDIR}
AC_CHECK_HEADER(Console.h, [],
             [AC_MSG_ERROR("cannot find header file ${JAGS_PKGINCLUDEDIR}/Console.h. Please use pkg-config to locate the JAGS library. See the INSTALL file for details.")])
CPPFLAGS=${OCPPFLAGS}
AC_LANG_POP()

# There is no consensus on where 64- versus 32-bit libraries should be
# installed.  We take a hint from R, assuming that if R is installed
# in, say, /usr/LIBnn/R then we can find the jags library in
# ${JAGS_PREFIX}/LIBnn.
#
# This is a reasonable assumption if both JAGS and R have been
# properly packaged, but may not be true if either has been manually
# installed. The user can set the environment variable LIBnn to
# override this.

if test -z "${LIBnn}"; then
   LIBnn=`"${R_HOME}/bin/R" CMD config LIBnn`
fi

# Set optional library path with option --with-jags-libdir or
# environment variable JAGS_LIBDIR. By default this is determined by
# JAGS_PREFIX and LIBnn.

AC_ARG_WITH([jags-libdir],
            AS_HELP_STRING([--with-jags-libdir=LIB_PATH],
                           [the value of libdir used when JAGS was installed]),
            [jags_lib_path=$withval])
if test [ -n "$jags_lib_path" ] ; then
   JAGS_LIBDIR=$jags_lib_path
else
   if test [ -z "${JAGS_LIBDIR}" ] ; then
      JAGS_LIBDIR=${JAGS_PREFIX}/${LIBnn}
   fi
fi

JAGS_CFLAGS="-I${JAGS_PKGINCLUDEDIR}"
JAGS_LIBS="-L${JAGS_LIBDIR} -ljags"
AC_MSG_NOTICE([Compile flags are $JAGS_CFLAGS])
AC_MSG_NOTICE([Link flags are $JAGS_LIBS])

fi
####################### END legacy configuration ##########################

# Support for private installations of JAGS. If the jags library is not
# on the linker path then it helps to hard-code the location of JAGS into
# the DSO for rjags. Suggested by Doug Bates, April 2011.

AC_ARG_ENABLE([rpath],
            AS_HELP_STRING([--enable-rpath],
                           [Hard link the location of the JAGS library]))
if test [ "x$enable_rpath" = "xyes"]; then
   JAGS_RPATH="-Wl,-rpath ${JAGS_LIBDIR}"
fi

# Sanity checks.
OLDFLAGS=${LDFLAGS}
OCFLAGS=${CFLAGS}
AC_LANG_PUSH(C)
CFLAGS=${JAGS_CFLAGS}
# Can we link to the library?
LDFLAGS="-L${JAGS_LIBDIR}"
AC_CHECK_LIB(jags, jags_version, [],
             [AC_MSG_ERROR("cannot link to JAGS library in ${JAGS_LIBDIR}.")])
AC_LANG_POP(C)
CFLAGS=${OCFLAGS}
# Does the library have the correct major version number?
AC_LANG_PUSH(C++)
OCXXFLAGS=${CXXFLAGS}
CXXFLAGS=${JAGS_CFLAGS}
LDFLAGS="${JAGS_LIBS} ${JAGS_RPATH}"
AC_MSG_CHECKING([version of JAGS library])
# Need to run this trivial test program first to detect runtime linkage problems
AC_RUN_IFELSE(
   [AC_LANG_SOURCE([
      #include <version.h>
      int main() { return 0; }])],
   [],
   [AC_MSG_RESULT([link error])
    AC_MSG_ERROR([Runtime link error: See the INSTALL file for help])])
# Now we can test
AC_LANG_POP()
LDFLAGS=${OLDFLAGS}
CXXFLAGS=${OCXXFLAGS}


AC_SUBST(JAGS_LIBS)
AC_SUBST(JAGS_CFLAGS)
AC_SUBST(JAGS_RPATH)


AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT()

# Initialize Autoconf
AC_INIT([cobrar], [0.2.1], [s.waschina@nutrinf.uni-kiel.de])

# Initialize variables
LIBSBML_LIBS=""
GLPK_LIBS=""
GMP_LIBS=""
LIBSBML_CFLAGS=""
GLPK_CFLAGS=""
GMP_CFLAGS=""

# Check for pkg-config
PKG_CHECK_EXISTS([pkg-config], [HAVE_PKG_CONFIG=yes], [HAVE_PKG_CONFIG=no])

# Check for user-provided libsbml path
AC_ARG_WITH([libsbml],
    AS_HELP_STRING([--with-libsbml=PATH], [Specify libsbml installation path]),
    [
        if test -d "$withval"; then
            LIBSBML_LIBS="-L$withval/lib -lsbml"
            LIBSBML_CFLAGS="-I$withval/include"
        else
            AC_MSG_ERROR([Directory $withval does not exist])
        fi
    ],
    [
        # Check for libsbml
        if test "$HAVE_PKG_CONFIG" = "yes"; then
          # Check for libsbml using pkg-config
    	    PKG_CHECK_MODULES([LIBSBML], [libsbml], , AC_MSG_NOTICE([libsbml not found using pkg-config]))
        fi

        # Last try (no user input, no pkg-config result)
        if test "$LIBSBML_LIBS" == ""; then
            AC_CHECK_LIB([sbml], [main], [
                LIBSBML_LIBS="-lsbml"
            ], [
                AC_MSG_ERROR([libsbml not found])
            ])
        fi
    ]
)


# Check for user-provided glpk path
AC_ARG_WITH([glpk],
    AS_HELP_STRING([--with-glpk=PATH], [Specify glpk installation path]),
    [
        if test -d "$withval"; then
            GLPK_LIBS="-L$withval/lib -lglpk"
            GLPK_CFLAGS="$withval/include"
        else
            AC_MSG_ERROR([Directory $withval does not exist])
        fi
    ],
    [
        # Check for glpk in default and Homebrew paths as before
        AC_CHECK_LIB([glpk], [glp_version], [
            GLPK_LIBS="-lglpk"
        ], [
            if test "$(uname)" = "Darwin"; then
                if test -n "$(command -v brew)"; then
                    BREW_PREFIX=$(brew --prefix)
                    GLPK_LIBS="-L$(find $BREW_PREFIX/Cellar/glpk -name lib) -lglpk"
                    GLPK_CFLAGS="-I$(dirname `find $BREW_PREFIX/Cellar/GLPK -name glpk.h`)"
                else
                    AC_MSG_ERROR([glpk not found. Consider installing it with homebrew via 'brew install glpk'.])
                fi
            fi
        ])
    ]
)

# Check for user-provided gmp path
AC_ARG_WITH([gmp],
    AS_HELP_STRING([--with-gmp=PATH], [Specify gmp (GNU multiple precision arithmetic library) installation path]),
    [
        if test -d "$withval"; then
            GMP_LIBS="-L$withval/lib -lgmp"
            GMP_CFLAGS="$withval/include"
        else
            AC_MSG_ERROR([Directory $withval does not exist])
        fi
    ],
    [
        # Check for gmp in default and Homebrew paths as before
        AC_CHECK_LIB([gmp], [__gmpz_init], [
            GMP_LIBS="-lgmp"
        ], [
            if test "$(uname)" = "Darwin"; then
                if test -n "$(command -v brew)"; then
                    BREW_PREFIX=$(brew --prefix)
                    GMP_LIBS="-L$(find $BREW_PREFIX/Cellar/gmp -name lib) -lgmp"
                    GMP_CFLAGS="-I$(dirname `find $BREW_PREFIX/Cellar/GMP -name gmp.h`)"
                else
                    AC_MSG_RESULT([gmp (optional) not found. Consider installing it with homebrew via 'brew install gmp'.])
                fi
            fi
        ])
    ]
)

# Output variables
AC_SUBST([LIBSBML_LIBS])
AC_SUBST([GLPK_LIBS])
AC_SUBST([GMP_LIBS])
AC_SUBST([LIBSBML_CFLAGS])
AC_SUBST([GLPK_CFLAGS])
AC_SUBST([GMP_CFLAGS])

# output Makevars
AC_CONFIG_FILES([src/Makevars])

# Output summary
AC_MSG_RESULT([libsbml found: $LIBSBML_CFLAGS $LIBSBML_LIBS])
AC_MSG_RESULT([glpk found: $GLPK_CFLAGS $GLPK_LIBS])
AC_MSG_RESULT([gmp found: $GMP_CFLAGS $GMP_LIBS])

# Output messages
AC_OUTPUT

# Initialize Autoconf
AC_INIT([cobrar], [0.1.0], [s.waschina@nutrinf.uni-kiel.de])

# Initialize variables
LIBSBML_LIBS=""
GLPK_LIBS=""
LIBSBML_CFLAGS=""
GLPK_CFLAGS=""

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
        # Check for libsbml in default and Homebrew paths as before
        if test "$HAVE_PKG_CONFIG" = "yes"; then
            # Check for libsbml using pkg-config
    	    PKG_CHECK_MODULES([LIBSBML], [libsbml])
    	    if test "$LIBSBML_LIBS" == ""; then
                AC_MSG_NOTICE([libsbml not found using pkg-config])
            fi
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
                    AC_MSG_ERROR([libsbml not found. Consider installing it with homebrew via 'brew install glpk'.])
                fi
            fi
        ])
    ]
)

# Output variables
AC_SUBST([LIBSBML_LIBS])
AC_SUBST([GLPK_LIBS])
AC_SUBST([LIBSBML_CFLAGS])
AC_SUBST([GLPK_CFLAGS])

# output Makevars
AC_CONFIG_FILES([src/Makevars])

# Output summary
AC_MSG_RESULT([libsbml found: $LIBSBML_CFLAGS $LIBSBML_LIBS])
AC_MSG_RESULT([glpk found: $GLPK_CFLAGS $GLPK_LIBS])

# Output messages
AC_OUTPUT

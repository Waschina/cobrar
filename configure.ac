# Initialize Autoconf
AC_INIT([cobrar], [0.1.0], [s.waschina@nutrinf.uni-kiel.de])

# Initialize variables
LIBSBML_LIBS=""
GLPK_LIBS=""
LIBSBML_INCLUDE_PATH=""
GLPK_INCLUDE_PATH=""

# Check for user-provided libsbml path
AC_ARG_WITH([libsbml],
    AS_HELP_STRING([--with-libsbml=<DIR>], [Specify libsbml installation path]),
    [
        if test -d "$withval"; then
            LIBSBML_LIBS="-L$withval/lib -lsbml"
            LIBSBML_INCLUDE_PATH="$withval/include"
        else
            AC_MSG_ERROR([Directory $withval does not exist])
        fi
    ],
    [
        # Check for libsbml in default and Homebrew paths as before
        AC_CHECK_LIB([sbml], [main], [
            LIBSBML_LIBS="-lsbml"
            # LIBSBML_INCLUDE_PATH="/usr/local/include" # Default path on macOS
        ], [
            if test "$(uname)" = "Darwin"; then
                if test -n "$(command -v brew)"; then
                    BREW_PREFIX=$(brew --prefix)
                    LIBSBML_LIBS="-L$BREW_PREFIX/Cellar/libsbml/lib -lsbml" # CHECK!
                    LIBSBML_INCLUDE_PATH="I$BREW_PREFIX/Cellar/libsbml" # CHECK!
                else
                    AC_MSG_ERROR([libsbml not found. Install libsbml or ensure Homebrew is installed.])
                fi
            else
                AC_MSG_ERROR([libsbml not found. Install libsbml.])
            fi
        ])
    ]
)

# Check for user-provided glpk path
AC_ARG_WITH([glpk],
    AS_HELP_STRING([--with-glpk=<DIR>], [Specify glpk installation path]),
    [
        if test -d "$withval"; then
            GLPK_LIBS="-L$withval/lib -lglpk"
            GLPK_INCLUDE_PATH="$withval/include"
        else
            AC_MSG_ERROR([Directory $withval does not exist])
        fi
    ],
    [
        # Check for glpk in default and Homebrew paths as before
        AC_CHECK_LIB([glpk], [glp_version], [
            GLPK_LIBS="-lglpk"
            # GLPK_INCLUDE_PATH="/usr/local/include" # Default path on macOS
        ], [
            if test "$(uname)" = "Darwin"; then
                if test -n "$(command -v brew)"; then
                    BREW_PREFIX=$(brew --prefix)
                    GLPK_LIBS="-L$($BREW_PREFIX/Cellar/glpk -name lib) -lglpk"
                    GLPK_INCLUDE_PATH="-I$(dirname `find $BREW_PREFIX/Cellar/glpk -name glpk.h`)"
                else
                    AC_MSG_ERROR([glpk not found. Install glpk or ensure Homebrew is installed.])
                fi
            else
                AC_MSG_ERROR([glpk not found. Install glpk.])
            fi
        ])
    ]
)

# Output variables
AC_SUBST([LIBSBML_LIBS])
AC_SUBST([GLPK_LIBS])
AC_SUBST([LIBSBML_INCLUDE_PATH])
AC_SUBST([GLPK_INCLUDE_PATH])

# output Makevars
AC_CONFIG_FILES([src/Makevars])

# Output summary
AC_MSG_RESULT([libsbml found at $LIBSBML_INCLUDE_PATH with $LIBSBML_LIBS])
AC_MSG_RESULT([glpk found at $GLPK_INCLUDE_PATH with $GLPK_LIBS])

# Output messages
AC_OUTPUT

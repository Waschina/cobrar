# Initialize Autoconf
AC_INIT([YourPackage], [1.0], [youremail@example.com])

# Check for GLPK
AC_ARG_WITH([glpk],
    AS_HELP_STRING([--with-glpk=DIR], [Specify the location of GLPK]),
    [GLPK_DIR="$withval"],
    [GLPK_DIR=""])

if test "$GLPK_DIR" = ""; then
    AC_CHECK_LIB([glpk], [glp_create_prob], [
        GLPK_LIBS="-lglpk"
        AC_DEFINE([HAVE_GLPK], [1], [Have GLPK library])
    ], [
        AC_MSG_ERROR([GLPK library not found])
    ])
else
    GLPK_LIBS="-L$GLPK_DIR/lib -lglpk"
    AC_DEFINE([HAVE_GLPK], [1], [Have GLPK library])
fi

# Check for libSBML
AC_ARG_WITH([libsbml],
    AS_HELP_STRING([--with-libsbml=DIR], [Specify the location of libSBML]),
    [LIBSBML_DIR="$withval"],
    [LIBSBML_DIR=""])

if test "$LIBSBML_DIR" = ""; then
    AC_CHECK_LIB([sbml], [SBML_parseL3Formula], [
        LIBSBML_LIBS="-lsbml"
        AC_DEFINE([HAVE_LIBSBML], [1], [Have libSBML library])
    ], [
        AC_MSG_ERROR([libSBML library not found])
    ])
else
    LIBSBML_LIBS="-L$LIBSBML_DIR/lib -lsbml"
    AC_DEFINE([HAVE_LIBSBML], [1], [Have libSBML library])
fi

# Output files
AC_SUBST(LIBSBML_LIBS)
AC_SUBST(GLPK_LIBS)
AC_CONFIG_FILES([src/Makevars])

# Output messages
AC_OUTPUT

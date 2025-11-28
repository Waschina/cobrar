# cobrar installtation instructions to reconfigure R environment in Windows

## R Environment Configuration for Package Building (Windows, Rtools + R 4.4.3)

> [!IMPORTANT]
> This guide describes how to configure your R environment to build R packages successfully on Windows using **R 4.4.3** and **Rtools44** (UCRT64 toolchain). Has also been tested in *R 4.4.1*. If you your system is using a newer Rtools version, remember to replace all "Rtools44" with e.g. "Rtools45".

The changes involve:

-   Installation of the libraries *libSBML* (to read and write files in SBML which is the current model standard) and *GLPK* (linear problem optimization).

-   Modification of the R environment to indicate the location of the libraries so that they can be used for package compilation.

### 0. Prerequisite: Install Rtools44

Before making any configuration changes, ensure that **Rtools44** is installed on your system. **Rtools44** provides *glpk* and enables installation of *libSBML*

Download and install from: <https://cran.r-project.org/bin/windows/Rtools/>

> [!NOTE]
> During installation, make sure to install the **UCRT64 toolchain** and remember the installation path (usually `C:/rtools44`).

### 1. Use Rtools to install **libSBML** and **GLPK**

After installing Rtools, open an **Rtools Bash terminal** , you can launch this APP from the Windows start menu). Within the terminal, run the following command to install **libSBML** and **GLPK**

``` bash
pacman -S mingw-w64-ucrt-x86_64-libsbml mingw-w64-ucrt-x86_64-glpk
```

### 2. Create and/or Edit the `.Renviron` File

In your **R home folder**, create (or edit if it already exists) a text file named `.Renviron`.

**Location:** The R home folder is the folder where R starts by defaulf. It can be obtained after starting R by typing

``` r
getwd()
```

The typical location is

```         
C:/Users/USERNAME/Documents/.Renviron
```

(Make sure that you replace your USERNAME by the relevant information in your system)

**Contents:** The file should define a PATH variable with information on the location of the R installation and of **Rtools**. It should also contain the location of a *Makevars.win* file (that will be generated in the followint steps)

``` r
PATH=C:/Program Files/R/R-4.4.3/bin/x64;${PATH};C:/rtools44/ucrt64/bin
R_MAKEVARS_USER=C:/Users/USERNAME/Documents/.R/Makevars.win
```

Make sure to replace the R version and the USERNAME by the relevant information in your system.

💡 **Notes:**

> -   Always use **forward slashes (`/`)** instead of backslashes (`\`) in R path settings, even on Windows.
> -   The `PATH` variable ensures R and Rtools are both available in your environment.
> -   The `R_MAKEVARS_USER` variable tells R where to find your custom `Makevars.win` file.

### 3. Create or Edit the `Makevars.win` File

Inside your **R home folder** check for a folder called **`.R`** . If it is not present, create one. Inside the **`.R`** folder create or edit a text file named `Makevars.win`:

**Location:**

```         
C:/Users/UNSERNAME/Documents/.R/Makevars.win
```

(Make sure that you replace your USERNAME by the relevant information in your system)

**Contents:**

``` make
UCRT64 = C:/rtools44/ucrt64
PKG_CPPFLAGS += -I"$(UCRT64)/include"
PKG_LIBS     += -L"$(UCRT64)/lib"
PKG_CPPFLAGS += -I"$(UCRT64)/include/sbml/packages"
PKG_LIBS     += -L"$(UCRT64)/lib"
LIBRARY_PATH += -I"$(UCRT64)/include/sbml/packages"
```

> 💡 **Explanation:**
>
> -   `UCRT64` defines the path to the Rtools64 UCRT directory.
>
> -   `PKG_CPPFLAGS` adds include paths for headers (including the required SBML packages).
>
> -   `PKG_LIBS` adds library search paths.
>
> -   `LIBRARY_PATH` extends the linker search path to find additional libraries.
>

### 4. Purpose of These Files

-   **`.Renviron`** — defines environment variables for R at startup, making sure it knows where to find compilers and tools.
-   **`Makevars.win`** — controls compiler and linker settings when you build R packages from source. Together, they ensure R uses the correct build toolchain (Rtools44 / UCRT64) and library paths.

### 5. Restart R and Verify the Setup

Before installing the package, make sure to **restart R** to ensure that the updated environment variables are loaded.You can confirm your environment is correctly configured by running the following in R:

``` r
Sys.getenv("PATH")
Sys.getenv("R_MAKEVARS_USER")
```

The location of the rtools should be shown in the PATH (among many other locations) and the `R_MAKEVARS_USER` should show the location of the Makevars.win

### 6. Install the cobrar Package

Once your environment is configured and verified, you can install the **COBRA-R** package directly from GitHub using the `remotes` package.

**Note:** You may need to install the `remotes` package first:

``` r
install.packages("remotes")
```

Once the `remotes` packages is installed you can just get cobrar by running

``` r
remotes::install_github("Waschina/cobrar")
```

After the installation is completed you can just load the package using

``` r
library(cobrar)
```

If it loads without errors, your configuration is correct and ready for use.

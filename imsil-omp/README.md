IMSIL 2024-03
=============

This README describes version 2025-04 of IMSIL (IMplant and Sputter sImuLator).

Copyright (2025) TU Wien.  
**Use of IMSIL requires a valid license.**
**Contact the principal author to obtain a license.**

Principal author:

> Gerhard Hobler  
> Institute of Solid-State Electronics, TU Wien  
> Gusshausstrasse 25-25a  
> A-1040 Wien  
> AUSTRIA  
> gerhard.hobler@tuwien.ac.at  


Requirements
------------

  * Compilation of the IMSIL source code requires a Fortran 2008 compatible 
    compiler. IMSIL has been tested with gfortran 10.5.0. 

  * Start and postprocessing scripts are written in Python 3. While
    it is possible to run IMSIL without these scripts, they are highly 
    recommended and so is an installation of Python 3. 

  * To run the `compare` script in the `tests` directory, Matplotlib and - for 
    the 3-D test example - Mayavi are required. On our system, Mayavi is 
    currently broken. This is an issue to be solved in the future.
  
  * If you want to compile the manual yourself, a latex engine is required.
    The makefile provided in the distribution uses pdflatex. However, a pdf
    of the manual is included, so it is usually not necessary to compile the
    manual.
    
  * Makefiles have been tested with GNU make.

IMSIL has been tested on workstations running Ubuntu 24.04. It should be 
possible to use other Linux distributions and compilers. IMSIL should also 
run under Windows with minor changes, but this has not been tested recently.


Directories
-----------

The install directory (`IMSIL_xxx`, where `xxx` denotes the version) contains 
the following subdirectories:

  * `code` ... the source code (Fortran)
  * `manual` ... the manual (LaTeX and a pdf)
  * `scripts` ... scripts that simplify running IMSIL and provide functions
      for reading and plotting results (Python)
  * `tables` ... tables needed for running IMSIL if specified in the input
      file. The most prominent example of a table is the scattering 
      table file SCATTAB.
  * `tests` ... some tests for checking successful installation. No results.
  * `tests_ref` ... same as `tests`, but including results. Do not run the
      examples in this directory but keep them as a reference.

If you plan to use the IMSIL start script (`scripts/imsil.py`), there are
certain requirements on where the executable and the working directories
(= where input and output files are kept) may reside. A simple solution is
to generate the executable code in the same place as the source code and to
create the working directory in the installation directory (e.g., 
`IMSIL_xxx/work`). So the directory structure could look like

        ~/programs/IMSIL_xxx  
            code        # source + executable code  
            manual      # tex files + pdf
            scripts     # scripts
            tables      # tables
            tests       # tests with own results
            tests_ref   # tests with reference results
            work        # (and subdirectories) input and output data
   
An alternative would be

        ~/programs/IMSIL_xxx  
            code        # source code  
            manual      # tex files  
            scripts     # scripts
            tables      # tables
            tests       # tests without results
            tests_ref   # tests with reference results
  
        ~/work/IMSIL_xxx  
            code        # executable  
            manual      # pdf  
            scripts     # link to ~/programs/IMSIL_xxx/scripts  
            tables      # link to ~/programs/IMSIL_xxx/tables  
            tests       # copy of ~/programs/IMSIL_xxx/tests + own results
            tests_ref   # link to ~/programs/IMSIL_xxx/tests_ref  
            project1    # (and subdirectories) IMSIL input and output data  
            project2    # (and subdirectories) IMSIL input and output data  
            ...  

Installation
------------

To install IMSIL, the source code must be compiled. A `makefile` is provided 
in the `code` directory. To use gfortran (recommended), just run `make` from 
the code directory. To produce a parallel version, run `make OMP=omp` (OpenMP 
has to be installed on your computer). To use the Intel compiler, source the 
`setvars.sh` script and run `make FC=ifort`. To produce an OpenMP version, add 
`OMP=omp`. In all cases, to produce a debug version, add `DB=db` to the make 
command.

There is also a `makefile` in `manual` for building the manual from LaTeX 
files. Unless you change the manual, there is no need to do so since a pdf 
`manual.pdf` is also provided. Another `makefile` is in the tests directory,
which can be used to run the tests.

In the makefile and associated files, and in the `imsil.py` script of the 
`scripts` directory, the following paths have to be adjusted:

  * `code/dirs.mk` and `code/dirs_db.mk`: These files are included by the 
    makefile when the optimized and debug versions are made, respectively.
    They define the paths of the source code (`SRCDIR`) and where the 
    executable is built (`EXEDIR`). In the first example of a directory 
    structure given above, set
    
            SRCDIR=~/programs/IMSIL_xxx/code
            EXEDIR=~/programs/IMSIL_xxx/code

       In the second example, this would be
    
            SRCDIR=~/programs/IMSIL_xxx/code
            EXEDIR=~/work/IMSIL_xxx/code
    
  * `manual/makefile`: This `makefile` uses the pdflatex engine. If you want 
    to use a different engine, adjust. Moreover, `OUTDIR` must be set to where
    you want to have the pdf, e.g., in the first directory structure
    
            OUTDIR=~/programs/IMSIL_xxx/manual
    
    in the second one

            OUTDIR=~/work/IMSIL_xxx/manual
            
    Note: Because of a peculiarity of pdflatex you must replace `~` by your
    home directory here.

  * `scripts/imsil.py`: In the IMSIL start script, `installdir` must be set to
    the path of `IMSIL_xxx` under which the start script is started, i.e., to 

            PROGRAM = 'IMSIL_xxx'
            installdir = os.path.join( os.path.expanduser('~'), 'programs',
                                       PROGRAM )
    
    in the first example and 
    
            PROGRAM = 'IMSIL_xxx'
            installdir = os.path.join( os.path.expanduser('~'), 'work',
                                       PROGRAM )
    
    in the second example.
    
  * `tests/makefile`: Replace `imsil` and `mimsil` by 
    `python ../scripts/imsil.py` and `python ../scripts/mimsil.py`, 
    respectively, or define aliases for `imsil` and `mimsil` (recommended).


Running IMSIL
-------------

#### Serial vs. Parallel processing

The `makefile` in the `code` directory, without the `OMP` parameter, produces 
the IMSIL executable `imsil.exe`. It processes trajectories sequentially. With 
the parameter `OMP=omp`, `make` procudes the IMSIL executable `imsil.omp.exe`.
Before running `imsil.omp.exe`, the environment variable `OMP_NUM_THREADS` has 
to be set to the number of parallel threads. This will differ according to how 
the code is run, see below.

This can be done via a command line
argument (`-t <n>` where `<n>` denotes the number of threads desired) when using 
the `(m)imsil.py` script ().

#### Using the `exe` file directly

The most straightforward way to run IMSIL is to call the executable with the 
input file `INP` and any other required files (such as `SCATTAB`) or
links to them in the current working directory. IMSIL will create all output 
files in the current directory. As a bare minimum, the `OUT` file will be 
created. The `OUT` file contains a copy of the input parameters and some 
derived parameters as well as some calculated quantities like yields and 
ranges. Other output files are described in the manual. They all have uppercase 
names. Before the parallel version is run, the environment variable 
`OMP_NUM_THREADS` has to be set from the command line to the number of desired 
parallel threads.

#### Using the `imsil.py` script

A more convenient way to run IMSIL is by using the `imsil.py` script, which is
provided in the `scripts` directory. With `imsil.py` the input file can have 
an arbitrary name with extension `.inp`. It is invoked by

        python <path-to-imsil.py> [options] [input-file ...]

The output files then have the same name as the input file with lowercase 
extensions corresponding to the filenames generated by the IMSIL executable. 
Moreover, the path of the tables directory is determined automatically and 
communicated to IMSIL.

With `imsil.py`, IMSIL is run in a subdirectory of the current directory.
This has several advantages: Multiple simulations can be run in the same 
directory. When a simulation is repeated, the old output files remain
unaffected till the end of the simulation. And the existence of the
subdirectory indicates that the simulation is still running. Given more
than one input files, the simulations can be run in parallel using the `-p`
flag. Each simulation can be run with multiple threads using the `-t` option.
In this case, `imsil.py` will set the environment variable `OMP_NUM_THREADS`
for you. Note that the `OMP_NUM_THREADS` is overridden when the number of 
threads is specified in the input file.

#### Using the `mimsil.py` script

`mimsil.py`, which is also provided in the `scripts` directory, allows to run 
a set of input files that are derived from a template. It calls `imsil.py`,
therefore to run `mimsil.py` requires to set `installdir` in `imsil.py`.
`mimsil.py` is typically used to run simulations that differ only in a few
input parameters. It requires a template file (`fname.tem`) and a parameters 
file (`fname.par`). It is invoked by

        python <path-to-mimsil.py> [options] [template-file ...]

or

        python <path-to-mimsil.py> [options] [parameters-file ...]

where the options are the same as for the `imsil.py` script.

The template file is an input file where the values of some parameters have 
been replaced by variables. Variables are names enclosed in `#` characters, 
like `#parameter#`. The values of the variables are defined in the parameters 
file. It may have either of two formats. The simpler of the two is

        FILE: <input/output file name without extension>
        parameter1=value1
        parameter2=value2
        ...
        
        FILE: <input/output file name without extension>
        parameter1=value3
        parameter2=value4
        ...
        
        ...
   
See Section 6.1 of the the manual for the alternative format and further 
information about the `imsil.py` and `mimsil.py` scripts.

See Section 7 of the manual and the `tests` directory for examples of template 
files and both formats of the parameters file.


Postprocessing
--------------

The `scripts` directory also provides the `read_out.py` file, which defines
Python functions for reading the most important types of IMSIL output files.
The functions defined in `read_out.py` are meant to be called from user-defined
Python scripts, typically but not necessarily to plot the results.
`read_out.py` may also be used as a script 
(`python /path/to/read_out.py fname.xxx`) to do best-guess postprocessing of 
an IMSIL output file.

See Section 6.2 of the manual for details on the `read_out.py` module.


Testing IMSIL
-------------

The `tests` subdirectory contains several subdirectories with example input
or template and parameters files. These can be run either directly from the
subdirectories, or from the `tests` directory using `makefile`:

        make <subdirectory>

to run the tests in a particular subdirectory, or

        make all

to run all tests at once. To start with, it is recommended to `make rrobin`,
since the `rrobin` jobs take only a few seconds. You may want to use the 
`-j` option of `make` to run simulations in parallel or to pass the `t` 
parameter to `make` using `t=<n>` to use `n` threads (cores) in parallel 
for each simulation.

Note that due to the way `make` is set up, the tests are only run if the `.out`
files or `.log` files are older than the `.inp`, `.tem`, and `.par` files, and 
the imsil executable. Thus, if you just want to repeat a simulation with a 
different number of threads, you must touch one of the latter files. 
Alternatively, you may touch the `imsil.exe` file.

The results may be plotted with the `compare` script:

        ./compare <filename>
   
where `<filename>` is the path to a Python file in one of the subdirectories
(they contain the plot scripts). Alternatively, to run all Python scripts in
subdirectory `<directory>`:

        ./compare <directory>

or to run all Python scripts in all subdirectories

        ./compare all

If the `tests_ref` directory exists as a sibling of `tests` and if it 
contains matching simulation results, these are displayed in the plots as
well (therefore the script name `compare`). Note that if the compared 
results are exactly the same, they usually cannot be distinguished in the 
plots. To confirm that indeed two versions of each result 
are plotted, the console output may be checked.

In most cases the plots also include difference plots between the two results 
compared, where shaded areas indicate differences larger than the mean 
expected difference due to statistical fluctuations. There are some caveats 
with this, however. First, by definition the difference must sometimes be 
larger than its mean expectation value. Second, the difference plot may be 
distorted if the histogram bins are not the same between the two results 
compared. Furthermore, the estimate of the mean difference between point 
defect profiles are usually underestimates since point defect generation is 
not statistically independent within a single cascade.


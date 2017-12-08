Misc:
    Units:
        1. Internal units are specified in the .yml input files.
        2. Specified in initial conditions to control their units (to cgs) and
        the units for snapshots. If these are not given, then assumed to be the
        same as the input internal units.


To print anything from the code normally, use
    mesasge( . . . )
        just like printf().
    error( . . . )
        similarly but it will stop the code as well.


SWIFT tour:
    Code in src/

    Input files, initial conditions, and python scripts in examples/
        examples/swift is the executable (swift_mpi for >~1e8 particles)
        examples/main.c

    scr/equation_of_state
    src/hydro/
        Contains folders that have all the SPH information.
        When compiling

    *.yml (YAML language) input file (not init cond)
        See examples/parameter_example.yml and other examples . . .

    src/parser.h
        Struct that contains all the info from the .yml

        To add a new input, just add it (properly) to the .yml, don't need to
        change the parser.

    src/runner_doiact.h
        All the neighbour searching and interaction things.
        e.g. DOPAIR1()
            Interacts two cells' particles together, calls the interaction
            function IACT_NONSYM() which will connect to the function defined
            in the chosen hydro/ folder.

    Default values for optional .yml parameters are defined at the top of the
    relevant files that call the parser, not consolidated to any one place.


Compiling:
    $ ./configure --help (Copied below in full...)

    --disable-vec
        to compile on cosma4

    --with-hydro=<scheme>
        to choose which folder in src/hydro/ to use


Modules:
    swift
    swift/c4/gcc/intelmpi/5.1.3
        which will automatically load others like fftw
        (for cosma4; replace c4/ with c5/ for cosma5)


Git basics:
	git checkout HASHNAME <or> BRANCHNAME
		Switch over to that commit or branch
	git add FILENAME
		Stages changes to a file
	git status
		show results of the above
	git checkout -- FILENAME <or> .
		Returns to the Head commit (for when something goes wrong!)
	git diff (FILENAME)
		Compares unstaged and committed files
	git diff --cached (FILENAME)
		Compares staged and committed files
	git commit
		Opens in vim, or -m " . . . ", to write commit message
	git push origin BRANCHNAME
		Push the commit to the branch
	git pull upstram BRANCHNAME

    To sync recent changes others have to the main master SWIFT thing:
	git checkout master
		Go to my master
	git pull upstream master
		Pull their version to my local one
	git push origin master
		Push those changes to my repository
	git checkout giant_impacts
        Switch to the giant_impacts branch
	git merge master
		Including commits from master (from upstream) to local repository
	git push origin giant_impacts
		Push changes to my remote repository


$ ./configure --help
    `configure' configures SWIFT 0.6.0 to adapt to many kinds of systems.

    Usage: ./configure [OPTION]... [VAR=VALUE]...

    To assign environment variables (e.g., CC, CFLAGS...), specify them as
    VAR=VALUE.  See below for descriptions of some of the useful variables.

    Defaults for the options are specified in brackets.

    Configuration:
      -h, --help              display this help and exit
          --help=short        display options specific to this package
          --help=recursive    display the short help of all the included packages
      -V, --version           display version information and exit
      -q, --quiet, --silent   do not print `checking...' messages
          --cache-file=FILE   cache test results in FILE [disabled]
      -C, --config-cache      alias for `--cache-file=config.cache'
      -n, --no-create         do not create output files
          --srcdir=DIR        find the sources in DIR [configure dir or `..']

    Installation directories:
      --prefix=PREFIX         install architecture-independent files in PREFIX
                              [/usr/local]
      --exec-prefix=EPREFIX   install architecture-dependent files in EPREFIX
                              [PREFIX]

    By default, `make install' will install all the files in
    `/usr/local/bin', `/usr/local/lib' etc.  You can specify
    an installation prefix other than `/usr/local' using `--prefix',
    for instance `--prefix=$HOME'.

    For better control, use the options below.

    Fine tuning of the installation directories:
      --bindir=DIR            user executables [EPREFIX/bin]
      --sbindir=DIR           system admin executables [EPREFIX/sbin]
      --libexecdir=DIR        program executables [EPREFIX/libexec]
      --sysconfdir=DIR        read-only single-machine data [PREFIX/etc]
      --sharedstatedir=DIR    modifiable architecture-independent data [PREFIX/com]
      --localstatedir=DIR     modifiable single-machine data [PREFIX/var]
      --libdir=DIR            object code libraries [EPREFIX/lib]
      --includedir=DIR        C header files [PREFIX/include]
      --oldincludedir=DIR     C header files for non-gcc [/usr/include]
      --datarootdir=DIR       read-only arch.-independent data root [PREFIX/share]
      --datadir=DIR           read-only architecture-independent data [DATAROOTDIR]
      --infodir=DIR           info documentation [DATAROOTDIR/info]
      --localedir=DIR         locale-dependent data [DATAROOTDIR/locale]
      --mandir=DIR            man documentation [DATAROOTDIR/man]
      --docdir=DIR            documentation root [DATAROOTDIR/doc/swift]
      --htmldir=DIR           html documentation [DOCDIR]
      --dvidir=DIR            dvi documentation [DOCDIR]
      --pdfdir=DIR            pdf documentation [DOCDIR]
      --psdir=DIR             ps documentation [DOCDIR]

    Program names:
      --program-prefix=PREFIX            prepend PREFIX to installed program names
      --program-suffix=SUFFIX            append SUFFIX to installed program names
      --program-transform-name=PROGRAM   run sed PROGRAM on installed program names

    System types:
      --build=BUILD     configure for building on BUILD [guessed]
      --host=HOST       cross-compile to build programs to run on HOST [BUILD]

    Optional Features:
      --disable-option-checking  ignore unrecognized --enable/--with options
      --disable-FEATURE       do not include FEATURE (same as --enable-FEATURE=no)
      --enable-FEATURE[=ARG]  include FEATURE [ARG=yes]
      --enable-debug=[yes/info/profile/no]
                              compile with debugging
      --disable-dependency-tracking  speeds up one-time build
      --enable-dependency-tracking   do not reject slow dependency extractors
      --enable-ipo            Enable interprocedural optimization [no/yes]
      --enable-mpi            Compile with functionality for distributed-memory
                              parallelism using MPI [yes/no]
      --enable-shared[=PKGS]  build shared libraries [default=yes]
      --enable-static[=PKGS]  build static libraries [default=yes]
      --enable-fast-install[=PKGS]
                              optimize for fast installation [default=yes]
      --disable-libtool-lock  avoid locking (might break parallel builds)
      --enable-task-debugging Store task timing information and generate task dump
                              files [yes/no]
      --enable-threadpool-debugging
                              Store threadpool mapper timing information and
                              generate threadpool dump files [yes/no]
      --enable-timers         Activate the basic timers [yes/no]
      --enable-debugging-checks
                              Activate expensive consistency checks [yes/no]
      --enable-naive-interactions
                              Activate use of naive cell interaction functions
                              [yes/no]
      --enable-gravity-force-checks
                              Activate expensive brute-force gravity checks for a
                              fraction 1/N of all particles [N]
      --enable-no-gravity-below-id
                              Zeros the gravitational acceleration of all
                              particles with an ID smaller than [N]
      --enable-optimization   Enable compile time optimization flags for host
                              [yes/no]
      --disable-vec           Disable vectorization
      --enable-portable-binary
                              disable compiler optimizations that would produce
                              unportable binaries
      --enable-sanitizer      Enable memory error detection using address
                              sanitizer [no/yes]
      --enable-parallel-hdf5  Enable parallel HDF5 library MPI functions if
                              available. [yes/no]
      --enable-compiler-warnings
                              Enable compile time warning flags, if compiler is
                              known [error/no/yes)]
      --disable-doxygen-doc   don't generate any doxygen documentation
      --enable-doxygen-dot    generate graphics for doxygen documentation
      --enable-doxygen-man    generate doxygen manual pages
      --enable-doxygen-rtf    generate doxygen RTF documentation
      --enable-doxygen-xml    generate doxygen XML documentation
      --enable-doxygen-chm    generate doxygen compressed HTML help documentation
      --enable-doxygen-chi    generate doxygen seperate compressed HTML help index
                              file
      --disable-doxygen-html  don't generate doxygen plain HTML documentation
      --disable-doxygen-ps    don't generate doxygen PostScript documentation
      --disable-doxygen-pdf   don't generate doxygen PDF documentation

    Optional Packages:
      --with-PACKAGE[=ARG]    use PACKAGE [ARG=yes]
      --without-PACKAGE       do not use PACKAGE (same as --with-PACKAGE=no)
      --with-pic              try to use only PIC/non-PIC objects [default=use
                              both]
      --with-gnu-ld           assume the C compiler uses GNU ld [default=no]
      --with-gcc-arch=<arch>  use architecture <arch> for gcc -march/-mtune,
                              instead of guessing
      --with-metis=PATH       root directory where metis is installed [yes/no]
      --with-tcmalloc         use tcmalloc library or specify the directory with
                              lib [yes/no]
      --with-profiler         use cpu profiler library or specify the directory
                              with lib [yes/no]
      --with-jemalloc         use jemalloc library or specify the directory with
                              lib [yes/no]
      --with-hdf5=PATH        location of h5cc or h5pcc for HDF5 configuration
      --with-hydro=<scheme>   Hydro dynamics to use [gadget2, minimal, hopkins,
                              default, gizmo, shadowfax debug default: gadget2]
      --with-kernel=<kernel>  Kernel function to use [cubic-spline,
                              quartic-spline, quintic-spline, wendland-C2,
                              wendland-C4, wendland-C6 default: cubic-spline]
      --with-hydro-dimension=<dim>
                              dimensionality of problem [3/2/1 default: 3]
      --with-equation-of-state=<EoS>
                              equation of state [ideal-gas, isothermal-gas
                              default: ideal-gas]
      --with-adiabatic-index=<gamma>
                              adiabatic index [5/3, 7/5, 4/3, 2 default: 5/3]
      --with-riemann-solver=<solver>
                              riemann solver (gizmo-sph only) [none, exact, trrs,
                              hllc, default: none]
      --with-cooling=<function>
                              cooling function [none, const-du, const-lambda,
                              grackle default: none]
      --with-ext-potential=<pot>
                              external potential [none, point-mass, isothermal,
                              softened-isothermal, disc-patch, sine-wave default:
                              none]
      --with-multipole-order=<order>
                              order of the multipole and gravitational field
                              expansion [ default: 5]

    Some influential environment variables:
      CC          C compiler command
      CFLAGS      C compiler flags
      LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
                  nonstandard directory <lib dir>
      LIBS        libraries to pass to the linker, e.g. -l<library>
      CPPFLAGS    C/C++/Objective C preprocessor flags, e.g. -I<include dir> if
                  you have headers in a nonstandard directory <include dir>
      CPP         C preprocessor
      MPICC       MPI C compiler command
      MPIRUN      Path to the mpirun command if non-standard
      DOXYGEN_PAPER_SIZE
                  a4wide (default), a4, letter, legal or executive

    Use these variables to override the choices made by `configure' or to help
    it to find libraries and programs with nonstandard names/locations.

    Report bugs to <https://gitlab.cosma.dur.ac.uk/swift/swiftsim>.








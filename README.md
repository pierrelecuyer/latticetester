# Lattice Tester

A C++ library for measuring the uniformity of integral lattices in the real space

## What this software is about

_Lattice Tester_ is a C++ software library to compute theoretical measures
of uniformity (figures of merit) for lattices in the $t$-dimensional integer space $Z^t$
(all lattice points have integer coordinates).
This integral property is often obtained after rescaling an original lattice by an integer factor $m > 1$.
Such lattices are encountered in particular in quasi-Monte Carlo integration
by lattice rules and in the analysis of uniform random number generators
defined by linear recurrences modulo a large integer.
Measures of uniformity include the length of the shortest nonzero vector 
in the lattice or in its dual (the spectral test), the Beyer ratio, 
as well as figures of merit that take normalized versions
of these measures over projections of the lattice on subsets of the $t$ coordinates,
and then take a weighted sum or the worst-case over the class of considered projections.

_Lattice Tester_ was built primarily as a base library for the software packages
[LatNet Builder](https://github.com/umontreal-simul/latbuilder),
and [LatMRG](https://github.com/umontreal-simul/latmrg), designed to construct and analyze
lattice rules (for quasi-Monte-Carlo) and multiple recursive linear 
congruential random number generators, respectively. 
It is also intended to be used for other applications related to lattices in the integer space.
It was developed to work under Linux environments (mostly because of NTL).

This version is a significant overhaul compared to the previous official version given 
[**HERE**](http://umontreal-simul.github.io/latticetester/).  
Some polishing remains to be done; it should be completed by February 2025.


## Documentation

A detailed description of _Lattice Tester_, its underlying theory, its organization, and examples, can be found in the 
[**Lattice Tester Guide**](https://www-labs.iro.umontreal.ca/~lecuyer/myftp/papers/lattester-guide-2025.pdf).
That .pdf document is the core documentation. 

The interface is specified in the 
[**API documentation**](http://pierrelecuyer.github.io/latticetester/namespaces.html).
This part still needs work.  It should be completed by February 2025.

## Compiling and Building

### Software Dependencies

Compiling *Lattice Tester* requires the following software to be installed:

* [GMP](https://gmplib.org/) compatible version with your NTL installation
* [NTL](http://www.shoup.net/ntl/index.html) 11.5.1 or later
* [Git](http://git-scm.com/) *(optional for downloading the source code)*
* [Python](https://www.python.org/) *(Needed by waf to compile and build the library)*
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) *(optional for generating
  the API documentation)*

You also need a recent compiler compliant with (at least) the C++14 standard.

### Configuring the Build

*Lattice Tester* relies on the
[waf meta build system](https://code.google.com/p/waf/) for configuring and
compiling the software source. Waf is included in the *Lattice Tester* source 
tree, but it depends on [Python](http://python.org/download), which must be 
available on the system on which *Lattice Tester* is to be compiled.

The commands below should work verbatim under Linux and MacOS systems.
*Microsoft Windows* users should replace every instance of `./waf` 
with the path under which the Python executable
(`python.exe`) or simply with `python waf`
if the Python installation path is accessible from the system `%PATH%`
environment variable.

Change the current directory to the root directory of the package, for example:

    cd git/latticetester

if you obtained the source code via `git`.
If you obtained the source code from the ZIP archive, the directory should be
named `latticetester-master` instead of `latticetester`.
At the root of the source tree, the `wscript` file contains the waf script that manages the build process.

To see the various commands and options, try

	./waf --help

First, the project must be configured with:

   ./waf configure --flags

A simple 

   ./waf configure

command (with no flags) should be enough to configure `waf` for a minimal build. 
Otherwise, some of the commonly used options include:

- `--out /path/to/build/location` specifies in which directory the
  build process will operate. The default is `./build`. You need permission
  to write in that directory.
- `--prefix /path/to/installation/location` specifies in which 
  directory *Lattice Tester* will be installed after it's compilation.
  The default is `/usr/local` on Linux (waf's default). You need permission
  to write in that directory.
- `--gmp /path/to/gmp` specifies the location of the gmp installation. 
  This is needed only if waf does not find gmp automatically.
- `--ntl /path/to/NTL` specifies the location of the NTL installation. 
  This is needed only if waf does not find NTL automatically.
 - `--build-docs` waf will build the documentation if and only if this flag is specified 
  in the configure step. 
- `--link-static` if this flag is specified, the compiler will link all the 
  libraries statically to the executable programs (the examples). 
  This might be convenient if NTL is installed in non standard paths.

To build the html documentation, 
[Doxygen](http://www.stack.nl/~dimitri/doxygen/) must be available on the system
and you must first configure with

   ./waf configure --build-docs

The waf script in `doc/wscript` manages the documentation build.  
The source of the main page is in `doc/dox/main.dox`.
The Doxygen options are selected in the file `Doxyfile.in`. 
For example, one may select if members from the `.cc` files are included or not
by changing the `FILE_PATTERNS` option, select if the `#include` statements at the head of a file
are shown or not with the `SHOW_INCLUDE_FILES` option, etc. 
The built documentation is placed in `build/doc/html/` with `index.html` as its main entry.
To deploy the documentation on the `pierrelecuyer.github.io/latticetester/` GitHub pages,
the contents of this `html` directory must be moved to the `gh-pages` branch of `latticetester`. 
The `README.md` in that branch provides more details on how to proceed for this deployment.

If NTL, and GMP are not part of the standard system installation and were
manually installed under, say, the `/opt/ntl`, and `/opt/gmp` directories,
which means that `/opt/ntl` and `/opt/gmp` all contain subdirectories named
`include` and `lib`, you should do:

    ./waf configure --ntl /opt/ntl --gmp /opt/gmp

It is possible to set the `CXX` environment variable to the path to a specific
C++ compiler to be used to build Lattice Tester, before running the `waf configure` command.

In a UNIX shell, it is also possible to write and run a simple `configure.sh`
script with `./configure.sh` to avoid typing the configure command by hand. 
This can be useful if you have many flags to include.

### Building and Installing

Once everything is configured correctly, you can build the
*Lattice Tester* library via:

    ./waf build

If the build process completed without errors, *Lattice Tester* can be installed to the
directory specified with the `--prefix` option during the configuration step, with:

    ./waf install

## Authors

François Blouin, Erwan Bourceret, Anna Bragina, Ajmal Chaumun, 
Raymond Couture, Marco Jacques, David Munger, François Paradis, Marc-Antoine Savard, Richard Simard, 
Mamadou Thiongane, Josée Turgeon, and Christian Weiss
have contributed to various versions of this software since around 1986,
under the lead of Pierre L'Ecuyer.

## License

_Lattice Tester_ is free open source software, distributed under the Apache 2.0 License.

# Lattice Tester

A C++ library for measuring the uniformity of integral lattices in the real space

## What this software is about

_Lattice Tester_ is a C++ software library to compute theoretical measures
of uniformity (figures of merit) for lattices in the $t$-dimensional integer space $Z^t$
(all lattice points have integer coordinates).
This integral property is often obtained after rescaling an original lattice by an integer factor $m > 1$.
Such lattices are encountered in particular in quasi-Monte Carlo integration
by lattice rules and in the analysis of uniform random number generators
defined by linear recurrences modulo a large integer.
Measures of uniformity include the length of the shortest nonzero vector 
in the lattice or in its dual (the spectral test), the Beyer ratio, 
as well as figures of merit that take normalized versions
of these measures over projections of the lattice on subsets of the $t$ coordinates,
and then take a weighted sum or the worst-case over the class of considered projections.

_Lattice Tester_ was built primarily as a base library for the software packages
[LatNet Builder](https://github.com/umontreal-simul/latbuilder),
and [LatMRG](https://github.com/umontreal-simul/latmrg), designed to construct and analyze
lattice rules (for quasi-Monte-Carlo) and multiple recursive linear 
congruential random number generators, respectively. 
It is also intended to be used for other applications related to lattices in the integer space.
It was developed to work under Linux environments (mostly because of NTL).

This version is a significant overhaul compared to the previous official version given 
[**HERE**](http://umontreal-simul.github.io/latticetester/).  
Some polishing remains to be done; it should be completed by February 2025.


## Documentation

A detailed description of _Lattice Tester_, its underlying theory, its organization, and examples, can be found in the 
[**Lattice Tester User's Guide**](https://www-labs.iro.umontreal.ca/~lecuyer/myftp/papers/lattester-guide-2025.pdf).
That .pdf document is the core documentation. 

The interface is specified in the 
[**API documentation**](http://pierrelecuyer.github.io/latticetester/namespaces.html).
This part still needs work.  It should be completed by February 2025.

## Compiling and Building

### Software Dependencies

Compiling *Lattice Tester* requires the following software to be installed:

* [GMP](https://gmplib.org/) compatible version with your NTL installation
* [NTL](http://www.shoup.net/ntl/index.html) 11.5.1 or later
* [Git](http://git-scm.com/) *(optional for downloading the source code)*
* [Python](https://www.python.org/) *(Needed by waf to compile and build the library)*
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) *(optional for generating
  the API documentation)*

You also need a recent compiler compliant with (at least) the C++14 standard.

### Configuring the Build

*Lattice Tester* relies on the
[waf meta build system](https://code.google.com/p/waf/) for configuring and
compiling the software source. Waf is included in the *Lattice Tester* source 
tree, but it depends on [Python](http://python.org/download), which must be 
available on the system on which *Lattice Tester* is to be compiled.

The commands below should work verbatim under Linux and MacOS systems.
*Microsoft Windows* users should replace every instance of `./waf` 
with the path under which the Python executable
(`python.exe`) or simply with `python waf`
if the Python installation path is accessible from the system `%PATH%`
environment variable.

Change the current directory to the root directory of the package, for example:

    cd git/latticetester

if you obtained the source code via `git`.
If you obtained the source code from the ZIP archive, the directory should be
named `latticetester-master` instead of `latticetester`.
At the root of the source tree, the `wscript` file contains the waf script that manages the build process.

To see the various commands and options, try

	./waf --help

First, the project must be configured with:

   ./waf configure --flags

A simple 

   ./waf configure

command (with no flags) should be enough to configure `waf` for a minimal build. 
Otherwise, some of the commonly used options include:

- `--out /path/to/build/location` specifies in which directory the
  build process will operate. The default is `./build`. You need permission
  to write in that directory.
- `--prefix /path/to/installation/location` specifies in which 
  directory *Lattice Tester* will be installed after it's compilation.
  The default is `/usr/local` on Linux (waf's default). You need permission
  to write in that directory.
- `--gmp /path/to/gmp` specifies the location of the gmp installation. 
  This is needed only if waf does not find gmp automatically.
- `--ntl /path/to/NTL` specifies the location of the NTL installation. 
  This is needed only if waf does not find NTL automatically.
 - `--build-docs` waf will build the documentation if and only if this flag is specified 
  in the configure step. 
- `--link-static` if this flag is specified, the compiler will link all the 
  libraries statically to the executable programs (the examples). 
  This might be convenient if NTL is installed in non standard paths.

To build the html documentation, 
[Doxygen](http://www.stack.nl/~dimitri/doxygen/) must be available on the system
and you must first configure with

   ./waf configure --build-docs

The waf script in `doc/wscript` manages the documentation build. The source of the main page is in `doc/dox/main.dox`.
The Doxygen options are selected in the file `doc/Doxyfile.in`. For example, one may select if members from the `.cc` files are included or not
by changing the `FILE_PATTERNS` option, select if the `#include` statements at the head of a file
are shown or not with the `SHOW_INCLUDE_FILES` option, etc. 
The built documentation is placed in `build/doc/html/` with `index.html` as its main entry.
To deploy the documentation on the `pierrelecuyer.github.io/latticetester/` GitHub pages,
the contents of this `html` directory must be moved to the `gh-pages` branch of `latticetester`. 
The `README.md` in that branch provides more details on how to proceed for this deployment.

If NTL and GMP are not part of the standard system installation and were
manually installed under, say, the `/opt/ntl`, and `/opt/gmp` directories,
which means that `/opt/ntl` and `/opt/gmp` all contain subdirectories named
`include` and `lib`, you should do:

    ./waf configure --ntl /opt/ntl --gmp /opt/gmp

It is possible to set the `CXX` environment variable to the path to a specific
C++ compiler to be used to build Lattice Tester, before running the `waf configure` command.

In a UNIX shell, it is also possible to write and run a simple `configure.sh`
script with `./configure.sh` to avoid typing the configure command by hand. 
This can be useful if you have many flags to include.

### Building and Installing

Once everything is configured correctly, you can build the
*Lattice Tester* library via:

    ./waf build

If the build process completed without errors, *Lattice Tester* can be installed to the
directory specified with the `--prefix` option during the configuration step, with:

    ./waf install

## Authors

François Blouin, Erwan Bourceret, Anna Bragina, Ajmal Chaumun, 
Raymond Couture, Marco Jacques, David Munger, François Paradis, Marc-Antoine Savard, Richard Simard, 
Mamadou Thiongane, Josée Turgeon, and Christian Weiss
have contributed to various versions of this software since around 1986,
under the lead of Pierre L'Ecuyer.

## License

_Lattice Tester_ is free open source software, distributed under the Apache 2.0 License.


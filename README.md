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


## Documentation

More details on _Lattice Tester_, its underlying theory, its organization, and examples, can be found in the 
[**Lattice Tester User's Guide** (in .pdf)](https://www-labs.iro.umontreal.ca/~lecuyer/guides/lattester-guide.pdf).
This is really the core documentation. 

The interface is specified in the 
[**API documentation**](http://pierrelecuyer.github.io/latticetester/namespaces.html).

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

    cd git/gilatticetester

if you obtained the source code via `git`.
If you obtained the source code from the ZIP archive, the directory should be
named `latticetester-master` instead of `latticetester`.
At the root of the source tree, the `wscript` file contains the waf script that manages the build process.

Try:

	./waf --help

to see the various commands and options.
Some options you might want to use:
- `--out /path/to/build/location` allows you to specify in which directory the
  build process will operate. The default is `./build`. You will need permission
  to write in that directory.
- `--prefix /path/to/installation/location` allows you to specify in which 
  directory you would like to install *Lattice Tester* after it's compilation.
  The default is `/usr/local` on Linux (waf's default). You will need permission
  to write in that directory.
- `--gmp /path/to/gmp` allows you to specify the location of your gmp
  installation. You will only need this flag if waf doesn't find your gmp
  installation automatically.
- `--ntl /path/to/NTL` allows you to specify the location of your NTL 
  installation. You will only need this flag if waf does not find your NTL
  installation automatically.
- `--build-docs` waf will build the documentation if this flag is specified and 
  will not build it if it is omitted.
- `--link-static` if this flag is specified, the compiler will link all the 
  libraries statically to the executable programs (the examples). 
  This might be convenient if you installed NTL in non standard paths.

First, the project must be configured with:

	./waf configure --needed --flags

A simple 

   ./waf configure

command should be enough to configure `waf` for a minimal build. 

If you also want to build the html documentation, 
[Doxygen](http://www.stack.nl/~dimitri/doxygen/) must be available on the system
and you should first configure with

   ./waf configure --build-docs

The waf script in `doc/wscript` manages the documentation build.  
The Doxygen options are selected in the file `Doxyfile.in` and the main page is in `doc/dox/main.dox`.
The built documentation is placed in `build/doc/html/` with `index.html` as its main entry.

If NTL, and GMP are not part of the standard system installation and were
manually installed under, say, the `/opt/ntl`, and `/opt/gmp` directories,
which means that `/opt/ntl` and `/opt/gmp` all contain subdirectories named
`include` and `lib`, you should do:

    ./waf configure --ntl /opt/ntl --gmp /opt/gmp

It is possible to set the `CXX` environment variable to the path to a specific
C++ compiler to be used to build Lattice Tester, before running the `waf configure` command.

In a UNIX shell, it is also possible to run the simple `configure.sh`
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


# MONGRAIL

This is Mongrail, a program for identifying hybrids and back crosses using genomic data.

Mongrail is a free software, you can redistribute it and/or modify it under the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be redistributed in proprietary programs.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Availability

The current version of Mongrail is available from ?

## Manual

The program manual is available at the [Mongrail wiki](https://github.com/mongrail/mongrail/wiki).

## Quickstart

### To Download the program

```
wget ? ; tar -xvzf mongrail-*.tar.gz ; cd mongrail-* ; make all

```

### To run the program

```
./run_mongrail.sh -A examples/small_owl/BO_phased_small.vcf.gz -B examples/small_owl/SO_phased_small.vcf.gz -i examples/small_owl/hybrids_small.vcf.gz -r 1.5

```

## Reporting Bugs

A list of known bugs can be found in the BUGS file. Details of compilation problems can be found in the INSTALL file.

If you find a bug which is not listed in these files please report it to snechakraborty@ucdavis.edu.

All bug reports should include:

```
   The version number of Mongrail, and where you obtained it.
   The hardware and operating system.
   The compiler used, including version number and compilation options.
   A description of the bug behaviour.
   A short program which reproducibly exercises the bug.

```

Any errors or omissions in the manual can also be reported to the same address.

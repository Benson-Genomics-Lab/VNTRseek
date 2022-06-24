# Important Notes

- Up-to-date installation documentation can be found on the project's
[GitHub wiki](https://github.com/yzhernand/VNTRseek/wiki)
- Only CentOS 6 and Ubuntu 12.10 and up have been tested.
- Only a Linux 64-bit version of our software is available at this time.
- Currently, the default install path is valid only on UNIX-like platforms.

# Requirements

Latest versions preferred.

- samtools
- bedtools

## Hardware

You will need a fast computer with plenty of RAM and disk space. This pipeline
is very CPU and IO intensive, and will require plenty of system memory and
space for output.

## Software

### Build requirements
Installation requires cmake (http://www.cmake.org/), minimum version 3.0

Additionally you will need GCC version 4.1.2 on Mac/Linux/CYGWIN or
a compatible compiler (only clang has been tested).

### Runtime requirements
The following programs are required for the pipeline to run:

- SQLite (>= 3.0.0)
- Perl (>= 5.14.0)

The Perl `DBI` and `DBD::SQLite` modules are required for interacting
with a MySQL server.

TRF is also required, but is downloaded during installation. If for some
reason the download fails, you can download it manually from
[the website](http://tandem.bu.edu/trf/trf409.linux64.download.html)
and save it as `trf409-ngs.linux.exe` in the build directory (see [Installation](#installation)
below).

To install all requirements at once (assuming you want to use gcc), copy the correct
command for your system:

### Ubuntu/Debian

Note: Replace `libmysqlclient-dev` with `libmariadbclient-dev` if using MariaDB.

```
sudo apt install cmake build-essential zlib1g-dev libdbi-perl libdbd-mysql-perl libmysqlclient-dev
```

### CentOS/Red Hat/Fedora

Note: Replace with `mysql-devel` with `mariadb-devel` if using MariaDB.

```
sudo yum install cmake gcc perl-DBI perl-DBD-MySQL mysql-devel
```

### Arch Linux (and derivatives)

Note: Arch Linux only uses MariaDB, unless the AUR is enabled.

```
sudo pacman -Sy cmake gcc zlib perl-dbi perl-dbd-mysql libmariadbclient
```

The clang compiler also works. If you have this installed, see below for how
to use clang over gcc if it is not your default compiler.

#Installation

To build and install to the default directory, simply run the following commands:

    tar xzvf vntrseekN.NN.tar.gz
    cd vntrseekN.NNsrc
	mkdir build
	cd build
	cmake ..
	make install # or sudo make install, if needed

By default, this will install the pipeline to `/usr/local/vntrseekN.NN` (eg,
`/usr/local/vntrseek1.09`).

If you would like to choose a different installation prefix, simply run:

	cmake -DCMAKE_INSTALL_PREFIX=<full path> ..

For example, to install to your home directory, `${HOME}/vntrseekN.NN`, use:

	cmake -DCMAKE_INSTALL_PREFIX=${HOME} ..

You can change your default compiler using `-DCMAKE_C_COMPILER=yourcompiler`.
For example, to use clang:

    cmake -DCMAKE_C_COMPILER=clang ..

**If you installed this pipeline as root, and are creating an INDIST
file** you may need to run it as root unless you give your user
permission to write to the installation directory.

If you installed to a non-standard location, you may need to add
`/path/to/prefix/bin` to your `PATH` variable (eg, if your prefix
was `/opt`, you will need to have `/opt/bin` in your `PATH`).

# Uninstalling

On UNIX, simply run:

	xargs rm < install_manifest.txt # or sudo xargs rm < install_manifest.txt

from the build directory you created above. The directory will
remain, however, so you will not lose any reference files.

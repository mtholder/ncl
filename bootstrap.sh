#!/bin/sh
# This script should be run  whenever the inputs to autoconf/automake change
#	(namely the configure.ac file or the m4 macros in config/m4)
# It appears that the autogenerated Makefiles enough contain dependency
#		information to know that it needs to run automake and config.status
#		whenever an Makefile.am changes.
set -x
aclocal -I config || exit
#libtoolize --copy --force || exit;
autoheader || exit
automake --add-missing || exit
autoconf

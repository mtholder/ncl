NCL - the NEXUS Class Library
==============================

This is version 2.1.21-dev

This git repo was forked from https://github.com/mtholder/ncl
so that all of the open tree of life developers would have write
permissions.
Mark is adding several example programs for dealing with trees from
the Open Tree project.
They are being added to this repo just because it is easier to distribute
NCL-dependent code as a part of the library than it is to document the 
procedure of installing C++ libraries and linking to them.
If we keep using these tools, we can fairly easily pull out these examples
as free standing client repo(s).

If we find any changes need to be made to the existing files in the `ncl`
subdir, we should try to remember to pull those back into the mtholder
fork of the repo so that other users of NCL can benefit.


Most documentation for this C++ class library is in the form of HTML
files in the html directory. Please begin by viewing the html/index.html
file in your web browser.

See the INSTALL file for important information about building and installing
the NCL and example programs, and incorporating the NCL into your own
applications.

As of March 09, 2012, NCL is available under a Simplified BSD license (see
BSDLicense.txt) in addition to the GPL license.

ACKNOWLEDGEMENTS
================

Many of the files used for testing were provided by Arlin Stoltzfus (see
http://www.molevol.org/camel/projects/nexus/ for more information), the Mesquite
package, and from TreeBase (thanks, Bill Piel!).

The SWIG bindings for version 2.1 were inspired by the work of David Suarez
Pascal in the 2007 Google Summer of Code.  They were simplified using new
features of NCL 2.1.  Thanks to David for blazing the way on the old version,
Google for funding, and NESCent (in particular Hilmar Lapp) for getting the
NESCent GSoC program going.

The 2010 GSoC effort also led to enhancements in terms of annotation storage and
xml parsing which are currently on:

 https://ncl.svn.sourceforge.net/svnroot/ncl/branches/xml

Thanks to NESCent, Google, and Michael Elliot for that support.


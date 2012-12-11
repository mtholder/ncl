================================================================================
First steps
================================================================================
To make the XCode project as close to the cross-platform build, I am using 
config.h created by running ${NCL_TOP}/configure from inside the 
${NCL_TOP}/xcode/ncl-xcode  directory.  The resulting config.h is referred to
from the NCL code (or at least could be).

You do not have to run "make" before using the xcode project (in fact you can 
delete all configure products other than config.h, if you like).


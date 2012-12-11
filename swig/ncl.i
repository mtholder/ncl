// minimal wrapper around NCL
/* Setting modulename to "ncl" with directors enabled */
%module(directors="1") nclwrapper
%{
 #include "ncl.h"
 #include "nxstreestream.h"
%}

#if defined(SWIGPERL) /* Perl-specific customizations because the SWIG (language-generic) exception handling wasn't working */

%exception {
	try {
		$action
	} catch (NxsException & e) {
		croak(PyExc_ValueError, e.nxs_what());
	}
	catch (...) {
		croak("unknown exception");
	}
}
#elif defined(SWIGPYTHON) /* Python-specific customizations because the SWIG (language-generic) exception handling wasn't working */

%exception {
	try {
		$action
	} catch (NxsException & e) {
		PyErr_SetString(PyExc_ValueError, e.nxs_what());
		return NULL;
	}
	catch (...) {
		PyErr_SetString(PyExc_ValueError, "unknown exception");
		return NULL;
	}
}

#endif // defined (SWIGPYTHON)

%feature("director") NxsTreeStream;


/* Allowing synonymous employment of certain STL classes */
typedef std::string string;
typedef std::ostream ostream;
typedef std::istream istream;
typedef long file_pos;

/* Enabling default typemaps for std::strings */
%include std_string.i

/* Ignored members of NCL */
%ignore NxsBlock::CloneBlock;
%ignore	NxsCharactersBlock::GetStateSymbolIndex;


/* Conditional processing based on language of wrappers being produced */


%include "std_vector.i"

namespace std {
   %template(vectori) vector< int >;
   %template(vectorvi) vector< vector<int> >;
   %template(vectord) vector< double >;
   %template(vectorst) vector< string >;
   %template(vectoru) vector< unsigned >;

class exception {
	public:
		exception();
		exception(const exception& rhs);
		virtual ~exception();
		virtual const char *what(void);
};

};


%include nxsexception.h
%include nxstoken.h
%include nxsblock.h
%include nxsreader.h
%include nxstaxablock.h
%include nxstreesblock.h
%include nxscharactersblock.h
%include nxsassumptionsblock.h
%include nxsdatablock.h
%include nxsdistancesblock.h
%include nxssetreader.h
%include nxspublicblocks.h
%include nxsmultiformat.h


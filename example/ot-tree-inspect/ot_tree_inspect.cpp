//	Copyright (C) 2007-2008 Mark T. Holder
//
//	This file is part of NCL (Nexus Class Library).
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc.,
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

/*******************************************************************************
 * This file contains code for 4 executables:
 *		NEXUSnormalizer, NEXUSvalidator, NEXUSinspector, and NEX_us2ml
 *	with conditional compilation used to determine the behavior.
 *
 *		* NEXUSnormalizer - writes a NEXUS version of the file with consistent
 *			ordering of blocks and commands. Ideally 2 equivalent files will
 *			produce the same normalized output. This version of tthe program is
 *			less ambitious. The goal is to be able to run (for any valid NEXUS
 *			in.nex file):
 *				$ NEXUSnormalizer in.nex > outOrig.nex
 *				$ NEXUSnormalizer outOrig.nex > outSecond.nex
 *				$ diff outOrig.nex outSecond.nex
 *			and find no differences.
 *		* NEXUSvalidator - reports errors and warnings to stderr. Invalid files
 *			cause exit with a non-zero return code
 *		* NEXUSinspector - writes a brief report of every block parsed
 *		* NEXUS_us2ml - writes a nexml version of the input (partially
 *			implemented, note that the code to write nexml is in us2ml.cpp).
 * See the processFilepath() function for an example of how to deal with NCL
 *	to read a file using the new MultiFormatReader class. When the file
 *	is correctly read, the processContent() function is called.
 *
 * All other code has to do with reading command line arguments and other
 * 	user-interface concerns.
 */
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include "ncl/nxsdefs.h"
#include "ncl/nxstreesblock.h"

using namespace std;
long gStrictLevel = 2;
long gInterleaveLen = -1;
bool gVerbose = false;
static long gBurnin = 0;
void processContent(PublicNexusReader & nexusReader, ostream *out);

enum Commands {
	POLYTOMY_COUNT = 1,
	TIP_LABEL_LIST = 2,
	NODE_LABEL_LIST = 3,
	PYDICT = 4
};
enum Commands gCommand = POLYTOMY_COUNT;

bool polytomyCountHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);
bool tipLabelsListHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);
bool nodeLabelsListHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);
bool pydictHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);
void listNodeLabels(const NxsSimpleTree &nst, bool tipsOnly, const NxsTaxaBlock &taxa);


bool polytomyCountHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	static unsigned long gTreeCount = 0;
	gTreeCount++;
	if (gVerbose)
		std::cerr << "Read tree " <<  gTreeCount<< '\n';
	std::vector<unsigned> outDegreeCounts;
	NxsSimpleTree nst(ftd, 0.0, 0, true);
	nst.RerootAt(0);
	std::vector<const NxsSimpleNode *> nodes =  nst.GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin();
		 nIt != nodes.end(); ++nIt) {
		const unsigned outDegree = (*nIt)->GetChildren().size();
		if (outDegree >= outDegreeCounts.size()) {
			outDegreeCounts.resize(outDegree + 1, 0);
		}
		outDegreeCounts[outDegree] += 1;
	}
	int o = 0;
	for (std::vector<unsigned>::const_iterator cIt = outDegreeCounts.begin();
		 cIt != outDegreeCounts.end(); ++cIt, ++o) {
		if (*cIt > 0) {
			std::cout << *cIt << " nodes had out-degree = " << o << '\n';
		}
	}
	std::cout << std::endl;
	return false;
}
bool nodeLabelsListHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	const NxsTaxaBlock * taxa  = dynamic_cast<const NxsTaxaBlock *>(treesB->GetTaxaBlockPtr(NULL));
	if (!taxa) {
		std::cerr << "could not get the out label collection!\n";
		throw new exception();
	}
	static unsigned long gTreeCount = 0;
	gTreeCount++;
	if (gVerbose)
		std::cerr << "Read tree " <<  gTreeCount<< '\n';
	std::vector<unsigned> outDegreeCounts;
	NxsSimpleTree nst(ftd, 0.0, 0, true);
	nst.RerootAt(0);
	listNodeLabels(nst, false, *taxa);
	return false;
}
bool tipLabelsListHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	const NxsTaxaBlock * taxa  = dynamic_cast<const NxsTaxaBlock *>(treesB->GetTaxaBlockPtr(NULL));
	if (!taxa) {
		std::cerr << "could not get the out label collection!\n";
		throw new exception();
	}
	static unsigned long gTreeCount = 0;
	gTreeCount++;
	if (gVerbose)
		std::cerr << "Read tree " <<  gTreeCount<< '\n';
	std::vector<unsigned> outDegreeCounts;
	NxsSimpleTree nst(ftd, 0.0, 0, true);
	nst.RerootAt(0);
	listNodeLabels(nst, true, *taxa);
	return false;
}

long extractOTTIDFromName(const std::string &n) {
	//std::cerr << "n = " << n << std::endl;
	const unsigned labelLen = n.length();
	if (labelLen < 3) {
		std::cerr << "label \""<< n << "\" found that does not end in ott##### pattern\n";
		throw exception();
	}
	unsigned i = labelLen - 1;
	for (;;) {
		if (n[i] == ' ' || n[i] == '_') {
			break;
		}
		if (i == 0) {
			std::cerr << "label \""<< n << "\" found that does not end in ott##### pattern\n";
			throw exception();
		}
		--i;
	}
	if (labelLen - i < 4 || n[i + 1] != 'o' || n[i + 2] != 't' || n[i + 3] != 't') {
		std::cerr << "label \""<< n << "\" found that does not end in ott##### pattern\n";
		throw exception();
	}
	long x = 0;
	std::string num = n.substr(i + 4);
	if (!NxsString::to_long(num.c_str(), &x)) {
		std::cerr << "label \""<< n << "\" found that does not end in ott##### pattern\n";
		throw exception();
	}
	return x;
}

bool pydictHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	const NxsTaxaBlock * taxa  = dynamic_cast<const NxsTaxaBlock *>(treesB->GetTaxaBlockPtr(NULL));
	if (!taxa) {
		std::cerr << "could not get the out label collection!\n";
		throw new exception();
	}
	static unsigned long gTreeCount = 0;
	gTreeCount++;
	if (gVerbose)
		std::cerr << "Read tree " <<  gTreeCount<< '\n';
	std::vector<unsigned> outDegreeCounts;
	NxsSimpleTree nst(ftd, 0.0, 0, true);
	nst.RerootAt(0);
	long unNamedID = -1;
	std::map<const NxsSimpleNode *, long> ndToIDForPyDict;
	std::cout << "id2parent_id = {\n";
	std::vector<const NxsSimpleNode *> nodes =  nst.GetPreorderTraversal();
	const std::string fn = "IDTabNameNewlineFileForPyDict.tsv";
	std::cerr << "The file \"" << fn << "\" will be overwritten with a tab-delimited content of two columns: the ID<TAB>Name\nfor each named node" << std::endl;
	std::ofstream namesFile;
	namesFile.open(fn.c_str(), std::ofstream::out|std::ios::binary);
	try {
			for (std::vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
			const NxsSimpleNode * nd = *nIt;
			long ottID = -1;
			if (nd->IsTip()) {
				const unsigned ind = nd->GetTaxonIndex();
				const std::string & n = taxa->GetTaxonLabel(ind);
				ottID = extractOTTIDFromName(n);
				namesFile << ottID << '\t' << n << '\n';
			} else {
				const std::string & n = nd->GetName();
				if (! n.empty()) {
					ottID = extractOTTIDFromName(n);
					namesFile << ottID << '\t' << n << '\n';
				} else {
					ottID = --unNamedID;
				}
			}
			ndToIDForPyDict[nd] = ottID;
			const NxsSimpleNode * par = 0L;
			const NxsSimpleEdge & edge = nd->GetEdgeToParentRef();
			par = nd->GetParent();
			if (par == 0) {
				std::cout << ottID <<": None,\n";
			} else {
				std::map<const NxsSimpleNode *, long>::const_iterator pIt = ndToIDForPyDict.find(par);
				if (pIt == ndToIDForPyDict.end()) {
					std::cerr << "parent node not found!\n";
					throw exception();
				}
				const long & pId = pIt->second;
				std::cout << ottID <<": " << pId <<" ,\n";
			}
		}
		std::cout << "}\n";
	} catch (...) {
		namesFile.close();
		throw;
	}
	namesFile.close();
	return false;
}

void listNodeLabels(const NxsSimpleTree &nst, bool tipsOnly, const NxsTaxaBlock &taxa)
{
	std::vector<const NxsSimpleNode *> nodes =  nst.GetPreorderTraversal();
	if (tipsOnly) {
	for (std::vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin();
			 nIt != nodes.end(); ++nIt) {
			if ((*nIt)->IsTip()) {
				const unsigned ind = (*nIt)->GetTaxonIndex();
				const std::string & n = taxa.GetTaxonLabel(ind);
				std::cout << n << '\n';
			}
		}
	} else {
		for (std::vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin();
			 nIt != nodes.end(); ++nIt) {
			if ((*nIt)->IsTip()) {
				const unsigned ind = (*nIt)->GetTaxonIndex();
				const std::string & n = taxa.GetTaxonLabel(ind);
				std::cout << n << std::endl;
			} else {
				std::cout << (*nIt)->GetName() << '\n';
			}
		}
	}
	std::cout << std::endl;
}
////////////////////////////////////////////////////////////////////////////////
// Takes NxsReader that has successfully read a file, and processes the
//	information stored in the reader.
//
// The caller is responsibel for calling DeleteBlocksFromFactories() to clean
//	up (if the reader uses the factory API).
////////////////////////////////////////////////////////////////////////////////
void processContent(PublicNexusReader & nexusReader, ostream *out)
{
	if (!out)
		return;

}

////////////////////////////////////////////////////////////////////////////////
// Creates a NxsReader, and tries to read the file `filename`.  If the
//	read succeeds, then processContent will be called.
////////////////////////////////////////////////////////////////////////////////
void processFilepath(
	const char * filename, // name of the file to be read
	ostream *out, // output stream to use (NULL for no output). Not that cerr is used to report errors.
	MultiFormatReader::DataFormatType fmt) // enum indicating the file format to expect.
	{
	assert(filename);
	try
		{
		MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
		nexusReader.SetWarningOutputLevel(NxsReader::AMBIGUOUS_CONTENT_WARNING);
		if (gStrictLevel != 2)
			nexusReader.SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
		NxsCharactersBlock * charsB = nexusReader.GetCharactersBlockTemplate();
		NxsDataBlock * dataB = nexusReader.GetDataBlockTemplate();
		charsB->SetAllowAugmentingOfSequenceSymbols(true);
		dataB->SetAllowAugmentingOfSequenceSymbols(true);
		if (gInterleaveLen > 0)
			{
			assert(charsB);
			charsB->SetWriteInterleaveLen(gInterleaveLen);
			dataB->SetWriteInterleaveLen(gInterleaveLen);
			}
		NxsTreesBlock * treesB = nexusReader.GetTreesBlockTemplate();
		assert(treesB);
		if (gStrictLevel < 2)
			treesB->SetAllowImplicitNames(true);
		if (gCommand == POLYTOMY_COUNT)
			treesB->setValidationCallbacks(polytomyCountHook, 0L);
		else if (gCommand == TIP_LABEL_LIST)
			treesB->setValidationCallbacks(tipLabelsListHook, 0L);
		else if (gCommand == NODE_LABEL_LIST)
			treesB->setValidationCallbacks(nodeLabelsListHook, 0L);
		else if (gCommand == PYDICT)
			treesB->setValidationCallbacks(pydictHook, 0L);
		if (gStrictLevel < 2)
			{
			NxsStoreTokensBlockReader *storerB =  nexusReader.GetUnknownBlockTemplate();
			assert(storerB);
			storerB->SetTolerateEOFInBlock(true);
			}
		try {
			nexusReader.ReadFilepath(filename, fmt);
			processContent(nexusReader, out);
			}
		catch(...)
			{
			nexusReader.DeleteBlocksFromFactories();
			throw;
			}
		nexusReader.DeleteBlocksFromFactories();
		}
	catch (const NxsException &x)
		{
		cerr << "Error:\n " << x.msg << endl;
		if (x.line >=0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		exit(2);
		}
	}

void readFilepathAsNEXUS(const char *filename, MultiFormatReader::DataFormatType fmt)
	{
	try {
		ostream * outStream = 0L;
		outStream = &cout;
		processFilepath(filename, outStream, fmt);
		}
	catch (...)
		{
		cerr << "Parsing of " << filename << " failed (with an exception)" << endl;
		exit(1);
		}
	}

void readFilesListedIsFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt)
	{
	ifstream masterStream(masterFilepath, ios::binary);
	if (masterStream.bad())
		{
		cerr << "Could not open " << masterFilepath << "." << endl;
		exit(3);
		}
	char filename[1024];
	while ((!masterStream.eof())  && masterStream.good())
		{
		masterStream.getline(filename, 1024);
		if (strlen(filename) > 0 && filename[0] != '#')
			readFilepathAsNEXUS(filename, fmt);
		}
	}

void printHelp(ostream & out)
	{
	out << "otTreeInspect.\n";
	out << "\nThe most common usage is simply:\n    NEXUStosplits <path to NEXUS file>\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -h on the command line shows this help message\n\n";
	out << "    -v verbose output\n\n";
	out << "    -c<command> to specify an action. <command> should be:\n";
	out << "        poly  for a polytomy count\n";
	out << "        tips  for tip labels (one per line)\n";
	out << "        labels  for all labels (one per line)\n";
	out << "        pydict  for a python dict of ID-> parent ID using OTT IDs when present and arbitrary (unstable) IDs for unnamed nodes.\n";
	out << "    -s<non-negative integer> controls the NEXUS strictness level.\n";
	out << "        the default level is equivalent to -s2 invoking the program with \n";
	out << "        -s3 or a higher number will convert some warnings into fatal errors.\n";
	out << "        Running with -s1 will cause the parser to accept dangerous constructs,\n";
	out << "        and running with -s0 will cause the parser make every attempt to finish\n";
	out << "        parsing the file (warning about very serious errors).\n\n";
	out << "        Note that when -s0 strictness level is used, and the parser fails to\n";
	out << "        finish, it will often be the result of an earlier error than the \n";
	out << "        error that is reported in the last message.\n";
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		//passs
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		//passs
#	else
		out << "    -i<format> specifies the input file format expected:\n";
#	endif
	out << "    -f<format> specifies the input file format expected:\n";
	out << "            -frelaxedphyliptree     newick tree (this is also the default)\n";
	out << "            -faafasta   Amino acid data in fasta\n";
	out << "            -fdnafasta  DNA data in fasta\n";
	out << "            -frnafasta  RNA data in fasta\n";
	out << "        The complete list of format names that can follow the -f flag is:\n";
	std::vector<std::string> fmtNames =  MultiFormatReader::getFormatNames();
	for (std::vector<std::string>::const_iterator n = fmtNames.begin(); n != fmtNames.end(); ++n)
		{
		out << "            "<< *n << "\n";
		}
	}

int main(int argc, char *argv[])
	{
	MultiFormatReader::DataFormatType f(MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT);

	bool readfile = false;
	bool el = false;
	bool depth = false;
	bool brief = false;
	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h')
			{
			printHelp(cout);
			return 0;
			}
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'v')
			gVerbose = true;
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 's')
			{
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gStrictLevel)))
				{
				cerr << "Expecting an integer after -s\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'f')
			{
			f = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2)
				{
				std::string fmtName(filepath + 2, slen - 2);
				f =  MultiFormatReader::formatNameToCode(fmtName);
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting a format after after -f\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'c')
			{
			std::string commandName;
			if (slen > 2)
				{
				commandName = string(filepath + 2, slen - 2);

				}
			if (commandName == "poly") {
				gCommand = POLYTOMY_COUNT;
			} else if (commandName == "tips") {
				gCommand = TIP_LABEL_LIST;
			} else if (commandName == "labels") {
				gCommand = NODE_LABEL_LIST;
			} else if (commandName == "pydict") {
				gCommand = PYDICT;
			}
			else
				{
				cerr << "Expecting a command after after -c\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		else
			{
			readfile = true;
			readFilepathAsNEXUS(filepath, f);
			}
		}
	if (!readfile)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
		}
	return 0;
	}


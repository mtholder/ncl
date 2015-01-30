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

#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"
#include <cassert>
#include "ncl/nxsdefs.h"
#include "ncl/nxstreesblock.h"

using namespace std;
long gStrictLevel = 2;
bool gVerbose = false;
void processContent(PublicNexusReader & nexusReader, ostream *out);


bool newTreeHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);

void describeUnnamedNode(const NxsTaxaBlockAPI* taxa, const NxsSimpleNode &, ostream & out);

void fillAncMRCA(const NxsSimpleNode * nd, map<const NxsSimpleNode *, set<long> > &n2m) {
	const vector<NxsSimpleNode *> children = nd->GetChildren();
	set<long> & mrca = n2m[nd];
	for (vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
		const NxsSimpleNode * c = *cIt;
		if (n2m.find(c) != n2m.end()) {
			set<long> & csl = n2m[c];
			mrca.insert(csl.begin(), csl.end());
		}
	}
}

/* use some globals, because I'm being lazy... */
string gCurrentFilename;
string gCurrTmpFilepath;

inline long ottIDFromName(const string & n) {
	//cout << "name \"" << n << "\"\n";
	if (n.empty()) {
		return -1;
	}
	const unsigned lastInd = n.length() - 1;
	unsigned currInd = lastInd;
	const char * c = n.c_str();
	if (strchr("0123456789", c[currInd]) == 0) {
		return -2;
	}
	while (currInd > 1) {
		--currInd;
		if (strchr("0123456789", c[currInd]) == 0) {
			++currInd;
			break;
		}
	}
	long conv = -2;
	NxsString::to_long(c + currInd, &conv);
	return conv;
}

inline long getOTTIndex(const NxsTaxaBlockAPI * taxa, const NxsSimpleNode & nd) {
	const string & name = nd.GetName();
	if (name.empty()) {
		const unsigned ind = nd.GetTaxonIndex();
		if (ind < taxa->GetNumTaxonLabels()) {
			const string tn = taxa->GetTaxonLabel(ind);
			return ottIDFromName(tn);
		}
		return -1;
	}
	return ottIDFromName(name);
}

set<long> gMRCADesignatorSet;
int gExitCode = 0;



void writeNewickSubtree(ostream & out, const NxsSimpleNode * sr, map<const NxsSimpleNode *, string > & leafNode2name) {
	assert(sr != 0);
	if (sr->IsTip()) {
		out << NxsString::GetEscaped(leafNode2name[sr]);
	} else {
		out << '(';
		bool first = true;
		const std::vector<NxsSimpleNode *> children = sr->GetChildren();
		for (std::vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
			const NxsSimpleNode * child = *cIt;
			if (first) {
				first = false;
			} else {
				out << ',';
			}
			writeNewickSubtree(out, child, leafNode2name);
		}
		assert(!first);
		out << ')';
	}

}

bool processRefTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	map<const NxsSimpleNode *, set<long> > refNdp2mrca;
	map<const NxsSimpleNode *, string > leafNode2name;
	const unsigned numMRCADesignators = gMRCADesignatorSet.size();
	assert(numMRCADesignators > 1);
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		if (nd->IsTip()) {
			long ottID = getOTTIndex(tb, *nd);
			assert(ottID >= 0);
			const unsigned ind = nd->GetTaxonIndex();
			assert(ind < tb->GetNumTaxonLabels());
			const string tn = tb->GetTaxonLabel(ind);
			leafNode2name[nd] = tn;
			if (gMRCADesignatorSet.find(ottID) != gMRCADesignatorSet.end()) {
				gMRCADesignatorSet.erase(ottID);
				refNdp2mrca[nd].insert(ottID);
			}
		} else {
			fillAncMRCA(nd, refNdp2mrca);
			if (refNdp2mrca[nd].size() == numMRCADesignators) {
				writeNewickSubtree(cout, nd, leafNode2name);
				cout << ";\n";
				return true;
			}
		}
	}
	if (gMRCADesignatorSet.empty()) {
		std::cerr << "Very odd: no node found that is an ancestor of all MRCA designators, but all designators found.\n";
	} else {
		std::cerr << "There following MRCA designator(s) not found (they all have to be leaf nodes):\n";
		for (set<long>::const_iterator mIt = gMRCADesignatorSet.begin(); mIt != gMRCADesignatorSet.end(); ++mIt) {
			std::cerr << *mIt << "\n";
		}
	}
	return false;
}

bool newTreeHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	static unsigned long gTreeCount = 0;
	const NxsTaxaBlockAPI * taxa = treesB->GetTaxaBlockPtr();
	gTreeCount++;
	if (gVerbose)
		cerr << "Read tree " <<  gTreeCount<< '\n';
	unsigned int nUnlabeledOutDegOne = 0;
	unsigned int nLabeledOutDegOne = 0;
	vector<string> parNames;
	NxsSimpleTree nst = NxsSimpleTree(ftd, 0.0, 0, true);
	if (!processRefTree(taxa, &nst)) {
		gExitCode = 1;
	}
	return false;
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
		if (gStrictLevel != 2)
			nexusReader.SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
		NxsCharactersBlock * charsB = nexusReader.GetCharactersBlockTemplate();
		NxsDataBlock * dataB = nexusReader.GetDataBlockTemplate();
		charsB->SetAllowAugmentingOfSequenceSymbols(true);
		dataB->SetAllowAugmentingOfSequenceSymbols(true);
		NxsTreesBlock * treesB = nexusReader.GetTreesBlockTemplate();
		assert(treesB);
		if (gStrictLevel < 2)
			treesB->SetAllowImplicitNames(true);
		treesB->setValidationCallbacks(newTreeHook, 0L);
		if (gStrictLevel < 2)
			{
			NxsStoreTokensBlockReader *storerB =  nexusReader.GetUnknownBlockTemplate();
			assert(storerB);
			storerB->SetTolerateEOFInBlock(true);
			}
		cerr << "Executing" <<endl;
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
	cerr << "[Reading " << filename << "	 ]" << endl;
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
	out << "NEXUStosplits takes reads a file and writes trees for each split that occurs in any tree in the file.\n";
	out << "\nThe most common usage is simply:\n    NEXUStosplits <path to NEXUS file>\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -h on the command line shows this help message\n\n";
	out << "    -v verbose output\n\n";
	out << "    -f<format> specifies the input file format expected:\n";
	out << "            -fnexus     NEXUS (this is also the default)\n";
	out << "            -frelaxedphyliptree  newick(this is also the default)\n";
	out << "        The complete list of format names that can follow the -f flag is:\n";
	vector<string> fmtNames =  MultiFormatReader::getFormatNames();
	for (vector<string>::const_iterator n = fmtNames.begin(); n != fmtNames.end(); ++n) {
		out << "            "<< *n << "\n";
	}
}

int main(int argc, char *argv[]) {
	MultiFormatReader::DataFormatType f(MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT);
	bool readfile = false;
	bool el = false;
	bool depth = false;
	bool brief = false;
	for (int i = 1; i < argc; ++i) {
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h')
			printHelp(cout);
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'v')
			gVerbose = true;
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 's') {
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gStrictLevel))) {
				cerr << "Expecting an integer after -s\n" << endl;
				printHelp(cerr);
				return 2;
			}
		} else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'f') {
			f = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2) {
				string fmtName(filepath + 2, slen - 2);
				f =  MultiFormatReader::formatNameToCode(fmtName);
			}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT) {
				cerr << "Expecting a format after after -f\n" << endl;
				printHelp(cerr);
				return 2;
			}
		} else {
			readfile = true;
			const string filepathstr(filepath);
			const size_t sp = filepathstr.find_last_of('/');
			if (sp == string::npos) {
				gCurrentFilename = filepathstr;
			} else {
				gCurrentFilename = filepathstr.substr(sp + 1);
			}
			try {
				for (int j = i + 1; j < argc; ++j) {
					long mott;
					if (!NxsString::to_long(argv[j], &mott) || mott < 1) {
						cerr << "Expecting positive integer for an OTT ID as a MRCA designators.\n";
						return 3;
					}
					gMRCADesignatorSet.insert(mott);
				}
				if (gMRCADesignatorSet.size() < 2) {
					cerr << "must specify at least 2 ottIDs as MRCA designators.\n";
					return 2;
				}
				readFilepathAsNEXUS(filepath, f);
				return gExitCode;
			}
			catch (...) {
				throw;
			}
		}
	}
	if (!readfile) {
		cerr << "Expecting the path to NEXUS file and a series of numeric, leaf OTT IDs as the only command line arguments!\n" << endl;
		printHelp(cerr);
		return 1;
	}
	return gExitCode;
}

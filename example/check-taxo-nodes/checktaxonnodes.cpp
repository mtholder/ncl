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

void describeUnnamedNode(const NxsTaxaBlockAPI* taxa, const NxsSimpleNode &, std::ostream & out);

void fillAncMRCA(const NxsSimpleNode * nd, std::map<const NxsSimpleNode *, std::set<long> > &n2m) {
	const std::vector<NxsSimpleNode *> children = nd->GetChildren();
	std::set<long> & mrca = n2m[nd];
	for (std::vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
		const NxsSimpleNode * c = *cIt;
		assert(n2m.find(c) != n2m.end());
		std::set<long> & csl = n2m[c];
		mrca.insert(csl.begin(), csl.end());
	}
}

/* use some globals, because I'm being lazy... */
NxsSimpleTree * gRefTree = 0;
NxsSimpleTree * gTaxonTree = 0;
std::map<long, const NxsSimpleNode *> gOttID2RefNode;
std::map<const NxsSimpleNode *, std::string> gRefTipToName;
std::map<long, const NxsSimpleNode *> gOttID2TaxNode;
std::map<const NxsSimpleNode *, long> gTaxNode2ottID;
std::set<const NxsSimpleNode *> gSupportedNodes;
std::string gCurrentFilename;
std::string gCurrTmpFilepath;
//std::ostream * gCurrTmpOstream = 0L;


std::string getLeftmostDesName(const NxsSimpleNode *nd) {
	const std::string & name = nd->GetName();
	if (!name.empty()) {
		return name;
	}
	const unsigned outDegree = nd->GetChildren().size();
	if (outDegree == 0) {
		return gRefTipToName[nd];
	}
	return getLeftmostDesName(nd->GetChildren()[0]);
}

std::string getRightmostDesName(const NxsSimpleNode *nd) {
	const std::string & name = nd->GetName();
	if (!name.empty()) {
		return name;
	}
	const unsigned outDegree = nd->GetChildren().size();
	if (outDegree == 0) {
		return gRefTipToName[nd];
	}
	const unsigned lastInd = outDegree - 1;
	return getRightmostDesName(nd->GetChildren()[lastInd]);
}

void describeUnnamedNode(const NxsSimpleNode *nd, std::ostream & out, unsigned int anc) {
	if (nd->GetName().length() > 0) {
		out << "ancestor " << anc << " node(s) before \"" << nd->GetName() << "\"" << std::endl;
	}
	std::vector<NxsSimpleNode *> children = nd->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree == 1U) {
		describeUnnamedNode(children[0], out, anc + 1);
	} else {
		std::string left = getLeftmostDesName(children[0]);
		std::string right = getRightmostDesName(children[outDegree - 1]);
		out << "ancestor " << anc << " node(s) before MRCA of \"" << left << "\" and ";
		out << "\"" << right <<'\"' << std::endl;
	}
}


void extendSupportedToRedundantNodes(const NxsSimpleTree * tree, std::set<const NxsSimpleNode *> & gSupportedNodes) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		std::vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree == 1 && gSupportedNodes.find(children[0]) != gSupportedNodes.end()) {
			gSupportedNodes.insert(nd);
		}
	}
}

bool singleDesSupportedOrNamed(const NxsSimpleNode *nd) {
	if (gSupportedNodes.find(nd) != gSupportedNodes.end()) {
		return true;
	}
	std::vector<NxsSimpleNode *> children = nd->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree == 1) {
		if (!nd->GetName().empty()) {
			return true;
		} else {
			return singleDesSupportedOrNamed(children[0]);
		}
	}
	return false;
}

void describeUnnamedUnsupported(std::ostream &out, const NxsSimpleTree * tree, const std::set<const NxsSimpleNode *> & gSupportedNodes) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	std::vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin();
	++nIt; //skip the root
	for (;nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		std::vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree > 0 && gSupportedNodes.find(nd) == gSupportedNodes.end()) {
			if (outDegree == 1) {
				if (singleDesSupportedOrNamed(nd)) {
					continue;
				}
			}
			if (nd->GetName().length() == 0) { //assume that it is from the taxonomy
				out << "Unsupported node ";
				describeUnnamedNode(nd, out, 0);
			}
		}
	}
}


inline long ottIDFromName(const std::string & n) {
	//std::cout << "name \"" << n << "\"\n";
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
	const std::string & name = nd.GetName();
	if (name.empty()) {
		const unsigned ind = nd.GetTaxonIndex();
		if (ind < taxa->GetNumTaxonLabels()) {
			const std::string tn = taxa->GetTaxonLabel(ind);
			return ottIDFromName(tn);
		}
		return -1;
	}
	return ottIDFromName(name);
}


std::map<const NxsSimpleNode *, std::set<long> > gRefNdp2mrca;
std::map<const NxsSimpleNode *, std::set<long> > gTaxNdp2mrca;
set<long> gRefLeafSet;
set<long> gTaxLeafSet;
map<const NxsSimpleNode *, long> gRefNamedNodes;


void processRefTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, *nd);
		if (nd->IsTip()) {
			const unsigned ind = nd->GetTaxonIndex();
			assert(ind < tb->GetNumTaxonLabels());
			const std::string tn = tb->GetTaxonLabel(ind);
			gRefTipToName[nd] = tn;
			assert(ottID >= 0);
			gRefNdp2mrca[nd].insert(ottID);
			gRefLeafSet.insert(ottID);
		} else {
			fillAncMRCA(nd, gRefNdp2mrca);
			if (ottID > 0) {
				gRefNamedNodes[nd] = ottID;
			}
		}
		if (ottID >= 0) {
			assert(gOttID2RefNode.find(ottID) == gOttID2RefNode.end());
			gOttID2RefNode[ottID] = nd;
		}
	}
}

void processTaxonomyTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, **nIt);
		assert(ottID >= 0);
		assert(gOttID2TaxNode.find(ottID) == gOttID2TaxNode.end());
		if (nd->IsTip()) {
			assert(gOttID2RefNode.find(ottID) != gOttID2RefNode.end());
			gTaxLeafSet.insert(ottID);
			gTaxNdp2mrca[nd].insert(ottID);
		} else {
			fillAncMRCA(nd, gTaxNdp2mrca);
		}
		gOttID2TaxNode[ottID] = *nIt;
		gTaxNode2ottID[*nIt] = ottID;
	}
	for (std::map<long, const NxsSimpleNode *>::const_iterator nit = gOttID2RefNode.begin(); nit != gOttID2RefNode.end(); ++nit) {
		assert(gOttID2TaxNode.find(nit->first) != gOttID2TaxNode.end());
	}
}


void writeSetDiff(std::ostream & out, const char *indent, const set<long> &fir, const char *firN, const set<long> & sec, const char *secN) {
		for (set<long>::const_iterator rIt = fir.begin(); rIt != fir.end(); ++rIt) {
			if (sec.find(*rIt) == sec.end()) {
				out << indent << "ott" << *rIt << " is in " << firN << " but not " << secN << "\n";
			}
		}
		for (set<long>::const_iterator rIt = sec.begin(); rIt != sec.end(); ++rIt) {
			if (fir.find(*rIt) == fir.end()) {
				out << indent << "ott" << *rIt << " is in " << secN << " but not " << firN << "\n";
			}
		}

}

bool isProperSubset(const set<long> & small, const set<long> & big) {
	if (big.size() <= small.size()) {
		return false;
	}
	for (set<long>::const_iterator rIt = small.begin(); rIt != small.end(); ++rIt) {
		if (big.find(*rIt) == big.end()) {
			return false;
		}
	}
	return true;
}

bool doCheckEquivalent(std::ostream &out, long ottID, const NxsSimpleNode * snode, std::map<const NxsSimpleNode *, std::set<long> > & srcLookup,
										  const NxsSimpleNode * tnode, std::map<const NxsSimpleNode *, std::set<long> > & taxLookup,
										  bool topLevel, bool climbSynth, bool climbTax) {
	std::map<const NxsSimpleNode *, std::set<long> >::const_iterator streeLSIt = srcLookup.find(snode);
	assert(streeLSIt != srcLookup.end());
	std::map<const NxsSimpleNode *, std::set<long> >::const_iterator taxtreeLSIt = taxLookup.find(tnode);
	assert(taxtreeLSIt != taxLookup.end());
	const std::set<long> & streeMRCA = streeLSIt->second;
	const std::set<long> & taxtreeMRCA = taxtreeLSIt->second;
	if (streeMRCA != taxtreeMRCA) {
		if (topLevel) {
			out << "ottID " << ottID << " incorrect:\n";
			writeSetDiff(out, "    ", streeMRCA, "synth", taxtreeMRCA, "taxonomy");
		}
		if (climbSynth && isProperSubset(streeMRCA, taxtreeMRCA)) {
			return doCheckEquivalent(out, ottID, snode->GetEdgeToParent().GetParent(), srcLookup, tnode, taxLookup, false, true, false);
		} else if (climbTax && isProperSubset(taxtreeMRCA, streeMRCA)) {
			return doCheckEquivalent(out, ottID, snode, srcLookup, tnode->GetEdgeToParent().GetParent(), taxLookup, false, false, true);
		} else {
			return false;
		}
	} else if (!topLevel) {
		out << "        Found identical leaf sets for the synthetic tree \"" << snode->GetName() << "\" and the taxonomic node \"" << tnode->GetName() << "\".\n";
	}
	return true;
}

void summarize(std::ostream & out) {
	for (map<const NxsSimpleNode *, long>::const_iterator rnit = gRefNamedNodes.begin(); rnit != gRefNamedNodes.end(); ++rnit) {
		const NxsSimpleNode * nd = rnit->first;
		const long ottID = rnit->second;
		std::map<long, const NxsSimpleNode *>::const_iterator tID2nd = gOttID2TaxNode.find(ottID);
		assert(tID2nd != gOttID2TaxNode.end());
		const NxsSimpleNode *taxNd = tID2nd->second;
		if (!doCheckEquivalent(out, ottID, nd, gRefNdp2mrca, taxNd, gTaxNdp2mrca, true, true, true)) {
			out << "        Could not find this set of leaves in the synth \"" << nd->GetName() <<"\" in any taxonomic node.\n";
		}
	}
	if (gTaxLeafSet != gRefLeafSet) {
		writeSetDiff(out, "", gRefLeafSet, "synth", gTaxLeafSet, "taxonomy");
	}
}

bool newTreeHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	static unsigned long gTreeCount = 0;
	const NxsTaxaBlockAPI * taxa = treesB->GetTaxaBlockPtr();
	gTreeCount++;
	if (gVerbose)
		std::cerr << "Read tree " <<  gTreeCount<< '\n';
	unsigned int nUnlabeledOutDegOne = 0;
	unsigned int nLabeledOutDegOne = 0;
	std::vector<std::string> parNames;
	NxsSimpleTree * nst = new NxsSimpleTree(ftd, 0.0, 0, true);
	//gExpanded.clear();
	//gTabooLeaf.clear();
	if (gRefTree == 0) {
		gRefTree = nst;
		processRefTree(taxa, nst);
	} else if (gTaxonTree == 0) {
		gTaxonTree = nst;
		processTaxonomyTree(taxa, nst);
	} else {
		std::cerr << "Exepting only 2 files: the synthetic tree, and then the taxonomy\n";
		std::exit(1);
	}
	if (gRefTree != nst && gTaxonTree != nst) {
		delete nst;
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
	out << "            -fnexus     NEXUS (this is also the default)\n";
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
	MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);

	bool readfile = false;
	bool el = false;
	bool depth = false;
	bool brief = false;
	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h')
			printHelp(cout);
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
		else
			{
			readfile = true;
			const std::string filepathstr(filepath);
			const size_t sp = filepathstr.find_last_of('/');
			if (sp == std::string::npos)
				{
				gCurrentFilename = filepathstr;
				}
			else
				{
				gCurrentFilename = filepathstr.substr(sp + 1);
				}
			gCurrTmpFilepath = std::string("tmp/") + gCurrentFilename;
			//std::cout << "gCurrTmpFilepath = " << gCurrTmpFilepath << '\n';
			//std::ofstream tostream(gCurrTmpFilepath);
			//gCurrTmpOstream = &tostream;
			try {
				readFilepathAsNEXUS(filepath, f);
				//tostream.close();
				}
			catch (...)
				{
				//tostream.close();
				throw;
				}
			//gCurrTmpOstream = 0L;
			}
		}
	if (!readfile)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
		}
	summarize(std::cout);
	return 0;
	}


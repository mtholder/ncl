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
long gInterleaveLen = -1;
bool gVerbose = false;
static long gBurnin = 0;
void processContent(PublicNexusReader & nexusReader, ostream *out);




bool newTreeHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);

void describeUnnamedNode(const NxsTaxaBlockAPI* taxa, const NxsSimpleNode &, std::ostream & out);

std::string getLeftmostDesName(const NxsTaxaBlockAPI* taxa, const NxsSimpleNode &nd) {
	const unsigned ind = nd.GetTaxonIndex();
	if (ind < taxa->GetNumTaxonLabels()) {
		return taxa->GetTaxonLabel(ind);
	}
	const std::string & name = nd.GetName();
	if (!name.empty()) {
		return name;
	}
	const unsigned outDegree = nd.GetChildren().size();
	if (outDegree == 0) {
		assert(false);
		return name;
	}
	return getLeftmostDesName(taxa, *nd.GetChildren()[0]);
}

std::string getRightmostDesName(const NxsTaxaBlockAPI* taxa, const NxsSimpleNode &nd) {
	const unsigned ind = nd.GetTaxonIndex();
	if (ind < taxa->GetNumTaxonLabels()) {
		return taxa->GetTaxonLabel(ind);
	}
	const std::string & name = nd.GetName();
	if (!name.empty()) {
		return name;
	}
	const unsigned outDegree = nd.GetChildren().size();
	if (outDegree == 0) {
		assert(false);
		return name;
	}
	const unsigned lastInd = outDegree - 1;
	return getRightmostDesName(taxa, *nd.GetChildren()[lastInd]);
}
void describeUnnamedNode(const NxsTaxaBlockAPI *taxa, const NxsSimpleNode &nd, std::ostream & out, unsigned int anc) {
	const NxsSimpleNode *c = nd.GetChildren()[0];
	assert(c);
	const std::string & cname = c->GetName();
	if (!cname.empty()) {
		out << "ancestor " << 1 + anc << " node(s) before \"" << cname << "\"\n";
		out.flush();
	} else {
		const unsigned coutDegree = c->GetChildren().size();
		if (coutDegree == 1U) {
			describeUnnamedNode(taxa, *c, out, anc + 1);
		} else {
			std::string left = getLeftmostDesName(taxa, *c);
			std::string right = getRightmostDesName(taxa, *c);
			out << "ancestor " << 1 + anc << " node(s) before MRCA of \"" << left << "\" and";
			out << "\"" << right <<'\"' << std::endl;
		}
	}
}

NxsSimpleTree * refTree = 0;

NxsSimpleTree * taxonomyTree = 0;
std::map<long, const NxsSimpleNode *> ottID2RefNode;
std::map<long, const NxsSimpleNode *> ottID2TaxNode;
std::map<const NxsSimpleNode *, long> taxNode2ottID;


std::set<const NxsSimpleNode *> gSupportedNodes;
void extendSupportedToRedundantNodes(const NxsSimpleTree * tree, std::set<const NxsSimpleNode *> & gSupportedNodes) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		std::vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree == 1 && gSupportedNodes.find(nd) != gSupportedNodes.end()) {
			gSupportedNodes.insert(children[0]);
		}
	}
}


void summarize(std::ostream & out) {
	extendSupportedToRedundantNodes(refTree, gSupportedNodes);
	out << gSupportedNodes.size() << " supported nodes.\n";
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

void processRefTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
		long ottID = getOTTIndex(tb, **nIt);
		if (ottID >= 0) {
			assert(ottID2RefNode.find(ottID) == ottID2RefNode.end());
			ottID2RefNode[ottID] = *nIt;
		}
	}
}

void processTaxonomyTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
		long ottID = getOTTIndex(tb, **nIt);
		assert(ottID >= 0);
		assert(ottID2TaxNode.find(ottID) == ottID2TaxNode.end());
		ottID2TaxNode[ottID] = *nIt;
		taxNode2ottID[*nIt] = ottID;
	}
	for (std::map<long, const NxsSimpleNode *>::const_iterator nit = ottID2RefNode.begin(); nit != ottID2RefNode.end(); ++nit) {
		assert(ottID2TaxNode.find(nit->first) != ottID2TaxNode.end());
	}
}

//@recursive!
void fillTipOTTIDs(const std::map<long, const NxsSimpleNode *> &taxonomy, long ottID, std::set<long> & tipOTTIDs) {
	std::map<long, const NxsSimpleNode *>::const_iterator tnIt = taxonomy.find(ottID);
	assert(tnIt != taxonomy.end());
	const NxsSimpleNode * tn = tnIt->second;
	std::vector<NxsSimpleNode *> children = tn->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree == 0) {
		tipOTTIDs.insert(taxNode2ottID[tn]);
	} else {
		for (std::vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
			long ct = taxNode2ottID[*cIt];
			fillTipOTTIDs(taxonomy, ct, tipOTTIDs);
		}
	}
}

const NxsSimpleNode * findMRCA(const std::map<long, const NxsSimpleNode *> & ref, const std::map<long, const NxsSimpleNode *> &taxonomy, long ottID) {
	std::set<long> tipOTTIDs;
	fillTipOTTIDs(taxonomy, ottID, tipOTTIDs);
	const unsigned nTips = tipOTTIDs.size();
	if (nTips < 2) {
		std::cerr << "findMRCA called on " << ottID << '\n';
		assert(false);
	}
	std::map<const NxsSimpleNode *, unsigned int> n2c;
	long shortestPathLen = -1;
	const NxsSimpleNode * shortestPathNode = 0;
	for (std::set<long>::const_iterator toIt = tipOTTIDs.begin(); toIt != tipOTTIDs.end(); ++toIt) {
		const NxsSimpleNode * nd = 0;
		std::map<long, const NxsSimpleNode *>::const_iterator rIt = ref.find(*toIt);
		if (rIt == ref.end()) {
			std::cerr << "tip " << *toIt << " a descendant of " << ottID << " not found.\n";
			assert(false);
		}
		nd = rIt->second;
		long currPathLen = 0;
		while (nd != 0) {
			n2c[nd] += 1;
			currPathLen += 1;
			nd = nd->GetEdgeToParentRef().GetParent();
		}
		if (shortestPathLen < 0 || currPathLen < shortestPathLen) {
			shortestPathLen = currPathLen;
			shortestPathNode = rIt->second;
		}
	}
	const NxsSimpleNode * cn = shortestPathNode;
	while (cn != 0) {
		if (n2c[cn] == nTips) {
			return cn;
		}
		cn = cn->GetEdgeToParentRef().GetParent();
	}
	assert(false);
	return 0L;
}

void markPathToRoot(std::map<const NxsSimpleNode *, std::set<long> > &n2m, long ottID) {
	const NxsSimpleNode * nd = 0;
	std::map<long, const NxsSimpleNode *>::const_iterator fIt = ottID2RefNode.find(ottID);
	if (fIt != ottID2RefNode.end()) {
		nd = fIt->second;
	} else {
		nd = findMRCA(ottID2RefNode, ottID2TaxNode, ottID);
	}
	assert(nd != 0);
	while (nd != 0) {
		n2m[nd].insert(ottID);
		nd = nd->GetEdgeToParentRef().GetParent();
	}
}

bool mrcaInThisLeafSet(const NxsSimpleNode * nd, const std::map<const NxsSimpleNode *, std::set<long> > & refNdp2mrca, const std::set<long> & leafSet) {
	std::vector<NxsSimpleNode *> children = nd->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree < 2) {
		return false;
	}
	bool foundFirstInf = false;
	for (std::vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
		const NxsSimpleNode * c = *cIt;
		if (refNdp2mrca.find(c) != refNdp2mrca.end()) {
			if (foundFirstInf) {
				return true;
			}
			foundFirstInf = true;
		}
	}
	return false;
}

void recordSupportedNodes(const std::map<const NxsSimpleNode *, std::set<long> > & refNdp2mrca,
						  std::set<std::set<long> > & sourceClades,
						  std::set<const NxsSimpleNode *> & supportedNodes,
						  const std::set<long> & leafSet) {
	for (std::map<const NxsSimpleNode *, std::set<long> >::const_iterator nsIt = refNdp2mrca.begin(); nsIt != refNdp2mrca.end(); ++nsIt) {
		const NxsSimpleNode * nd = nsIt->first;
		const NxsSimpleNode * par = nd->GetEdgeToParentRef().GetParent();
		if (par != 0) {
			const std::set<long> & nm = nsIt->second;
			if (mrcaInThisLeafSet(nd, refNdp2mrca, leafSet) && sourceClades.find(nm) != sourceClades.end()) {
				supportedNodes.insert(nd);
			}
		}
	}
}

void processSourceTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	std::map<const NxsSimpleNode *, std::set<long> > ndp2mrca;
	std::map<const NxsSimpleNode *, std::set<long> > refNdp2mrca;
	std::set<long> leafSet;
	std::vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode & nd = **nIt;
		std::vector<NxsSimpleNode *> children = nd.GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree > 0) {
			std::set<long> & mrca = ndp2mrca[&nd];
			for (std::vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
				const NxsSimpleNode * c = *cIt;
				assert(ndp2mrca.find(c) != ndp2mrca.end());
				std::set<long> & csl = ndp2mrca[c];
				mrca.insert(csl.begin(), csl.end());
			}
		} else if (outDegree == 0) {
			long ottID = getOTTIndex(tb, **nIt);
			assert(ottID >= 0);
			assert(ottID2TaxNode.find(ottID) != ottID2TaxNode.end());
			ndp2mrca[&nd].insert(ottID);
			markPathToRoot(refNdp2mrca, ottID);
			leafSet.insert(ottID);
		}
	}
	std::set<std::set<long> > sourceClades;
	for (std::map<const NxsSimpleNode *, std::set<long> >::const_iterator nsIt = ndp2mrca.begin(); nsIt != ndp2mrca.end(); ++nsIt) {
		const NxsSimpleNode * nd = nsIt->first;
		const NxsSimpleNode * par = nd->GetEdgeToParentRef().GetParent();
		if (par != 0) {
			const std::set<long> & nm = nsIt->second;
			sourceClades.insert(nm);
		}
	}
	recordSupportedNodes(refNdp2mrca, sourceClades, gSupportedNodes, leafSet);
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
	if (refTree == 0) {
		refTree = nst;
		processRefTree(taxa, nst);
	} else if (taxonomyTree == 0) {
		taxonomyTree = nst;
		processTaxonomyTree(taxa, nst);
	} else {
		processSourceTree(taxa, nst);
	}
	if (refTree != nst && taxonomyTree != nst) {
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
			readFilepathAsNEXUS(filepath, f);
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


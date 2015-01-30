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
const NxsSimpleNode * findMRCAFromIDSet(const map<long, const NxsSimpleNode *> & ref, const set<long> & idSet, long trigger);
template<typename T>
const std::string & getStrOrThrow(const T & nd, const std::map<T, std::string> & tipNameMap);


const NxsSimpleNode * findNextSignificantNode(const NxsSimpleNode * node, const map<const NxsSimpleNode *, set<long> > & ndp2mrca) {
	const NxsSimpleNode * currNode = node;
	for (;;) {
		map<const NxsSimpleNode *, set<long> >::const_iterator mIt = ndp2mrca.find(currNode);
		assert(mIt != ndp2mrca.end());
		const set<long> & oset = mIt->second;
		vector<NxsSimpleNode *> children = currNode->GetChildren();
		const NxsSimpleNode * sc = 0L;
		for (vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
			const NxsSimpleNode * child = *cIt;
			const map<const NxsSimpleNode *, set<long> >::const_iterator nIt = ndp2mrca.find(child);
			if (nIt != ndp2mrca.end()) {
				if (sc == 0L) {
					sc = child;
				} else {
					return currNode; // more than one child is in ndp2mrca, so this node is significant
				}
			}
		}
		if (sc == 0L) {
			cerr << "Failing. Node found with ottIDs marked, but no children with ottIDs marked:\n";
			for (set<long>::const_iterator oIt = oset.begin(); oIt != oset.end(); ++oIt) {
				if (oIt != oset.begin()) {
					cerr << ", ";
				}
				cerr << *oIt;
			}
			cerr << endl;
			assert(false);
		}
		const set<long> & dset = ndp2mrca.find(sc)->second;
		if (dset != oset) {
			cerr << "Failing. Internal node found with an ottID assignment. At this point the ottID should map to leaves. Par ottIDs:\n";
			for (set<long>::const_iterator oIt = oset.begin(); oIt != oset.end(); ++oIt) {
				if (oIt != oset.begin()) {
					cerr << ", ";
				}
				cerr << *oIt;
			}
			cerr << endl;
			cerr << "child ottIDs:\n";
			for (set<long>::const_iterator oIt = dset.begin(); oIt != dset.end(); ++oIt) {
				if (oIt != dset.begin()) {
					cerr << ", ";
				}
				cerr << *oIt;
			}
			cerr << endl;
			assert(false);
		}
		currNode = sc;
	}

}


void writeSubtreeNewickOTTIDs(ostream &out, const NxsSimpleNode * node, const map<const NxsSimpleNode *, set<long> > & ndp2mrca) {
	map<const NxsSimpleNode *, set<long> >::const_iterator nIt = ndp2mrca.find(node);
	assert(nIt != ndp2mrca.end());
	const set<long> & ottIDSet = nIt->second;
	if (nIt->second.size() == 1) {
		const long ottID = *ottIDSet.begin();
		out << "ott" << ottID;
	} else {
		assert(nIt->second.size() > 1);
		const NxsSimpleNode * nsn = findNextSignificantNode(node, ndp2mrca);
		out << '(';
		unsigned numcwritten = 0;
		vector<NxsSimpleNode *> children = nsn->GetChildren();
		for (vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
			const NxsSimpleNode * child = *cIt;
			nIt = ndp2mrca.find(child);
			if (nIt != ndp2mrca.end()) {
				if (numcwritten > 0) {
					out << ',';
				}
				writeSubtreeNewickOTTIDs(out, child, ndp2mrca);
				++numcwritten;
			}
		}
		assert(numcwritten > 1);
		out << ')';
	}
}

const NxsSimpleNode * findMRCAFromIDSet(const map<long, const NxsSimpleNode *> & ref,
	                           			const set<long> & idSet, long trigger) {
	map<const NxsSimpleNode *, unsigned int> n2c;
	long shortestPathLen = -1;
	const NxsSimpleNode * shortestPathNode = 0;
	for (set<long>::const_iterator toIt = idSet.begin(); toIt != idSet.end(); ++toIt) {
		const NxsSimpleNode * nd = 0;
		map<long, const NxsSimpleNode *>::const_iterator rIt = ref.find(*toIt);
		if (rIt == ref.end()) {
			cerr << "tip " << *toIt << " a descendant of " << trigger << " not found.\n";
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
	const unsigned nTips = idSet.size();
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

void writeSet(std::ostream & out, const char *indent, const set<long> &fir, const char * sep) {
	for (set<long>::const_iterator rIt = fir.begin(); rIt != fir.end(); ++rIt) {
		if (rIt != fir.begin()) {
			out << sep;
		}
		out << indent << "ott" << *rIt;
	}
}

template<typename T>
inline const std::string & getStrOrThrow(const T & nd, const std::map<T, std::string> & tipNameMap) {
	typename std::map<T, string>::const_iterator tnIt = tipNameMap.find(nd);
	if (tnIt == tipNameMap.end()) {
		throw NxsException("AssertionError: Key not found");
	}
	return tnIt->second;
}

const string & getLeftmostDesName(const NxsSimpleNode *nd, const map<const NxsSimpleNode *, string> & tipNameMap, bool useNdNames) {
	if (useNdNames) {
		const string & name = nd->GetName();
		if (!name.empty()) {
			return name;
		}
	}
	if (nd->IsTip()) {
		return getStrOrThrow(nd, tipNameMap);
	}
	return getLeftmostDesName(nd->GetFirstChild(), tipNameMap, useNdNames);
}

const string &  getRightmostDesName(const NxsSimpleNode *nd, const map<const NxsSimpleNode *, string> & tipNameMap, bool useNdNames) {
	if (useNdNames) {
		const string & name = nd->GetName();
		if (!name.empty()) {
			return name;
		}
	}
	if (nd->IsTip()) {
		return getStrOrThrow(nd, tipNameMap);
	}
	return getRightmostDesName(nd->GetLastChild(), tipNameMap, useNdNames);
}

void describeUnnamedNode(const NxsSimpleNode *nd, ostream & out, unsigned int anc, const map<const NxsSimpleNode *, string> & tipNameMap, bool useNdNames) {
	if (useNdNames && nd->GetName().length() > 0) {
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before \"" << nd->GetName() << "\"" << endl;
		} else {
			out << "the node \"" << nd->GetName() << "\"" << endl;
		}
		return;
	}
	vector<NxsSimpleNode *> children = nd->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree == 0) {
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before the leaf \"" << getStrOrThrow(nd, tipNameMap)  << "\"" << endl;
		} else {
			out << "the leaf \"" << getStrOrThrow(nd, tipNameMap)  << "\"" << endl;
		}
	} else if (outDegree == 1U) {
		describeUnnamedNode(children[0], out, anc + 1, tipNameMap, useNdNames);
	} else {
		string left = getLeftmostDesName(children[0], tipNameMap, useNdNames);
		string right = getRightmostDesName(children[outDegree - 1], tipNameMap, useNdNames);
		if (anc > 0) {
			out << "ancestor " << anc << " node(s) before MRCA of \"" << left << "\" and " << "\"" << right <<'\"' << endl;
		} else {
			out <<  "MRCA of \"" << left << "\" and " << "\"" << right <<'\"' << endl;
		}
	}
}


/* use some globals, because I'm being lazy... */
NxsSimpleTree * gRefTree = 0;
NxsSimpleTree * gTaxonTree = 0;
map<long, const NxsSimpleNode *> gOttID2RefNode;
map<const NxsSimpleNode *, string> gRefTipToName;
map<long, const NxsSimpleNode *> gOttID2TaxNode;
map<const NxsSimpleNode *, long> gTaxNode2ottID;
set<const NxsSimpleNode *> gSupportedNodes;
string gCurrentFilename;
string gCurrTmpFilepath;
ostream * gCurrTmpOstream = 0L;
bool gReadingTxtFile = false;
map<long, set<long> > gNonMono;
const bool gTrustNamedNodes = true;
map<const NxsSimpleNode *, long> gExpanded;
map<long, const NxsSimpleNode *> gTabooLeaf;
set<long> gTaxLeafOTTIDs;
map<const NxsSimpleNode *, set<long> > gAPrioriProblemNode;
bool gNoAprioriTests = true;
int gRefTreeNumNamedInternalsNodes = 0;
int gExitCode = 0;


void extendSupportedToRedundantNodes(const NxsSimpleTree * tree, set<const NxsSimpleNode *> & gSupportedNodes) {
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree == 1) {
			if (gSupportedNodes.find(children[0]) != gSupportedNodes.end() || children[0]->GetName().length() > 0) {
				if (gAPrioriProblemNode.find(nd) != gAPrioriProblemNode.end()) {
					assert(false); // shouldn't get out-degree one nodes w/ our designators
				}
				gSupportedNodes.insert(nd);
			}
		}
	}
}

bool singleDesSupportedOrNamed(const NxsSimpleNode *nd) {
	if (gSupportedNodes.find(nd) != gSupportedNodes.end()) {
		return true;
	}
	if (nd->GetOutDegree() == 1) {
		if (!nd->GetName().empty()) {
			return true;
		} else {
			return singleDesSupportedOrNamed(nd->GetFirstChild());
		}
	}
	return false;
}
bool IsRedundantNodeAroundTip(const NxsSimpleNode * nd) {
	if (nd->IsTip()) {
		return true;
	}
	if (nd->GetOutDegree() == 1) {
		return IsRedundantNodeAroundTip(nd->GetFirstChild());
	}
	return false;
}
int describeUnnamedUnsupported(ostream &out, const NxsSimpleTree * tree, const set<const NxsSimpleNode *> & supported) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin();
	int numUnsupported = 0;
	++nIt; //skip the root
	for (;nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree > 0 && supported.find(nd) == supported.end()) {
			if (IsRedundantNodeAroundTip(nd)) {
				//pass
			} else if (outDegree == 1 && gTrustNamedNodes && singleDesSupportedOrNamed(nd)) {
				//pass
			} else if (nd->GetName().length() == 0) { //assume that it is from the taxonomy
				if (gNoAprioriTests) {
					out << "Unsupported node ";
				} else {
					map<const NxsSimpleNode *, set<long> >::const_iterator gaIt = gAPrioriProblemNode.find(nd);
					if (gaIt == gAPrioriProblemNode.end()) {
						out << "Novel unsupported node ";
					} else {
						out << "Confirmation of unsupported node (designators =";
						writeSet(out, "", gaIt->second, " ");
						out << ") ";
					}
				}
				describeUnnamedNode(nd, out, 0, gRefTipToName, false);
				numUnsupported += 1;
			}
		}
	}
	return numUnsupported;
}


void summarize(ostream & out) {
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	extendSupportedToRedundantNodes(gRefTree, gSupportedNodes);
	int numUnsupported = describeUnnamedUnsupported(out, gRefTree, gSupportedNodes);
	for (map<const NxsSimpleNode *, set<long> >::const_iterator gaIt = gAPrioriProblemNode.begin(); gaIt != gAPrioriProblemNode.end(); ++gaIt) {
		if (gSupportedNodes.find(gaIt->first) != gSupportedNodes.end()) {
			out << "Claim of unsupported apparently refuted for designators: ";
			writeSet(out, "", gaIt->second, " ");
			out << ". See standard error stream for details.\n";
		}
	} 
	int supNonNamed = 0;
	int numSupportedInternals = 0;
	for (set<const NxsSimpleNode *>::const_iterator rIt = gSupportedNodes.begin(); rIt != gSupportedNodes.end(); ++rIt) {
		if (!(*rIt)->IsTip()) {
			numSupportedInternals += 1;
			if ((*rIt)->GetName().length() == 0) {
				supNonNamed += 1;
			}
		}
	}
	out << "\n\nFinal summary:\n";
	out << gRefTreeNumNamedInternalsNodes << " internal nodes were named in the reference tree. These were not rigorously checked against the taxonomy. They may not be detected as errors.\n";
	out << numSupportedInternals << " internal nodes where flagged as being supported by an input (including taxonomy).\n";
	int supNamed = numSupportedInternals - supNonNamed;
	out << "    " << supNamed << " of these were named (some of the support could just be the taxonomic expansion of tips).\n";
	out << "    " << supNonNamed << " of these were unnamed.\n";
	out << numUnsupported << " unsupported nodes.\n";
	out << endl;
	if (gExitCode < 0) {
		gExitCode -= numUnsupported;
	} else {
		gExitCode = numUnsupported;
	}
}



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
	const NxsSimpleNode * ndp = &nd;
	map<const NxsSimpleNode *, long>::const_iterator expIt = gExpanded.find(ndp);
	if (expIt != gExpanded.end()) {
		//cerr << "  shortcircuit returning " << expIt->second << endl;
		return expIt->second;
	}
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

void processRefTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, *nd);
		if (nd->IsTip()) {
			const unsigned ind = nd->GetTaxonIndex();
			assert(ind < tb->GetNumTaxonLabels());
			const string tn = tb->GetTaxonLabel(ind);
			gRefTipToName[nd] = tn;
		} else if (nd->GetName().length() > 0) {
			gRefTreeNumNamedInternalsNodes += 1;
		}
		if (ottID >= 0) {
			assert(gOttID2RefNode.find(ottID) == gOttID2RefNode.end());
			gOttID2RefNode[ottID] = nd;
		}
	}
}

void processTaxonomyTree(const NxsTaxaBlockAPI * tb, const NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_iterator nIt = nodes.begin(); nIt != nodes.end(); ++nIt) {
		const NxsSimpleNode * nd = *nIt;
		long ottID = getOTTIndex(tb, **nIt);
		assert(ottID >= 0);
		assert(gOttID2TaxNode.find(ottID) == gOttID2TaxNode.end());
		if (nd->IsTip()) {
			assert(gOttID2RefNode.find(ottID) != gOttID2RefNode.end());
			gTaxLeafOTTIDs.insert(ottID);
		}
		gOttID2TaxNode[ottID] = *nIt;
		gTaxNode2ottID[*nIt] = ottID;
	}
	for (map<long, const NxsSimpleNode *>::const_iterator nit = gOttID2RefNode.begin(); nit != gOttID2RefNode.end(); ++nit) {
		assert(gOttID2TaxNode.find(nit->first) != gOttID2TaxNode.end());
	}
}

//@recursive!
void fillTipOTTIDs(const map<long, const NxsSimpleNode *> &taxonomy, long ottID, set<long> & tipOTTIDs) {
	map<long, const NxsSimpleNode *>::const_iterator tnIt = taxonomy.find(ottID);
	assert(tnIt != taxonomy.end());
	const NxsSimpleNode * tn = tnIt->second;
	vector<NxsSimpleNode *> children = tn->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree == 0) {
		tipOTTIDs.insert(gTaxNode2ottID[tn]);
	} else {
		for (vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
			long ct = gTaxNode2ottID[*cIt];
			fillTipOTTIDs(taxonomy, ct, tipOTTIDs);
		}
	}
}

const NxsSimpleNode * findMRCA(const map<long, const NxsSimpleNode *> & ref,
	                           const map<long, const NxsSimpleNode *> &taxonomy, long ottID) {
	set<long> tipOTTIDs;
	fillTipOTTIDs(taxonomy, ottID, tipOTTIDs);
	const unsigned nTips = tipOTTIDs.size();
	if (nTips < 2) {
		cerr << "findMRCA called on " << ottID << '\n';
		assert(false);
	}
	gNonMono[ottID] = tipOTTIDs;
	return findMRCAFromIDSet(ref, tipOTTIDs, ottID);
}


void markPathToRoot(map<const NxsSimpleNode *, set<long> > &n2m, long ottID) {
	const NxsSimpleNode * nd = 0;
	map<long, const NxsSimpleNode *>::const_iterator fIt = gOttID2RefNode.find(ottID);
	if (fIt != gOttID2RefNode.end()) {
		nd = fIt->second;
	} else {
		assert(false);
		nd = findMRCA(gOttID2RefNode, gOttID2TaxNode, ottID);
	}
	assert(nd != 0);
	while (nd != 0) {
		n2m[nd].insert(ottID);
		nd = nd->GetEdgeToParentRef().GetParent();
	}
}

bool mrcaInThisLeafSet(const NxsSimpleNode * nd,
					   const map<const NxsSimpleNode *, set<long> > & refNdp2mrcaThisLeafSet) {
	vector<NxsSimpleNode *> children = nd->GetChildren();
	const unsigned outDegree = children.size();
	if (outDegree < 2) {
		return false;
	}
	bool foundFirstInf = false;
	for (vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
		const NxsSimpleNode * c = *cIt;
		const map<const NxsSimpleNode *, set<long> >::const_iterator rmIt = refNdp2mrcaThisLeafSet.find(c);
		if (rmIt != refNdp2mrcaThisLeafSet.end()) {
			if (foundFirstInf) {
				return true;
			}
			foundFirstInf = true;
		}
	}
	return false;
}

// find most recent anc of nd with out-degree > 1
const NxsSimpleNode * findFirstBranchingAnc(const NxsSimpleNode *nd) {
	const NxsSimpleNode * anc = nd->GetEdgeToParentRef().GetParent();
	assert(anc);
	while (anc->GetOutDegree() == 1) {
		anc = anc->GetEdgeToParentRef().GetParent();
	}
	return anc;
}

void recordSupportedNodes(const map<const NxsSimpleNode *, set<long> > & refNdp2mrca,
						  const map<set<long>, const NxsSimpleNode *> & sourceClades,
						  set<const NxsSimpleNode *> & supportedNodes,
						  const set<long> & leafSet,
						  const map<const NxsSimpleNode *, set<long> > & srcNdp2mrca) {
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	if (false) {//debugging
		cerr << "sourceClades:\n";
		for (map<set<long>, const NxsSimpleNode *>::const_iterator scit = sourceClades.begin(); scit != sourceClades.end(); ++scit) {
			cerr << "  clade: ";
			writeSet(cerr, "", scit->first, " ");
			cerr << "\n";
		}
	}

	for (map<const NxsSimpleNode *, set<long> >::const_iterator nsIt = refNdp2mrca.begin(); nsIt != refNdp2mrca.end(); ++nsIt) {
		const NxsSimpleNode * nd = nsIt->first;
		const bool printDB = false; //gAPrioriProblemNode.find(nd) != gAPrioriProblemNode.end();
		const NxsSimpleNode * par = nd->GetEdgeToParentRef().GetParent();
		if (par != 0) {
			const set<long> & nm = nsIt->second;
			const NxsSimpleNode * firstBranchingAnc = findFirstBranchingAnc(nd);
			map<const NxsSimpleNode *, set<long> >::const_iterator ancIt = refNdp2mrca.find(firstBranchingAnc);
			assert(ancIt != refNdp2mrca.end());
			const set<long> & anm = ancIt->second;
			if (printDB) { //debugging
				cerr << "DEBUGGING refTree node " << (long)nd << ": ";
				writeSet(cerr, "", nm, " ");
				cerr << "\n";
				cerr << "refTree par " << (long)firstBranchingAnc << ": ";
				writeSet(cerr, "", anm, " ");
				cerr << "\n";
			}
			if (mrcaInThisLeafSet(nd, refNdp2mrca)) {
				if (anm != nm) {
					map<set<long>, const NxsSimpleNode *>::const_iterator scIt = sourceClades.find(nm);
					if (scIt != sourceClades.end()) {
						if (gAPrioriProblemNode.find(nd) != gAPrioriProblemNode.end()) {
							ostream & out = cerr;
							map<const NxsSimpleNode *, set<long> >::const_iterator apIt = gAPrioriProblemNode.find(nd);
							out << "ERROR!: a priori unsupported node found. Designators were ";
							writeSet(out, "", apIt->second, " ");
							out << ". A node was found, which (when pruned to the leaf set of an input tree) contained:\n";
							writeSet(out, "    ", nm, " ");
							out << "\nThe subtree from the source was: ";
							const NxsSimpleNode * srcNd = scIt->second;
							writeSubtreeNewickOTTIDs(out, srcNd, srcNdp2mrca);
							gExitCode = -1;
						}
						supportedNodes.insert(nd);
						if (printDB) { //debugging
							cerr <<"SUPPORTED!\n";
						}
					} else {
						if (printDB) { //debugging
							cerr <<"UNsupported. node mrca set not in source tree\n";
						}
					} 
				} else {
					if (printDB) { //debugging
						cerr <<"UNsupported. node mrca set == first branching ancestors mrca set on this leaf set\n";
					}
				}
			} else {
				if (printDB) { //debugging
					cerr <<"UNsupported. not a mrcaInThisLeafSet\n";
				}
			}
		}
	}
}

void writeNewickOTTIDs(ostream &out, const NxsSimpleTree * tree, const map<const NxsSimpleNode *, set<long> > & ndp2mrca) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	const NxsSimpleNode * root = nodes[0];
	writeSubtreeNewickOTTIDs(out, root, ndp2mrca);
	out << ";\n";
}

void expandOTTInternalsWhichAreLeaves(const NxsTaxaBlockAPI * tb, NxsSimpleTree * tree) {
	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	map<NxsSimpleNode *, set<long> > replaceNodes;
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode & nd = **nIt;
		const unsigned outDegree = nd.GetOutDegree();
		if (outDegree == 0) {
			long ottID = getOTTIndex(tb, **nIt);
			assert(ottID >= 0);
			assert(gOttID2TaxNode.find(ottID) != gOttID2TaxNode.end());
			if (gTaxLeafOTTIDs.find(ottID) == gTaxLeafOTTIDs.end()) {
				set<long> leafSet;
				fillTipOTTIDs(gOttID2TaxNode, ottID, leafSet);
				replaceNodes[const_cast<NxsSimpleNode *>(&nd)] = leafSet;
			}
		}
	}
	for (map<NxsSimpleNode *, set<long> >::const_iterator rIt = replaceNodes.begin(); rIt != replaceNodes.end(); ++rIt) {
		NxsSimpleNode * oldNode = rIt->first;
		const set<long> & leafSet = rIt->second;
		assert(leafSet.size() > 0);
		oldNode->SetTaxonIndex(UINT_MAX); // make this no longer appear to be a tip
		for (set<long>::const_iterator lsIt = leafSet.begin(); lsIt != leafSet.end(); ++lsIt) {
			NxsSimpleNode *newNode =  tree->AllocNewNode(oldNode);
			oldNode->AddChild(newNode);
			gExpanded[newNode] = *lsIt;
			assert(gTabooLeaf.find(*lsIt) == gTabooLeaf.end());
			gTabooLeaf[*lsIt] = newNode;
		}
	}
}

void processSourceTree(const NxsTaxaBlockAPI * tb, NxsSimpleTree * tree) {
	expandOTTInternalsWhichAreLeaves(tb, tree);
	map<const NxsSimpleNode *, set<long> > ndp2mrca;
	map<const NxsSimpleNode *, set<long> > refNdp2mrca;
	set<long> leafSet;

	vector<const NxsSimpleNode *> nodes =  tree->GetPreorderTraversal();
	for (vector<const NxsSimpleNode *>::const_reverse_iterator nIt = nodes.rbegin(); nIt != nodes.rend(); ++nIt) {
		const NxsSimpleNode & nd = **nIt;
		vector<NxsSimpleNode *> children = nd.GetChildren();
		const unsigned outDegree = children.size();
		if (outDegree > 0) {
			set<long> & mrca = ndp2mrca[&nd];
			for (vector<NxsSimpleNode *>::const_iterator cIt = children.begin(); cIt != children.end(); ++cIt) {
				const NxsSimpleNode * c = *cIt;
				assert(ndp2mrca.find(c) != ndp2mrca.end());
				set<long> & csl = ndp2mrca[c];
				mrca.insert(csl.begin(), csl.end());
			}
		} else if (outDegree == 0) {
			long ottID = getOTTIndex(tb, **nIt);
			map<long, const NxsSimpleNode *>::const_iterator tlIt = gTabooLeaf.find(ottID);
			if (tlIt != gTabooLeaf.end()) {
				assert(tlIt->second == *nIt);
			}
			assert(ottID >= 0);
			assert(gOttID2TaxNode.find(ottID) != gOttID2TaxNode.end());
			ndp2mrca[&nd].insert(ottID);
			markPathToRoot(refNdp2mrca, ottID);
			leafSet.insert(ottID);
		}
	}
	map<set<long>, const NxsSimpleNode *> sourceClades;
	for (map<const NxsSimpleNode *, set<long> >::const_iterator nsIt = ndp2mrca.begin(); nsIt != ndp2mrca.end(); ++nsIt) {
		const NxsSimpleNode * nd = nsIt->first;
		const NxsSimpleNode * par = nd->GetEdgeToParentRef().GetParent();
		if (par != 0 && !nd->IsTip()) {
			const set<long> & nm = nsIt->second;
			sourceClades[nm] = nsIt->first;
		}
	}
	recordSupportedNodes(refNdp2mrca, sourceClades, gSupportedNodes, leafSet, ndp2mrca);
	if (gCurrTmpOstream != 0) {
		*gCurrTmpOstream << "#NEXUS\nBEGIN TREES;\n";

		*gCurrTmpOstream << "   Tree pruned_synth = [&R] ";
		writeNewickOTTIDs(*gCurrTmpOstream, gRefTree, refNdp2mrca);
		
		*gCurrTmpOstream << "   Tree input = [&R] ";
		writeNewickOTTIDs(*gCurrTmpOstream, tree, ndp2mrca);
		
		*gCurrTmpOstream << "END;\n";
	}
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
	NxsSimpleTree * nst = new NxsSimpleTree(ftd, 0.0, 0, true);
	gExpanded.clear();
	gTabooLeaf.clear();
	if (gRefTree == 0) {
		gRefTree = nst;
		processRefTree(taxa, nst);
	} else if (gTaxonTree == 0) {
		gTaxonTree = nst;
		processTaxonomyTree(taxa, nst);
	} else {
		processSourceTree(taxa, nst);
	}
	if (gRefTree != nst && gTaxonTree != nst) {
		delete nst;
	}
	return false;
}

void markSuspectNode(const set<long> & designators) {
	const NxsSimpleNode * mrca = findMRCAFromIDSet(gOttID2RefNode, designators, -1);
	assert(mrca->GetName().length() == 0);
	gAPrioriProblemNode[mrca] = designators;
}


void parseAndProcessMRCADesignatorsFile(string filepath) {
	if (gRefTree == 0L || gTaxonTree != 0L) {
		cerr << "gRefTree" << (long )gRefTree << "\n";
		cerr << "gTaxonTree" << (long )gTaxonTree << "\n";
		cerr << "\nDesignators file must come after the full tree estimate, but before the taxonomy in the argument length\n";
		throw exception();
	}
	ifstream inpf(filepath.c_str());
	assert(inpf.good());
	string line;
	while (getline(inpf, line)) {
		string stripped = NxsString::strip_surrounding_whitespace(line);
		if (!stripped.empty()) {
			list<string> words;
			NxsString::split(stripped, &words);
			if (words.size() < 2) {
				cerr << "Expecting >1 designator. Found: " << line << "\n";
				throw exception();
			}
			set<long> designators;
			for (list<string>::const_iterator dIt = words.begin(); dIt != words.end(); ++dIt) {
				long d;
				if (!NxsString::to_long(dIt->c_str(), &d)) {
					cerr << "Expecting numeric designator. Found: " << line << "\n";
					throw exception();
				}
				designators.insert(d);
			}
			markSuspectNode(designators);
		}
	}
	assert(gAPrioriProblemNode.empty() == gNoAprioriTests);
	inpf.close();
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
	vector<string> fmtNames =  MultiFormatReader::getFormatNames();
	for (vector<string>::const_iterator n = fmtNames.begin(); n != fmtNames.end(); ++n)
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
				string fmtName(filepath + 2, slen - 2);
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
			const string filepathstr(filepath);
			const size_t sp = filepathstr.find_last_of('/');
			if (sp == string::npos)
				{
				gCurrentFilename = filepathstr;
				}
			else
				{
				gCurrentFilename = filepathstr.substr(sp + 1);
				}
			gCurrTmpFilepath = string("tmp/") + gCurrentFilename;
			cerr << "gCurrTmpFilepath = " << gCurrTmpFilepath << '\n';
			ofstream tostream(gCurrTmpFilepath.c_str());
			gCurrTmpOstream = &tostream;
			const size_t fnl = gCurrentFilename.length();
			if (gCurrentFilename.substr(fnl - 4) == string(".txt")) {
				gReadingTxtFile = true;
			} else {
				gReadingTxtFile = false;
			}
			try {
				if (gReadingTxtFile) {
					gNoAprioriTests = false;
					cerr << "found text file. Treating each line as a series of MRCA designators. Nodes defined by these designators are asserted to be unsupported. If support is found, execution will stop at the supporting tree.\n";
					parseAndProcessMRCADesignatorsFile(filepath);
				} else {
					readFilepathAsNEXUS(filepath, f);
				}
				tostream.close();
				}
			catch (...)
				{
				tostream.close();
				cerr << "\nExiting due to and exception.\n";
				return 1;
				}
			gCurrTmpOstream = 0L;
			}
		}
	if (!readfile)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
		}
	summarize(cout);
	return gExitCode;
	}


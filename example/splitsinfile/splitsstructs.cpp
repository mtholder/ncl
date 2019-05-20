#include "splitsstructs.h"
#include <algorithm>
bool TreesToSplits::gTrackTrivial = false;
bool TreesToSplits::gTreatAsRooted = false;
bool TreesToSplits::gTrackFreq = false;
bool TreesToSplits::gTrackOccurrence = false;
bool TreesToSplits::gTrackEdgeLen = false;
bool TreesToSplits::gTrackEdgeLenSummary = false;
bool TreesToSplits::gTrackHeight = false;
bool TreesToSplits::gTrackHeightSummary = false;


const WORD_INT_TYPE gLastMask[] = {1, 3, 7, 0x0F, 0x1F, 0x3F, 0x7F, 0x7F };
const WORD_INT_TYPE gBit[] = {1, 2, 4, 8, 0x10, 0x20, 0x40, 0x80 };
Split::Split(unsigned nTax)
	:splitRep(1 + ((nTax - 1)/NBITS_IN_WORD), 0),
	lastMask(gLastMask[(nTax-1) % NBITS_IN_WORD])
	{
	}

void Split::setIndex(unsigned ind)
{
	const unsigned elIndex = ind/NBITS_IN_WORD;
	assert(elIndex < splitRep.size());
	WORD_INT_TYPE & el = splitRep[elIndex];
	el |= gBit[ind % NBITS_IN_WORD ];
}


unsigned gCounter = 0;

void TreesToSplits::recordTree(const NxsFullTreeDescription & ftd, NxsTaxaBlockAPI *taxB) {
	TaxaBlockToSplitsMap::iterator tsmIt = taxaBlockToSplitsMap.find(taxB);
	if (tsmIt == taxaBlockToSplitsMap.end()) {
		taxaBlockToSplitsMap[taxB] = NTreesSplitsMap(0, SplitsMap());
		tsmIt = taxaBlockToSplitsMap.find(taxB);
	}
	NTreesSplitsMap & ntsm = tsmIt->second;
	SplitsMap & sm = ntsm.second;
	unsigned nextNTax = taxB->GetNTax();
	if (nextNTax != nTax)
		splits.clear();
	nTax = nextNTax;
	recordTreeToSplitsMap(ftd, sm, ntsm.first);
	ntsm.first += 1;
}

void TreesToSplits::recordTreeToSplitsMap(const NxsFullTreeDescription & ftd, TreesToSplits::SplitsMap & sm, unsigned treeIndex) {
	nst.Initialize(ftd);
	if (!treatAsRooted)
		nst.RerootAt(0);
	std::vector<const NxsSimpleNode *> nodes =  nst.GetPreorderTraversal();
	if (nodes.empty())
		return;
	std::vector<const NxsSimpleNode *>::reverse_iterator last =  nodes.rend();
	--last; /* root is a trivial split */
	while (splits.size() < nodes.size())
		splits.push_back(SplitAndScratch(Split(nTax), NULL));

	const bool doHeights = (trackHeight || trackHeightSummary);
	if (doHeights) {
		while (dblScratch.size() < nodes.size())
			dblScratch.push_back(0.0);
	}
	std::vector<double>::iterator hIt = dblScratch.begin();

	std::vector<SplitAndScratch>::iterator spIt = splits.begin();
	for (std::vector<const NxsSimpleNode *>::reverse_iterator ndIt = nodes.rbegin(); ndIt != last; ++ndIt, ++spIt) {
		const NxsSimpleNode * nd = *ndIt;
		//cout << "visiting node " << (long) nd << '\n';
		assert(nd);
		assert(spIt != splits.end());
		SplitAndScratch & ss = *spIt;
		Split & split = ss.first;
		split.reset();
		const NxsSimpleNode * c = nd->GetFirstChild();
		nd->scratch = (void *)&ss;
		double * h = 0L;
		if (doHeights) {
			h = &(*hIt);
			ss.second = (void *) h;
			++hIt;
			*h = 0.0;
		}
		if (c) {
			while (c) {
				const SplitAndScratch * css = (SplitAndScratch *)(c->scratch);
				const Split * cSplit = &(css->first);
				assert(cSplit);
				split.setToUnion(*cSplit);
				if (doHeights) {
					NxsSimpleEdge nse = c->GetEdgeToParent();
					double n = nse.GetDblEdgeLen();
					double * ch = (double*)(css->second);
					n += *ch;
					*h = std::max(*h, n);
				}
				c = c->GetNextSib();
			}
			if (!treatAsRooted)
				split.invertIfNeeded();
			recordSplit(split, nd, sm, treeIndex);
		}
		else {
			split.setIndex(nd->GetTaxonIndex());
			if (!treatAsRooted)
				split.invertIfNeeded();
			if (trackTrivial)
				recordSplit(split, nd, sm, treeIndex);
		}
	}
}

void TreesToSplits::recordSplit(const Split &split, const NxsSimpleNode *nd, TreesToSplits::SplitsMap & sm, unsigned treeIndex) {
	gCounter++;
	TreesToSplits::SplitsMap::iterator smIt = sm.find(split);
	if (smIt == sm.end()) {
		sm[split] = SplitInfo();
		smIt = sm.find(split);
		//std::cout << "New Split " << gCounter << ": " << split.getNewick() << '\n';
	}
	else {
		//std::cout << "Repeated Split " << gCounter << ": " << split.getNewick() << '\n';
	}

	SplitInfo & info = smIt->second;
	if (trackFreq || trackEdgeLenSummary || trackEdgeLen)
		info.nTimes += 1;
	if (trackOccurrence)
		info.inclusion.insert(treeIndex);

	if (trackEdgeLenSummary) {
		NxsSimpleEdge nse = nd->GetEdgeToParent();
		const double el = nse.GetDblEdgeLen();
		info.edgeLenSum += el;
		info.edgeLenSumSq += el*el;
	}
	else if (trackEdgeLen){
		NxsSimpleEdge nse = nd->GetEdgeToParent();
		const double el = nse.GetDblEdgeLen();
		info.edgeLen.push_back(el);
	}

	if (trackHeightSummary)
		{
		const SplitAndScratch * css = (SplitAndScratch *)(nd->scratch);
		double *h = (double *)css->second;
		info.heightSum += *h;
		}
	else if (trackHeight)
		{}
}

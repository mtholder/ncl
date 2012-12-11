#if ! defined(SPLITSSTRUCTS_HPP)
#define SPLITSSTRUCTS_HPP
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>

#include "ncl/nxsdefs.h"
#include "ncl/nxstreesblock.h"
class NxsTaxaBlockAPI;
class SplitInfo
{
	public:
		SplitInfo()
			:edgeLenSum(0.0),
			edgeLenSumSq(0.0),
			heightSum(0.0),
			heightSumSq(0.0),
			nTimes(0)
			{}
		void reset() {
			edgeLen.clear();
			heights.clear();
			inclusion.clear();
			edgeLenSum = 0.0;
			edgeLenSumSq = 0.0;
			heightSum = 0.0;
			heightSumSq = 0.0;
			nTimes = 0;
		}
		std::vector<double> edgeLen;
		std::vector<double> heights;
		NxsUnsignedSet inclusion;
		double edgeLenSum;
		double edgeLenSumSq;
		double heightSum;
		double heightSumSq;
		unsigned nTimes;
};

#define WORD_INT_TYPE unsigned char
#define NBITS_IN_WORD (8*sizeof(WORD_INT_TYPE))
class Split
{
	public:
		Split(unsigned nTax);
		void reset() {
			unsigned size = splitRep.size();
			splitRep.assign(size, 0);
		}
		void invertIfNeeded() {
			if (1 & (splitRep[0])) {
				//cout << "Inverting: " << getNewick(true) << '\n';
				std::vector<WORD_INT_TYPE>::iterator sIt = splitRep.begin();
				for (; sIt != splitRep.end(); ++sIt) {
					const WORD_INT_TYPE comp = ~(*sIt);
					*sIt = comp;
				}
				WORD_INT_TYPE & lastWord = *(splitRep.rbegin());
				lastWord &= lastMask;
				//cout << "Done inverting: " << getNewick(true) << '\n';
			}
		}
		void setIndex(unsigned ind);
		void setToUnion(const Split & other) {
			assert(splitRep.size() == other.splitRep.size());
			std::vector<WORD_INT_TYPE>::iterator sIt = splitRep.begin();
			std::vector<WORD_INT_TYPE>::const_iterator oIt = other.splitRep.begin();
			for (; sIt != splitRep.end(); ++sIt, ++oIt) {
				*sIt |= *oIt;
			}
		}
		std::string getNewick(bool evenIfTrivial=true) const {
			std::ostringstream out;
			if (writeNewick(out, 0.0, false, evenIfTrivial))
				return out.str();
			return std::string();
		}
		bool writeNewick(std::ostream &out, bool evenIfTrivial) const {
			return writeNewick(out, 0.0, false, evenIfTrivial);
		}
		bool writeNewick(std::ostream &out, double edgeLen, bool evenIfTrivial) const {
			return writeNewick(out, edgeLen, true, evenIfTrivial);
		}
		void dumpIndices(NxsUnsignedSet * inc,  NxsUnsignedSet *exc) const {
			WORD_INT_TYPE bit;
			unsigned word = 0;
			unsigned index = 0;
			const unsigned nWords = splitRep.size();
			if (nWords == 0)
				return;
			if (nWords > 1) {
				for (;word < nWords - 1; ++word) {
					const WORD_INT_TYPE &uc =  splitRep[word];
					bit = 1;
					for (unsigned i = 0; i < NBITS_IN_WORD; ++i, ++index) {
						if (bit & uc) {
							if (inc)
								inc->insert(index);
						}
						else if (exc)
							exc->insert(index);
						bit <<= 1;
					}
				}
			}
			const WORD_INT_TYPE &lastc =  splitRep[nWords - 1];
			bit = 1;
			for (unsigned i = 0; i < NBITS_IN_WORD; ++i, ++index) {
				if (bit >lastMask)
					return;
				if (bit & lastc) {
					if (inc)
						inc->insert(index);
				}
				else if (exc)
					exc->insert(index);
				bit <<= 1;
			}
		}
		bool writeNewick(std::ostream &out, double edgeLen, bool writeEdgeLen, bool evenIfTrivial)  const {
			NxsUnsignedSet inc;
			NxsUnsignedSet exc;
			dumpIndices(&inc, &exc);
			if (!evenIfTrivial && (inc.size() == 1 || exc.size() == 1))
				return false;
			out << '(';
			NxsUnsignedSet::const_iterator sIt = inc.begin();
			if (!inc.empty()) {
				out << '(' << (1 + *sIt);
				for (++sIt; sIt != inc.end(); ++sIt)
					out << ',' << (1 + *sIt);
				out << ')';
				if (writeEdgeLen)
					out << ':' << edgeLen;
			}

			if (!exc.empty()) {
				for (sIt = exc.begin(); sIt != exc.end(); ++sIt)
					out << ',' << (1 + *sIt);
				out << ')';
			}
			return true;
		}

		friend bool operator<(const Split & one, const Split & another);

	private:
		std::vector<WORD_INT_TYPE> splitRep;
		WORD_INT_TYPE lastMask;
};

inline bool operator<(const Split & one, const Split & another)
{
	const unsigned osize = one.splitRep.size();
	const unsigned asize = another.splitRep.size();
	if (osize == asize) {
		std::vector<WORD_INT_TYPE>::const_iterator oIt = one.splitRep.begin();
		std::vector<WORD_INT_TYPE>::const_iterator aIt = another.splitRep.begin();
		for (; oIt != one.splitRep.end(); ++oIt, ++aIt) {
			if (*aIt != *oIt) {
				return *oIt < *aIt;
			}
		}
		return false;
	}
	return (osize < asize);
}

class TreesToSplits
{
	public:
		typedef std::map<Split, SplitInfo> SplitsMap;
		typedef std::pair<unsigned, SplitsMap> NTreesSplitsMap;
		typedef std::pair<Split, void *> SplitAndScratch;
		typedef std::map<NxsTaxaBlockAPI *,  NTreesSplitsMap > TaxaBlockToSplitsMap;
		static bool gTrackTrivial;
		static bool gTreatAsRooted;
		static bool gTrackFreq;
		static bool gTrackOccurrence;
		static bool gTrackEdgeLen;
		static bool gTrackEdgeLenSummary;
		static bool gTrackHeight;
		static bool gTrackHeightSummary;

		TreesToSplits()
		  :trackTrivial(gTrackTrivial),
		  treatAsRooted(gTreatAsRooted),
		  trackFreq(gTrackFreq),
		  trackOccurrence(gTrackOccurrence),
		  trackEdgeLen(gTrackEdgeLen),
		  trackEdgeLenSummary(gTrackEdgeLenSummary),
		  trackHeight(gTrackHeight),
		  trackHeightSummary(gTrackHeightSummary),
		  nst(0, 0.0),
		  nTax(0)
		   	{
		   	}

		void recordTree(const NxsFullTreeDescription & ftd, NxsTaxaBlockAPI *taxB);

		bool trackTrivial;
		bool treatAsRooted;
		bool trackFreq;
		bool trackOccurrence;
		bool trackEdgeLen;
		bool trackEdgeLenSummary;
		bool trackHeight;
		bool trackHeightSummary;
		TaxaBlockToSplitsMap taxaBlockToSplitsMap;

	protected:
		void recordTreeToSplitsMap(const NxsFullTreeDescription & ftd, SplitsMap & sm, unsigned treeIndex);

		void recordSplit(const Split &split, const NxsSimpleNode *nd, SplitsMap & sm, unsigned treeIndex);


		NxsSimpleTree nst;
		std::vector<SplitAndScratch> splits;
		std::vector<double> dblScratch;
		unsigned nTax;
};

#endif

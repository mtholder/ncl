#ifndef NORMALIZER_H
#define NORMALIZER_H

#include <string>
#include <iostream>
class PublicNexusReader;

void writeAsNexus(PublicNexusReader & nexusReader, std::ostream & os);

class TranslatingConventions {
	public:
		TranslatingConventions()
			:idPrefix(),
			globalIncrementingIDs(false),
			currentOTUIndex(0),
			currentOTUsIndex(0),
			currentCharsIndex(0),
			currentCharIndex(0),
			currentStateIndex(0),
			currentStateSetIndex(0),
			currentTreesIndex(0),
			currentTreeIndex(0),
			currentNodeIndex(0),
			currentEdgeIndex(0),
			treatNodeLabelsAsStrings(false),
			emitTreesAndTaxaOnly(false) {
			}
	std::string idPrefix;
	bool globalIncrementingIDs;
	unsigned currentOTUIndex;
	unsigned currentOTUsIndex;
	unsigned currentCharsIndex;
	unsigned currentCharIndex;
	unsigned currentRowIndex;
	unsigned currentStateIndex;
	unsigned currentStateSetIndex;
	unsigned currentTreesIndex;
	unsigned currentTreeIndex;
	unsigned currentNodeIndex;
	unsigned currentEdgeIndex;
	bool treatNodeLabelsAsStrings;
	bool emitTreesAndTaxaOnly;
};

#endif
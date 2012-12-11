#include<string>
#include<vector>
#include<cassert>
#include "ncl/nxsmultiformat.h"

/*!
	Fills the `matrix` with a new encoding for the data.

	This demonstrates how you can access the internal representation of the class
		NxsCharactersBlock (the v2.1 character querying API)

	In this simple example, we reject any characters block that does not have
		a simple unadulterated.

	For simplicity's sake, we will also convert any missing data, and gaps
		to the binary code 1111 (0xF in hexadecimal) which indicates that all
		four bases could be found.  This is easy to do because the "fundamental"
		states will be the ones with state codes [0,nstates), in other words 0,1,2,3 will be
		the codes for A, C, G, and T in DNA.  missing and gaps will be the only valid
		states with negative codes, and all ambiguity codes will will have a value
		>= nstates.

	We'll recode the data using the following mapping in which each base represents
		a different bit:
		1 <-> A,
		2 <-> C
		4 <-> G
		8 <-> T

*/
void convertToBitsForMatrix(NCL_COULD_BE_CONST NxsCharactersBlock *charBlock, std::vector<std::vector<char> > * converted)
{
	if (converted == 0L)
		return;
	converted->clear();
	if (charBlock == 0L)
		return;
	const NxsDiscreteDatatypeMapper *mapper = charBlock->GetDatatypeMapperForChar(0);
	const bool hasGaps = charBlock->GetGapSymbol() != '\0';

	NxsDiscreteDatatypeMapper defaultMapper(NxsCharactersBlock::dna, hasGaps);
	if (!mapper->IsSemanticallyEquivalent(defaultMapper))
		{
		std::string errormsg = "The convertToBitsForMatrix example function only supports DNA datatypes without any new combinations of states";
		throw NxsException(errormsg);
		}

	/* First we'll create a vector that maps NCL's state codes to the appropriate
	 binary encoding.

	*/
	std::vector<char> stateCodeToBinaryMapping;
	const NxsDiscreteStateCell highestLegalStateCode = mapper->GetHighestStateCode();
	assert(highestLegalStateCode < 15);// there are 15 possible states and combinations of states for dna
	stateCodeToBinaryMapping.resize(highestLegalStateCode);

	/*
	Let's start by recording the "fundamental" (unambiguous) states which are
	 guaranteed by NCL's api to be the state codes from 0 up to number_of_states - 1
	*/
	NxsDiscreteStateCell nclCode = mapper->GetStateCodeStored('A');
	assert(nclCode >= 0 && nclCode < 4);
	std::cerr << "Adding fundamental state "  << nclCode << " 1\n";
	stateCodeToBinaryMapping[nclCode] = (char) 1; // ncl code for A goes to 1 in our binary coding...

	nclCode = mapper->GetStateCodeStored('C');
	assert(nclCode >= 0 && nclCode < 4);
	std::cerr << "Adding fundamental state "  << nclCode << " 2\n";
	stateCodeToBinaryMapping[nclCode] = (char) 2; // ncl code for C goes to 2 in our binary coding...

	nclCode = mapper->GetStateCodeStored('G');
	assert(nclCode >= 0 && nclCode < 4);
	std::cerr << "Adding fundamental state "  << nclCode << " 4\n";
	stateCodeToBinaryMapping[nclCode] = (char) 4; // ncl code for G goes to 4 in our binary coding...

	nclCode = mapper->GetStateCodeStored('T');
	assert(nclCode >= 0 && nclCode < 4);
	std::cerr << "Adding fundamental state "  << nclCode << " 8\n";
	stateCodeToBinaryMapping[nclCode] = (char) 8; // ncl code for T goes to 8 in our binary coding

	std::cerr << "Adding state sets up to  "  << highestLegalStateCode << "\n";
	/* Now we can deal with all of the ambiguity state codes by producing the binary
		code that is the bitwise union of all states.  We'll walk through all of the number_of_states .. highest state code
		and register the conversion of NCL code to our representation.
	*/
	for (NxsDiscreteStateCell currNCLCode = 4; currNCLCode <= highestLegalStateCode; ++currNCLCode)
		{
		const std::set<NxsDiscreteStateCell> & stateSet = mapper->GetStateSetForCode(currNCLCode);
		char unionOfStateBits = 0;
		std::set<NxsDiscreteStateCell>::const_iterator stateSetIt = stateSet.begin();
		for (; stateSetIt != stateSet.end(); ++stateSetIt)
			{
			const NxsDiscreteStateCell fundamentalStateCode = *stateSetIt;
			assert(fundamentalStateCode >= 0 && fundamentalStateCode < 4); // all of the states in the state set are the "fundamental" states -- never codes for sets of states
												 // hence we can assume that they will already be added to stateCodeToBinaryMapping
												 // when we added A,C,G,and T above.
			char bitForCurrentState = stateCodeToBinaryMapping[(int) fundamentalStateCode];
			unionOfStateBits |= bitForCurrentState;
			}
		stateCodeToBinaryMapping[currNCLCode] = unionOfStateBits;
		std::cerr << "Adding state set" << currNCLCode <<  ' ' << highestLegalStateCode << "\n";

		}

	std::cerr << "stateCodeToBinaryMapping is ";
	for (unsigned i = 0; i < highestLegalStateCode; ++i)
		std::cerr << i << ' ' << (int) stateCodeToBinaryMapping[i] << std::endl;

	/* Now we are ready to convert each row of the matrix to the binary encoding

		Because the matrices can be large, we'll try to avoid creating a vector of
		binary codes and pushing it onto the `converted` vector (that involves copying).

		Instead we will:
			1. resize `converted` to the correct size (so that it does not have to reallocate during the subsequent operations).
			2. get a reference to the vector that makes up each row.
			3. resize that row.
			4. add the codes to that row.
	*/
	const char allBasesBinaryCode = 15; // the union of the bits for A, C, G, and T is 15
	const unsigned nc = charBlock->GetNCharTotal();
	const unsigned nt = charBlock->GetNTaxTotal();
	converted->resize(nt); // create the nt rows in the converted matrix.
	for (unsigned taxIndex = 0; taxIndex < nt; ++taxIndex)
		{
		std::vector<char> & binaryRow = (*converted)[taxIndex];  // make sure to store the reference here, so we don't make a copy of the row.
		binaryRow.reserve(nc); // make the row alloc enough space
		const NxsDiscreteStateRow & nclRow = charBlock->GetDiscreteMatrixRow(taxIndex);
		NxsDiscreteStateRow::const_iterator nclCodeIt = nclRow.begin();
		for (; nclCodeIt != nclRow.end(); ++nclCodeIt)
			{
			const NxsDiscreteStateCell & currNCLStateCode = *nclCodeIt;
			if (currNCLStateCode < 0)
				{ // missing or gap. We could check which one by testing for whether
					// the code is NXS_MISSING_CODE or NXS_GAP_STATE_CODE
					// since this scrip maps both to the allBasesCode (15)
					// we won't bother checking here.
				binaryRow.push_back(allBasesBinaryCode);
				}
			else
				{
				assert(currNCLStateCode <= highestLegalStateCode); // this constraint is guaranteed by NCL's API
				const char currBinaryCode = stateCodeToBinaryMapping[currNCLStateCode];
				binaryRow.push_back(currBinaryCode);
				}
			}
		}

	/* That's it!
		We passed in the converted matrix as pointer so that we would not have to
		copy it on return.

		We created a mapping of NCL's state codes to our own.
		Becaus we know that NCL's codes will be small integers, then a vector
		is more efficient than a std::map<NxsDiscreteStateCell, char> (but a std::map would
		have worked, too).

		Specifically we know that:

		-2 for NXS_GAP_STATE_CODE,
		-1 for NXS_MISSING_CODE,
		0 up to (number_of_states - 1) for the "fundamental" states (and they will be in the
			order of the symbols string)
		number_of_states up to and including highestLegalStateCode for ambiguous states

		So we can write conversion code that exploits these facts
	*/



}


/*!
	writes the contents of the matrix out to `out' as NEXUS.


	This demonstrates the "old" (the v2.0 character querying API).

*/
void printMatrixWithoutLabels(std::ostream & out, NCL_COULD_BE_CONST  NxsCharactersBlock *cb)
{
	const char missingChar = cb->GetMissingSymbol();
	char gapChar = cb->GetGapSymbol();
	const unsigned nc = cb->GetNCharTotal();
	const unsigned nt = cb->GetNTaxTotal();
	for (unsigned t = 0; t < nt; ++t)
		{
		for (unsigned c = 0; c < nc; ++c)
			{
			if (cb->IsMissingState(t, c))
				out << missingChar;
			else if (cb->IsGapState(t, c))
				{
				assert (gapChar != '\0'); //matrices withouth gapSymbols should not have cells with gaps!
				out << '-';
				}
			else
				{
				const unsigned ns = cb->GetNumStates(t, c);
				if (ns == 1)
					{
					out << cb->GetState(t, c, 0);
					}
				else
					{
					const bool isPoly = cb->IsPolymorphic(t, c);
					if (isPoly)
						out << '(';
					else
						out << '{';
					for (unsigned s = 0; s < ns; ++s)
						out << cb->GetState(t, c, s);
					if (isPoly)
						out << ')';
					else
						out << '}';
					}
				}
			}
		out<< '\n';
		}
}



int main(int argc, char * argv[])
{
	/* Usually it is easiest to surround interactions with NCL in a try block
		that catches NxsException instances.
		Reading files can certainly generate these errors, but even queries
		after the parse can result in these exceptions.
	*/
	try {
		MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);

		/* See discussion in "Making NCL less strict" */
		NxsCharactersBlock * charsB = nexusReader.GetCharactersBlockTemplate();
		NxsDataBlock * dataB = nexusReader.GetDataBlockTemplate();
		charsB->SetAllowAugmentingOfSequenceSymbols(true);
		dataB->SetAllowAugmentingOfSequenceSymbols(true);

		NxsTreesBlock * treesB = nexusReader.GetTreesBlockTemplate();
		assert(treesB);
		treesB->SetAllowImplicitNames(true);

		nexusReader.cullIdenticalTaxaBlocks(true);

		/* End of code related to the section on making NCL less strict */


		for (int argn = 1; argn < argc; ++argn)
			{
			std::cerr << "Reading " << argv[argn] << "\n";
			try {
				nexusReader.DemoteBlocks();
				nexusReader.ReadFilepath(argv[argn], MultiFormatReader::NEXUS_FORMAT);
				}
			catch(const NxsException &x)
				{
				std::cerr << "Error:\n " << x.msg << std::endl;
				if (x.line > 0 || x.pos > 0)
					std::cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << std::endl;
				return 2;
				}
			catch(...)
				{
				nexusReader.DeleteBlocksFromFactories();
				std::cerr << "Exiting with an unknown error" << std::endl;
				return 1;
				}

			/* See discussion in "Getting information out of NCL" */

			int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
			std::cout << numTaxaBlocks << " TAXA block(s) read.\n";
			for (int i = 0; i < numTaxaBlocks; ++i)
				{
				NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
				std::string taxaBlockTitle = taxaBlock->GetTitle();
				std::cout << "Taxa block index " << i << " has the Title \"" << taxaBlockTitle << "\"\n";

				const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
				std::cout  <<  nCharBlocks << " CHARACTERS/DATA block(s) refer to this TAXA block\n";
				for (unsigned j = 0; j < nCharBlocks; ++j)
					{
					NCL_COULD_BE_CONST  NxsCharactersBlock * charBlock = nexusReader.GetCharactersBlock(taxaBlock, j);
					std::string charBlockTitle = charBlock->GetTitle();
					std::cout << "Char block index " << j << " has the Title \"" << charBlockTitle << "\"\n";

					/* call a function that uses the older 2.0 API to access data from the character block */

					std::cout << "Here is its matrix (without taxon labels):\n";
					printMatrixWithoutLabels(std::cout, charBlock);
					std::vector<const NxsDiscreteDatatypeMapper *> mappers = charBlock->GetAllDatatypeMappers();
					std::cout << "The charblock uses " << mappers.size() << " NxsDiscreteDatatypeMapper(s):\n";
					for (std::vector<const NxsDiscreteDatatypeMapper *>::const_iterator mIt = mappers.begin(); mIt != mappers.end(); ++mIt)
						{
						std::cout << "Mapper description:\n";
						(*mIt)->DebugWriteMapperFields(std::cout);
						}

					/* call a function that uses the new 2.1 API to convert the matrix to a binary encoding */
					std::vector<std::vector<char> > convertedMatrix;
					std::cout << "The matrix as binary codes:\n";
					convertToBitsForMatrix(charBlock, &convertedMatrix);
					for (unsigned taxIndex = 0; taxIndex < convertedMatrix.size(); ++taxIndex)
						{
						std::cout << taxaBlock->GetTaxonLabel(taxIndex) << " ==> ";
						const std::vector<char> & binaryCodedRow = convertedMatrix[taxIndex];
						std::vector<char>::const_iterator bcIt = binaryCodedRow.begin();
						for (; bcIt != binaryCodedRow.end(); ++bcIt)
							{
							std::cout << ' ' <<(int) *bcIt;
							}
						std::cout << '\n';
						}
					}


				const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
				std::cout  <<  nTreesBlocks << " Trees block(s) refer to this TAXA block\n";
				for (unsigned j = 0; j < nTreesBlocks; ++j)
					{
					const NxsTreesBlock * treeBlock = nexusReader.GetTreesBlock(taxaBlock, j);
					std::string treesBlockTitle = treeBlock->GetTitle();
					std::cout << "Taxa block index " << j << " has the Title \"" << treesBlockTitle << "\"\n";
					}
				}

				/* End of code related to "Getting information out of NCL" */

			}

		nexusReader.DeleteBlocksFromFactories();
	}
	catch(const NxsException &x)
		{
		std::cerr << "Error:\n " << x.msg << std::endl;
		if (x.line > 0 || x.pos > 0)
			std::cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << std::endl;
		return 2;
		}
	catch(...)
		{
		std::cerr << "Exiting with an unknown error" << std::endl;
		return 1;
		}
	return 0;
}

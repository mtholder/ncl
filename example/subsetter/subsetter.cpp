#include<string>
#include<vector>
#include<cassert>
#include "ncl/nxsmultiformat.h"







int main(int argc, char * argv[])
{
	/* Usually it is easiest to surround interactions with NCL in a try block
		that catches NxsException instances.
		Reading files can certainly generate these errors, but even queries
		after the parse can result in these exceptions.
	*/
	if (argc < 3)
		{
		std::cerr << "Expecting two arguments: a file name and a TaxSet/CharSetName\n";
		return 1;
		}
	if (argv[1][0] == '-' &&  argv[1][1] == 'h' && argv[1][2] == '\0' )
		{
		std::cerr << "Takes two arguments: a file name and a TaxSet/CharSetName\n Prints a NEXUS file that contains only the TaxSet/CharSet indicated.  Also suppresses taxa that do not have any unambiguous characters for the indicated charset.\n";
		return 0;
		}

	std::string filename(argv[1]);
	std::string setName(argv[2]);
	try {
		int blocksToRead =  (PublicNexusReader::NEXUS_TAXA_BLOCK_BIT
							| PublicNexusReader::NEXUS_CHARACTERS_BLOCK_BIT
							| PublicNexusReader::NEXUS_ASSUMPTIONS_BLOCK_BIT
							| PublicNexusReader::NEXUS_SETS_BLOCK_BIT);
		MultiFormatReader nexusReader(blocksToRead, NxsReader::WARNINGS_TO_STDERR);

		/* See discussion in "Making NCL less strict" */
		NxsCharactersBlock * charsB = nexusReader.GetCharactersBlockTemplate();
		NxsDataBlock * dataB = nexusReader.GetDataBlockTemplate();
		charsB->SetAllowAugmentingOfSequenceSymbols(true);
		dataB->SetAllowAugmentingOfSequenceSymbols(true);


		/* End of code related to the section on making NCL less strict */


		std::cerr << "Reading " << filename << "\n";
		try {
			nexusReader.DemoteBlocks();
			nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
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

		const unsigned numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
		if (numTaxaBlocks != 1)
			{
			std::cerr << "Expecting a file with exactly 1 TAXA block, but found " << numTaxaBlocks << " in the file " << filename << ".\n";
			return 2;
			}
		NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(0);
		const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
		if (nCharBlocks != 1)
			{
			std::cerr << "Expecting a file with exactly 1 CHARACTERS/DATA block, but found " << nCharBlocks << " in the file " << filename << ".\n";
			return 3;
			}
		NCL_COULD_BE_CONST  NxsCharactersBlock * charBlock = nexusReader.GetCharactersBlock(taxaBlock, 0);
		const unsigned totalNumTaxa = taxaBlock->GetNTax();
		const unsigned totalNumChars = charBlock->GetNCharTotal();
		NxsUnsignedSet taxaToInclude;
		NxsUnsignedSet charsToInclude;
		unsigned nTaxaInSet = taxaBlock->GetIndexSet(setName, &taxaToInclude);
		unsigned nCharInSet = charBlock->GetIndexSet(setName, &charsToInclude);
		if (nTaxaInSet + nCharInSet == 0)
			{
			std::cerr << "Expecting \"" << setName << "\" to refer to a TaxSet or CharSet.\n";
			return 3;
			}
		if (nTaxaInSet == 0)
			{
			for (unsigned i = 0; i < totalNumTaxa; ++i)
				taxaToInclude.insert(i);
			nTaxaInSet = taxaToInclude.size();
			}
		if (nCharInSet == 0)
			{
			for (unsigned i = 0; i < totalNumChars; ++i)
				charsToInclude.insert(i);
			nCharInSet = charsToInclude.size();
			}
		NxsUnsignedSet taxaWithUnambigData;
		NxsUnsignedSet::const_iterator taxIter = taxaToInclude.begin();
		std::vector<std::string> escapedNamesVec;
		size_t maxLabelLen = 0;
		std::vector<const NxsDiscreteDatatypeMapper *> mapperVec =  charBlock->GetAllDatatypeMappers();
		if (mapperVec.size() != 1)
			{
			std::cerr << "Expecting the characters block to contain one type of data (MIXED datatype is not supported).\n";
			return 3;
			}
		const NxsDiscreteDatatypeMapper *mapper = mapperVec[0];
		NxsDiscreteStateCell numStates = (NxsDiscreteStateCell)mapper->GetNumStates();
		for (; taxIter != taxaToInclude.end(); ++taxIter)
			{
			const NxsDiscreteStateRow & row = charBlock->GetDiscreteMatrixRow(*taxIter);
			bool hasData = false;
			for (NxsUnsignedSet::const_iterator cIt = charsToInclude.begin(); cIt != charsToInclude.end(); ++cIt)
				{
				NxsDiscreteStateCell c = row[*cIt];
				if (c >= 0 && c < numStates)
					{
					hasData = true;
					break;
					}
				}
			if (hasData)
				{
				taxaWithUnambigData.insert(*taxIter);
				NxsString lab = taxaBlock->GetTaxonLabel(*taxIter);
				std::string e = NxsString::GetEscaped(lab);
				size_t ll = e.length();
				if (ll > maxLabelLen)
					maxLabelLen = ll;
					escapedNamesVec.push_back(e);
				}
			}

		if (taxaWithUnambigData.empty())
			{
			std::cerr << "No selected taxa had unambiguous data for the selected characters!\n";
			return 5;
			}

		std::ostream & out = std::cout;

		out << "#NEXUS\nbegin TAXA;\n Dimensions ntax = " << taxaWithUnambigData.size() << "\n; TaxLabels\n ";
		std::vector<std::string>::const_iterator nIt = escapedNamesVec.begin();
		for (; nIt != escapedNamesVec.end(); ++nIt)
			{
			out << ' ' << *nIt;
			}
		out << ";\nEND;\nBegin CHARACTERS;\n Dimensions nchar = " << nCharInSet << ";\n";
		charBlock->WriteFormatCommand(out);
		out << "\n Matrix\n";

		taxIter = taxaWithUnambigData.begin();
		nIt = escapedNamesVec.begin();
		for (; nIt != escapedNamesVec.end(); ++nIt, ++taxIter)
			{
			const std::string & e = *nIt;
			std::string s(maxLabelLen + 2 - e.length(), ' ');
			out << e << s << ' ';
			const NxsDiscreteStateRow & row = charBlock->GetDiscreteMatrixRow(*taxIter);
			for (NxsUnsignedSet::const_iterator cIt = charsToInclude.begin(); cIt != charsToInclude.end(); ++cIt)
				{
				mapper->WriteStateCodeAsNexusString(out, row[*cIt]);
				}
			out << '\n';
			}
		out << ";\nEND;\n";
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

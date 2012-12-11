#include<string>
#include<vector>
#include<cassert>
#include "ncl/nxscxxdiscretematrix.h"
#include "ncl/nxsmultiformat.h"

bool gCullGapped = false;


int main(int argc, char * argv[])
{
	/* Usually it is easiest to surround interactions with NCL in a try block
		that catches NxsException instances.
		Reading files can certainly generate these errors, but even queries
		after the parse can result in these exceptions.
	*/
	std::string filename;
	for (int argi = 1; argi < argc; ++argi) 
	    {
	    if (strlen(argv[argi]) > 1 && argv[argi][0] == '-')
	        {
	        if (argv[argi][1] == 'h')
                {
                std::cerr << "Takes a arguments: The path to a NEXUS file with a single characters block.  Optional arguments: -g to cull positions with gaps or missing data\n";
                return 0;
                }
            else if (argv[argi][1] == 'g')
                gCullGapped = true;
            else
                {
                std::cerr << "Option flag -" << argv[argi][1] << " is not supported\n.";
                return 1;
                }
            }
        else if (filename.empty())
            filename.assign(argv[argi]);
        else
            {
    		std::cerr << "Expecting one file name as an argument.\n";
            return 1;
            }
        }
            

	if (filename.empty())
		{
		std::cerr << "Expecting one arguments: a file name\n";
		return 1;
		}

	try {
		int blocksToRead =  (PublicNexusReader::NEXUS_TAXA_BLOCK_BIT
							| PublicNexusReader::NEXUS_CHARACTERS_BLOCK_BIT
							| PublicNexusReader::NEXUS_ASSUMPTIONS_BLOCK_BIT
							| PublicNexusReader::NEXUS_SETS_BLOCK_BIT
							);
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
		
		std::vector<const NxsDiscreteDatatypeMapper *> mappers = charBlock->GetAllDatatypeMappers();
		if (mappers.size() != 1)
			{
			std::cerr << "Expecting an unmixed characters block, but found a matrix with datatype = mixed or a datatype with augmented symbols\n";
			return 4;
			}
		
		const NxsDiscreteDatatypeMapper * dm = mappers[0];
		
		ScopedTwoDMatrix<NxsCDiscreteStateSet> compressedMatrix;
        std::vector<unsigned> patternCounts;
        std::vector<double> patternWeights;
        bool hasWeights = true;
        bool hasIntWeights = true;
        
        
        std::vector<NxsCharacterPattern> compressedTransposedMatrix;
        if (true)
            {
            std::vector<int> originalIndexToCompressed;
            std::vector<std::set<unsigned> > compressedIndexToOriginal;
            if (true) 
                {
                NxsCXXDiscreteMatrix cxxMat(*charBlock, false, 0L, false);

                hasWeights = cxxMat.hasWeights();
                hasIntWeights = cxxMat.hasIntWeights();
                NxsCompressDiscreteMatrix(cxxMat, compressedTransposedMatrix, &originalIndexToCompressed, &compressedIndexToOriginal);
                }
           std::vector<double> * wtsPtr = (hasWeights ? &patternWeights : 0L);   
           NxsTransposeCompressedMatrix(compressedTransposedMatrix, compressedMatrix, &patternCounts, wtsPtr);
           }
		
		NxsCDiscreteStateSet ** matrixAlias = compressedMatrix.GetAlias();
		const unsigned ntaxTotal =  charBlock->GetNTaxTotal();
		const unsigned numPatterns = patternCounts.size();
		std::set<unsigned> culledSet;
		if (gCullGapped)
		    {
		    unsigned ind = 0;
		    for (std::vector<NxsCharacterPattern>::const_iterator ctmIt = compressedTransposedMatrix.begin();
		                                                          ctmIt != compressedTransposedMatrix.end(); 
		                                                          ++ctmIt, ++ind)
                {
                if (!ctmIt->StateCodesAreNonNegative())
                    culledSet.insert(ind);
                }
		    }
		const unsigned numCulled = culledSet.size();
		
		
		std::cout << "#NEXUS\nBEGIN DATA;\n\tDimensions ntax = " << ntaxTotal << " nchar = " << numPatterns - numCulled << ";\n\t";
		charBlock->WriteFormatCommand(std::cout);
		std::cout << "Matrix\n";
		const unsigned width = taxaBlock->GetMaxTaxonLabelLength();
        for (unsigned i = 0; i < ntaxTotal; i++)
			{
            const std::string currTaxonLabel = NxsString::GetEscaped(taxaBlock->GetTaxonLabel(i));
            std::cout << currTaxonLabel;
            unsigned currTaxonLabelLen = (unsigned)currTaxonLabel.size();
            unsigned diff = width - currTaxonLabelLen;
            for (unsigned k = 0; k < diff+5; k++)
                std::cout << ' ';
            NxsCDiscreteStateSet * matrixRow = matrixAlias[i];
            /*
            for (unsigned j = 0; j < numPatterns; ++j)
                std::cout << ' ' << (int) matrixRow[j] ;
            std::cout << '\n';
            */
            for (unsigned j = 0; j < numPatterns; ++j)
                {
                if (culledSet.find(j) == culledSet.end()) 
                    dm->WriteStateCodeAsNexusString(std::cout, matrixRow[j], true);
                }
            /*
            for (unsigned j = 0; j < numPatterns; ++j)
                std::cout << ' ' << (int) matrixRow[j];
            */
            std::cout << '\n';
            }
        const char * sp = (hasWeights ? " " : " * ");

        unsigned ind = 0;
        std::cout << ";\nEND;\nBEGIN ASSUMPTIONS;\n\tWTSET" << sp << " counts ( vector ) =";
        for (std::vector<unsigned>::const_iterator cIt = patternCounts.begin(); cIt != patternCounts.end(); ++cIt, ++ind)
            {
            if (culledSet.find(ind) == culledSet.end()) 
                std::cout << ' ' << *cIt;
            }
        std::cout << ";\n";
        if (hasWeights)
            {
            std::cout << "\tWTSET * sum_of_weights ( vector ) =";
            if (hasIntWeights)
                {
                for (std::vector<double>::const_iterator cIt = patternWeights.begin(); cIt != patternWeights.end(); ++cIt, ++ind)
                    {
                    int w = int(0.01 + *cIt);
                    if (culledSet.find(ind) == culledSet.end()) 
                        std::cout << ' ' << w;
                    }

                }
            else
                for (std::vector<double>::const_iterator cIt = patternWeights.begin(); cIt != patternWeights.end(); ++cIt, ++ind)
                    {
                    if (culledSet.find(ind) == culledSet.end()) 
                        std::cout << ' ' << *cIt;
                    }
            std::cout << ";\n";
            }
        std::cout << "END;\n";
       
        

		
		
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

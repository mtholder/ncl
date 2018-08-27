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

/*******************************************************************************
 * This file contains code for 4 executables:
 *		NEXUSnormalizer, NEXUSvalidator, NEXUSinspector, and NEX_us2ml
 *	with conditional compilation used to determine the behavior.
 *
 *		* NEXUSnormalizer - writes a NEXUS version of the file with consistent
 *			ordering of blocks and commands. Ideally 2 equivalent files will
 *			produce the same normalized output. This version of tthe program is
 *			less ambitious. The goal is to be able to run (for any valid NEXUS
 *			in.nex file):
 *				$ NEXUSnormalizer in.nex > outOrig.nex
 *				$ NEXUSnormalizer outOrig.nex > outSecond.nex
 *				$ diff outOrig.nex outSecond.nex
 *			and find no differences.
 *		* NEXUSvalidator - reports errors and warnings to stderr. Invalid files
 *			cause exit with a non-zero return code
 *		* NEXUSinspector - writes a brief report of every block parsed
 *		* NEXUS_us2ml - writes a nexml version of the input (partially
 *			implemented, note that the code to write nexml is in us2ml.cpp).
 * See the processFilepath() function for an example of how to deal with NCL
 *	to read a file using the new MultiFormatReader class. When the file
 *	is correctly read, the processContent() function is called.
 *
 * All other code has to do with reading command line arguments and other
 * 	user-interface concerns.
 */
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"
#include "normalizer.h"
#include <cassert>
using namespace std;
//#include "ncl/nxscxxdiscretematrix.h"

#if defined(TO_NEXML_CONVERTER) && TO_NEXML_CONVERTER
	void	writeAsNexml(PublicNexusReader & nexusReader, ostream & os, TranslatingConventions & guidTag);
#endif

bool gQuietMode = false;
std::ofstream gCommonFileStream;
std::ostream * gCommonOstream = 0L;
#if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
	enum ExportFormatEnum
		{
		NEXUS_EXPORT_FORMAT,
		PHYLIP_EXPORT_FORMAT,
		RELAXED_PHYLIP_EXPORT_FORMAT,
		FASTA_EXPORT_FORMAT,
		NEXML_EXPORT_FORMAT,
		TREE_NEXML_EXPORT_FORMAT,
		UNSUPPORTED_EXPORT_FORMAT
		};
	ExportFormatEnum gExportFormat = NEXML_EXPORT_FORMAT;
	std::string gExportPrefix("out");
	ExportFormatEnum readExportFormatName(const std::string &);
	void exportData(PublicNexusReader & nexusReader, MultiFormatReader::DataFormatType f, long interleavLen, std::string prefix, std::ostream *, TranslatingConventions &);


	ExportFormatEnum readExportFormatName(const std::string & s)
	{
		const char * gExportFormatNames[] = {   "nexus",
												"phylip",
												"relaxedphylip",
												"fasta",
												"nexml",
												"treenexml"
												};
		const unsigned gNumExportFormats = 6;
	
		NxsString l(s.c_str());
		NxsString::to_lower(l);
		int ind = NxsString::index_in_array(l, gExportFormatNames, gNumExportFormats);
		if (ind < 0)
			return UNSUPPORTED_EXPORT_FORMAT;
		return ExportFormatEnum(ind);
	}
	std::string gNEXUSSafeNamesToWrite;
	std::string gNEXUSSafeNamesToRead;
	bool gProduceEvenTrivalTranslation = false;
	
	void reverseTranslateNames(PublicNexusReader & reader, std::string filepath);
	void substituteSafeTaxaLabels(PublicNexusReader & reader, std::string filepath, bool evenTrivial);

	
#endif  // if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP

bool gAltNexus = false;


long gStrictLevel = 2;
bool gUnderscoresToSpaces = false;
bool gValidateInternals = true;
bool gTreesViaInMemoryStruct = true;
long gInterleaveLen = -1;
bool blocksReadInValidation = false;
bool gSuppressingNameTranslationFile = false;
bool gAllowNumericInterpretationOfTaxLabels = true;
TranslatingConventions gTranslatingConventions;

enum ProcessActionsEnum
	{
	REPORT_BLOCKS,
	OUTPUT_NORMALIZED_NEXUS,
	OUTPUT_ANY_FORMAT,
	OUTPUT_NEXML,
	VALIDATE_ONLY
	};


void processContent(PublicNexusReader & nexusReader, ostream *os, ProcessActionsEnum currentAction);
MultiFormatReader * instantiateReader();

#	if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
	MultiFormatReader * gNexusReader = NULL;
#	endif


void reportNexusStats(const PublicNexusReader & nexusReader, ostream *os)
{
	if (!os)
		return;

	const unsigned nTaxaBlocks = nexusReader.GetNumTaxaBlocks();
	*os <<  nTaxaBlocks << " taxa block(s) read.\n";
	for (unsigned t = 0; t < nTaxaBlocks; ++t)
		{
		NxsTaxaBlock * tb = nexusReader.GetTaxaBlock(t);
		*os << "Taxa block #" << t + 1 << ".\n";
		tb->Report(*os);
		const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(tb);
		*os <<  nCharBlocks << " Characters/Data block(s) read that link to this Taxa block.\n";
		for (unsigned i = 0; i < nCharBlocks; ++i)
			{
			NxsCharactersBlock * cb = nexusReader.GetCharactersBlock(tb, i);

			//NxsCXXDiscreteMatrix mat(*cb, true);

			*os << "Character block #" << i + 1 << " for this Taxa block.\n";
			cb->Report(*os);
			const unsigned nAssumpBlocks = nexusReader.GetNumAssumptionsBlocks(cb);
			*os <<  nAssumpBlocks << " Assumptions block(s) read that link to this Characters block.\n";
			for (unsigned j= 0; j < nAssumpBlocks; ++j)
				{
				NxsAssumptionsBlock * ab = nexusReader.GetAssumptionsBlock(cb, j);
				*os << "Assumptions block #" << j + 1 << " for this Characters block.\n";
				ab->Report(*os);
				}
			}
		const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(tb);
		*os <<  nTreesBlocks << " Trees/Data block(s) read that link to this Taxa block.\n";
		for (unsigned i = 0; i < nTreesBlocks; ++i)
			{
			NxsTreesBlock * cb = nexusReader.GetTreesBlock(tb, i);
			*os << "Trees block #" << i + 1 << " for this Taxa block.\n";
			cb->Report(*os);
			const unsigned nAssumpBlocks = nexusReader.GetNumAssumptionsBlocks(cb);
			*os <<  nAssumpBlocks << " Assumptions block(s) read that link to this Trees block.\n";
			for (unsigned j= 0; j < nAssumpBlocks; ++j)
				{
				NxsAssumptionsBlock * ab = nexusReader.GetAssumptionsBlock(cb, j);
				*os << "Assumptions block #" << j + 1 << " for this Trees block.\n";
				ab->Report(*os);
				}
			}
		const unsigned nAssumpBlocks = nexusReader.GetNumAssumptionsBlocks(tb);
		*os <<  nAssumpBlocks << " Assumptions block(s) read that link to this Taxa block.\n";
		for (unsigned j= 0; j < nAssumpBlocks; ++j)
			{
			NxsAssumptionsBlock * ab = nexusReader.GetAssumptionsBlock(tb, j);
			*os << "Assumptions block #" << j + 1 << " for this Taxa block.\n";
			ab->Report(*os);
			}
		const unsigned nDistancesBlocks = nexusReader.GetNumDistancesBlocks(tb);
		*os <<  nDistancesBlocks << " Distances block(s) read that link to this Taxa block.\n";
		for (unsigned j= 0; j < nDistancesBlocks; ++j)
			{
			NxsDistancesBlock * ab = nexusReader.GetDistancesBlock(tb, j);
			*os << "Distances block #" << j + 1 << " for this Taxa block.\n";
			ab->Report(*os);
			}
		const unsigned nUnalignedBlocks = nexusReader.GetNumUnalignedBlocks(tb);
		*os <<  nUnalignedBlocks << " Unaligned block(s) read that link to this Taxa block.\n";
		for (unsigned j= 0; j < nUnalignedBlocks; ++j)
			{
			NxsUnalignedBlock * ab = nexusReader.GetUnalignedBlock(tb, j);
			*os << "Unaligned block #" << j + 1 << " for this Taxa block.\n";
			ab->Report(*os);
			}
		*os << "\n\n";
		}
	const unsigned nUnknown = nexusReader.GetNumUnknownBlocks();
	*os <<  nUnknown << " private block(s) read.\n";
	for (unsigned t = 0; t < nUnknown; ++t)
		{
		NxsStoreTokensBlockReader * ub = nexusReader.GetUnknownBlock(t);
		*os << "Private block #" << t + 1 << " is a " << ub->GetID() << " block.\n";
		}
}


void writeAsNexus(PublicNexusReader & nexusReader, ostream & os)
{
	BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();
	os << "#NEXUS\n";
	for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt)
		{
		NxsBlock * b = *bIt;
		if (b)
			b->WriteAsNexus(os);
		}
}



////////////////////////////////////////////////////////////////////////////////
// Takes NxsReader that has successfully read a file, and processes the
//	information stored in the reader.
//
// The caller is responsibel for calling DeleteBlocksFromFactories() to clean
//	up (if the reader uses the factory API).
////////////////////////////////////////////////////////////////////////////////
void processContent(PublicNexusReader & nexusReader, ostream *os, ProcessActionsEnum currentAction)
	{
	BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();
	if (blocks.size() == 0) {
		cerr << "Error:\n No understandable content was found.\n";
		exit(1);
	}
	if (currentAction == REPORT_BLOCKS)
		reportNexusStats(nexusReader, os);
	else if (OUTPUT_NORMALIZED_NEXUS == currentAction && os)
		{
		writeAsNexus(nexusReader, *os);
		}
	else if (OUTPUT_NEXML == currentAction && os)
		{
#		if defined(TO_NEXML_CONVERTER) && TO_NEXML_CONVERTER
			writeAsNexml(nexusReader, *os, gTranslatingConventions);
#		else
			cerr << "Error nexml conversion not implemented\n";
			exit(1);
#		endif
		}
	else if (OUTPUT_ANY_FORMAT == currentAction && os)
		{
#		if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
			if (!gNEXUSSafeNamesToRead.empty()) {
				reverseTranslateNames(nexusReader, gNEXUSSafeNamesToRead);
			}
			else if (!gNEXUSSafeNamesToWrite.empty()) {
				substituteSafeTaxaLabels(nexusReader, gNEXUSSafeNamesToWrite, gProduceEvenTrivalTranslation);
			}


			std::string fullExportPrefix;
			MultiFormatReader::DataFormatType f = MultiFormatReader::NEXUS_FORMAT;
			if (gExportFormat == NEXUS_EXPORT_FORMAT)
				exportData(nexusReader, MultiFormatReader::NEXUS_FORMAT, gInterleaveLen, gExportPrefix, gCommonOstream, gTranslatingConventions);
			else if (gExportFormat == NEXML_EXPORT_FORMAT || gExportFormat == TREE_NEXML_EXPORT_FORMAT ) {
				if (gExportFormat == TREE_NEXML_EXPORT_FORMAT) {
					gTranslatingConventions.emitTreesAndTaxaOnly = true;
				}
				exportData(nexusReader, MultiFormatReader::NEXML_FORMAT, gInterleaveLen, gExportPrefix, gCommonOstream, gTranslatingConventions);
			}
			else if (gExportFormat == PHYLIP_EXPORT_FORMAT) {
				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".dna");
				f = (gInterleaveLen < 0 ? MultiFormatReader::PHYLIP_DNA_FORMAT : MultiFormatReader::INTERLEAVED_PHYLIP_DNA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".rna");
				f = (gInterleaveLen < 0 ? MultiFormatReader::PHYLIP_RNA_FORMAT : MultiFormatReader::INTERLEAVED_PHYLIP_RNA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".aa");
				f = (gInterleaveLen < 0 ? MultiFormatReader::PHYLIP_AA_FORMAT : MultiFormatReader::INTERLEAVED_PHYLIP_AA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				//fullExportPrefix = gExportPrefix;
				//fullExportPrefix.append(".discrete");
				//f = (gInterleaveLen < 0 ? MultiFormatReader::PHYLIP_DISC_FORMAT : MultiFormatReader::INTERLEAVED_PHYLIP_DISC_FORMAT);
				//exportData(nexusReader, f, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				fullExportPrefix = gExportPrefix;
				exportData(nexusReader, MultiFormatReader::PHYLIP_TREE_FORMAT, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);
			}
			else if (gExportFormat == RELAXED_PHYLIP_EXPORT_FORMAT) {
				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".dna");
				f = (gInterleaveLen < 0 ? MultiFormatReader::RELAXED_PHYLIP_DNA_FORMAT : MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".rna");
				f = (gInterleaveLen < 0 ? MultiFormatReader::RELAXED_PHYLIP_RNA_FORMAT : MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".aa");
				f = (gInterleaveLen < 0 ? MultiFormatReader::RELAXED_PHYLIP_AA_FORMAT : MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				//fullExportPrefix = gExportPrefix;
				//fullExportPrefix.append(".discrete");
				//f = (gInterleaveLen < 0 ? MultiFormatReader::RELAXED_PHYLIP_DISC_FORMAT : MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DISC_FORMAT);
				//exportData(nexusReader, f, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				fullExportPrefix = gExportPrefix;
				exportData(nexusReader, MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);
			}
			else if (gExportFormat == FASTA_EXPORT_FORMAT) {
				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".dna");
				exportData(nexusReader,  MultiFormatReader::FASTA_DNA_FORMAT, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".rna");
				exportData(nexusReader,  MultiFormatReader::FASTA_RNA_FORMAT, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".aa");
				exportData(nexusReader,  MultiFormatReader::FASTA_AA_FORMAT, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);

				fullExportPrefix = gExportPrefix;
				exportData(nexusReader, MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT, gInterleaveLen, fullExportPrefix, gCommonOstream, gTranslatingConventions);
			}
			else {
				cerr << "Unsupported export format requested.\n";
				exit(1);
			}
#		else
			cerr << "Exporting to any format is implemented by this executable.\n";
			exit(1);
#		endif
		}
	else if (VALIDATE_ONLY == currentAction)
		{
		if (!blocks.empty())
			blocksReadInValidation = true;
		}
	}


MultiFormatReader * instantiateReader()
{
	MultiFormatReader * nexusReader = new MultiFormatReader(-1, NxsReader::WARNINGS_TO_STDERR);
	if (gQuietMode)
		nexusReader->SetWarningOutputLevel(NxsReader::SKIPPING_CONTENT_WARNING);
	if (gStrictLevel != 2)
		nexusReader->SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
	if (gUnderscoresToSpaces) 
		nexusReader->SetCoerceUnderscoresToSpaces(true);
	NxsCharactersBlock * charsB = nexusReader->GetCharactersBlockTemplate();
	NxsDataBlock * dataB = nexusReader->GetDataBlockTemplate();
	charsB->SetAllowAugmentingOfSequenceSymbols(true);
	dataB->SetAllowAugmentingOfSequenceSymbols(true);
	if (gInterleaveLen > 0)
		{
		assert(charsB);
		charsB->SetWriteInterleaveLen(gInterleaveLen);
		dataB->SetWriteInterleaveLen(gInterleaveLen);
		}

	NxsTreesBlock * treesB = nexusReader->GetTreesBlockTemplate();
	assert(treesB);
	if (gStrictLevel < 2)
		treesB->SetAllowImplicitNames(true);
	treesB->SetWriteFromNodeEdgeDataStructure(gTreesViaInMemoryStruct);
	treesB->setValidateInternalNodeLabels(gValidateInternals);
	if (!gValidateInternals) {
		gTranslatingConventions.treatNodeLabelsAsStrings = true;
	}
	treesB->setAllowNumericInterpretationOfTaxLabels(gAllowNumericInterpretationOfTaxLabels);
	if (gAltNexus)
		treesB->setWriteTranslateTable(false);
	if (gStrictLevel < 2)
		{
		NxsStoreTokensBlockReader *storerB =  nexusReader->GetUnknownBlockTemplate();
		assert(storerB);
		storerB->SetTolerateEOFInBlock(true);
		}
	nexusReader->conversionOutputRecord.addNumbersToDisambiguateNames = true;
	
	if (gSuppressingNameTranslationFile)
		nexusReader->conversionOutputRecord.writeNameTranslationFile = false;
	return nexusReader;
}

////////////////////////////////////////////////////////////////////////////////
// Creates a NxsReader, and tries to read the file `filename`.  If the
//	read succeeds, then processContent will be called.
//	\returns 0 on success
////////////////////////////////////////////////////////////////////////////////
#if !defined(MULTIFILE_NEXUS_READER) ||  !MULTIFILE_NEXUS_READER
	int processFilepath(
		const char * filename, // name of the file to be read
		ostream * os, // output stream to use (NULL for no output). Not that cerr is used to report errors.
		MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
		ProcessActionsEnum currentAction) // enum that is passed on to processContent to indicate what should be done with the content of the file.
#else
	int processFilepath(
		const char * filename, // name of the file to be read
		ostream * , // output stream to use (NULL for no output). Not that cerr is used to report errors.
		MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
		ProcessActionsEnum ) // enum that is passed on to processContent to indicate what should be done with the content of the file.
#endif 
	{
	assert(filename);
	try
		{
		MultiFormatReader * nexusReader;
#		if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
			nexusReader = gNexusReader;
#		else
			nexusReader = instantiateReader();
#		endif

		if (!gQuietMode)
			cerr << "Executing" << endl;
		try {
			nexusReader->DemoteBlocks();
			nexusReader->ReadFilepath(filename, fmt);
#			if !defined(MULTIFILE_NEXUS_READER) ||  !MULTIFILE_NEXUS_READER
				processContent(*nexusReader, os, currentAction);
#			endif
			}
		catch(...)
			{
			nexusReader->DeleteBlocksFromFactories();
#			if ! defined(MULTIFILE_NEXUS_READER) || !MULTIFILE_NEXUS_READER
				delete nexusReader;
#			endif
			throw;
			}
#		if !defined(MULTIFILE_NEXUS_READER) ||  !MULTIFILE_NEXUS_READER
			nexusReader->DeleteBlocksFromFactories();
			delete nexusReader;
#		endif
		return 0;
		}
	catch (const NxsException &x)
		{
		cerr << "Error:\n " << x.msg << endl;
		if (x.line > 0 || x.pos > 0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		return 2;
		}
	}

/*! \returns 0 on success*/
int readFilepathAsNEXUS(const char *filename, MultiFormatReader::DataFormatType fmt, ProcessActionsEnum currentAction)
	{
	if (!gQuietMode)
		cerr << "[Reading " << filename << "	 ]" << endl;
	try {
		ostream * outStream = 0L;
		if (currentAction != VALIDATE_ONLY)
			outStream = &cout;
		return processFilepath(filename, outStream, fmt, currentAction);

		}
	catch (...)
		{
		cerr << "Normalizing of " << filename << " failed (with an exception)" << endl;
		return 1;
		}
	}

/*! \returns 0 on success*/
int readFilesListedIsFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt, ProcessActionsEnum currentAction)
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
			{
			int rc = readFilepathAsNEXUS(filename, fmt, currentAction);
			if (rc != 0)
				return rc;
			}
		}
	return 0;
	}

#if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
	const char * gExeName = "NEXUSvalidator";
#elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
	const char * gExeName = "NEXUSinspector";
#elif defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
	const char * gExeName = "NCLconverter";
#else
#	if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
		const char * gExeName = "NEXUSunion";
#	else
		const char * gExeName = "NEXUSnormalizer";
#	endif
# endif
	
void printHelp(ostream & out)
	{
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		out << "NEXUSvalidator takes reads a file and exits with a success (return code 0) if the file is valid.\n";
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		out << "NEXUSinspector takes reads a file and writes a report of the content to standard out.\n";
#	elif defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
		out << "NCLconverter takes reads a file and writes a report of the content to a file prefix (specified with the -o flag) in the chosen output format (specified with the -e flag).\n";
#	else
#		if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
			out << "NEXUSunion reads a series of NEXUS file and writes the union of all of their content to standard out (using the NEXUSnormalizer conventions of indentation and syntax).\n";
#		else
			out << "NEXUSnormalizer takes reads a file and rewrites the file to standard out with consistent indentation and syntax.\n";
#		endif
# 	endif
	out << "\nThe most common usage is simply:\n    " << gExeName << " <path to NEXUS file>\n";
	out << "\nCommand-line flags:\n\n";
#if !defined(JUST_VALIDATE_NEXUS) && !defined(JUST_REPORT_NEXUS) && !defined(TO_NEXML_CONVERTER)
	out << "    -a AltNexus output (no translation table in trees)\n\n";
#endif
#if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
	out << "    -d<fn> specifies the single output destination. Or you can use -d- to indicate that\n";
	out << "             output should be directed to standard output.Warning use of this option may result\n";
	out << "             in an invalid output due to concatenation of separate \"blocks\" of information\n";
	out << "             into a single file!  \n";
	out << "    -e<format> specifies the output file format expected:\n";
	out << "            -enexus  \"normalized\" NEXUS output\n";
	out << "            -efasta  Character data in fasta (could result in multiple output files)\n";
	out << "            -ephylip  Trees and character data in phylip (could result in multiple output files)\n";
	out << "            -erelaxedphylip  Trees and character data in relaxed phylip (could result in multiple output files)\n";
	out << "            -enexml  nexml output (this is also the default)\n";
	out << "            -etreenexml  nexml output of trees and taxa only\n";
#endif
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
#if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
	out << "    -g     Global ids (otu1, otu2 ...) in NeXML rather than IDs\n";
	out << "             which are composed by appending contained element to\n";
	out << "             their parents. Instead of g0n0 for the first tree in the first TREES block\n";
	out << "             and g1n0 for the first tree in the second trees block, with the\n";
	out << "             -g option these would be tree0, tree1... (sequentially numbered regardless\n";
	out << "             of the block.";
#endif
	out << "    -h help. on the command line shows this help message\n\n";
#if !defined(JUST_VALIDATE_NEXUS) && !defined(JUST_REPORT_NEXUS) && !defined(TO_NEXML_CONVERTER)
	out << "    -i<number> specifies the length of the interleaved pages to create\n";
#endif
#if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
	out << "    -j     Suppress the creation of a NameTranslationFile\n";
#endif
#if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
	out << "    -o<fn> specifies the output prefix.  An appropriate suffix and extension are added\n";
	out << "    -pe# the index of the first edge in global Id NeXML output mode\n";
	out << "    -pn# the index of the first node in global Id NeXML output mode\n";
	out << "    -po# the index of the first otu in global Id NeXML output mode\n";
	out << "    -pO# the index of the first otus in global Id NeXML output mode\n";
	out << "    -pt# the index of the first tree in global Id NeXML output mode\n";
	out << "    -pT# the index of the first trees in global Id NeXML output mode\n";
#endif
	out << "    -q quiet. suppress NCL status messages while reading files\n\n";
	out << "    -l<path> reads a file and treats each line of the file as a path to NEXUS file\n\n";
	out << "    -s<non-negative integer> controls the NEXUS strictness level.\n";
	out << "        the default level is equivalent to -s2 invoking the program with \n";
	out << "        -s3 or a higher number will convert some warnings into fatal errors.\n";
	out << "        Running with -s1 will cause the parser to accept dangerous constructs,\n";
	out << "        and running with -s0 will cause the parser make every attempt to finish\n";
	out << "        parsing the file (warning about very serious errors).\n\n";
	out << "        Note that when -s0 strictness level is used, and the parser fails to\n";
	out << "        finish, it will often be the result of an earlier error than the \n";
	out << "        error that is reported in the last message.\n";
#if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
	out << "    -t<tag>     tag used as a prefix for nexmlids.\n";
	out << "    -u     converts underscores to spaces in formats other than NEXUS.\n";
#endif
	out << "    -x do NOT validate internal labels in trees as taxa labels\n\n";
	out << "    -X do NOT treat numbers in trees as taxon numbers, treat them as arbitrary\n        labels (should not be used with NEXUS files).\n\n";
#if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
	out << "    -y<filename> translate to \"safe\" taxon names and store the new names as a NEXUS.\n";
	out << "             file called <filename> with a TaxaAssociation block. The first taxa block\n";
	out << "             in the association block will hold the original names, and the second will\n";
	out << "             hold the \"safe\" names\n";
	out << "    -Y<filename> behaves like -y, except with -Y a translation file will be produced even";
	out << "             if the original names were already \"safe\"\n";
	out << "    -z<filename> use the NEXUS-formatted file called <filename> with a TaxaAssociation block\n";
	out << "             to restore original names.  Assumes that the first taxa block in the TaxaAssociation\n";
	out << "             block holds the original name and the second is the current name. This function\n";
	out << "             is useful for \"undoing\" the effects of the -y option.\n";
#endif
	}

int do_main(int argc, char *argv[])
	{
	NxsReader::setNCLCatchesSignals(true);
	MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		ProcessActionsEnum currentAction = VALIDATE_ONLY;
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		ProcessActionsEnum currentAction = REPORT_BLOCKS;
#	elif defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
		ProcessActionsEnum currentAction = OUTPUT_ANY_FORMAT;
#	elif defined(TO_NEXML_CONVERTER) && TO_NEXML_CONVERTER
		ProcessActionsEnum currentAction = OUTPUT_NEXML;
#	else
		ProcessActionsEnum currentAction = OUTPUT_NORMALIZED_NEXUS;
# 	endif

	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (slen < 2 || filepath[0] != '-')
			continue;
		if (filepath[1] == 'h')
			{
			printHelp(cout);
			return 1;
			}
		else if (filepath[1] == 'q')
			gQuietMode = true;
		else if (filepath[1] == 'x')
			gValidateInternals = false;
		else if (filepath[1] == 'X')
			gAllowNumericInterpretationOfTaxLabels = false;
		else if (filepath[1] == 'u')
			gUnderscoresToSpaces = true;
		else if (filepath[1] == 'j')
			gSuppressingNameTranslationFile = true;
		else if (filepath[1] == 's')
			{
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gStrictLevel)))
				{
				cerr << "Expecting an integer after -s\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		//pass
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		//pass
#	elif defined(TO_NEXML_CONVERTER) && TO_NEXML_CONVERTER
		//pass
#	else
		else if (filepath[1] == 'i')
			{
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gInterleaveLen)) || gInterleaveLen < 1)
				{
				cerr << "Expecting a positive integer after -i\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		else if (filepath[1] == 'a')
			{
			if ((slen != 2))
				{
				cerr << "Not expecting a value after -a\n" << endl;
				printHelp(cerr);
				return 2;
				}
			gAltNexus = true;
			}
#	endif
#	if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
		else if (filepath[1] == 'e') {
			if (slen > 2)
				{
				std::string efmtName(filepath + 2, slen - 2);
				gExportFormat = readExportFormatName(efmtName);
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting a format after -e\n" << endl;
				printHelp(cerr);
				return 2;
				}
		}
		else if (filepath[1] == 'p') {
			long x = -1;
			if ((slen <= 3) || (!NxsString::to_long(filepath + 3, &x)) || x < 0)
				{
				cerr << "Expecting a positive integer after -p* flag (where * is a placeholder).\n" << endl;
				printHelp(cerr);
				return 2;
				}
			char code = filepath[2];
			if (code == 'e')
				gTranslatingConventions.currentEdgeIndex = (unsigned)x;
			else if (code == 'n')
				gTranslatingConventions.currentNodeIndex = (unsigned)x;
			else if (code == 'o')
				gTranslatingConventions.currentOTUIndex = (unsigned)x;
			else if (code == 'O')
				gTranslatingConventions.currentOTUsIndex = (unsigned)x;
			else if (code == 't')
				gTranslatingConventions.currentTreeIndex = (unsigned)x;
			else if (code == 'T')
				gTranslatingConventions.currentTreesIndex = (unsigned)x;
			else {
				cerr << "Expecting e, n, o, O, t, or T after -p.\n" << endl;
				printHelp(cerr);
				return 2;
				}
		}
		else if (filepath[1] == 'g') {
			gTranslatingConventions.globalIncrementingIDs = true;
		}
		else if (filepath[1] == 'o') {
			if (slen > 2)
				{
				std::string oname(filepath + 2, slen - 2);
				gExportPrefix = oname;
				}
			if (gExportPrefix.empty())
				{
				cerr << "Expecting an output file prefix after -o\n" << endl;
				printHelp(cerr);
				return 2;
				}
		}
		else if (filepath[1] == 'd') {
			if (slen > 2)
				{
				if (gCommonOstream != 0)
					{
					cerr << "Expecting only one -d flag per invocation.\n" << endl;
					return 2;
					}
				std::string dname(filepath + 2, slen - 2);
				if (dname == "-")
					gCommonOstream = &std::cout;
				else
					{
					gCommonFileStream.open(dname.c_str());
					if (!gCommonFileStream.good())
						{
						cerr << "Error opening " << dname << ".\n" << flush;
						return 2;
						}
					if (!gQuietMode) {
						cerr << "Directing all output to " << dname << ".\n" << flush;
					}
					gCommonOstream = & gCommonFileStream;
					}
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting an output file prefix after -o\n" << endl;
				printHelp(cerr);
				return 2;
				}
		}
		else if (filepath[1] == 'y' || filepath[1] == 'Y') {
			if (slen > 2)
				{
				std::string oname(filepath + 2, slen - 2);
				if (!gNEXUSSafeNamesToRead.empty()) 
					{
					cerr << "The -y and -z options cannot both be specified!\n" << endl;
					printHelp(cerr);
					return 2;
					}
				if (!gNEXUSSafeNamesToWrite.empty()) 
					{
					cerr << "The -y only be specified once per invocation!\n" << endl;
					printHelp(cerr);
					return 2;
					}
				gNEXUSSafeNamesToWrite = oname;
				if (gNEXUSSafeNamesToWrite.empty())
					{
					cerr << "Expecting an output file prefix after -o\n" << endl;
					printHelp(cerr);
					return 2;
					}
				}
				if (filepath[1] == 'Y')
					gProduceEvenTrivalTranslation = true;
		}
		else if (filepath[1] == 'z') {
			if (slen > 2)
				{
				std::string oname(filepath + 2, slen - 2);
				if (!gNEXUSSafeNamesToWrite.empty()) 
					{
					cerr << "The -y and -z options cannot both be specified!\n" << endl;
					printHelp(cerr);
					return 2;
					}
				if (!gNEXUSSafeNamesToRead.empty()) 
					{
					cerr << "The -z only be specified once per invocation!\n" << endl;
					printHelp(cerr);
					return 2;
					}
				gNEXUSSafeNamesToRead = oname;
				if (gNEXUSSafeNamesToRead.empty())
					{
					cerr << "Expecting gNEXUSSafeNamesToRead after -z\n" << endl;
					printHelp(cerr);
					return 2;
					}
				}
		}
		else if (filepath[1] == 't') {
			if (slen > 2)
				{
				std::string guidTag(filepath + 2, slen - 2);
				gTranslatingConventions.idPrefix = guidTag;
				if (gTranslatingConventions.idPrefix.empty())
					{
					cerr << "Expecting an ID prefix after -t\n" << endl;
					printHelp(cerr);
					return 2;
					}
				}
		}
#	endif
		else if (filepath[1] == 'f')
			{
			f = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2)
				{
				std::string fmtName(filepath + 2, slen - 2);
				f =  MultiFormatReader::formatNameToCode(fmtName);
				if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
					{
					cerr << "Unknown format \"" << fmtName << "\" after -f\n" << endl;
					printHelp(cerr);
					return 3;
					}
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting a format after -f\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		}
#	if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
		gNexusReader = instantiateReader();
		gNexusReader->cullIdenticalTaxaBlocks(true);
#	endif
	bool readfile = false;
	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (slen < 1)
			continue;
		if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 'l')
			{
			readfile = true;
			int rc = readFilesListedIsFile(filepath+2, f, currentAction);
			if (rc != 0)
				return rc;
			}
		else if (filepath[0] != '-')
			{
			readfile = true;
			int rc = readFilepathAsNEXUS(filepath, f, currentAction);
			if (rc != 0)
				return rc;
			}
		}
#	if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
		if (gNexusReader)
			{
			processContent(*gNexusReader, &std::cout, OUTPUT_NORMALIZED_NEXUS);
			gNexusReader->DeleteBlocksFromFactories();
			delete gNexusReader;
			}
#	endif

	if (!readfile)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
		}
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		if (blocksReadInValidation)
			return  0;
		std::cerr << "No blocks read\n";
		return 1;
#	else
		return 0;
#	endif
	}

int main(int argc, char *argv[])
	{
	int rc = do_main(argc, argv);
	if (gCommonOstream != 0L && gCommonOstream == &gCommonFileStream)
		gCommonFileStream.close();
	return rc;
	}

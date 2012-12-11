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
#include <cassert>

#include "splitsstructs.h"
using namespace std;
long gStrictLevel = 2;
long gInterleaveLen = -1;
bool gVerbose = false;
static long gBurnin = 0;
void processContent(PublicNexusReader & nexusReader, ostream *out);




bool newTreeHook(NxsFullTreeDescription &, void *, NxsTreesBlock *);



void writeTranslateCommand(NxsTaxaBlockAPI * taxa, std::ostream & out)
{
	if (!taxa)
		return;
	out << "    TRANSLATE" << "\n";
	const unsigned nt = taxa->GetNTaxTotal();
	for (unsigned i = 0; i < nt; ++i)
		{
		if (i > 0)
				out << ",\n";
		out << "        " << i + 1 << ' ' << NxsString::GetEscaped(taxa->GetTaxonLabel(i));
		}
	out << ";\n";
}

TreesToSplits * gTreesToSplitsB = NULL;
void writeStarTreeCommand(NxsTaxaBlockAPI * taxa, std::ostream & out)
	{
	out << "Tree star = " << " [&";
	out << (gTreesToSplitsB->treatAsRooted ? 'R' : 'U');
	out << "] (1";
	const unsigned nt = taxa->GetNTaxTotal();
	for (unsigned i = 1; i < nt; ++i)
		out << ',' << (1+i);
	out << ");\n";
	}

bool newTreeHook(NxsFullTreeDescription &ftd, void * arg, NxsTreesBlock *treesB)
{
	static unsigned long gTreeCount = 0;
	gTreeCount++;
	if (gTreeCount <= (unsigned long) gBurnin)
		{
		if (gVerbose)
			std::cerr << "Skipping tree " <<  gTreeCount<< '\n';
		return false;
		}
	if (gVerbose)
		std::cerr << "Read tree " <<  gTreeCount<< '\n';
	NxsTaxaBlockAPI *taxB = 0L;
	if (treesB)
		taxB = treesB->GetTaxaBlockPtr(NULL);
	if (!arg)
		return false;
	TreesToSplits * tts = (TreesToSplits *)(arg);
	tts->recordTree(ftd, taxB);
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
	BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();

	*out << "#NEXUS\n";
	for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt)
		{
		NxsBlock * b = *bIt;
		if (b->GetID() == "TAXA")
			b->WriteAsNexus(*out);
		}
	TreesToSplits::TaxaBlockToSplitsMap::const_iterator tbtsm = gTreesToSplitsB->taxaBlockToSplitsMap.begin();
	for (; tbtsm != gTreesToSplitsB->taxaBlockToSplitsMap.end(); ++tbtsm) {
		*out << "Begin Trees;\n";
		NxsTaxaBlockAPI * tb = tbtsm->first;
		if (tb && gTreesToSplitsB->taxaBlockToSplitsMap.size() > 1)
			*out << "Link Taxa = " << NxsString::GetEscaped(tb->GetTitle()) << ";\n";
		writeTranslateCommand(tb, *out);
		const TreesToSplits::NTreesSplitsMap & ntsm = tbtsm->second;
		const unsigned nt = ntsm.first;
		const TreesToSplits::SplitsMap & sm = ntsm.second;
		unsigned n = 1;
		writeStarTreeCommand(tb, *out);


		for (TreesToSplits::SplitsMap::const_iterator smIt = sm.begin(); smIt != sm.end(); ++smIt, ++n) {
			const SplitInfo & si = smIt->second;
			*out << "Tree split_" << n << " = " << " [&";
			*out << (gTreesToSplitsB->treatAsRooted ? 'R' : 'U');
			*out << "] [";
			if (gTreesToSplitsB->trackFreq || gTreesToSplitsB->trackEdgeLenSummary || gTreesToSplitsB->trackEdgeLen)
				*out << "&W " << ((double)si.nTimes)/((double)nt) << " ] [";
			if (gTreesToSplitsB->trackHeightSummary) {
				*out << "meanH =" << ((double)si.heightSum)/((double)si.nTimes);
			}
			*out << "] ";

			if (gTreesToSplitsB->trackEdgeLenSummary)
				smIt->first.writeNewick(*out, si.edgeLenSum/((double)si.nTimes), true);
			else
				smIt->first.writeNewick(*out, true);
			*out << ";\n";

		}
		*out << "End;\n";
	}

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
		if (!gTreesToSplitsB)
			gTreesToSplitsB = new TreesToSplits();
		treesB->setValidationCallbacks(newTreeHook, gTreesToSplitsB);
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
	out << "    -d report depth of splits\n\n";
	out << "    -e report edge length of splits\n\n";
	out << "    -b brief report of splits\n\n";
	out << "    -o report which trees a split occurs in.\n\n";
	out << "    -p report the proportion of trees that contain the split occurs in.\n\n";
	out << "    -r treat splits as rooted\n\n";
	out << "    -t report trivial splits.\n\n";
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
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'r')
			TreesToSplits::gTreatAsRooted = true;
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'v')
			gVerbose = true;
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 't')
			TreesToSplits::gTrackTrivial = true;
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'p')
			TreesToSplits::gTrackFreq = true;
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'o')
			TreesToSplits::gTrackOccurrence = true;
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'e')
			{
			el = true;
			}
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'd')
			{
			depth = true;
			}
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'b')
			{
			brief = true;
			}
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 's')
			{
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gStrictLevel)))
				{
				cerr << "Expecting an integer after -s\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'x')
			{
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gBurnin)))
				{
				cerr << "Expecting an integer after -x\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		//pass
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		//pass
#	else
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'i')
			{
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gInterleaveLen)) || gInterleaveLen < 1)
				{
				cerr << "Expecting a positive integer after -d\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
#	endif
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
			if (brief)
				{
				if (el)
					{
					TreesToSplits::gTrackEdgeLenSummary = true;
					TreesToSplits::gTrackFreq = true;
					}
				if (depth)
					{
					TreesToSplits::gTrackHeightSummary = true;
					TreesToSplits::gTrackFreq = true;
					}
				}
			else
				{
				if (el)
					TreesToSplits::gTrackEdgeLen = true;
				if (depth)
					TreesToSplits::gTrackHeight = true;
				}

			if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 'l')
				{
				readfile = true;
				readFilesListedIsFile(filepath+2, f);
				}
			else
				{
				readfile = true;
				readFilepathAsNEXUS(filepath, f);
				}
			}
		}
	if (!readfile)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
		}
	return 0;
	}


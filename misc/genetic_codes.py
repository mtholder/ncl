#!/usr/bin/env python
# Parsing of information from ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
# into code for NCL's (alphabetical) order of codons
import sys
from cStringIO import StringIO

codes = []

aaOrder = "ACDEFGHIKLMNPQRSTVWY*"
class Code(object):
    def __init__(self):
        self.short_name = ""
    def to_order(self, output_nuc_order):
        assert len(self.aa) == len(self.nuc_order[0])
        assert len(self.aa) == len(self.nuc_order[1])
        assert len(self.aa) == len(self.nuc_order[2])
        gen_code = {}
        for codon_n in xrange(len(self.aa)):
            aa = self.aa[codon_n]
            f = self.nuc_order[0][codon_n]
            s = self.nuc_order[1][codon_n]
            t = self.nuc_order[2][codon_n]
            codon = f + s + t
            gen_code[codon] = aa
        n = StringIO()
        for f in output_nuc_order:
            for s in output_nuc_order:
                for t in output_nuc_order:
                    codon = f + s + t
                    n.write(gen_code[codon])
        return n.getvalue()
    def to_order_ind(self, output_nuc_order):
        s = self.to_order(output_nuc_order)
        return [aaOrder.index(c) for c in s]



a = Code()
a.name = "Standard"
a.short_name = "SGC0"
a.id = 1
a.aa =         "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa =   "---M---------------M---------------M----------------------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
               "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
               "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Vertebrate Mitochondrial"
a.short_name = "SGC1"
a.id = 2
a.aa = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
a.start_aa = "--------------------------------MMMM---------------M------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Yeast Mitochondrial"
a.short_name = "SGC2"
a.id = 3
a.aa = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "----------------------------------MM----------------------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma"
a.short_name = "SGC3"
a.id = 4
a.aa = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "--MM---------------M------------MMMM---------------M------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Invertebrate Mitochondrial"
a.short_name = "SGC4"
a.id = 5
a.aa = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
a.start_aa = "---M----------------------------MMMM---------------M------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"
a.short_name = "SGC5"
a.id = 6
a.aa = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "-----------------------------------M----------------------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Echinoderm Mitochondrial; Flatworm Mitochondrial"
a.short_name = "SGC8"
a.id = 9
a.aa = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
a.start_aa = "-----------------------------------M---------------M------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Euplotid Nuclear"
a.short_name = "SGC9"
a.id = 10
a.aa = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "-----------------------------------M----------------------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Bacterial and Plant Plastid"
a.id = 11
a.aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "---M---------------M------------MMMM---------------M------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Alternative Yeast Nuclear"
a.id = 12
a.aa = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "-------------------M---------------M----------------------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Ascidian Mitochondrial"
a.id = 13
a.aa = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"
a.start_aa = "---M------------------------------MM---------------M------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Alternative Flatworm Mitochondrial"
a.id = 14
a.aa = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
a.start_aa = "-----------------------------------M----------------------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Blepharisma Macronuclear"
a.id = 15
a.aa = "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "-----------------------------------M----------------------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Chlorophycean Mitochondrial"
a.id = 16
a.aa = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "-----------------------------------M----------------------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Trematode Mitochondrial"
a.id = 21
a.aa = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
a.start_aa = "-----------------------------------M---------------M------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Scenedesmus obliquus Mitochondrial"
a.id = 22
a.aa = "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "-----------------------------------M----------------------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)

a = Code()
a.name = "Thraustochytrium Mitochondrial"
a.id = 23
a.aa = "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
a.start_aa = "--------------------------------M--M---------------M------------"
a.nuc_order = ["TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
              "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
              "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"]
codes.append(a)


maxId = 0
for code in codes:
    maxId = max(maxId, code.id)
sys.stdout.write("std::vector<std::string> cn(%d);\n" % (maxId))
sys.stdout.write("std::vector<std::string> code(%d);\n" % (maxId))

desired_order = "ACGT"
sys.stdout.write("/*\n")
for code in codes:
    assert code.nuc_order == a.nuc_order

    sys.stdout.write('  code index %d => "%s"\n'  % (code.id - 1, code.name))
sys.stdout.write("*/\n\n" )


for code in codes:
    n = code.name
    n = n.upper()
    sys.stdout.write(' NXS_GCODE_%s = %d,\n'  % (n, code.id - 1))
sys.stdout.write("\n\n" )

for code in codes:
    sys.stdout.write('code[%d] = "%s";\n'  % (code.id - 1, code.to_order(desired_order)))

sys.stdout.write('\n\n')
for n, code in enumerate(codes):
    l = code.to_order_ind(desired_order)
    if n == 0:
        x = l
        for j in xrange(64):
            sys.stdout.write('\taaInd[%d] = %d;\n' % (j, l[j]))
    else:
        sys.stdout.write('\tif (codeIndex = %s) {\n' % code.name)
        for j in xrange(64):
            if l[j] != x[j]:
                sys.stdout.write('\t\taaInd[%d] = %d;\n' % (j, l[j]))
        sys.stdout.write('\t}\n')

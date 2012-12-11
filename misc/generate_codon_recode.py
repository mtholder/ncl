#!/usr/bin/env python
import sys
f = open(sys.argv[1], "rU")
bases = "ACGT"
aaInd = "ACDEFGHIKLMNPQRSTVWY"
def getCodonStr(n):
    f = (n//16)
    s = ((n&15)//4)
    t = n&3
    return "%s%s%s" %(bases[f], bases[s], bases[t])

for line in f:
    enumName, aas = line.strip().split()
    el = enumName.lower()
    if False:
        #sys.stdout.write("%s %d \n" % (aas, len(aas)))
        sys.stdout.write("\tif(gCode == %s) {\n" % enumName)
        sys.stdout.write("\t\tconst int ccitac%s[] = {" % el)
        comCodedInt = 0
        for allCodonsInd, aaLetter in enumerate(aas):
            if aaLetter != "*":
                if comCodedInt != 0:
                    sys.stdout.write(", ")
                sys.stdout.write("%d" % allCodonsInd)
                comCodedInt += 1
        sys.stdout.write("};\n")
        sys.stdout.write("\t\tn = %d;\n" % (comCodedInt))
        sys.stdout.write("\t\tconst int caaind%s[] = {" % el)
        comCodedInt = 0
        for allCodonsInd, aaLetter in enumerate(aas):
            if aaLetter != "*":
                if comCodedInt != 0:
                    sys.stdout.write(", ")
                sys.stdout.write("%d" %  aaInd.index(aaLetter))
                comCodedInt += 1
        sys.stdout.write("};\n")

        sys.stdout.write("\t\tconst char * ccodstr%s[] = {" % el)
        comCodedInt = 0
        for allCodonsInd, aaLetter in enumerate(aas):
            if aaLetter != "*":
                if comCodedInt != 0:
                    sys.stdout.write(", ")
                sys.stdout.write('"%s"' % getCodonStr(allCodonsInd))
                comCodedInt += 1
        sys.stdout.write("};\n");
        sys.stdout.write("\t\tstd::copy(ccitac%s, ccitac%s + n, back_inserter(c.compressedCodonIndToAllCodonsInd));\n" % (el, el))
        sys.stdout.write("\t\tstd::copy(caaind%s, caaind%s + n, back_inserter(c.aaInd));\n"% (el, el))
        sys.stdout.write("\t\tstd::copy(ccodstr%s, ccodstr%s + n, back_inserter(c.codonStrings));\n" % (el, el))
        sys.stdout.write("\t\treturn c;\n\t}\n")
    else:
            #sys.stdout.write("%s %d \n" % (aas, len(aas)))
        sys.stdout.write("\tif(gCode == %s) {\n" % enumName)
        sys.stdout.write("\t\tconst int tr%s[] = {" % el)
        comCodedInt = 0
        for allCodonsInd, aaLetter in enumerate(aas):
            if comCodedInt != 0:
                sys.stdout.write(", ")
            if aaLetter != "*":
                sys.stdout.write("%d" % comCodedInt)
                comCodedInt += 1
            else:
                sys.stdout.write("-1")
        sys.stdout.write("};\n")
        sys.stdout.write("\t\tstd::copy(tr%s, tr%s + 64, back_inserter(v));\n" % (el, el))
        sys.stdout.write("\t\treturn v;\n\t}\n")

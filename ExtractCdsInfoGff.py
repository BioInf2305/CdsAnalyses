##this script extract cds coordinates and associated protein and mRNA id from NCBI gff file

import sys
import re
from collections import OrderedDict

def ExtractCdsInfo(gff_file,out_file):
    dest=open(out_file,"w")
    dest1=open("non_protein_coding.cds","w")
    with open(gff_file) as gff:
        mrna=''
        protein=''
        for line in gff:
            if not line.startswith("#"):
                a=line.strip().split("\t")
                pattern=re.compile('Name=(.[^;]*)')
                match=re.findall(pattern,line)
                if a[2]=="mRNA":
                    mrna=match[0]
                elif a[2]=="CDS":
                    #print(line)
                    if len(match)==0:
                        pattern=re.compile('GeneID:([0-9]+)')
                        match=re.findall(pattern,line)
                        dest1.write(a[0]+"\t"+a[3]+"\t"+a[4]+"\t"+a[6]+"\t"+mrna+"\t"+match[0])
                        dest1.write("\n")
                    else:
                        protein=match[0]
                        if mrna=='':
                            print("ERROR: no mrna information before cds")
                            sys.exit(0)
                        else:
                            dest.write(a[0]+"\t"+a[3]+"\t"+a[4]+"\t"+a[6]+"\t"+mrna+"\t"+protein)
                            dest.write("\n")
    dest.close()
inp=sys.argv[1]
outp=sys.argv[2]
ExtractCdsInfo(inp,outp)

import sys
import re
from collections import OrderedDict

def ExtractLongIso(inFile,outPrefix):
    dest=open(outPrefix+"_longIso.fa","w")
    dest1=open(outPrefix+"_longIso_gene_transcript_prot_id.txt","w")
    geneSeqDict=OrderedDict()
    geneProtDict=OrderedDict()
    geneTranscrDict=OrderedDict()
    geneRepeat=0
    seq=""
    with open(inFile) as inputFasta:
        for line in inputFasta:
            if line.startswith(">"):
                if seq!="":
                    if geneRepeat==1:
                        geneRepeat=0
                        if len(geneSeqDict[gene])<len(seq):
                            geneSeqDict[gene]=seq
                            geneProtDict[gene]=prot
                            geneTranscrDict[gene]=transcript
                    else:
                        geneSeqDict[gene]=seq
                        geneProtDict[gene]=prot
                        geneTranscrDict[gene]=transcript
                a=line.split()
                prot=a[0]
                pattern=re.compile(r'gene:([^\s]+)\stranscript:([^\s]+)')
                match=re.findall(pattern,line)
                transcript=match[0][1]
                gene=match[0][0]
                seq=""
                if match[0][0] in geneSeqDict:
                    geneRepeat=1
                else:geneRepeat=0
            else:
                seq+=line.rstrip()
    if geneRepeat==1:
        if len(geneSeqDict[gene])<len(seq):
            geneSeqDict[gene]=seq
            geneProtDict[gene]=prot
            geneTranscrDict[gene]=transcript
    for gene in geneSeqDict:
        protSeq=geneSeqDict[gene]
        if (len(protSeq)>=50) and ("*" not in protSeq):
            dest.write(">"+gene+"\n")
            dest.write("\n".join([protSeq[i:i+60] for i in range(0,len(protSeq),60)]))
            dest.write("\n")
            dest1.write(gene+"\t"+geneProtDict[gene]+"\t"+geneTranscrDict[gene])
            dest1.write("\n")
    dest.close()
    dest1.close()

if __name__=="__main__":
    ExtractLongIso(sys.argv[1],sys.argv[2])

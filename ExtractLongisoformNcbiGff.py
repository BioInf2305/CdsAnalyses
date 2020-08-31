import sys
import re
from collections import OrderedDict

'given NCBI gff and protein fasta file, this script can extract the protein sequence of the longest isoform of the gene'

def MainFunction(gff,prot,protOut):
    gffDict=gffToDict(gff)
    protDict=protToDict(prot)
    dest=open(protOut,"w")
    for gene in gffDict:
        if len(gffDict[gene])==1:
            dest.write(">"+gffDict[gene][0])
            dest.write("\n")
            tmpSeq=protDict[gffDict[gene][0]]
        else:
            tmpProtLength=[]
            tmpProtList=gffDict[gene]
            for proteins in tmpProtList:
                tmpProtLength.append(len(protDict[proteins]))
            dest.write(">"+tmpProtList[tmpProtLength.index(max(tmpProtLength))])
            dest.write("\n")
            tmpSeq=protDict[tmpProtList[tmpProtLength.index(max(tmpProtLength))]]
        tmpSeqList=[tmpSeq[i:i+60] for i in range(0,len(tmpSeq),60)]
        dest.write("\n".join(tmpSeqList))
        dest.write("\n")
    dest.close()

def gffToDict(gff):
    gffDict=OrderedDict()
    with open(gff) as source:
        for line in source:
            if not line.startswith("#") and "protein_id" in line:
                a=line.split()
                if a[2]=="CDS":
                    pattern1=re.compile(r'(GeneID)\:([0-9]+)')
                    tmpGeneId=pattern1.findall(line)
                    pattern2=re.compile(r'(Name)\=([A-Za-z0-9\.\_]+)')
                    tmpProtId=pattern2.findall(line)
                    if len(tmpGeneId[0])!=2 or len(tmpProtId[0])!=2:
                        print("Error at this line "+line)
                        sys.exit(1)
                    if tmpGeneId[0][1] not in gffDict:
                        gffDict[tmpGeneId[0][1]]=[]
                    if tmpProtId[0][1] not in gffDict[tmpGeneId[0][1]]:
                        gffDict[tmpGeneId[0][1]].append(tmpProtId[0][1])
    return gffDict

def protToDict(prot):
    protDict=OrderedDict()
    seq=""
    with open(prot) as source:
        for line in source:
            if line.startswith(">"):
                if seq!="":
                    protDict[nameSeq]=seq
                nameSeq=line.split()[0][1:]
                seq=""
            else:
                seq+=line.rstrip()
    protDict[nameSeq]=seq
    return protDict

if __name__=="__main__":
    MainFunction(sys.argv[1],sys.argv[2],sys.argv[3])

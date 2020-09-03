import sys
import os.path
import subprocess
import argparse
from multiprocessing import Pool
from collections import OrderedDict

class RunPrankGblockAlignment:

    def __init__(self,orthoInput,fileInfo,cores,outDir):
        self.orthoInput=orthoInput
        self.file1=fileInfo
        self.cores=int(cores)
        self.outDir=outDir
        self.speciesFileDict=OrderedDict() ##dict of species a key and a fasta file a value
        self.allSeqDict=OrderedDict()
        self.allSeqDictList=[]

    def ReadInputFile(self):
        with open(self.file1) as source1:
            for lin in source1:
                b=lin.split()
                self.speciesFileDict[b[0]]=b[1]

    def ReadFasta(self):
        self.ReadInputFile()
        for species in self.speciesFileDict:
            self.allSeqDict[species]=OrderedDict()
            with open(self.speciesFileDict[species]) as source2:
                seq=""
                head=""
                for line in source2:
                    line=line.rstrip()
                    a=line.split()
                    if line.startswith(">"):
                        if seq!="":
                            self.allSeqDict[species][head]=[]
                            self.allSeqDict[species][head].append(seq)
                        head=a[0][1:]
                        seq=""
                    else:
                        seq+=line
                if seq!="":
                    self.allSeqDict[species][head]=[]
                    self.allSeqDict[species][head].append(seq)


    def MakeSeparateDict(self):
        self.ReadFasta()
        speciesList=list(self.speciesFileDict.keys())
        os.mkdir(self.outDir)
        with open(self.orthoInput) as source:
            for line in source:
                line=line.rstrip()
                a=line.split()
                seqIds=a[1:]
                tmpSeqString=OrderedDict()
                tmpSeqString[a[1]]=OrderedDict()
                for i in range(len(seqIds)):
                    if self.allSeqDict[speciesList[i]].get(seqIds[i],"-1")!="-1":
                        tmpSeqString[a[1]][speciesList[i]]=self.allSeqDict[speciesList[i]][seqIds[i]]
                if len(tmpSeqString[a[1]])==len(speciesList):
                    self.allSeqDictList.append(tmpSeqString)
                else:
                    print(len(tmpSeqString[a[1]]))
                    print(tmpSeqString[a[1]])
                    print("Errror orthologous sequence missing on this line for one or more species "+line)
                    sys.exit(1)

    def RunPrankGblock(self,seqDict):
        seqName=list(seqDict.keys())[0]
        speciesSeq=seqDict[seqName]
        directory="./"+self.outDir+"/"+seqName
        os.mkdir(directory)
        logFile="logFile.txt"
        logFilePath=os.path.join(directory,logFile)
        fileName=seqName+".fasta"
        prankOut=seqName
        gblockOut=seqName+".best.fas"
        filePath=os.path.join(directory,fileName)
        filePrankOut=os.path.join(directory,prankOut)
        fileGblockOut=os.path.join(directory,gblockOut)
        fil=open(filePath,"w")
        logF=open(logFilePath,"w")
        for species in speciesSeq:
            fil.write(">"+species)
            fil.write("\n")
            seqString=speciesSeq[species][0]
            seqStringList=[seqString[i:i+60] for i in range(0,len(seqString),60)]
            fil.write("\n".join(seqStringList))
            fil.write("\n")
        fil.close()
        subprocess.call(["/home/maulik/software/orthologs/prank/prank-msa/src/prank","-codon","-F",\
        "-d="+filePath,"-o="+filePrankOut],stdout=logF,stderr=logF)
        cmdGblock="o"+"\n"+fileGblockOut+"\n"+"t"+"\n"+"g"+"\n"+"q"
        subprocess.call(["/home/maulik/software/orthologs/Gblocks_0.91b/Gblocks << here\n"+cmdGblock+"\nhere\n"],stdout=logF,stderr=logF,shell=True)
        logF.close()

    def RunParallel(self):
        self.MakeSeparateDict()
        for i in range(0,len(self.allSeqDictList),int(self.cores)):
            processList=self.allSeqDictList[i:i+int(self.cores)]
            print(processList)
            with Pool(processes=len(processList)) as pool:
                pool.map(self.RunPrankGblock,self.allSeqDictList)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="This python script will run the prank and gblock in parallel given the orthologs coding sequences of several species")
    parser.add_argument("-orthoF","--orthoFile",metavar="File",help="File containing information about orthologs('OrthoInfo1.txt')",required=True)
    parser.add_argument("-listF","--listFiles",metavar="File",help="File containing names of species and associated fasta files", required=True)
    parser.add_argument("-coreN","--numCores",metavar="Int",type=int,help="Number of cores to use",required=True)
    parser.add_argument("-outD","--outDir",metavar="Str",help="Name of the output dir (the script will make it)",required=True)
    args=parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        runPrankGblock=RunPrankGblockAlignment(args.orthoFile,args.listFiles,args.numCores,args.outDir)
        runPrankGblock.RunParallel()

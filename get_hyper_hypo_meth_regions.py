"""
Input: bedfile, methylome in wig format
Output: 2 files: 1 contains bed regions with average methyl level >=0.5, the other <0.5

command: python function.py bedfile wigfile
"""
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
from sys import argv


def calcaveragemeth(region,methlome):
    chrom=region[0]
    start=region[1]
    end=region[2]
    mlevels=[]
    for i in range(start,end):
        try:
            methlevel=methlome[chrom][i]
            #print methlevel,i
            mlevels+=[methlevel]
        except KeyError:
            pass
    
    try:
        return sum(mlevels)/len(mlevels)
    except ZeroDivisionError:
        return 0.0

def readmethylome(wigfile):
    seq=open(wigfile).read()
    seq=seq.strip().split('chrom')
    seq=seq[1:]
    for i in range(len(seq)):
        seq[i]=seq[i].strip().split('\n')
    metdata={}
    for i in range(len(seq)):
        chrom=seq[i][0][1:]
        print chrom
        metdata[chrom]={}
        print len(seq[i]),"positions"
        for j in range(1,len(seq[i])):
            tmp=seq[i][j].strip().split("\t")
            try:
                coor=int(tmp[0])
                beta=float(tmp[1])
                metdata[chrom][coor]=beta
            except:
                #print seq[i][j] #this is at the end of each chromosome
                continue
    return metdata

def readbedfile(bedfile):
    regions=[]
    lines=open(bedfile).read().strip().split('\n')
    for l in lines:
        tmp=l.strip().split('\t')
        region=[tmp[0],int(tmp[1]),int(tmp[2])]
        regions+=[region]
    return regions

def main():
    bedfile=argv[1]
    methfile=argv[2]
    outfile=argv[3]
    print "Reading bed regions..."
    DHS=readbedfile(bedfile)
    print "Reading Methylation data..."
    methdata=readmethylome(methfile)
    print "calculating average methylation levels per region..."
    avgscores=[]
    out1=open(outfile+".hyper.bed",'w')
    out2=open(outfile+".hypo.bed",'w')
    for region in DHS:
    	score=calcaveragemeth(region,methdata)
        avgscores+=[score]
        if score>=0.5: #if hyper
        	out1.write(region[0]+'\t'+str(region[1])+'\t'+str(region[2])+'\n')
        else: #if hypo
        	out2.write(region[0]+'\t'+str(region[1])+'\t'+str(region[2])+'\n')
        if len(avgscores)%10000==0:
            print len(avgscores)
    out1.close()
    out2.close()
    #print "plotting histogram..."
    #plt.hist(avgscores,bins=100);
    #plt.plot([1,2,3])
    #plt.savefig(outfile)
    return

if __name__=="__main__":
	main()

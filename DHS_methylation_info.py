'''This program gives information about DHS and their methylation levels

Input: bedfile, wig files for methylation
Output: 
- bed file with: 
	+ length of region
	+ average methylation levels per region
	+ number of methylation sites
- distro. of methyl. levels in DHS as histogram in png format, 2 modes:
	+ average methylation level of DHS regions with at least 2 methylated sites
	+ ---- no lower limit of methylated sites
- distro. of methyl. levels in non-DHS as histogram 
- distro. of methyl. levels in the genome as histogram

'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sys import argv


def calcaveragemeth(region,methlome):
    chrom=region[0]
    start=region[1]
    end=region[2]
    mlevels=[]
    for i in range(start,end,2): #step size is 2 because both C and G are marked as methylated
        try:
            methlevel=methlome[chrom][i]
            #print methlevel,i
            mlevels+=[methlevel]
        except KeyError:
            pass
    
    '''try:
        return sum(mlevels)/len(mlevels)
    except ZeroDivisionError:
        return 0.0'''
    numberofmods=len(mlevels)
    if numberofmods!=0:
    	return sum(mlevels)/len(mlevels),numberofmods
    else:
    	return 0.0,numberofmods

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
        region=[tmp[0],int(tmp[1]),int(tmp[2])] #(chrom,start,end)
        regions+=[region]
    return regions

def main():
    bedfile = argv[1]
    methfile = argv[2]
    tag = argv[3]
    print "Reading bed regions..."
    DHS = readbedfile(bedfile)
    print "Reading Methylation data..."
    methdata = readmethylome(methfile)
    print "calculating average methylation levels per region..."
    avgscores = []
    numberofmods = []
    regionlens = []
    for region in DHS:
    	#regionlen = region[2]-region[1]
        average,nummods = calcaveragemeth(region,methdata)
        avgscores+=[average]
        numberofmods+=[nummods]
        #regionlens+=[regionlen]
        if len(avgscores)%10000==0:
            print len(avgscores)
    

    print "Outputing info file"
    outfilebed=tag+'.methinfo.bed'
    target=open(outfilebed,'w')
    for i in range(len(DHS)):
    	region=DHS[i]
    	line = '\t'.join(map(str,region))+'\t'+str(region[2]-region[1])+'\t'+str(avgscores[i])+'\t'+str(numberofmods[i])
    	target.write(line+'\n')
    target.close()

        
    '''    
    print "plotting histogram..."
    plt.hist(avgscores,bins=100);
    #plt.plot([1,2,3])
    plt.savefig(outfile)
    '''
    return


if __name__=="__main__":
	main()

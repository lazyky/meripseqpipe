import os
import sys
designfile = sys.argv[1]
tophat2_skip = sys.argv[2]
hisat2_skip = sys.argv[3]
bwa_skip = sys.argv[4]
star_skip = sys.argv[5]
with open(designfile) as name_list:
    name_list.readline()
    for line in name_list:
        data = line.replace('\n','').replace('\r','').split(',')
        if tophat2_skip == "false":
            oldname = data[0]+'_tophat2.bam'
            newname = data[0]+'_'+data[1]+'_'+data[2]+'_'+data[3]+'_tophat2.bam' 
            os.rename(oldname,newname)
        if hisat2_skip == "false":
            oldname = data[0]+'_hisat2.bam'
            newname = data[0]+'_'+data[1]+'_'+data[2]+'_'+data[3]+'_hisat2.bam'
            os.rename(oldname,newname)
        if bwa_skip == "false":
            oldname = data[0]+'_bwa.bam'
            newname = data[0]+'_'+data[1]+'_'+data[2]+'_'+data[3]+'_bwa.bam'
            os.rename(oldname,newname)
        if star_skip == "false":
            oldname = data[0]+'_star.bam'
            newname = data[0]+'_'+data[1]+'_'+data[2]+'_'+data[3]+'_star.bam'
            os.rename(oldname,newname)
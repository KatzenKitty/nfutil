#!/usr/bin/env python
"""
# Created: Tue, 18 Oct 2011 13:35:27 +1000

%prog attempts to resolve scaffolding errors in a assembled contigs using read mapping. 

TODO: Detailed description
"""
# IMPORTS
import sys, os, traceback, optparse
import time
import re
import pysam
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from numpy import *

# MetaData
__author__ = "Nabil-Fareed Alikhan"
__licence__ = "GPLv3"
__version__ = "0.3"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"

# Globals config vars
contigsFile = '/media/backup/assembly/n1/rerun/ins100/contigs.fa'
GBKident = 'n1'
GBKdesc = ''
outDir = 'tmp/'
readLocation = '/media/backup/assembly/n1/in_reads/n1_proc.fastq'
numThreads = 4
numSDEV = 1
allowMisMatch = 2

def main ():

    global options, args
    contigsFile = options.filename
    readLocation = options.reads
    outDir = options.out
    if options.annotate: print 'Running annotation-only mode (-a flag), script will concatenate multi-fasta, create new fasta file & create GenBank file (annotating contig boundaries & Ns)'
    else: print 'Starting assembly patching...'
    # Create biopython seq-annotation struct - masterSeq
    nullseq = Seq("", generic_dna)
    masterSeq = SeqRecord(nullseq, id=GBKident, name=GBKident, description=GBKdesc)
    # Load contigs fasta sequence 
    offset = 0
    if options.verbose:
        print "contig_name\tcontig_length\tN_start\tN_end\tN_length\tString"
    for seq_record in SeqIO.parse(contigsFile, "fasta", generic_dna):
        masterSeq.seq += seq_record.seq
    # When N is detected:
        for match in re.finditer( '[Nn]+', str(seq_record.seq) ) :
            # Mark all locations of Ns, start location but save number of Ns as a note
            sf = SeqFeature(FeatureLocation(offset + match.start() , offset + match.end()), type="misc_feature", qualifiers = { "label" : [ str(match.end() - match.start()) + " Scaffolding Ns" ]})
            masterSeq.features.append(sf)
            # If running verbose, output location of N runs. 
            if options.verbose:
                print seq_record.id + "\t" + str(len(seq_record.seq)) + "\t" + str(match.start()) + "\t" + str(match.end()) + "\t" + str(match.end() - match.start() ) + "\t" + match.group()
        ### Closed loop for N search ### 
        # When contig boundary is detected:
        contStart = offset
        offset += len (seq_record.seq)

        # Save contig boundary as feature
        sf = SeqFeature(FeatureLocation(contStart , offset ), type="fasta_record", qualifiers = { "label" : [ seq_record.id]})
        masterSeq.features.append(sf)
    # Annotate mode ends here
    if options.annotate:
        print 'Writing annotated gbk & fna to ' + outDir + "/" + GBKident
        SeqIO.write(masterSeq, outDir + "/" + GBKident + ".gbk", "genbank")
        SeqIO.write(masterSeq, outDir + "/" + GBKident + ".fna", "fasta")
    else:
        # If running with verbose, create a temporary snapshow of sequences for debugging purposes
        if options.verbose:
            print 'Writing debug gbk to ' + outDir + "/" + GBKident + "Debug.gbk" 
            SeqIO.write(masterSeq, outDir + "/" + GBKident + "Debug.gbk", "genbank")
            SeqIO.write(masterSeq, outDir + "/" + GBKident + "Debug.fna", "fasta")
        # Shift subsequent annotation by cumulative number of Ns
        for feat in masterSeq.features:
            # How many Ns before the start ? 
            nBeforeStart = masterSeq.seq.count('N', 0, int(feat.location.nofuzzy_start) )
            # How many Ns before the end ? 
            nBeforeEnd = masterSeq.seq.count('N', 0, int(feat.location.nofuzzy_end))
            # Ammend the feature
            feat.location = FeatureLocation( feat.location.nofuzzy_start - nBeforeStart , feat.location.nofuzzy_end - nBeforeEnd)
        # Delete all Ns    
        noNs = Seq(str(masterSeq.seq).replace('N','') , generic_dna) 
        masterSeq.seq = noNs
        # START BWA sub-routine
        if os.path.exists(outDir + "/" + GBKident + ".sorted.bam") == False:
            # Create new sequence 
            bwaFasta = outDir + "/tmpBWA.fna"
            SeqIO.write(masterSeq, bwaFasta, "fasta")
            # Use BWA for read mapping onto new sequence
                # Index reference
            bwaIndex = 'bwa index ' + bwaFasta    
            print 'Executing ' + bwaIndex 
            os.system(bwaIndex)
                # Align with aln, samse
            bwaAln = 'bwa aln -t ' + str(numThreads) + " -n " + str(allowMisMatch) + " " + bwaFasta + ' ' + readLocation + ' > ' + outDir +'/tmpBWA.sai'    
            bwaSamse = 'bwa samse ' + bwaFasta + ' ' + outDir + '/tmpBWA.sai ' + readLocation + ' > ' + outDir + '/tmpBWA.sam'
            print 'Executing ' + bwaAln
            os.system(bwaAln)
            print 'Executing ' + bwaSamse
            os.system(bwaSamse)
                # Convert output into BAM
            samBam =  'samtools view -bS ' + outDir + '/tmpBWA.sam > ' + outDir + "/tmpBWA.bam"   
            print 'Executing ' + samBam
            os.system(samBam)
                # Sort BAM
            bamSort = 'samtools sort ' + outDir + "/tmpBWA.bam " + outDir + "/" + GBKident + ".sorted"
            print 'Executing ' + bamSort 
            os.system(bamSort)
            bamIndex = 'samtools index ' + outDir + "/" + GBKident + ".sorted.bam"
            print 'Executing ' + bamIndex
            os.system(bamIndex)
            # Tidy up unneeded files from bwa/sam output
        else:
            print 'Existing BAM detected: ' + outDir + "/" + GBKident + ".sorted.bam. This is the right one, no? "
        # Scout N locations from BAM file and whether reads have spanned this region. Coverage? 
        # We assume that coverage should not deviate from average contig coverage by 1 standard deviation. 

        print 'Reading samfile...'
        samfile = pysam.Samfile( outDir + "/" + GBKident + ".sorted.bam", "rb" )
        curContig = ''
        contigTable = { }
        if options.verbose: print 'Rating\tstart\tstop\tCov\tCov_cutoff\tContigMean\tContigDev'
        for feat in masterSeq.features:
            # what is current contig?
            if feat.type == 'fasta_record':
                curContig = feat    
            # Misc_feature i.e scaffolding, coverage should not be  less than one standard deviation of contig mean. 
            if feat.type == 'misc_feature':
                # Look up contig mean coverage and std dev. If not found, calculate and add to store values.
                # print contigTable[ curContig.qualifiers['label'][0] ] 
                if curContig.qualifiers['label'][0] not in contigTable:
                    coverArray = [ ] 
                    for pileupcolumn in samfile.pileup( 'n1', curContig.location.nofuzzy_start , curContig.location.nofuzzy_end): 
                        if pileupcolumn.pos > (curContig.location.nofuzzy_start -1 ) and pileupcolumn.pos < (curContig.location.nofuzzy_end ):
                            coverArray.append(int(pileupcolumn.n) )
                    coverArray.sort()
   
                    contigTable[curContig.qualifiers['label'][0] ]  = [ array(coverArray).mean(), array(coverArray).std() , coverArray[len(coverArray)/2], coverArray[len(coverArray)/4] ]  
                # Find coverage for scaffolding point 
                scafCov = -1
                covcov = [] 
                for pileupcolumn in samfile.pileup( 'n1', feat.location.nofuzzy_start -2 , feat.location.nofuzzy_end+2    ):
                    if pileupcolumn.pos >= (feat.location.nofuzzy_start -2 ) and pileupcolumn.pos < (feat.location.nofuzzy_end+2   ):
                        # print "S " + str(pileupcolumn.n) + " "  + str(pileupcolumn.pos)
                        covcov.append(pileupcolumn.n)
                while len(covcov) < 4:
                    covcov.append(0)
                covcov.sort()
#                print covcov
                scafCov = array(covcov).min()        
                # If coverage < (mean - std dev), Scaffolding is good, repair, print stats to support
                contCovMean = contigTable[curContig.qualifiers['label'][0]][0]
                contCovDev = contigTable[curContig.qualifiers['label'][0]] [1]
                contMedian = contigTable[curContig.qualifiers['label'][0]] [2]
                contquart = contigTable[curContig.qualifiers['label'][0]] [3]
                contCovCut = contCovMean - (contCovDev * numSDEV)
                if scafCov <= contCovCut:
                    print ('GOOD'+"\t" + str(feat.location.nofuzzy_start) +"\t" +  str(feat.location.nofuzzy_end ) + "\t" + str(scafCov) + "\t"+  str(contCovCut) + "\t" + str(contCovMean) + "\t" +str(contCovDev) + "\t" +str(contMedian) + "\t" +str(contquart))
                else:
                # If not, leave deleted.
                    print 'BAD'+'\t' + str(feat.location.nofuzzy_start) +"\t" +  str(feat.location.nofuzzy_end ) + "\t" + str(scafCov) + "\t"+  str(contCovCut) + "\t" + str(contCovMean) + "\t" +str(contCovDev) + "\t" +str(contMedian) + "\t" +str(contquart)
                    # Re-insert Ns
        # Output new fasta & Genbank of new file
        print 'Writing genbank & fasta to ' + outDir + "/" + GBKident + ".gbk"
        SeqIO.write(masterSeq, outDir + "/" + GBKident + ".gbk", "genbank")
        # Re-insert contig boundaries, if flagged to do so (--contigs)
        if options.contigs:
            # TODO:
            print 'Re-inserting contig boundaries'
        SeqIO.write(masterSeq, outDir + "/" + GBKident + ".fna", "fasta")



if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[-2]
        parser = optparse.OptionParser(epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option ('-a', '--annotate', action='store_true', default=False,help='only annotate contig boundaries and N characters.')
        parser.add_option ('-c', '--contigs', action='store_true', default=False,help='Final FASTA file should be a multi-fasta of contigs.')
        parser.add_option("-f", "--fasta", action="store", type="string", dest="filename", help='location of contigs.fa')
        parser.add_option("-r", "--reads", action="store", type="string", dest="reads", help='location of reads to map (fastq)')
        parser.add_option("-o", "--out", action="store", type="string", dest="out", help='location of output')
        (options, args) = parser.parse_args()
        if options.verbose: print "Executing @ " + time.asctime()
        main()
        if options.verbose: print "Ended @ " + time.asctime()
        if options.verbose: print 'total time in minutes:',
        if options.verbose: print (time.time() - start_time) / 60.0
        sys.exit(0)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)


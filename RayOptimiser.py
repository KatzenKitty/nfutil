#!/usr/bin/env python
"""
# Created: Thu, 15 Dec 2011 15:46:54 +1000

Will determine optimal assembly if given a directory full of Ray assemblies

TODO: Detailed description
"""
import sys, os, traceback, optparse
import time
import re
import numpy
import shutil
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from pexpect import run, spawn
__author__ = "Nabil-Fareed Alikhan"
__licence__ = "GPLv3"
__version__ = "0.3"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"

def main ():

    global options, args
    # INPUT: Directory of Ray assemblies of different kmers
    # ###### Output directory to copy optimal assembly
    # To find optimal kmer. 
    contigstats = [] 
    scafstats = []
    # For each folder in specified folder
    path = options.dir
    alldirList = os.listdir(path)
    # Filter folders that aren't Ray assemblies
    dirList = []
    if options.verbose: print "###" + options.dir + "###"
    for fname in alldirList:
        check = 0
        if os.path.isfile(os.path.abspath(path) +  "/" + fname + "/LibraryStatistics.txt"):
            check += 1 
        if os.path.isfile(os.path.abspath(path) +  "/" + fname + "/Scaffolds.fasta"):
            check += 1
        if os.path.isfile(os.path.abspath(path) +  "/" + fname + "/Contigs.fasta"):
            check += 1
        if check == 3: dirList.append( fname )    
    peak0Mean = []
    peak0Sdev = []
    if options.verbose: print 'Folder\tPeak\tMean\tSdev'
    for fname in dirList:
        lib = os.path.abspath(path) +  "/" + fname + "/LibraryStatistics.txt"
        # Find Peak0 mean & Sdev. Set ave. scaffold length limit = mean + 2 sdev
        f = open(lib, 'r')
        peakZ = 0
        stats = fname + '\t'
        for line in f:
            if re.search( "Peak" , line ):
                if re.search( "Peak 0" , line ):
                    peakZ = 0
                else:
                    peakZ += 1
                stats += str(peakZ) + '\t' 
                # Find insert size peaks
                # Get peak 0 details
                # If more than 1, add to twinpeaks array 
            if re.search ("AverageOuterDistance: \d+", line):
                m = re.match("\s+AverageOuterDistance: (\d+)", line)
                stats += m.group(1) + "\t"
                if peakZ == 0:
                    peak0Mean.append(int( m.group(1) ))
            if re.search ("StandardDeviation: \d+", line):
                m = re.match("\s+StandardDeviation: (\d+)",line)
                stats += m.group(1) + "\t"
                if peakZ == 0:
                    peak0Sdev.append(int ( m.group(1) ) )
        f.close()
        if options.verbose: print stats
    scafCut = numpy.median(peak0Mean) + 3 * numpy.median(peak0Sdev)        
    if options.verbose: print "Scaffolding cutoff is: " + str (scafCut)
    # Screen according to scaffolding length limit and not more than 20k Ns
    scafFilt = []
    if options.verbose: print 'Folder\tNo.Ns\tN-gaps\tMin\tMax\tMean'
    numScaff = {}
    for fname in dirList:
        cont = os.path.abspath(path) +  "/" + fname + "/Scaffolds.fasta"
        allNs = []
        for seq_record in SeqIO.parse( cont , "fasta"):
            if fname in numScaff.keys():
                numScaff[fname] += 1 
            else: 
                numScaff[fname] = 0 
            contigNs = re.findall('[Nn]+', str(seq_record.seq)) 
            for co in contigNs:
                allNs.append(co)
         # Screen kmer here:
        Nlen = []
        for run in allNs:
            Nlen.append( len(run) )
        Nmean = numpy.average(Nlen)
        NoNs = sum(Nlen)
        if NoNs < 20000 and  Nmean < scafCut:
            # Ignore nmer 17
            if fname != '17':
                scafFilt.append(fname)
        # Print Scaffolding stats
        if options.verbose: print fname +'\t' + str(sum(Nlen)) +'\t'+ str(len(allNs)) +'\t'+ str(numpy.min(Nlen)) +'\t'+ str(numpy.max(Nlen)) +'\t'+ str(Nmean)
    if options.verbose: print "Choosing from kmers: " + str(scafFilt)
    # Determine highest contigs n50 as optimal kmer, from remaining kmers.
    if options.verbose: print 'Folder\tn50\tn75\tn90\tNumContigs'
    bestFold = '' 
    bestN50 = 0 
    for fname in dirList:
        cont = os.path.abspath(path) +  "/" + fname + "/Contigs.fasta"
        ContigLen = []
        for seq_record in SeqIO.parse( cont , "fasta"):
            noN = len(str(seq_record.seq))
            ContigLen.append( noN )
        ContigLen.sort(reverse=True)
        n50 = 0 
        n75 = 0
        n90 = 0
        sumit = 0
        totalbase = sum(ContigLen)
        for j in ContigLen:
            sumit += j
            if sumit >= totalbase / 2 and n50  == 0:
                n50 = j
            elif sumit >= totalbase *0.75 and n75 == 0:
                n75 = j
            elif sumit >= totalbase * 0.9 and n90 == 0:
                n90 = j
        # Print Contigs stats
        if options.verbose: print fname + '\t'+str(n50) + '\t' +str(n75)+ '\t'+  str(n90)+'\t'+str( len(ContigLen) )
        isValid = False
        for scf in scafFilt:
            if scf == fname: isValid = True
        if n50 > bestN50 and isValid:
            bestFold = fname
            bestN50 = n50 
    if bestN50 == 0:
        print options.dir + '\t-1'
    else:
        print options.dir + '\t' +bestFold + "\t" + str(bestN50)+ "\t" + str(numScaff[bestFold])
    if options.verbose: print ' '    
    # COPY files across to output dir
    if options.out != None and bestFold != '':
        bestdir = os.path.abspath(path) +  "/" + bestFold
        copyList = os.listdir(bestdir)
        for fname in copyList:        
            shutil.copy(bestdir + '/' + fname, options.out + '/' + fname)

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[-2]
        parser = optparse.OptionParser(epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option ("-i", "--dir", action="store", type="string", dest="dir", help='Directory of Ray assemblies')
        parser.add_option ('-o','--out', action="store",type="string",dest="out", help='Directory to copy best Ray assembly')
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


#!/usr/bin/env python
"""
# Created: Mon, 31 Oct 2011 13:19:03 +1000

%prog bins paired in reads in BAM file based off insert size

WARNING: be sure to run with python 2.6 (load module python on barrine)
"""
import sys
sys.path.append('/home/uqnalikh/lib/pysam/')
import os
import traceback, optparse
import time
import re
import pysam 

#from pexpect import run, spawn
__author__ = "Nabil-Fareed Alikhan"
__licence__ = "GPLv3"
__version__ = "0.3"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"

Bins = [ [120, 180], [181,270], [271,330]  ] 
output = {} 

def main ():
    inFile = options.input
    outDir = options.out
    name = options.prefix 
    samfile = pysam.Samfile( inFile , "rb" )
    iter = samfile.fetch()
    # There will be len(Bins) + 3 special bins, unmapped & mapped 
    # but didnt fit a bin & neg insert pairs reads. 
    cleanOutfiles(outDir, name )
    noGood = 0 
    unmapCount = 0
#    overCount = 0 
    binCount = 0
    unbinCount = 0
    totalReads = 0 

    unmap = (outDir +"/" + name + "-unmapped.fq")
#    overmap = (outDir +"/" + name + "-negisize.fq")
    nogoodmap = (outDir + "/" + name + "-nogood.fq")
    unbinmap = (outDir+ "/" + name + "-unbinned.fq")
    binmap = [ ]
    output[outDir +"/" + name + "-unmapped.fq"] =  []
 #   output[outDir +"/" + name + "-negisize.fq"] = []
    output[outDir + "/" + name + "-nogood.fq"] = []
    output[outDir+ "/" + name + "-unbinned.fq"] = []
    for bin in Bins:
  #      binner = open(outDir +"/" + name + "-"+ str(bin[0]) + "." + str(bin[1])+".fq",'w')
        binner = outDir +"/" + name + "-"+ str(bin[0]) + "." + str(bin[1])+".fq"
        output[binner] = []
    mates = {} 
    if options.verbose: print 'Precaching mates' 
    for read in iter:
        totalReads += 1 
        if options.verbose and totalReads % 1000000 ==  0: print str(totalReads) +'\t' + time.asctime()
        if read.is_read2:
            mates[read.qname] = [read.seq, read.qual] 
    totalReads = 0
    # Mates are loaded. Now read reads.
    if options.verbose: print 'Binning reads'
    samfile = pysam.Samfile( inFile , "rb" )
    iter = samfile.fetch()
    for read in iter:
        totalReads += 1 
        if options.verbose and totalReads % 1000000 ==  0:
            print str(totalReads) +'\t' + time.asctime()
            cleanbinCache()
        # Place in unmapped bin if applicable, either read and pair are unmapped
        # i.e insert cannot be determined 
        if read.is_read1:
            if read.is_unmapped == True:
                 matemateToBin(unmap, read,  mates[read.qname] , False )
                 unmapCount = unmapCount + 1 
       
       # Place in negative insert size bin if applicable 
      #      elif read.isize < 0:
       #         matemateToBin(overmap, read,  mates[read.qname], False  )
        #        overCount = overCount + 1 
                #print read  
       # Iterate through bins, place if it fits in one
            elif read.is_proper_pair:  
                binned = 0 
                box = 0 
                for bin in Bins:
                    start = bin[0]
                    stop = bin[1]
                    if abs(read.isize) >= start and abs(read.isize) <= stop:
                        binbin =  outDir +"/" + name + "-"+ str(bin[0]) + "." + str(bin[1])+".fq"

                        matemateToBin(binbin, read,  mates[read.qname] , False )
                        binned = binned + 1 
                        binCount = binCount + 1 
                    box = box  + 1
                if binned == 0 :
                    # else, dump in mapped but didnt fit bin. 
                #    writeToBin( unbinmap, read,   samfile.mate(read) )
                    matemateToBin(unbinmap , read,  mates[read.qname], False  )
                    unbinCount = unbinCount + 1 
                
            else: 
                noGood = noGood +1
                matemateToBin(nogoodmap, read,  mates[read.qname],False   )
    cleanbinCache();
    print 'Reads & pair errored : ' + str( noGood )
    print 'Reads & pair unmapped: ' + str( unmapCount ) 
#    print 'Reads & pair negative insert size: ' + str( overCount ) 
    print 'Reads & pair binned : ' + str( binCount )
    print 'Reads & pair mapped but unbinned: ' + str( unbinCount )
    print 'Total reads in BAM: ' + str( totalReads )

def cleanOutfiles(outDir, name ):
    unmapped =  outDir +"/" + name + "-unmapped.fq" 
    unbinned =  outDir +"/" + name + "-unbinned.fq"
    overlapped = outDir +"/" + name + "-negisize.fq"
    unbinmap = outDir+ "/" + name + "-unbinned.fq"

    if os.path.exists( unmapped  ):
        os.remove( unmapped )
    if os.path.exists( unbinned  ):
        os.remove(unbinned )
    if os.path.exists( overlapped ):
        os.remove(overlapped)
    if os.path.exists( unbinmap):
        os.remove(unbinmap) 
    for bin in Bins:
        fil = outDir +"/" + name + "-"+ str(bin[0]) + "." + str(bin[1])+" .fq"
        if os.path.exists(fil):
            os.remove(fil)

def cleanbinCache() :
    if options.verbose: print 'Writing reads to file'
    for file in output:
        wr = open(file,'a')
        wr.write(''.join(output[file]) )
        output[file] = []


def matemateToBin(file, read, matearray, swap ):
    dump = ''
    if(swap == True and matearray != '' ):
        dump += ('@' + read.qname + '/1\n')
        dump += (matearray[0] + '\n')
        dump += ("+" + read.qname + '/1\n')
        dump += (matearray[1] + '\n' )
        dump = ('@' + read.qname + "/2\n")
        dump += (read.seq + '\n')
        dump += ("+" + read.qname + '/2\n')
        dump += (read.qual + '\n')
    else:
        dump += ('@' + read.qname + "/1\n")
        dump += (read.seq + '\n')
        dump += ("+" + read.qname + '/1\n')
        dump += (read.qual + '\n')
        if matearray != '':
            dump += ('@' + read.qname + '/2\n')
            dump += (matearray[0] + '\n')
            dump += ("+" + read.qname + '/2\n')
            dump += (matearray[1] + '\n' )
    output[file].append(dump) 

def writeToBin(file, read, pair ):
    
    dump = ('@' + read.qname + "/1\n")
    dump += (read.seq + '\n')
    dump += ("+" + read.qname + '/1\n')
    dump += (read.qual + '\n')
    if pair != '':
        dump += ('@' + pair.qname + '/2\n')
        dump += (pair.seq + '\n')
        dump += ("+" + pair.qname + '/2\n')
        dump += (pair.qual + '\n' )
    output[file].append(dump) 
    if( len(output[file]) > 1000000):
        if options.verbose:
            print 'Writing to ' + file 
        wr = open(file,'a')
        wr.write(''.join(output[file]) )
        output[file] = [] 

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[-2]
        parser = optparse.OptionParser(epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option("-i", "--input", action="store", type="string", dest="input", help='input BAM file')
        parser.add_option("-o", "--output", action="store", type="string", dest="out", help='output folder')
        parser.add_option("-p", "--prefix", action="store", type="string", dest="prefix", help='output prefix')
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


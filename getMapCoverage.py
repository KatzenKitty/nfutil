#!/usr/bin/env python
"""
# Created: Mon, 31 Oct 2011 13:19:03 +1000

%prog fetches depth of coverage from a BAM file 

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

def main ():
    inFile = args[0]
    f = open(inFile,'r')
    lines = f.readlines()
    fout = open(inFile + '.tab','w')
    dat = []
    head = [] 
    for line in lines:
        samfile = pysam.Samfile( line.strip() , "rb" )
        head.append( os.path.basename(line) ) 
        col = [  ] 
        for pileupcolumn in samfile.pileup():
            if pileupcolumn.pos >= len(col):
                zer = [0] * ( 1 + pileupcolumn.pos - len(col) )
                col +=  zer
            col[ pileupcolumn.pos ] = pileupcolumn.n
#        col.insert(0, os.path.basename(line.strip()) ) 
        samfile.close()
        dat.append(col)
    print 'transposing' 
    transpose = [list(c) for c in zip(*dat)]
    d = 0 
    print 'writing to ' + inFile + '.tab'
    for z in head:
        fout.write(str(z).strip() + '\t')
    fout.write('\n' )
    for x in transpose: 
        fout.write(str(d))
        for t in x:
            fout.write( '\t' + str(t))
        fout.write('\n' )
        d += 1
    print 'meep ' + str(d)

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[-2]
        parser = optparse.OptionParser(epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
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


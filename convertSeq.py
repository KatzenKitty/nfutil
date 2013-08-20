#!/usr/bin/env python
"""
# Created: Mon, 11 Mar 2013 16:11:10 +1000

Converts files between GBK, EMBL, FASTA using BioPython

A simple script that converts between different genome file formats. 
This script requires Biopython (v1.52+) installed to run.
See <http://biopython.org>

### CHANGE LOG ### 
2013-03-11 Nabil-Fareed Alikhan <n.alikhan@.uq.edu.au>
    * Formatted for NF argparse python template 
"""
import sys, os, traceback, argparse
import time
from Bio import SeqIO

__author__ = "Nabil-Fareed Alikhan"
__licence__ = "GPLv3"
__version__ = "0.3"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__+  " by " + __author__ + " <" + __email__ + ">"

def convertFile(fileLoc, intype, outfile, outtype):

    if outfile != None:
        count = SeqIO.convert(fileLoc, intype, outfile, outtype)
        print "Converted %i records" % count
    else:
       seq_record = SeqIO.read(fileLoc, intype)
       print seq_record.format(outtype)

def main ():
    #TODO: Auto checking of file format, so it doesn't have to be Specified 
    global args
    convertFile(args.input, args.inFormat, args.output, args.outFormat)

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[1].strip()
        parser = argparse.ArgumentParser(description=desc,epilog=epi)
        parser.add_argument ('-v', '--verbose', action='store_true', \
                default=False, help='verbose output')
        parser.add_argument('--version', action='version', \
                version='%(prog)s ' + __version__)

        # VARIABLE ARGUMENTS 
        parser.add_argument('-o','--output',action='store', \
                help='output file location [Default: print to STDOUT]')

        # POSITIONAL ARGUMENT
        parser.add_argument ('inFormat', action='store', \
                choices=['fasta', 'genbank', 'embl'],    \
                help='Input file format either fasta, genbank or embl')
        parser.add_argument ('input', action='store', help='Input File')
        parser.add_argument ('outFormat', action='store',\
                choices=['fasta', 'genbank', 'embl'],    \
                help='Output Format either fasta, genbank or embl')

        args = parser.parse_args()
        if args.verbose: print "Executing @ " + time.asctime()
        main()
        if args.verbose: print "Ended @ " + time.asctime()
        if args.verbose: print 'total time in minutes:',
        if args.verbose: print (time.time() - start_time) / 60.0
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


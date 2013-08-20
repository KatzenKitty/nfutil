#!/usr/bin/env python
"""
# Created: Wed, 14 Dec 2011 14:15:06 +1000

TODO: Single line description. 

TODO: Detailed description
"""
import sys, os, traceback, optparse
import time
import re
#from pexpect import run, spawn
__author__ = "Nabil-Fareed Alikhan"
__licence__ = "GPLv3"
__version__ = "0.3"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"

def main ():

    global options, args
    path = options.dir
    # For each folder in specified folder
    dirList = os.listdir(path) 
    twinpeaks = {}
    for fname in dirList:
        lib = os.path.abspath(path) +  "/" + fname + "/21/LibraryStatistics.txt"
        if os.path.exists( lib ):
            # Open kmer 21 folder /21  
            # Open LibraryStatistics.txt 
            stats = fname + '\t'
            f = open(lib, 'r')
            peakZ = 0
            peak0Mean = 0
            peak0Sdev = 0
            curPeak = [ ]
            for line in f:
#                print line
                if re.search( "Peak" , line ):
                    if re.search( "Peak 0" , line ):    
                        peakZ = 0
                        stats += str(peakZ) + '\t'
                    else:
                        peakZ += 1
                        stats += str(peakZ) + "\t"
                # Find insert size peaks
                # Get peak 0 details
                # If more than 1, add to twinpeaks array 
                if re.search ("AverageOuterDistance: \d+", line):
                    m = re.match("\s+AverageOuterDistance: (\d+)", line)
                    stats += m.group(1) + "\t" 
                    curPeak.append( m.group(1))
                    if peakZ == 0:
                        peak0Mean = m.group(1)
                        
                if re.search ("StandardDeviation: \d+", line):
                    m = re.match("\s+StandardDeviation: (\d+)",line)
                    stats += m.group(1) + "\t"
                    curPeak.append( m.group(1))
                    if peakZ == 0:
                        peak0Sdev = m.group(1)
                    if peakZ == 1:
                        twinpeaks[fname] = curPeak

            f.close()
            print stats                    
    # Twinpeaks holds assemblies and insert sizes with more than one peak
    # modify to run with strict insert size
    for key in twinpeaks.iterkeys():
        pbspath = os.path.abspath( options.pbs ) +  "/ray_" + key + ".pbs"
        newpbs = ''
        val = twinpeaks[key] 
        # modify to run with strict insert size
        if os.path.exists( pbspath ):
            f = open( pbspath, 'r')
            for line in f:
                if re.match('mpirun', line):
                    newpbs += line.replace('-k', str(val[0]) + " " + str(val[1]) + " -k")
                else:     
                    newpbs += line
            #output to folder 
            o = open(options.out + "/ray_" + key + ".pbs", 'w')
            o.write(newpbs)
            if options.runjob:
                os.system("qsub " + options.out + "/ray_" + key + ".pbs")

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[-2]
        parser = optparse.OptionParser(epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        parser.add_option ("-i", "--dir", action="store", type="string", dest="dir", help='Directory of Ray assemblies') 
        parser.add_option ("-p", "--pbs", action="store", type="string", dest="pbs", help='Directory of PBS scripts to modify')
        parser.add_option ("-o", "--out", action="store", type="string", dest="out", help='Directory of modified PBS scripts')
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option ('-r', '--runjobs', action="store_true", default=False, dest="runjob", help='automatically launch job')        
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


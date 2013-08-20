#!/usr/bin/env python
"""
# Created: Mon, 21 Nov 2011 13:45:58 +1000

Ray Assembler Pipeline script

TODO: Detailed description
"""
import sys, os, traceback, optparse
import time
import re
#from pexpect import run, spawn

#METADATA
__author__ = "Mitchell Stanton Cook & Nabil-Fareed Alikhan"
__licence__ = "?"
__version__ = "0.3"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"

def main ():
    #PBS DEFAULT PARAMETERS
    PWD         =  os.getcwd()
    MEM         = '22'
    NODE        = 'NodeType=medium'
    TIME        = '48' 
    OUTD        = '$PBS_O_WORKDIR'
    READD       = ''
    # READ FILTERING PARAMETERS
    ILL2VEL     = '/work2/UQ/Science/SCMB/Beatson/scripts/ill2vel.py'
    ILL2VEL33   = '/work2/UQ/Science/SCMB/Beatson/scripts/ill_new_2vel.py'
    MINLENGTH   = '50'
    MINQUAL     = '20'
    ENC         = 'Q64' 
    # RAY PARAMETERS
    KMERSTART   = '17'
    KMERSTOP    = '35' 

    global options, args
    if (options.filelist == None or not os.path.isfile(options.filelist)):
        print "Please specify a valid filelist"
        sys.exit(1)
    # CONFIGURE GLOBAL OPTIONS
    if (options.startkmer != None):
        KMERSTART = str(options.startkmer)
    if (options.endkmer != None):
        KMERSTOP = str(options.endkmer)
    if (options.walltime != None):
        TIME = str(options.walltime)
    if (options.memory != None):
        MEM = str(options.memory)
        if int(MEM) > 22:
            NODE = 'NodeType=large'
    if (options.encode != None):            
        ENC = str(options.encode)
    if (options.minqual != None):
        MINQUAL = str(options.minqual)
    if (options.minlength != None):
        MINLEN = str(options.minlength)
    if (options.output != None and os.path.isdir(options.output) ):
        OUTD = os.path.abspath(options.output) + '/'
    if (options.path != None and os.path.isdir(options.path) ):
        READD = os.path.abspath(options.path) + '/'
    flist   = options.filelist
    indat   = open(flist).readlines()

    if options.verbose: 
        print 'Running file ' + flist
        print 'Options: KMERSTART ' + KMERSTART + ', KMERSTOP ' + KMERSTOP + ', WALLTIME ' + TIME + ', MEMORY ' + MEM + ', NODE ' + NODE + ', ENCODING ' + ENC + ', MINQUALITY ' + MINQUAL + ', MINLENGTH ' + MINLENGTH +', OUTPUTDIR ' + OUTD + ', READDIR ' + READD
    i = 0
    while i < len(indat):
        vars = indat[i].split('_') 
        F1 = READD+str(indat[i]).strip('\n')
        F2 = READD+str(indat[i+1]).strip('\n')

        ##Write the submission script
        fout = open('ray_'+vars[1]+'.pbs', 'w')
        fout.write('#!/bin/bash\n')
        fout.write('#PBS -A uq-Beatson\n')
        # Could modify mpiprocs here for larger CPU grab
        fout.write('#PBS -l select=1:ncpus=8:mpiprocs=8:mem='+MEM+'gb:'+NODE+'\n')
        fout.write('#PBS -N rayass_'+vars[1]+'\n')
        fout.write('#PBS -l walltime='+TIME+':00:00\n')
        fout.write('\n')
        fout.write('#LOADING REQUIRED MODULES \n')
        fout.write('module load OpenMPI/1.5.3\n')
        fout.write('module load compiler/gcc-4.5.2\n')
        fout.write('\n')

        ##Setup preprocessing
        fout.write('#PRE-PROCESSING \n')
        fout.write('cp ' + F1 + ' $TMPDIR\n')
        fout.write('cp ' + F2 + ' $TMPDIR\n') 
        fout.write('cd $TMPDIR \n')
        # HANDLE COMPRESSION
        OFF = 0
        if F1[len(F1)-2:] == "z2":
            fout.write('bunzip2 '+F1.split("/")[-1]+'\n')
            fout.write('bunzip2 '+F2.split("/")[-1]+'\n')
            OFF = -4
        elif F1[len(F1)-2:] == "ip":
            fout.write('unzip '+F1.split("/")[-1]+'\n')
            fout.write('unzip '+F2.split("/")[-1]+'\n')
            OFF = -4
        elif F1[len(F1)-2:] == "gz":
            if F1[len(F1)-2:] == "tar.gz":
                fout.write('tar -zxvf '+F1.split("/")[-1]+'\n')
                fout.write('tar -zxvf '+F2.split("/")[-1]+'\n')
                OFF = -6
            else:
                fout.write('gunzip '+F1.split("/")[-1]+'\n')
                fout.write('gunzip '+F2.split("/")[-1]+'\n')
                OFF = -3
        else:
            print 'Assuming uncompressed'

        if ENC == 'Q33':
            fout.write('python '+ILL2VEL33+' '+F1.split("/")[-1][:OFF]+' '+F2.split("/")[-1][:OFF]+' '+vars[1]+'_proc.fastq'+' '+MINLENGTH+' '+MINQUAL+'\n')
        else:
            fout.write('python '+ILL2VEL+' '+F1.split("/")[-1][:OFF]+' '+F2.split("/")[-1][:OFF]+' '+vars[1]+'_proc.fastq'+' '+MINLENGTH+' '+MINQUAL+'\n')
        fout.write('\n')
    
        ##Carry out computations
        fout.write('#LAUNCH THE CALCULATION \n')
        fout.write('NP=`cat $PBS_NODEFILE | wc -l`\n')
        kmercount = int(KMERSTART)
        kmerdirs = [] 
        while ( kmercount <= int(KMERSTOP) ): 
            fout.write('mpirun -np $NP Ray -i '+vars[1]+'_proc.fastq -k '+ str(kmercount)+ '  -o ' +str(kmercount) +' -show-distance-summary\n')
            kmerdirs.append(str(kmercount))
            kmercount += 2 
        fout.write('\n\n')

        ##Do postprocessing
        fout.write('#POST-PROCESSING \n')
        #outl = F1.split('/')
        fout.write('cd '+ OUTD +'\n')
        fout.write('mkdir RAYassembly \n')
        fout.write('cd RAYassembly \n')
        fout.write('mkdir '+vars[1]+'\n')
        fout.write('cd '+vars[1]+'\n')
        ##In reads done
        sumString = ''
        for k in kmerdirs:
            fout.write('cp -r $TMPDIR/' + k + '  . \n')
            sumString += k +'/OutputNumbers.txt '
        # SUMMARY STATS
        fout.write( 'python /home/uqnalikh/Beatson_shared/scripts/extract_ray_results.py ' + sumString + ' > SUMMARY.txt\n') 
        fout.write('cp -r $TMPDIR/' + vars[1]+'_proc.fastq . \n')
        fout.close()
        if options.runjob:
            os.system("qsub ray_"+vars[1]+'.pbs')
        i=i+2

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[-2]
        parser = optparse.OptionParser(epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option ("-f", "--filelist", action="store", type="string", dest="filelist", help='filelist of reads [Required]' )
        parser.add_option ('-s', '--kmerstart', action="store", type="int",dest="startkmer", help='start of kmer range')
        parser.add_option ('-e', '--endstart', action="store", type="int",dest="endkmer", help='end of kmer range')
        parser.add_option ('-o', '--output', action="store", type="string",dest="output", help='path to output directory [optional]')
        parser.add_option ('-p', '--path', action="store", type="string",dest="path", help='optional path to append to filenames in filelist')
        pbsgroup = optparse.OptionGroup(parser, "PBS Options", "")
        pbsgroup.add_option ('-w', '--walltime', action="store", type="int",dest="walltime", help='walltime (hours)')
        pbsgroup.add_option ('-m', '--memory', action="store", type="int",dest="memory", help='allocated memory (G)')
        pbsgroup.add_option ('-r', '--runjobs', action="store_true", default=False, dest="runjob", help='automatically launch job')
        parser.add_option_group(pbsgroup)
        filtgroup = optparse.OptionGroup(parser, "Read filtering Options", "")
        filtgroup.add_option ('-q', '--quality', action="store", type="string",dest="encode", help='read quality encoding')
        filtgroup.add_option ('-u', '--minqual', action="store", type="int",dest="minqual", help='minimum read quality')
        filtgroup.add_option ('-l', '--minlength', action="store", type="int",dest="minlength", help='minimum read length')
        parser.add_option_group(filtgroup)
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


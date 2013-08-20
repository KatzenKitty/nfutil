#!/usr/bin/env python
"""
# Created: Tue, 25 Oct 2011 13:19:03 +1000

%prog calculates a range of assembly statistics for a velvet assembly. 

TODO: Detailed description
"""
import sys
sys.path.append('/home/uqnalikh/lib/matlibplot/')
import os
import matplotlib as matp
matp.use('Agg')
import matplotlib.pyplot as plt
import shutil
import matplotlib.mlab as mlab
import traceback, optparse
import time
import re
import numpy as np
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

#from pexpect import run, spawn
__author__ = "Nabil-Fareed Alikhan"
__licence__ = "GPLv3"
__version__ = "0.7"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"

insertCutoff = 700
coverageCutoff = 600 

### CHANGELOG ###
# v0.3: Intial build
# v0.4: Added Support for RAY
# v0.5: Improved RAY support
# v0.6: Contigs stats mode for RAY assembly
# v0.7: Fixed bug with number of N-runs
### CHANGELOG ###

def main ():
#Final assembly parameters for all
#    Read location
#    All contigs
#    All protein
#    Number of contigs
#    N50 N90
#   For each assembly:
#       #Access to contigs.fa
#       #Velvet parameters
#       #Number of contigs, Number of contigs > 1kb :> contigs.fa
#       #N50, N75 N90 :> contigs.fa
#       #Average contig size, largest contig :> contigs.fa
#       #Number of Ns & distribution of Ns. :> contigs.fa
#       #Read coverage for contigs :>  stats.txt
#       #Observed insert size :> Sequence & Graph2 
#       Graph A50 + coverage depth
#       #Nmer graph :> SampleAssemblies.txt
#       &Mauve statistics (?) :> 
    global options, args
    velvetDir = options.dir
    outputDir = options.out
    if options.name != None:
        name = options.name
    else:
        if options.dir.split("/")[-1] == None or options.dir.split("/")[-1] == '':
            name = options.dir.split("/")[-2]
        if options.dir.split("/")[-2] == '.':
            name = options.dir.split("/")[-1]
        else:
            name = options.dir.split("/")[-2] + "-" + options.dir.split("/")[-1]
    assemble = 'Velvet'
    if options.ray:
        assemble = 'Ray'
    finalOut = outputDir  + "/" + name
    if not os.path.exists(finalOut):
        os.makedirs(finalOut)
    html = open(finalOut +  '/index.html', 'w')   
    dat = open(finalOut + '/stats.dat','w')
    SCAFNAME = 'Contigs'
    head = "name\tnmer\tvelvet version\tnumber of contigs\tcontigs > 1kb\tave. contig length\t    largest contig\tn50\tn75\tn90\tnumber of Ns\tnumber of N runs\ttotal length of contigs with N\tnumber of contigs with Ns\tmedian coverage depth\tLog:nodes\tLog:N50\tNumber of filtered reads\tNumber of submitted reads\tMedian observed insert\tMode observed insert\tobserved insert sdev\ttotal bases assembled"
    allstats = [0] *  (len(head.split("\t")  ) +1 ) 
    allstats[0] = name 
    # Write Html header
    html.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n')
    html.write('<html>\n')
    html.write('\t<head>\n')
    html.write('\t\t<link rel="stylesheet" href="../styles.css" />')
    html.write('\t</head>\n')
    html.write('\t<body>\n')
    html.write('\t\t<div id="conteneur">\n')
    # The main header
    html.write('\t\t\t<div id="header">' + assemble + ' assembly information for ' + name + '</div>\n')
    # Configure the html menu here
    html.write('\t\t\t<div id="haut">\n')
    html.write('\t\t\t\t<ul class="menuhaut">\n')
    html.write('\t\t\t\t\t<li><a href="../">Back</a></li>\n')
    html.write('\t\t\t\t</ul>\n')
    html.write('\t\t\t</div>\n')
    # Begin HTML Content/results
    html.write('\t\t<div id="centre">\n')
    html.write('\t\t<b>Multi-FASTA file of scaffolds, produced by ' + assemble )
    kmer = 0 
    MedCov = '' 
    vers = ''
    finalGraph = '' 
    if options.ray: 
        html.write(' <a href="Scaffolds.fasta">[Scaffolds.fasta]</a></b><br>\n')
        shutil.copyfile(velvetDir +"/Scaffolds.fasta" , finalOut+ "/Scaffolds.fasta")
        #Ray Parameters
        f = open(velvetDir + "/RayCommand.txt", 'r+')
        param = '' 
        for line in f:
            param += line.replace('\\',' ').strip() + " "
        if param.find('-p') != -1 or param.find('-i') != -1:
            SCAFNAME = 'Scaffolds'
        html.write('<h2>Ray parameters</h2>\n')
        html.write( param+ '<br>\n' )
        f = open(velvetDir + "/RayVersion.txt", 'r+')
        for line in f:
            vers +=  line.strip()
        allstats[2] = vers
        html.write( vers + '<br>\n' )   
        f = open(velvetDir + "/CoverageDistributionAnalysis.txt", 'r+')
        for line in f:
            if line.find('k-mer length:') != -1 :
                kmer = int(re.match('k-mer length:\s+(\d+)', line ).group(1))
                allstats[1] = kmer
            if line.find('PeakCoverage:') != -1 :
                MedCov = re.match('PeakCoverage:\s+(\d+)', line ).group(1)
                allstats[14] = MedCov
        f = open(velvetDir + "/NumberOfSequences.txt", 'r+')
        for line in f: 
            reads = '0' 
            if line.find('\s+NumberOfSequences:') != -1 :
                reads = re.match('\s+NumberOfSequences:\s+(\d+)', line ).group(1)
            allstats[15] = '0' 
            allstats[16] = '0' 
            allstats[17] = reads
            allstats[18] = reads
    else:
        html.write(' <a href="contigs.fa">[contigs.fa]</a></b><br>\n')
        shutil.copyfile(velvetDir +"/contigs.fa" , finalOut+ "/contigs.fa") 
        shutil.copyfile(velvetDir +"/stats.txt" , finalOut+ "/stats.txt")
        # Velvet Parameters:
        f = open(velvetDir + "/Log", 'r+')
        html.write('<h2>Velvet parameters</h2>')
        velh = ''
        velg = ''
        vers = '' 
        for line in f:
            if (line.find("velveth") != -1 ):
                velh = line[line.find("velveth"):len(line)]
                kmer = int(line.split()[2])
                allstats[1] = kmer
            if (line.find("velvetg") != -1 ):    
                velg = line[line.find("velvetg"):len(line)]    
            if line.find("Version") != -1:
                vers = line    
                allstats[2] = vers.strip() 
            if line.find('Median coverage depth') != -1 :
                MedCov = re.match('Median coverage depth = (\d+\.\d+)', line ).group(1)
                allstats[14] = MedCov
            if (line.find('Final graph has ') != -1 ):
                finalGraph = line
                graphmatch = re.match("Final graph has (\d+) nodes and n50 of (\d+).+using (\d+)/(\d+) reads"  , finalGraph)
                allstats[15] = graphmatch.group(1) 
                allstats[16] = graphmatch.group(2)
                allstats[17] = graphmatch.group(3)
                allstats[18] = graphmatch.group(4)
                
        html.write(velh + "<br>")        
        html.write(velg + "<br>")
        html.write(vers + "<br>")
    # GRAPH 1: Contig size, number of ns
    cont = ''
    if (os.path.isfile(velvetDir + "/contigs.fa") == True):
        cont = "contigs.fa" 
    elif ( os.path.isfile(velvetDir + "/"+ name +  "_contigs.fa") == True):
        cont = name + "_contigs.fa"
    elif ( os.path.isfile(velvetDir + "/contigs/contigs_"+ name +  ".fa") == True):
        cont = "/contigs/contigs_" + name + ".fa"
    elif os.path.isfile(velvetDir +"/Scaffolds.fasta") == True:
        cont = 'Scaffolds.fasta'
    if (len(cont) > 0 ):
        ContigLen = []
        contiggreater = 0
        largestContig = 0
        largeName = ''
        contigNLen = []
        ncount = []
        nRuns  = []
        for seq_record in SeqIO.parse(velvetDir +"/"+ cont , "fasta"):
            noN = len(str(seq_record.seq))
            ContigLen.append( noN )
            nRuns.append( re.findall('[Nn]+', str(seq_record.seq) ) )
            if(noN > largestContig):
                largestContig = noN
                largeName = seq_record.id
            Ns = str(seq_record.seq).count("N")
            if noN >= 1000:
                contiggreater = contiggreater + 1
            if ( Ns != 0 ):
                ncount.append(Ns)
                contigNLen.append(noN)
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        ax.set_xlabel(SCAFNAME + ' size (bp)')
        ax.set_ylabel('Number of placeholder characters in the scaffold')
        plt.title( name + ": "+SCAFNAME+" size vs placeholder characters")
        ax.plot(contigNLen, ncount, 'bx')
        plt.savefig(finalOut + "/Contiggraph-" + name +  ".svg")
        totalbase = sum(ContigLen)
        html.write('<h1>' + SCAFNAME +' information</h1>')
        html.write('\nTotal bases assembled: ' + str(totalbase))
        html.write('\n<br>Total Number of '+SCAFNAME+': ' + str(len(ContigLen)))
        html.write('\n<br>Total Number of '+SCAFNAME+' greater than 1kb: ' + str(contiggreater))
        html.write('\n<br>Average '+SCAFNAME+' size: ' + str(sum(ContigLen) / len(ContigLen)) )
        html.write('\n<br>Largest '+SCAFNAME+' size: ' + str(largestContig )  + " [" + largeName +"]")
        allstats[22] = totalbase
        allstats[3] = len(ContigLen)
        allstats[4] = contiggreater
        allstats[5] = sum(ContigLen) / len(ContigLen)
        allstats[6] = largestContig
        ContigLen.sort(reverse=True)
        n50 = 0
        n75 = 0
        n90 = 0
        sumit = 0 
        for j in ContigLen:
            sumit += j
            if sumit >= totalbase / 2 and n50  == 0:
                n50 = j
            elif sumit >= totalbase *0.75 and n75 == 0:
                n75 = j
            elif sumit >= totalbase * 0.9 and n90 == 0:
                n90 = j
        if options.ray == True:
            f = open(velvetDir + "/NumberOfSequences.txt", 'r+')
            numRead = ''
            for line in f:
                if re.match('\s+NumberOfSequences: \d+',line):
                    allstats[18] = re.match('\s+NumberOfSequences: (\d+)',line).group(1)
                    numRead = allstats[18]
                    allstats[17] = '0'
            f.close()
            html.write('<br>Number of reads: '+ numRead)         
        html.write('<br><br>N50: ' + str(n50) )
        html.write('<br>N75: ' + str(n75) )
        html.write('<br>N90: ' + str(n90) )
        html.write('<br><br>Number of Ns: ' + str(sum(ncount) ) + " in " + str(len(nRuns)) + " runs")
        html.write('\n<br>Total length of '+SCAFNAME+' with Ns: ' + str( sum(contigNLen)) )
        html.write('\n<br>Total Number of '+SCAFNAME+' that have Ns: ' + str( len(contigNLen)))
        allstats[7] = n50
        allstats[8] = n75
        allstats[9] = n90
        allstats[10] = sum(ncount)
        allstats[11] = len(nRuns)
        allstats[12] = sum(contigNLen)
        allstats[13] = len(contigNLen)
        html.write('<br>Median coverage depth: ' + MedCov)
        html.write('<br>' + finalGraph + "<br>")
        html.write('\t\t<br><img src="Contiggraph-' + name +'.svg"  width="500"  /><br>\n')
    # GRAPH 2: Insert size graph 
    html.write('\t\t<h1>Insert size information</h1>')
    foundIn = False
    size = [] 
    obs = [] 
    if (os.path.isfile(velvetDir + "/Graph2") == True):
        if ( os.path.isfile(velvetDir + "/insertGraph.csv") == False or options.redo == True):
            os.system("observed-insert-length.pl " + velvetDir + " > " + velvetDir + "/insertGraph.csv")
        f  = open(velvetDir + "/insertGraph.csv", 'r+')
        foundIn = True
        for line in  f:
            try:
                if (int(line.split()[0]) <= insertCutoff):
                    size.append( int(line.split()[0]))
                    obs.append( int (line.split()[2]))
            except ValueError:
                if re.match('Observed', line ) != None:
                    html.write(line + "<br>") 
                    if re.match('Observed median.+length: (\d+)', line ) != None:
                        allstats[19] = re.match('Observed median.+length: (\d+)', line ).group(1)        
                    if re.match('Observed mode.+length: (\d+)', line ) != None:
                        allstats[20] = re.match('Observed mode.+length: (\d+)', line ).group(1)
                    if re.match('Observed .+deviation: (\d+)', line ) != None:
                        allstats[21] = re.match('Observed .+deviation: (\d+)', line ).group(1)
        shutil.copyfile(velvetDir +"/insertGraph.csv" , finalOut+ "/insertGraph.csv")
    elif (os.path.isfile(velvetDir + "/Library0.txt") == True):
        f =  open(velvetDir + "/Library0.txt", 'r+')
        foundIn = True
        for line in f: 
            if (int(line.split()[0]) <= insertCutoff):
                size.append( int(line.split()[0]))
                obs.append( int (line.split()[1]))
        f.close()        
        shutil.copyfile(velvetDir +"/Library0.txt" , finalOut+ "/insertGraph.csv")
    if foundIn:
       fig2 = plt.figure()
       ax = fig2.add_subplot(111)
       plt.title( name + ": Observed insert size")
       ax.set_xlabel('insert size (bp)')
       ax.set_ylabel('Frequency')
       ax.plot(size,obs)
       plt.savefig(finalOut + "/Insertgraph-" + name +  ".svg")
       html.write('\t\t<img src="Insertgraph-' + name +'.svg"  width="500"  /><br>\n')
       html.write('N.B x-axis on this graph has been limited to ' + str(insertCutoff) + "<br>")
       html.write('\t\t<b>Raw data:<a href="insertGraph.csv">[insertGraph.csv]</a></b><br>\n')
    # GRAPH 3: n50 graph
    sample = False
    n50 = [] 
    nmer = [] 
    numContigs = [] 
    aveContigs = []
    datloc = '' 
    if (os.path.isfile(velvetDir + "/sampleallassemblies.txt") == True):
        f = open (velvetDir  + "/sampleallassemblies.txt")
        lines = f.readlines()
        sample = True
        f.close()
        lines.sort()
        datloc = '<b>Raw data:<a href="sampleassemblies.csv">[sampleassemblies.csv]</a></b><br>\n'
        shutil.copyfile(velvetDir +"/sampleallassemblies.txt" , finalOut+ "/sampleassemblies.csv")
        for line in lines:
            if (line.split()[0] != "nmer"):
                nmer.append(float(line.split()[0]))
                n50.append(float(line.split()[5]) )
                numContigs.append(int(line.split()[1]))
                aveContigs.append(int(line.split()[4]))
    elif os.path.isfile(velvetDir + "/SUMMARY.txt") == True:
        f = open(velvetDir + "/SUMMARY.txt")
        lines = f.readlines()
        f.close()
        lines.sort()
        sample = True
        shadow = 17
        datloc = '<b>Raw data:<a href="SUMMARY.csv">[SUMMARY.csv]</a></b><br>\n'
        shutil.copyfile(velvetDir +"/SUMMARY.txt" , finalOut+ "/SUMMARY.csv")        
        for line in lines:
            if (line.split()[0] != "ID"):
                if (line.split()[0] == 'OutputNumbers.txt'):
                    nmer.append(float(shadow))
                else:
                    nmer.append(float(line.split()[0]))
                if options.cont:
                    n50.append(float(line.split()[4]) )
                    numContigs.append(int(line.split()[1]) )
                    aveContigs.append(int(line.split()[3] ) )
                else:
                    n50.append(float(line.split()[14]) )
                    numContigs.append(int(line.split()[11]))
                    aveContigs.append(int(line.split()[13]))
                shadow += 2 
    if sample:
        figA = plt.figure()
        ax = figA.add_subplot(111)
        ax.set_xlabel('nmer'  )
        ax.set_ylabel('number of contigs')
        plt.title( name + ": nmer vs number of contigs")
        plt.grid(True)
        ax.plot(nmer,numContigs)
        plt.savefig(finalOut + "/numContiggraph-" + name +  ".svg")
        html.write('<h1>Nmer and effect on assembly quality</h1>\n')
        html.write("Final assembly used nmer of " + str(kmer) +"<br>")
        html.write('\t\t<img src="numContiggraph-' + name +'.svg"  width="500"  /><br>\n ')
        html.write(datloc)
        figb = plt.figure()
        ax = figb.add_subplot(111)
        ax.set_xlabel('nmer'  )
        ax.set_ylabel('average contig size (bp)')
        plt.title( name + ": nmer vs average contig size")
        plt.grid(True)
        ax.plot(nmer,aveContigs)
        plt.savefig(finalOut + "/aveContiggraph-" + name +  ".svg")
        html.write('\t\t<img src="aveContiggraph-' + name +'.svg"  width="500"  /><br>\n ')
        html.write(datloc)
        fig3 = plt.figure()
        ax = fig3.add_subplot(111)
        ax.set_xlabel('nmer'  )
        ax.set_ylabel('n50 (bp)')
        plt.title( name + ": nmer vs n50")
        plt.grid(True)
        ax.text(0.05, 0.05, "Final assembly used nmer of " + str(kmer), transform=ax.transAxes, fontsize=10, verticalalignment='top')
        ax.plot(nmer,n50)
        plt.savefig(finalOut + "/n50graph-" + name +  ".svg")
        html.write('\t\t<img src="n50graph-' + name +'.svg"  width="500"  /><br>\n ')
        html.write(datloc)
    # GRAPH 4: Coverage cutoff
    if os.path.isfile(velvetDir + "/stats.txt") == True :
        if ( os.path.isfile(velvetDir + "/estcov.csv") == False or options.redo == True):
            os.system("velvet-estimate-exp_cov.pl " + velvetDir + "/stats.txt > " + velvetDir + "/estcov.csv" )
        f = open(velvetDir + "/estcov.csv", 'r+') 
        lines = f.readlines()
        f.close()
        lines.sort()
        freq = []  
        coverage =  []  
        for line in lines:
            if (line.split()[0] != "Predicted" and line.split()[0] != "velvetg" and float(line.split()[0]) <= coverageCutoff):
                coverage.append(float(line.split()[0]))
                freq.append(float(line.split()[2]) )
        fig4 = plt.figure()
        ax = fig4.add_subplot(111)
        ax.set_xlabel('Coverage')
        ax.set_ylabel('Frequency (Number of contigs)')
        plt.title( name + ": Average read coverage for contigs")
        plt.grid(True)
        ax.plot(coverage,freq)
        plt.savefig(finalOut + "/Coveragegraph-" + name +  ".svg")
        html.write('<h1>Read coverage</h1>\n')
        html.write('\t\t<img src="Coveragegraph-' + name +'.svg"  width="500"  /><br>\n')
        shutil.copyfile(velvetDir +"/estcov.csv" , finalOut+ "/estcov.csv")
        html.write('N.B x-axis on this graph has been limited to ' + str(coverageCutoff) + "<br>")
        html.write('\t\t<b>Raw data:<a href="estcov.csv">[estcov.csv]</a></b><br>\n')
        html.write('\t\t<b>Also see:<a href="stats.txt">[stats.txt]</a></b><br>\n')
    html.write('\t</div>\n\t<div id="pied">Nabil-Fareed Alikhan. 2011. n.alikhan@uq.edu.au<br>\n')
    html.write('CSS Design by Nicolas Fafchamps</div></div>\n')
    html.write('\t</body>\n</html>\n')
    # Find number of Reads for RAY
    if options.ray == True:
        f = open(velvetDir + "/NumberOfSequences.txt", 'r+')
        for line in f:
            if re.match('\s+NumberOfSequences: \d+',line):
                allstats[18] = re.match('\s+NumberOfSequences: (\d+)',line).group(1)
                allstats[17] = '0'
        f.close()
    html.close()
    dat.close()
    outstr = ''
    if os.path.isfile( outputDir + '/summary.csv' )  == False or options.clean == True:
        summary = open(outputDir + '/summary.csv','w')
        summary.write( head + "\n")
        outstr += head + "\n"
        summary.close()
    summary = open(outputDir + '/summary.csv','a')
    for el in allstats:
        summary.write( str(el) +'\t'   );
        outstr += str(el) + '\t'
    print outstr + '\n'
    summary.write('\n');
    summary.close()

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[-2]
        parser = optparse.OptionParser(epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option("-f", "--folder", action="store", type="string", dest="dir", help='Location of Velvet assembly (or Ray with -y flag)')
        parser.add_option("-o", "--output", action="store", type="string", dest="out", help='HTML output folder')
        parser.add_option("-r","--redo", action='store_true',default=False, dest="redo", help='Re-run insert size & exp coverage script')
        parser.add_option("-y","--ray", action='store_true',default=False, dest="ray", help='Folder is a Ray Assembly')
        parser.add_option("-t","--contigs",action='store_true',default=False,dest="cont",help='Use Contigs assembly stats')
        parser.add_option("-c","--clean",action='store_true',default=False, dest="clean", help='clean output directory')
        parser.add_option("-n","--name",action='store', type="string", help='name of assembly')
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


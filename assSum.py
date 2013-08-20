#!/usr/bin/env python
"""
# Created: Thu, 27 Oct 2011 15:19:32 +1000

TODO: Single line description. 

TODO: Detailed description
"""
import sys, os, traceback, optparse
import time
import re
import csv
import sys
sys.path.append('/home/uqnalikh/matplotlib/build/lib.linux-x86_64-2.7/')
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


### CHANGELOG ###
# v0.3: Intial build
# v0.4: Fixed titles, Added runs from N chart 
### CHANGELOG ###

__author__ = "Nabil-Fareed Alikhan"
__licence__ = "GPLv3"
__version__ = "0.4"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"

def main ():

    global options, args
    wid = len(readData('name') ) * 0.7
    bargraph("N50",readData('name'), readData('n50'), 'Strain', 'n50',  options.output + "/n50" , 8,wid , [] ,  '', '' )
    bargraph("Number of reads", readData('name'),readData('Number of filtered reads'),'Strain', 'Number of reads',  options.output + "/filtReads", 8,wid, readData('Number of submitted reads'), 'Filtered reads', 'Total reads'   ) 
    bargraph("Number of contigs",readData('name'), readData('number of contigs'), 'Strain', 'Number of contigs',  options.output +     "/contigs" , 8,wid , readData("contigs > 1kb") ,  'Total contigs', 'Contigs > 1kb' )    
    bargraph("Average contig length",readData('name'), readData('ave. contig length'), 'Strain', 'Average contig length',  options.output +     "/avecontigs" , 8,wid , [] ,  '', '' )
    bargraph("Number of scaffolding characters",readData('name'), readData('number of Ns'), 'Strain', 'Number of scaffolding characters',  options.output +     "/NsNo" , 8,wid , [] ,  '', '' )  
    bargraph("Number of N-gaps",readData('name'), readData('number of N runs'), 'Strain', 'Number of N-gaps',  options.output +     "/NsRuns" , 8,wid , [] ,  '', '' )
    html = open(options.output +  '/index.html', 'w')
    html.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n')
    html.write('<html>\n')
    html.write('\t<head>\n')
    html.write('\t\t<link rel="stylesheet" href="./styles.css" />')
    html.write('\t</head>\n')
    html.write('\t<body>\n')
    html.write('\t\t<div id="conteneur">\n')
    # The main header
    html.write('\t\t\t<div id="header">Assembly summary page</div>\n')
    # Configure the html menu here
    html.write('\t\t\t<div id="haut">\n')
    html.write('\t\t\t\t<ul class="menuhaut">\n')
    html.write('\t\t\t\t\t<li><a href="./summary.csv">All statistics (csv)</a></li>\n')
    html.write('\t\t\t\t</ul>\n')
    html.write('\t\t\t</div>\n')
    # Begin HTML Content/results
    html.write('\t\t<div id="centre">')
    html.write('<h1>Links to individual summary pages</h1>')
    quotes=open( options.input , "rb" )
    csvReader= csv.DictReader( quotes ,  delimiter='\t')
    c = 1 
    for data in csvReader:
        html.write(  '<a href="' +  data['name'] + '">' + data['name']+ '</a>')
        if ( c % 10 == 0):
            html.write('<br>\n')
        else:
            html.write('&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;')
        c = c+1 
    html.write('<h1>N50</h1>')
    html.write('\t\t<br><a href="n50.svg"><img src="n50.svg"  width="700"/></a><br>\n')
    html.write('<h1>Average contig length</h1>')
    html.write('\t\t<br><a href="avecontigs.svg"><img src="avecontigs.svg"  width="700"/></a><br>\n')
    html.write('<h1>Number of contigs</h1>')
    html.write('\t\t<br><a href="contigs.svg"><img src="contigs.svg"  width="700"/></a><br>\n')
    html.write('<h1>Number of Scaffolding placeholder characters (number of Ns)</h1>')
    html.write('\t\t<br><a href="NsNo.svg"><img src="NsNo.svg"  width="700"/></a><br>\n')
    html.write('<h1>Number of N-gaps (runs of Ns)</h1>')
    html.write('\t\t<br><a href="NsRuns.svg"><img src="NsRuns.svg"  width="700"/></a><br>\n')
    html.write('<h1>Number of reads</h1>')
    html.write('\t\t<br><a href="filtReads.svg"><img src="filtReads.svg"  width="700"/></a><br>\n')
    html.write('\t</div>\n\t<div id="pied">Nabil-Fareed Alikhan. 2011. n.alikhan@uq.edu.au<br>\n')
    html.write('CSS Design by Nicolas Fafchamps</div></div>\n')
    html.write('\t</body>\n</html>\n')
    html.close()

def readData(col):
    quotes=open( options.input , "rb" )
    csvReader= csv.DictReader( quotes ,  delimiter='\t')
    doop = []
    for data in csvReader:
        doop.append( data[col]  )
    return doop 

def bargraph(title, xcol, ycol,xaxis, yaxis, outFile , height, width, ycol2, leg1,leg2):
    fig1 = plt.figure( figsize=(width,height) )
    ax = fig1.add_subplot(111, autoscale_on=True)
    ax.set_adjustable("datalim")
    ax.set_xlabel( xaxis )
    barw = 0.55
    if len(ycol2) > 0: 
        barw = 0.35
    ind = np.arange(len(xcol) ) 
    ax.set_ylabel( yaxis )
    ax.set_title( title )
    bar1 = ax.bar(ind ,  map(int, ycol), barw, color='r' )
    if len(ycol2) > 0: 
        bar2 = ax.bar(ind + barw ,  map(int, ycol2), barw, color='y' )
        ax.legend( (bar1[0], bar2[0]), ( leg1, leg2) )
    plt.xticks( ind , xcol , size ='x-small')
    plt.savefig(outFile + ".svg")

def stackedbar(title, xlabel,  xcol, ycol,xaxis, yaxis, outFile , height, width):
    print ''
    
    
    
    
if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[-2]
        parser = optparse.OptionParser(epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option ('-i', '--input', action='store', type='string', help='input csv')
        parser.add_option ('-o', '--output', action='store', type='string', help='output folder')
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


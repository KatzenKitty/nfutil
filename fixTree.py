#!/usr/bin/env python
"""
# Created: Thu, 21 Jun 2012 14:03:49 +1000
TODO: Single line description. 

TODO: Detailed description

### CHANGE LOG ### 
YYYY-MM-DD <name> <address>
    * Initial build
"""
import sys, os, traceback, optparse
import time
from lxml import etree

__author__ = "Nabil-Fareed Alikhan"
__licence__ = "GPLv3"
__version__ = "0.3"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"
USAGE = "%prog [options] svg_file Template_file"

def main ():

    global options, args
    root = etree.parse(args[0]).getroot()
    out = open('edit' + args[0],'w')
    f = open(args[1],'r')
    convert = {}
    for line in f.readlines():
        convert[line.split('|')[0].strip()] = line.split('|')[1:]
    for elem in root.iter():
        if str(elem.tag).endswith('text'):
            if convert.has_key(elem.text):
                
                ent = convert[elem.text]
                doo = elem.text
                elem.text = ent[0].strip()
                print  doo + '\t' + elem.text
                if len(ent) >= 2:
                    elem.attrib['fill'] = ent[1].strip()
    out.write(etree.tostring(root, pretty_print=True))

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[1].strip()
        parser = optparse.OptionParser(usage=USAGE, epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        # EXAMPLE OF BOOLEAN FLAG: verbose
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        # EXAMPLE OF command line variable: 'output'
        parser.add_option('-o','--output',action='store',type='string',help='output prefix')
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


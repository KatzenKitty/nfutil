#!/usr/bin/env python
"""
# Created: Tue, 17 Apr 2012 13:47:47 +1000
TODO: Single line description. 

TODO: Detailed description

### CHANGE LOG ### 
YYYY-MM-DD <name> <address>
    * Initial build
"""
import sys, os, traceback, optparse
import time
import markdown
import tarfile

__author__ = "Nabil-Fareed Alikhan"
__licence__ = "GPLv3"
__version__ = "0.3"
__email__ = "n.alikhan@uq.edu.au"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"
USAGE = "%prog [options] arg1 args2"

def main ():

    global options, args
    if len(args) < 1:
        print 'Valid directory location required'
        exit(1)
    if not os.path.isdir(args[0]):
        print 'Not a valid directory. Exiting'
        exit(1)
    DIR = args[0]
    # Iterate through all directories:
    for path, dirs, files in os.walk(DIR):
    # If directory has file README
        for f in files:
            if f == "README":
                R =  os.path.join(path, f)
                # Create HTML file. Formating text file
                print 'formatting ' + R 
                readfile = open(R,'r')
                outfile = open(os.path.join(path,'index.html'),'w')
                README = ''
                DESC = {}
                TAR = []
                tabCount = 0 
                for line in readfile.readlines():
                    if line.startswith('%TARGZ'):
                        TAR.append(line.split('|')[1].strip())
                    elif line.startswith('%'):
                        DESC[line.split('|')[0].replace('%','')] = line.split('|')[1].strip()
                    else:
                        README +=  line
                # list directory files as relative links.
                html = markdown.markdown(README, output_format='html4')
                outfile.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"\n \
                        "http://www.w3.org/TR/html4/strict.dtd">\n')
                outfile.write(tabCount*'\t' +'<html>\n')
                tabCount += 1
                outfile.write(tabCount*'\t' +'<head>\n')
                outfile.write(tabCount*'\t' +'</head>\n')
                outfile.write(tabCount*'\t' + '<body>\n')
                tabCount += 1 
                outfile.write(tabCount*'\t' +'Quick links: <a href="#dir">Directory listing</a>\n')
                outfile.write(tabCount*'\t' +'<a href="#file">File listing</a>\n')
                if len(TAR) > 0:
                    outfile.write(tabCount*'\t' +'<a href="#tar">Archive listing</a>\n')
                outfile.write(tabCount*'\t'+html + '\n')
                # Print Directories
                outfile.write(tabCount*'\t' +'<a name="dir"><h2>Directory listing of ' + path + '</h2></a>\n')
                tabCount += 1 
                outfile.write(tabCount*'\t' +'<table border="0" cellpadding="5">\n')
                tabCount += 1 
                outfile.write(tabCount*'\t' +'<tr><td><b>Name</b></td><td><b>Last Modified</b></td><td><b>Description</b></td>\n')
                if path != DIR:
                    outfile.write(tabCount*'\t' +'<tr><td><a href="..">Parent Directory</a></td><td><td/><td></td></tr>\n')
                dirs.sort()
                for g in dirs:
                    fileStrin = '' 
                    fileStrin += '<tr><td><a href="'+g +'">' + g + '</a></td>'
                    fileStrin += '<td>' + str(time.ctime(os.path.getmtime(os.path.join(path, g)))) + '</td>'
                    if DESC.has_key(g):
                        fileStrin += '<td>' + DESC[g] + '</td>'
                    else:
                        fileStrin += '<td></td>'
                    fileStrin += '</tr>\n'
                    outfile.write(tabCount*'\t' +fileStrin)
                tabCount -= 1
                outfile.write(tabCount*'\t' +'</table>\n')
                tabCount -= 1
                # PRINT CREATE TARS
                if len(TAR) > 0:
                    TAR.sort()
                    outfile.write(tabCount*'\t' +'<a name="tar"><h2>Archive listing of ' + path + '</h2></a>\n')
                    tabCount += 1 
                    outfile.write(tabCount*'\t' +'<table border="0" cellpadding="5">\n')
                    tabCount += 1 
                    outfile.write(tabCount*'\t' +'<tr><td><b>Name</b></td><td><b>Last Modified</b></td><td><b>Size</b></td><td><b>Description</b></td>\n')
                    for t in TAR:
                        tarName = "ALL" + t + ".tar.gz"
                        tarLoc = os.path.join(path,"ALL" + t + ".tar.gz")
                        if not os.path.exists(tarLoc) or options.forcetar:
                            print 'creating tar ' + tarLoc
                            with tarfile.open(tarLoc, "w:gz") as tar:
                                for fi in files:
                                    if fi.endswith(t):
                                        tar.add(os.path.join(path,fi) )
                            tar.close();
                        fileStrin = '' 
                        fileStrin += '<tr><td><a href="'+ tarName +'">' + tarName + '</a></td>'
                        gloc = os.path.join(path, g)
                        fileStrin += '<td>' + str(time.ctime(os.path.getmtime(tarLoc))) + '</td>'
                        fileStrin += '<td>' + str(convert_bytes(os.path.getsize(tarLoc)))+ '</td>'
                        fileStrin += '<td>Archive of all .' + t + ' files</td>'
                        fileStrin += '</tr>\n'
                    outfile.write(tabCount*'\t' +fileStrin+'\n')
                    tabCount -= 1 
                    outfile.write(tabCount*'\t' + '</table>\n')
                # Print files
                outfile.write(tabCount*'\t' +'<a name="file"><h2>File listing of ' + path + '</h2></a>\n')
                tabCount += 1 
                outfile.write(tabCount*'\t' +'<table border="0" cellpadding="5">\n')
                tabCount += 1 
                outfile.write(tabCount*'\t' +'<tr><td><b>Name</b></td><td><b>Last Modified</b></td><td><b>Size</b></td><td><b>Description</b></td>\n')
                files.sort()
                for g in files:
                    fileStrin = '' 
                    fileStrin += '<tr><td><a href="'+g +'">' + g + '</a></td>'
                    gloc = os.path.join(path, g)
                    fileStrin += '<td>' + str(time.ctime(os.path.getmtime(gloc))) + '</td>'
                    fileStrin += '<td>' + str(convert_bytes(os.path.getsize(gloc)))+ '</td>'
                    
                    if DESC.has_key(g):
                        fileStrin += '<td>' + DESC[g] + '</td>'
                    elif g.endswith('.gbk') or g.endswith('.embl'):
                        gopen = open(gloc,'r')
                        for line in gopen.readlines():
                            if line.startswith('DE'):
                                fileStrin += '<td>' + line[line.find(' '):].strip() + '</td>'
                                break
                        gopen.close()
                    else:
                        fileStrin += '<td></td>'
                    fileStrin += '</tr>\n'
                    outfile.write(tabCount*'\t' +fileStrin)
                tabCount -= 1 
                outfile.write(tabCount*'\t' +'</table>\n')
                tabCount -= 2 
                outfile.write(tabCount*'\t' +'</body>\n')
                tabCount -= 1 
                outfile.write(tabCount*'\t' +'</html>\n')

def convert_bytes(bytes):
        bytes = float(bytes)
        if bytes >= 1099511627776:
            terabytes = bytes / 1099511627776
            size = '%.2fT' % terabytes
        elif bytes >= 1073741824:
            gigabytes = bytes / 1073741824
            size = '%.2fG' % gigabytes
        elif bytes >= 1048576:
            megabytes = bytes / 1048576
            size = '%.2fM' % megabytes
        elif bytes >= 1024:
            kilobytes = bytes / 1024
            size = '%.2fK' % kilobytes
        else:
            size = '%.2fb' % bytes
        return size

if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[1].strip()
        parser = optparse.OptionParser(usage=USAGE, epilog = epi, formatter=optparse.IndentedHelpFormatter(), description=desc,  version='%prog v' + __version__)
        # EXAMPLE OF BOOLEAN FLAG: verbose
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option ('-t', '--forcetar', action='store_true', default=False, help='Force overwrite of tars')
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


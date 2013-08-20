#!/usr/bin/env python
"""
# Created: Mon, 11 Mar 2013 15:35:30 +1000

TODO: Single line description. 

TODO: Detailed description

### CHANGE LOG ### 
2013-03-11 Nabil-Fareed Alikhan <n.alikhan@uq.edu.au>
    * Formatted to use NF python template & argparse
"""
import sys, os, traceback, argparse
import time
import subprocess
import string
import urllib, urllib2

__author__ = "Mitchell Sullivan"
__licence__ = "??"
__version__ = "0.1"
__email__ = "??"
epi = "Licence: "+ __licence__ +  " by " + __author__ + " <" + __email__ + ">"

# install a custom handler to prevent following of redirects automatically.
class SmartRedirectHandler(urllib2.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers
opener = urllib2.build_opener(SmartRedirectHandler())
urllib2.install_opener(opener);

class contig:
    def __init__(self, name, coverage, forseq, revseq):
        self.name = name
        self.coverage = coverage
        self.to = []
        self.fr = []
        self.forseq = forseq
        self.revseq = revseq
        self.blast = []
        self.placed = False
    def length(self):
        return len(self.forseq)
    def repeat(self):
        if len(self.to) > 1 or self.fr > 1:
            return True


def getShortestPath(contigDict, startContig, endContig, count, wayin, direct, direct2):
    count -= 1
    wayin = wayin + (startContig,)
    if startContig == endContig:
        if direct == direct2:
            return [wayin]
    elif count == -1:
        return None
    else:
        outlist = []
        if direct:
            for i in contigDict[startContig].to:
                if i[0][0] == '-':
                    x = getShortestPath(contigDict, i[0][1:], endContig, count, wayin, False, direct2)
                else:
                    x = getShortestPath(contigDict, i[0], endContig, count, wayin, True, direct2)
                if x != None:
                    for j in x:
                        outlist.append(j)
        else:
            for i in contigDict[startContig].fr:
                if i[0][0] == '-':
                    x = getShortestPath(contigDict, i[0][1:], endContig, count, wayin, False, direct2)
                else:
                    x = getShortestPath(contigDict, i[0], endContig, count, wayin, True, direct2)
                if x != None:
                    for j in x:
                        outlist.append(j)
        if outlist != []:
            return outlist
        else:
            return None

def reverse_seq(sequence):
    transtab = string.maketrans('atgcATGC', 'tacgTACG')
    seq = sequence[::-1]
    return seq.translate(transtab)

def translate_dna(sequence):

    #dictionary with the genetic code
   gencode = {
   'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
   'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
   'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
   'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
   'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
   'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
   'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
   'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
   'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
   'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
   'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
   'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
   'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
   'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
   'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
   'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
   }

   proteinseq = ''
   sequence = string.upper(sequence)
   #loop to read DNA sequence in codons, 3 nucleotides at a time
   for n in range(0,len(sequence),3):
      #checking to see if the dictionary has the key
      if gencode.has_key(sequence[n:n+3]) == True:
         proteinseq += gencode[sequence[n:n+3]]
      else:
         proteinseq += 'X'
   #return protein sequence
   return proteinseq

# takes a reference file and makes a FA and FAA file
def makeRefFiles(filename):
    genbank = open(filename)
    gettrans = 0
    prots = []
    getseq = False
    for line in genbank:
        if line.startswith('     CDS             '):
            gettrans = 1
        elif gettrans == 1 and line.startswith('                     /translation="'):
            prot = line.rstrip()[35:]
            if line.rstrip()[-1] == '"':
                gettrans = 0
                prots.append(prot[:-1])
            else:
                gettrans = 2
        elif gettrans == 2:
            prot += line.rstrip()[21:]
            if prot[-1] == '"':
                prot = prot[:-1]
                gettrans = 0
                prots.append(prot)
        elif line.startswith('ORIGIN'):
            getseq = True
            seq = ''
        elif line.startswith('//'):
            getseq = False
        elif getseq:
            seq += ''.join(line.split()[1:])
    out = open('pipe_temp/ref.faa', 'w')
    count = 0
    for i in prots:
        count += 1
        out.write('>' + str(count) + '\n')
        out.write(i + '\n')
    out.close()
    out = open('pipe_temp/ref.fa', 'w')
    out.write('>ref_fa\n' + seq + '\n')
    out.close()

def orderContigs(filename, graphfile, tag, outDir):
    subprocess.Popen('makeblastdb -dbtype nucl -in pipe_temp/ref.fa -out pipe_temp/ref',
                     shell=True).wait()
    subprocess.Popen('blastn -db pipe_temp/ref -outfmt 6 -evalue 0.001 -num_threads 8\
                     -out pipe_temp/contigs_ref.out -query ' + filename, shell=True).wait()
    blast = open('pipe_temp/contigs_ref.out')
    lastquery = ''
    first = True
    contigOrderList = []
    for line in blast:
        query, subject, ident, length, mismatch, indel, qStart, qEnd, rStart, rEnd, eVal, bitScore = line.split()
        ident = float(ident)
        length = int(length)
        mismatch = int(mismatch)
        indel = int(indel)
        qStart = int(qStart)
        qEnd = int(qEnd)
        rStart = int(rStart)
        rEnd = int(rEnd)
        eVal = float(eVal)
        if query != lastquery:
            if first:
                first = False
            else:
                if len(blasthit) == 1:
                    contigOrderList.append((min(blasthit[0]), blasthit[0][0] < blasthit[0][1], lastquery.split('_')[1]))
            lastquery = query
            blasthit = []
            if length > 200:
                blasthit.append((rStart, rEnd, bitScore))
        else:
            if length > 200 and blasthit != [] and bitScore == blasthit[0][2]:
                blasthit.append((rStart, rEnd, bitScore))
    if len(blasthit) == 1:
        contigOrderList.append((min(blasthit[0]), blasthit[0][0] < blasthit[0][1], lastquery.split('_')[1]))
    contigOrderList.sort()
    contigDict = {}
    graph = open(graphfile)
    getnode = 0
    for line in graph:
        if line.startswith('NODE'):
            node = line.split()[1]
            coverage = int(line.split()[3]) * 1.0 / int(line.split()[2])
            getnode = 1
        elif getnode == 1:
            forseq = line.rstrip()
            getnode = 2
        elif getnode == 2:
            revseq = line.rstrip()
            getnode = 0
            aninstance = contig(node, coverage, forseq, revseq)
            contigDict[node] = aninstance
        elif line.startswith('ARC'):
            arcstart, arcend, mult = line.split()[1:]
            mult = int(mult)
            if arcstart[0] == '-':
                try:
                    contigDict[arcstart[1:]].fr.append((arcend, mult))
                except:
                    pass
            else:
                try:
                    contigDict[arcstart].to.append((arcend, mult))
                except:
                    pass
            if arcstart[0] == '-':
                arcstart = arcstart[1:]
            else:
                arcstart = '-' + arcstart
            if arcend[0] == '-':
                try:
                    contigDict[arcend[1:]].to.append((arcstart, mult))
                except:
                    pass
            else:
                try:
                    contigDict[arcend].fr.append((arcstart, mult))
                except:
                    pass
    graph.close()
    coverageList = []
    for i in contigDict:
        if len(contigDict[i].fr) == 1 and len(contigDict[i].to) == 1:
            coverageList.append(contigDict[i].coverage)
      #      print 'ding'
    covAverage = sum(coverageList) * 1.0 / len(coverageList)
    sdSum = 0
    for i in coverageList:
        sdSum += (i - covAverage) ** 2
    covSD = (sdSum / len(coverageList)) ** 0.5
    contigOrderListTemp = []
    print covAverage, covSD
    for i in contigOrderList:
        if contigDict[i[2]].coverage < covAverage + 3 * covSD:
            contigOrderListTemp.append(i)
    finalOrder = [contigOrderList[0]]
    count = 0
    contigOrderList = contigOrderListTemp
    for i in range(0, len(contigOrderList) -1 ):
        startContig = contigOrderList[i][2]
        endContig = contigOrderList[i+1][2]
        dir = contigOrderList[i][1]
        dir2 = contigOrderList[i + 1][1]
        paths = getShortestPath(contigDict, startContig, endContig, 20, tuple(), dir, dir2)
        maxpath = 100000000
        if paths != None:
          #  print paths
            for q in paths:
                length = 0
                for j in q[1:-1]:
                    length += contigDict[j].length()
                if length < maxpath:
                    path = q
                    maxpath = length
            if maxpath < 10000:
                
                if dir:
                    wayin = True
                else:
                    wayin = False
                for q in range(0, len(path) -1):
                    if wayin:
                        toit = []
                        for j in contigDict[path[q]].to:
                            toit.append(j[0])
                        if path[q + 1] in toit:
                            finalOrder.append((None, True, path[q+1]))
                            wayin = True
                        else:
                            finalOrder.append((None, False, path[q+1]))
                            wayin = False
                    else:
                        frit = []
                        for j in contigDict[path[q]].fr:
                            frit.append(j[0])
                        if path[q + 1] in frit:
                            finalOrder.append((None, True, path[q + 1]))
                            wayin = True
                        else:
                            finalOrder.append((None, False, path[q + 1]))
                            wayin = False
            else:
                finalOrder += ['n', contigOrderList[i + 1]]
        else:
            finalOrder += ['n', contigOrderList[i + 1]]
    finalFNA = ''
    lenlist = []
    gotset = set()
   # print finalOrder
    for i in finalOrder:
        if i == 'n':
            finalFNA += 'n' * 100
            lenlist.append((100, 'n'))
        else:
            if i[1]:
                finalFNA += contigDict[i[2]].forseq
                lenlist.append((len(contigDict[i[2]].forseq), i[2]))
                gotset.add(i[2])
            else:
                finalFNA += contigDict[i[2]].revseq
                lenlist.append((len(contigDict[i[2]].revseq), i[2]))
                gotset.add(i[2])
    fasta = open(filename)
    first = True
    for line in fasta:
        if line.startswith('>'):
            if not first:
                if addseq:
                    lenlist.append((thelen, lastname))
            else:
                first = False
            lastname = line.split('_')[1]
            if line.split('_')[1] in gotset:
                addseq = False
            else:
                addseq = True
                finalFNA += 'n' * 100
                lenlist.append((100, 'n'))
                thelen = 0
        else:
            if addseq:
                finalFNA += line.rstrip()
                thelen += len(line.rstrip())
    if addseq:
        lenlist.append((thelen, lastname))
    out = open(outDir + "/" + tag + '.fna', 'w')
    out.write('>' + tag + '\n')
    for i in range(0, len(finalFNA), 60):
        out.write(finalFNA[i:i+60] + '\n')
    #out = open(tag + '.gb', 'w')
    #pos = 0
    #for i in lenlist:
    #    out.write('     gene            ' + str(pos + 1) + '..' + str(pos + i[0]) + '\n')
    #    out.write('                     /gene=' + i[1] + '\n')
    #    pos += i[0]
    return lenlist



def callGenes(tag, database, outDir):
    subprocess.Popen('build-icm pipe_temp/temp.icm <  ' + database, shell=True).wait()
    subprocess.Popen('glimmer3 ' + outDir +"/"  + tag + '.fna pipe_temp/temp.icm pipe_temp/temp_glim', shell=True).wait()
    glim = open('pipe_temp/temp_glim.predict')
    genelist = []
    for line in glim:
        if not line.startswith('>'):
            genelist.append((int(line.split()[1]), int(line.split()[2])))
    glim.close()
    fastaSeq = ''
    fnaFile = open( outDir +"/"  +tag + '.fna')
    for line in fnaFile:
        if not line.startswith('>'):
            fastaSeq += line.rstrip()
    fnaFile.close()
    out = open('pipe_temp/genes.faa', 'w')
    count = 1
    for i in genelist:
        out.write('>' + str(count) + '\n')
        count += 1
        if i[0] < i[1]:
            out.write(translate_dna(fastaSeq[i[0] - 1:i[1]])[:-1] + '\n')
        else:
            out.write(translate_dna(reverse_seq(fastaSeq[i[1] - 1:i[0]]))[:-1] + '\n')
    out.close()
    return genelist, fastaSeq

def hmmertime(db, seq, stype):
    #print db, seq
    if db == 'pfam':
        parameters = {
                      'hmmdb':db,
                      'seq':seq
                      }
    else:
        parameters = {
                      'seqdb':db,
                      'seq':seq
                      }
    enc_params = urllib.urlencode(parameters);
    #post the seqrch request to the server
    request = urllib2.Request('http://hmmer.janelia.org/search/' + stype, enc_params)
    #get the url where the results can be fetched from
    #print urllib2.urlopen(request).values()
    results_url = urllib2.urlopen(request).values()[7]
    res_params = {
                  'output':'text',
                  'range':'1,1'
                 }
    
    # add the parameters to your request for the results
    enc_res_params = urllib.urlencode(res_params)
    modified_res_url = results_url + '?' + enc_res_params
    # send a GET request to the server
    results_request = urllib2.Request(modified_res_url)
    data = urllib2.urlopen(results_request)
    
    # print out the results
    x = data.read()
    getit = False
    for line in x.split('\n'):
       # print line
        if line.startswith('============================================================================================================'):
            if getit:
                getit = False
            else:
                getit = True
        elif line.startswith('-------------------------------------------------------------------------------'):
            pass
        elif getit:
            splitline = line.split('\t')
            if len(splitline) == 1 or len(splitline) == 0:
                return None
                getit = False
            else:
                return splitline
 #               if float(splitline) < 0.01:
  #                  return splitline
            

def getAnn(tag):
    faa = open('pipe_temp/genes.faa')
    genes = []
    anns = []
    for line in faa:
        if not line.startswith('>'):
            genes.append(line.rstrip())
    for i in range(0, len(genes)):
        db = None
        ann = hmmertime('swissprot', '>test\n' + genes[i], 'phmmer')
        if ann != None and float(ann[3]) > 0.01:
            ann = None
        elif ann != None:
            db = 'swissprot'
        if ann == None:
            ann = hmmertime('nr', 'tset\n' + genes[i], 'phmmer')
        if ann != None and float(ann[3]) > 0.01:
            ann = None
        elif ann != None and db == None:
            db = 'nr'
        if ann == None:
            ann = hmmertime('pfam', 'test\n' + genes[i], 'hmmscan')
        if ann != None and float(ann[3]) > 0.01:
            ann = None
        elif ann != None and db == None:
            db = 'pfam'
        if ann != None:
            ann = 'Matches ' + ann[6].strip() + ' in ' + ann[5] + '. Match found in ' + db + ' database'
            if db == 'pfam':
                ann += ' using hmmscan.'
            else:
                ann += ' using phmmer.'
        else:
            ann = 'No significant matches found.'
        anns.append(ann)
    return anns

def writeGenbank(genes, fastaSeq, anns, tag, conLengths, tag2, outDir):
    out = open(outDir +"/" + tag + '.gbk', 'w')
    out.write('FEATURES             Location/Qualifiers\n')
    out.write('     source          1..' + str(len(fastaSeq)) + '\n')
    pos = 0
    annList = []
    count = 0
    for i in range(0, len(genes)):
        count += 5
        annList.append((tag2 + str(count).zfill(6), anns[i], genes[i][0], genes[i][1]))
    for i in conLengths:
        out.write('     misc_feature            ' + str(pos + 1) + '..' + str(pos + i[0]) + '\n')
        out.write('                     /note="node_' + i[1] + '"\n')
        pos += i[0]
    for i in annList:
        if i[2] < i[3]:
            out.write('     gene            ' + str(i[2]) + '..' + str(i[3]) + '\n')
        else:
            out.write('     gene            complement(' + str(i[3]) + '..' + str(i[2]) + ')\n')
        out.write('                     /locus_tag="' + i[0] + '"\n')
        if i[2] < i[3]:
            out.write('     CDS             ' + str(i[2]) + '..' + str(i[3]) + '\n')
        else:
            out.write('     CDS             complement(' + str(i[3]) + '..' + str(i[2]) + ')\n')
        out.write('                     /locus_tag="' + i[0] + '"\n')
        if i[1] != 'No significant matches found.':
            temp2 = 0
            temp1 = 0
            tempstring = i[1]
            while temp2 != -1:
                cutitoff = temp1
                temp1 = temp2
                temp2 = tempstring.find(' in ', temp1 + 1)
            prodstring = tempstring[8:cutitoff]
            outstring = '/product="' + prodstring + '"'
            outlist = outstring.split()
            count = 21
            ostring = '                     '
            for j in outlist:
                count += len(j) + 1
                if count > 80:
                    ostring = ostring[:-1] + '\n                     '
                    count = 21 + len(j) + 1
                ostring += j + ' '
            ostring += '\n'
            out.write(ostring)
        outstring = '/note="' + i[1] + '"'
        outlist = outstring.split()
        count = 21
        ostring = '                     '
        for j in outlist:
            count += len(j) + 1
            if count > 80:
                ostring = ostring[:-1] + '\n                     '
                count = 21 + len(j) + 1
            ostring += j + ' '
        ostring += '\n'
        out.write(ostring)
        if i[2] < i[3]:
            aastring = translate_dna(fastaSeq[i[2] - 1:i[3]])[:-1]
        else:
            aastring = translate_dna(reverse_seq(fastaSeq[i[3] - 1:i[2]]))[:-1]
        aastring = '/translation="' + aastring + '"'
        for j in range(0, len(aastring), 59):
            out.write('                     ' + aastring[j:j+59] + '\n')
    out.write('ORIGIN\n')
    for i in range(0, len(fastaSeq), 60):
        out.write(' ' * (9 - len(str(i + 1))) + str(i + 1))
        for j in range(i, i + 60, 10):
            out.write(' ' + fastaSeq[j:j+10])
        out.write('\n')
    out.write('//')


    os.mkdir('pipe_temp')

def main ():

    global args
    print 'Hello world!'
    if args.output != None:
        print 'Output: ' + args.output
    print args
        # USAGE: Reference genbank, Contig.fa, LastGraph, Output prefix, Locus tag, Glimmer database for predictions, Output dir
    prefix = ''
    locus = 'GENE'
    outputDir = './'
    if args.prefix != None:
        prefix = args.prefix
    if args.locus != None:
        locus = args.locus
    if args.output != None:
        outputDir = args.output
    makeRefFiles(args.refGBK)
    orderout = orderContigs(args.contigsFa, args.LastGraph, locus, outputDir)
    orderout = []
    genes, fastaFNA = callGenes(locus, args.GlimmerDb, outputDir )
    anns = getAnn(locus)
    writeGenbank(genes, fastaFNA, anns, prefix, orderout, locus, outputDir)


if __name__ == '__main__':
    try:
        start_time = time.time()
        desc = __doc__.split('\n\n')[1].strip()
        parser = argparse.ArgumentParser(description=desc,epilog=epi)
        # EXAMPLE OF BOOLEAN FLAG: verbose
        parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
        parser.add_argument('-p','--prefix',action='store',help='output file prefix')
        parser.add_argument('-o','--output',action='store',help='output directory')
        parser.add_argument('-l','--locus',action='store',help='locus_tag')
        parser.add_argument ('refGBK', action='store', help='Reference Genbank')
        parser.add_argument ('contigsFa', action='store', help='Contigs multi-FASTA file')
        parser.add_argument ('LastGraph', action='store', help='Velvet LastGraph file')
        parser.add_argument ('GlimmerDb', action='store', help='Glimmer database for predictions')
        # USAGE: Reference genbank, Contig.fa, LastGraph, Output prefix, Locus tag, Glimmer database for predictions, Output dir
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


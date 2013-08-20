/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package signalsort;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 *
 * @author nabil
 */
public class Signal implements Runnable {

    static PrintStream OUT;
    final List INPUTFILES;
    final boolean VERBOSE;
    final boolean SIGNAL;
    final String[] LOCALISATION;
    final String[] FOLDERS;
    final String BLAST;

    static void threadMessage(String message) {
        String threadName = Thread.currentThread().getName();
        OUT.print("\n[" + threadName + "] " + message);
    }

    static void threadMessage(String message, boolean tick) {
        if (tick) {
            OUT.print(message);
        }
    }

    public Signal(List inInputFiles, OutputStream inOut, boolean inverbose, String[] inFolders, boolean inSig, String[] inlocal, String inblast) {
        INPUTFILES = inInputFiles;
        OUT = (PrintStream) inOut;
        VERBOSE = inverbose;
        LOCALISATION = inlocal;
        SIGNAL = inSig;
        FOLDERS = inFolders;
        BLAST = inblast;
    }

    @Override
    public void run() {
        // Load Psort results into datastructure
        try {
            threadMessage("Loading pSort results");
            // Create dat folder
            File dat = new File("dat/");
            if (!dat.exists()) {
                dat.mkdir();
            }
            double Ident = 80.0;
            double Len = 80.0;
            double XLen = 80.0;
            double XIdent = 40.0;
            
            pStruct REF = null;
            //  BUILD OR UPDATE PSORTB TABLE DATABASE
            for (int i = 0; i < INPUTFILES.size(); i++) {
                REF = new pStruct(new File(INPUTFILES.get(i).toString()));
            }
            checkPsortDb(REF, dat);

            // --- preprocess folders, create all results --- //
            //  FOR EACH FOLDER CHECK FOR GENBANK OR [TODO: FASTA-NUC FILE ]
            if (FOLDERS != null) {
                for (int i = 0; i < FOLDERS.length; i++) {
                    File folder = new File(FOLDERS[i]);
                    for (File foo : folder.listFiles()) {
                        File fasta = foo;
                        if (foo.getName().endsWith(".gbk")) {
                            threadMessage("Converting " + foo.getName());
                            fasta = GenToFasta(fasta, dat);
                        } else {
                            fasta = new File(dat.getCanonicalPath() + File.separator + foo.getName());
                            copyTo(foo, fasta, false);
                        }
                        //  CHECK EACH AND LOAD FILE POINTER AS ELEMENT IN FOLDER NAME SET

                        // CHECK BLAST RESULTS AND RUN BLAST IF REQUIRED
                        File blastResults = null;
                        if (BLAST.compareTo("blastx") == 0 ) {
                            blastResults = new File(dat + File.separator + REF.fnafasta.getName() + "Vs" + folder.getName() + "-" + fasta.getName() + ".blastx.tab");
                        } else if (BLAST.compareTo("tblastn") == 0 ) {
                            blastResults = new File(dat + File.separator + REF.fnafasta.getName() + "Vs" + folder.getName() + "-" + fasta.getName() + ".tblastn.tab");
                        } else {
                            blastResults = new File(dat + File.separator + folder.getName() + "-" + fasta.getName() + "Vs" + REF.fnafasta.getName() + ".tab");

                        }
                        if (!blastResults.exists()) {
                             if (BLAST.compareTo("blastx") == 0 ) {
                                execThis("makeblastdb -in " + REF.faafasta.getCanonicalPath(), OUT);
                                String exec = "blastx   -evalue 0.00005 num_threads 8  -outfmt 6  -query " + fasta.getCanonicalPath()
                                        + " -db " + REF.faafasta.getCanonicalPath() + " -out " + blastResults.getCanonicalPath();
                                execThis(exec, OUT);
                            } else if (BLAST.compareTo("tblastn") == 0 ) {
                                execThis("makeblastdb -dbtype nucl -in " + fasta.getCanonicalPath(), OUT);
                                String exec = "tblastn -num_threads 8 -evalue 0.00005 -outfmt 6 -query " + REF.faafasta.getCanonicalPath()
                                        + " -db " + fasta.getCanonicalPath() + " -out " + blastResults.getCanonicalPath();
                                execThis(exec, OUT);
                            } else {
                                execThis("makeblastdb -dbtype nucl -in " + fasta.getCanonicalPath(), OUT);
                                String exec = "blastn -num_threads 8 -dust no -evalue 0.00005 -outfmt 6 -query " + REF.fnafasta.getCanonicalPath()
                                        + " -db " + fasta.getCanonicalPath() + " -out " + blastResults.getCanonicalPath();
                                execThis(exec, OUT);
                            }
                        }

                        //  CHECK PRESENCE OF  MATCH FOR EACH FILE. CHECK FILE then pStruct
                        BufferedReader in = new BufferedReader(new FileReader(blastResults));
                        String line = "";

                        while ((line = in.readLine()) != null) {
                            String[] lineArray = line.split("\\s+");
                            if (lineArray.length >= 2) {
                                //  identity > Ident
                                 if (BLAST.compareTo("blastx") == 0 ) {
                                    String blastacc = lineArray[1].split("\\|")[3];
                                    if (blastacc.contains(".")){
                                        blastacc = blastacc.substring(0, blastacc.indexOf("."));
                                    }
                                     if (Double.parseDouble(lineArray[2]) >= XIdent
                                             && (((double) REF.entries.get(blastacc).length * 3.00) / Double.parseDouble(lineArray[3])) * 100 >= XLen) {
                                         pEntry temp = REF.entries.get(blastacc);
                                         temp.addBlastMatch(foo.getCanonicalPath(), line);
                                         REF.entries.put(temp.acc, temp);
                                     }

                                 } else {
                                     String blastacc = "";
                                     if (BLAST.compareTo("tblastn") == 0) {
                                         blastacc = lineArray[0].split("\\|")[3];
                                     } else {
                                         blastacc = lineArray[0].split("\\|")[1];
                                     }
                                     if (blastacc.contains(".")) {
                                         blastacc = blastacc.substring(0, blastacc.indexOf("."));
                                     }
                                     if (Double.parseDouble(lineArray[2]) >= Ident
                                             && ((double) REF.entries.get(blastacc).length / Double.parseDouble(lineArray[3])) * 100 >= Len) {
                                         pEntry temp = REF.entries.get(blastacc);
                                        temp.addBlastMatch(foo.getCanonicalPath(), line);
                                        REF.entries.put(temp.acc, temp);
                                    }
                                }
                            }
                           
                            //  CHECK PRESENCE OF PROTEIN MATCH FOR EACH FILE

                        }
                         in.close();
                    }
                }
            }
        

        //  RUN SIGNALP FOR EACH PROTEIN & GRAB RESULTS
            REF = signalP(REF);
            //  FORMAT RESULTS
            this.writePstruct(REF, new File(REF.location.getCanonicalPath() +".FINAL"));




            threadMessage("Done.\n");
            // Done.
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static void copyTo(File inputFile, File outputFile, boolean binary) throws FileNotFoundException, IOException {
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(inputFile));
        BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(outputFile));
        if( binary){
            bis = new BufferedInputStream(new FileInputStream(inputFile), 4096);
            bos = new BufferedOutputStream(new FileOutputStream(outputFile), 4096);
        }
        int c;
        while ((c = bis.read()) != -1) {
            bos.write(c);
        }
        bis.close();
        bos.close();
    }
    public pStruct signalP(pStruct ref) {

        threadMessage("Running SignalP");
        /*
        File masterSignal = new File(masterFile.getCanonicalPath() + ".signalp");
        masterSignal.createNewFile();

        HashMap signal = new HashMap();
        in = new BufferedReader(new FileReader(masterSignal));
        line = "";
        while ((line = in.readLine()) != null) {
            String[] lineArray = line.split("\\s+");
            signal.put(lineArray[0], lineArray[1]);
        }
        in.close();
        for (pStruct cur : psorts) {
            int count = 1;
            for (pEntry tem : cur.entries) {
                // if exists check if it contains all proteins described in psort files.
                // Fill fields if possible
                if (signal.get(tem.acc) != null) {
                    String sig = signal.get(tem.acc).toString();
                    if (sig.compareTo("YES") == 0) {
                        tem.signalp = true;
                    } else {
                        tem.signalp = false;
                    }

                } else {
                    // if not run signalp.
                    // add result to data structure
                    File temP = new File(dat + File.separator + "temp-SignalP");
                    out = new BufferedWriter(new FileWriter(temP));
                    out.write(">" + tem.acc);
                    out.newLine();
                    out.write(tem.fasta);
                    out.newLine();
                    out.close();
                    String exec = "signalp-3.0/signalp -t gram- -f summary " + temP.getCanonicalPath();
                    Process q = Runtime.getRuntime().exec(exec);
                    InputStream istrm = q.getInputStream();
                    InputStreamReader istrmrdr = new InputStreamReader(istrm);
                    String data;
                    BufferedReader buffrdr = new BufferedReader(istrmrdr);
                    boolean NN = false;
                    int HMMcoun = 0;
                    while ((data = buffrdr.readLine()) != null) {
                        if (data != null) {
                            if (data.contains("Prediction: Signal peptide")) {
                                NN = true;
                            }
                            String[] linArray = data.split("\\s+");
                            if (!data.contains("#")) {
                                if (linArray.length == 6 && linArray[5].compareTo("YES") == 0) {
                                    HMMcoun++;
                                }
                                if (linArray.length == 7 && linArray[6].compareTo("YES") == 0) {
                                    HMMcoun++;
                                }
                            }
                        }
                    }
                }
            }
        }*/
        return ref;
    }

    public void writePstruct(pStruct ref, File output) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(output));
        String Head = "ACCESSION\tLOCUS_TAG\tDescription\tGene-length\tLOCALISATION\tLOCALISATION-SCORE\tKnown-surface-protein\tSIGNAL-PEPTIDE(Sigp-HMM)\tProbability\tSIGNAL-PEPTIDE(Sigp-NN)\tConsensus\tSIGNAL-PEPTIDE(pSort)\tMotif\tTransmembrane-result\tTransmembrane-note\t";
        for (String fold : FOLDERS) {
            Head += new File(fold).getName() + "\tDetails\t";
        }
        out.write(Head);

        out.newLine();
        for (String temp : ref.entries.keySet()) {
            boolean loc = false;
            pEntry tem = ref.entries.get(temp);

            if (LOCALISATION != null) {
                for (String local : LOCALISATION) {
                    if (tem.desc.split("\t")[31].compareToIgnoreCase(local) == 0) {
                        loc = true;
                    }
                }
            }
            if (LOCALISATION == null || loc) {
                if ((SIGNAL && tem.desc.split("\t")[7].compareTo("Signal peptide detected") == 0) || SIGNAL == false) {
                    out.write(tem.acc + "\t");
                    out.write(tem.locus + "\t");
                    // Description
                    out.write(tem.desc.split("\t")[1]+ "\t");
                    // Gene length
                    out.write(tem.length + "\t");
                    // Localisation and score from psort
                    out.write(tem.desc.split("\t")[31] + "\t");
                    out.write(tem.desc.split("\t")[35] + "\t");
                    // Known surface protein
                    out.write(tem.desc.split("\t")[13] + "\t");
                    // Signalp HMM results
                    out.write(tem.signalpHMMResult + "\t");
                    out.write(tem.signalpHMMDesc + "\t");
                    // SignalP NN results
                    out.write(tem.signalpNNResult + "\t");
                    out.write(tem.signalpNNDesc + "\t");
                    // Signal peptide (psort)
                    out.write(tem.desc.split("\t")[7] + "\t");
                    // Motif location
                    out.write(tem.desc.split("\t")[19] + "\t");
                    //Transmembrane Result & note
                    out.write(tem.desc.split("\t")[22] + "\t");
                    out.write(tem.desc.split("\t")[23] + "\t");
                    // ADD BLAST RESULTS

                    for(String fold : FOLDERS ){
                        int count = 0;
                        File folder = new File(fold);
                        String blsDetails = "" ;
                        for (File foo : folder.listFiles()) {
                            if( tem.blastMatch.containsKey(foo.getCanonicalPath())){
                                String bls = tem.blastMatch.get(foo.getCanonicalPath() ).toString() ;
                                count++;
                                 if (BLAST.compareTo("blastx") == 0 ) {
                                blsDetails += bls.split("\t")[0] + ":" + bls.split("\t")[6]+".."+bls.split("\t")[7] + ",";
                                }else{
                                blsDetails += bls.split("\t")[1] + ":" + bls.split("\t")[8]+".."+bls.split("\t")[9] + ",";
                                    
                                }
                            }
                        }
                        out.write(count +  "/" + folder.listFiles().length +"\t"+blsDetails+"\t");
                    }
                    out.newLine();
                }
            }
        }
        /*   for (String temp : ref.entries.keySet()) {
            pEntry tem = ref.entries.get(temp);
            boolean printIt = true;
            if (SIGNAL) {
                if (tem.signalpDesc == false) {
                    printIt = false;
                }
            }
            if (LOCALISATION != null) {
                int hit = 0;
                for (String l : LOCALISATION) {
                    if (tem.desc.split("\t")[31].compareToIgnoreCase(l) == 0) {
                        hit++;
                    }
                }
                if (hit == 0) {
                    printIt = false;
                }
            }
            if (printIt) {
                out.write(tem.acc + "\t");
                out.write(tem.locus + "\t");
                if (tem.signalp) {
                    out.write("YES\t");
                } else {
                    out.write("NO\t");
                }
                if (tem.desc.split("\t")[33].contains("multiple localization sites")) {
                    out.write(tem.desc.split("\t")[31] + "*\t");
                } else {
                    out.write(tem.desc.split("\t")[31] + "\t");
                }
                out.write(tem.desc.split("\t")[1] + "\t");
                out.write(tem.acc + "\t");
                out.newLine();
            }
        }*/
        out.close();

    }

    public static void execThis(String exec, PrintStream str) {
        try {
            threadMessage("Running " + exec);
            Process q = Runtime.getRuntime().exec(exec);
            InputStream istrm = q.getErrorStream();
            InputStreamReader istrmrdr = new InputStreamReader(istrm);
            String data;
            BufferedReader buffrdr = new BufferedReader(istrmrdr);
            while ((data = buffrdr.readLine()) != null) {
                threadMessage(data);
            }
            q.destroy();
        } catch (Exception e) {
            e.printStackTrace();
            threadMessage(e.getMessage());
        }
    }

    public File GenToFasta(File gbk, File wd) throws FileNotFoundException, BioException, IOException {
        BufferedReader in = new BufferedReader(new FileReader(gbk));
        File outFasta = new File(wd + File.separator + gbk.getName().replaceAll(".gbk", "") + ".fna");
        if (!outFasta.exists()) {
            OutputStream faaout = new FileOutputStream(outFasta);
            RichSequenceIterator s = RichSequence.IOTools.readGenbankDNA(in, null);
            while (s.hasNext()) {
                RichSequence zz = s.nextRichSequence();
                RichSequence.IOTools.writeFasta(faaout, zz, null);
            }
        }
        return outFasta;
    }

    public void checkPsortDb(pStruct pSort, File wd) throws Exception {
        pSort.faafasta = new File(wd.getAbsolutePath() + File.separator + pSort.location.getName() + ".faa");
        pSort.fnafasta = new File(wd.getAbsolutePath() + File.separator + pSort.location.getName() + ".fna");
        if (!pSort.faafasta.exists()) {
            // Check for master protein fasta file in dat folder
            try {
                pSort.faafasta.createNewFile();
                pSort.fnafasta.createNewFile();
                OutputStream faaout = new FileOutputStream(pSort.faafasta);
                OutputStream fnaout = new FileOutputStream(pSort.fnafasta);
                threadMessage("Downloading sequences for " + pSort.faafasta.getName());
                String proNames = "";
                int proCount = 0;
                boolean first = true;
                File gbk = null;
                BufferedReader gbkget = null;
                RichSequence zz = null;
                Pattern p = Pattern.compile("\\W*([a-zA-Z]{2}_.+):(\\d+)..(\\d+)");
                for (String gotAcc : pSort.entries.keySet()) {
                    pEntry cur = pSort.entries.get(gotAcc);
                    proNames += cur.acc + ",";
                    // Bundle 50 at a time and fetch from intewebs
                    proCount++;
                    if (proCount % 50 == 0 || proCount == pSort.entries.size()) {
                        threadMessage("Getting " + proNames);
                        String ooo = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=" + "refseqp" + "&id=" + proNames + "&style=raw";
                        URLConnection connection = new URL(ooo).openConnection();
                        BufferedReader cont = new BufferedReader(new InputStreamReader(connection.getInputStream()));
                        RichSequenceIterator ss = RichSequence.IOTools.readGenbankProtein(cont, new SimpleNamespace("ref"));
                        while (ss.hasNext()) {
                            RichSequence zs = ss.nextRichSequence();
                            for (Feature teemp : zs.getFeatureSet()) {
                                if (teemp.getAnnotation().containsProperty("coded_by")) {
                                    String coded = (teemp.getAnnotation().getProperty("coded_by")).toString();
                                    threadMessage(coded);

                                    Matcher m = p.matcher(coded);
                                    if (m.find()) {
                                        String acc = m.group(1);
                                        if (first) {
                                            gbk = GetGBK(acc, wd);
                                            gbkget = new BufferedReader(new FileReader(gbk));
                                            first = false;
                                            RichSequenceIterator s = RichSequence.IOTools.readGenbankDNA(gbkget, null);
                                            while (s.hasNext()) {
                                                zz = s.nextRichSequence();
                                            }
                                        }
                                        int start = Integer.parseInt(m.group(2));
                                        int stop = Integer.parseInt(m.group(3));
                                        Sequence fn = DNATools.createDNASequence(zz.subStr(start, stop),
                                                zs.getAccession() + "|" + zs.getDescription().replaceAll("\n", " "));
                                        RichSequence.IOTools.writeFasta(fnaout, fn, new SimpleNamespace("ref"));
                                    }
                                }
                            }
                            RichSequence.IOTools.writeFasta(faaout, zs, new SimpleNamespace("ref"));
                        }
                        proNames = "";
                        threadMessage("Fetching..." + proCount + "/" + pSort.entries.size());
                    }
                }
                gbkget.close();
                fnaout.close();
                faaout.close();
            } catch (Exception e) {
                threadMessage(e.getMessage());
                e.printStackTrace();
            }
        }

            LinkedHashMap<String, DNASequence> fastaRef = FastaReaderHelper.readFastaDNASequence(pSort.fnafasta);
            for (String gotAcc : pSort.entries.keySet()) {
                pEntry cur = pSort.entries.get(gotAcc);
                cur.length = fastaRef.get(gotAcc).getLength();
                pSort.entries.put(gotAcc, cur);
        }

    }
    /*    String datfile = "";
    ArrayList<pStruct> psorts = new ArrayList<pStruct>();
    for (int i = 0; i < INPUTFILES.size(); i++) {
    datfile += new File(INPUTFILES.get(i).toString()).getName();
    psorts.add(new pStruct(new File(INPUTFILES.get(i).toString())));
    threadMessage(".", true);
    }
    threadMessage("OK!", true);

    // Create dat folder
    File dat = new File("dat/");
    dat.mkdir();

    File masterFile = new File(dat.getAbsolutePath() + File.separator + "dat-" +datfile.hashCode());
    // Check for master protein fasta file in dat folder
    // Name should be of String.hashcode of input files.
    File masterFasta = new File(masterFile.getCanonicalPath() + ".faa");
    masterFasta.createNewFile();
    threadMessage("Your database file is: " + masterFasta.getName());
    // if exists check if it contains all proteins described in psort files.
    // Fill fields if possible
    ArrayList<String> allPro = new ArrayList<String>();
    for (pStruct cur : psorts) {
    for (pEntry tem : cur.entries) {
    allPro.add(tem.acc);
    }
    }
    int olSize = allPro.size();
    BufferedReader in = new BufferedReader(new FileReader(masterFasta));
    String line = "";
    ArrayList<String> removeList = new ArrayList<String>();
    while ((line = in.readLine()) != null) {
    if (line.startsWith(">") && !line.isEmpty()) {
    String ac = line.substring(1).split("\\s+")[0].split("\\|")[0].trim();
    removeList.add(ac);
    }
    }
    allPro.removeAll(removeList);
    in.close();
    threadMessage("Searching for new proteins: " + allPro.size() + "/" + olSize);
    // fetch and ammend missing proteins.
    // Header MUST be >acessionNo|locustag
    String proNames = "";
    int proCount = 0;
    BufferedWriter out = new BufferedWriter(new FileWriter(masterFasta, true));
    for (String cur : allPro) {
    proNames += cur + ",";
    // Bundle 50 at a time and fetch from intewebs
    proCount++;
    if (proCount % 50 == 0 || proCount == allPro.size()) {
    String fas = fetchFasta(proNames, "refseqp");
    if (fas.contains("ERROR")) {
    threadMessage(fas);
    }
    out.write(fas);
    out.newLine();
    proNames = "";
    threadMessage("Fetching..." + proCount + "/" + allPro.size());
    }
    }

    out.close();
    // Create blast database.
    threadMessage("Building database");
    pThread.execThis("makeblastdb -in " + masterFasta.getAbsolutePath(), OUT);

    // Filter psort results according to input
    threadMessage("Building query sequences");
    // Fetch fasta sequence for each entry from database.faa (one pass)
    LinkedHashMap<String, ProteinSequence> a = FastaReaderHelper.readFastaProteinSequence(masterFasta);
    //FastaReaderHelper.readFastaDNASequence for DNA sequences
    for (pStruct cur : psorts) {
    for (pEntry tem : cur.entries) {
    //    if (tem.acc.compareTo(entry.getValue().getOriginalHeader().split("\\s+")[0].split("\\|")[0].trim()) == 0) {
    //    tem.locus = entry.getValue().getOriginalHeader().split("\\s+")[0].split("\\|")[1].trim();
    //      cur.geneLocus = tem.locus.split("_")[0];
    //tem.locus = entry.getValue().getOriginalHeader().split("\\s+")[0].split("\\|")[1].trim();
    tem.locus = a.get(tem.acc).getOriginalHeader().split("\\s+")[0].split("\\|")[1].trim();
    tem.fasta = a.get(tem.acc).toString();
    //}
    }
    }

    // Seperate psort datastructure into number of threads
    threadMessage("Running BLAST");
    // Load BLAST datafile as a hash table.
    // check if BLAST result in BLAST hash.
    File masterBlast = new File(masterFile.getCanonicalPath() + ".Blast");
    masterBlast.createNewFile();

    HashMap blast = new HashMap();
    in = new BufferedReader(new FileReader(masterBlast));
    line = "";
    while ((line = in.readLine()) != null) {
    String[] lineArray = line.split("\\s+");
    if(lineArray.length >= 2){
    blast.put(lineArray[0], lineArray[1]);
    }
    }
    in.close();
    // BLAST against database, modify blast result for each entry
    for (pStruct cur : psorts) {
    int count = 1;
    for (pEntry tem : cur.entries) {

    if (count % 500 == 0 || count == cur.entries.size()) {
    threadMessage("BLASTp: [" + cur.location.getName() + "] " + count + "/" + cur.entries.size() + " complete.");
    }
    count++;
    if (blast.get(tem.acc) != null) {
    tem.blast = blast.get(tem.acc ).toString();

    } else {
    File tempB = new File(dat + File.separator + "temp-BLAST");
    out = new BufferedWriter(new FileWriter(tempB));
    out.write(">" + tem.acc);
    out.newLine();
    out.write(tem.fasta);
    out.newLine();
    out.close();
    String exec = "blastp -evalue 0.000000005 -outfmt 6 -query " + tempB.getCanonicalPath()
    + " -db " + masterFasta.getCanonicalPath();
    Process q = Runtime.getRuntime().exec(exec);
    InputStream istrm = q.getInputStream();
    InputStreamReader istrmrdr = new InputStreamReader(istrm);
    String data;
    BufferedReader buffrdr = new BufferedReader(istrmrdr);
    while ((data = buffrdr.readLine()) != null) {
    if (Double.parseDouble(data.split("\\s+")[2]) > 90.0
    && (Integer.parseInt(data.split("\\s+")[3]) > (0.90 * tem.fasta.length()))) {
    if (tem.acc.compareTo(data.split("\\s+")[1].split("\\|")[0]) != 0) {
    tem.blast += (data.split("\\s+")[1]) + ",";
    }
    }
    }
    q.destroy();
    }
    }
    }
    out = new BufferedWriter(new FileWriter(masterBlast));
    for (pStruct cur : psorts) {
    for (pEntry tem : cur.entries) {
    out.write(tem.acc + "\t" + tem.blast );
    out.newLine();

    }
    }
    out.close();
    // Load signalp datafile as a hash table.
    // check if signalp result in signalp hash.
    File masterSignal = new File(masterFile.getCanonicalPath() + ".signalp");
    masterSignal.createNewFile();

    HashMap signal = new HashMap();
    in = new BufferedReader(new FileReader(masterSignal));
    line = "";
    while ((line = in.readLine()) != null) {
    String[] lineArray = line.split("\\s+");
    signal.put(lineArray[0], lineArray[1]);
    }
    in.close();

    threadMessage("Running SignalP");

    for (pStruct cur : psorts) {
    int count = 1;
    for (pEntry tem : cur.entries) {
    // if exists check if it contains all proteins described in psort files.
    // Fill fields if possible
    if (signal.get(tem.acc) != null) {
    String sig = signal.get(tem.acc).toString();
    if (sig.compareTo("YES") == 0) {
    tem.signalp = true;
    } else {
    tem.signalp = false;
    }

    } else {
    // if not run signalp.
    // add result to data structure
    File temP = new File(dat + File.separator + "temp-SignalP");
    out = new BufferedWriter(new FileWriter(temP));
    out.write(">" + tem.acc);
    out.newLine();
    out.write(tem.fasta);
    out.newLine();
    out.close();
    String exec = "signalp-3.0/signalp -t gram- -f summary " + temP.getCanonicalPath();
    Process q = Runtime.getRuntime().exec(exec);
    InputStream istrm = q.getInputStream();
    InputStreamReader istrmrdr = new InputStreamReader(istrm);
    String data;
    BufferedReader buffrdr = new BufferedReader(istrmrdr);
    boolean NN = false;
    int HMMcoun = 0;
    while ((data = buffrdr.readLine()) != null) {
    if (data != null) {
    if (data.contains("Prediction: Signal peptide")) {
    NN = true;
    }
    String[] linArray = data.split("\\s+");
    if (!data.contains("#")) {
    if (linArray.length == 6 && linArray[5].compareTo("YES") == 0) {
    HMMcoun++;
    }
    if (linArray.length == 7 && linArray[6].compareTo("YES") == 0) {
    HMMcoun++;
    }
    }
    }
    }
    q.destroy();
    if (NN && (HMMcoun == 5)) {
    tem.signalp = true;
    } else {
    tem.signalp = false;
    }
    }
    if (count % 500 == 0 || count == cur.entries.size()) {
    threadMessage("SignalP: [" + cur.location.getName() + "] " + count + "/" + cur.entries.size() + " complete.");
    }
    count++;
    }
    }
    // Generate output signalp datafile (for next time)
    out = new BufferedWriter(new FileWriter(masterSignal));
    for (pStruct cur : psorts) {
    for (pEntry tem : cur.entries) {
    out.write(tem.acc + "\t");
    if (tem.signalp) {
    out.write("YES");
    } else {
    out.write("NO");
    }
    out.newLine();

    }
    }
    threadMessage("Signalp complete.");
    out.close();
    //       final boolean VERBOSE;
    //     final boolean CONSERVED;
    //     final boolean SIGNAL;
    //     final String[] LOCALISATION;

    // construct HASHMAPS for searching
    HashMap ha = new HashMap();
    if (CONSERVED) {
    for (pStruct cur : psorts) {
    for (pEntry te : cur.entries) {
    ha.put(te.acc, cur.location);
    }
    }
    }
    // create table
    for (pStruct cur : psorts) {
    // Sort psort results
    //  Collections.sort( cur.entries);

    out = new BufferedWriter(new FileWriter(cur.location.getCanonicalPath() + ".FINAL"));
    out.write("ACCESSION\tLOCUS_TAG\tSIGNAL-PEPTIDE\tLOCALISATION\tDESCRIPTION\tBLAST MATCHES");
    out.newLine();
    for (pEntry tem : cur.entries) {
    boolean printIt = true;
    if (SIGNAL) {
    if (tem.signalp == false) {
    printIt = false;
    }
    }
    if (LOCALISATION != null) {
    int hit = 0 ;
    for (String l : LOCALISATION){
    if( tem.desc.split("\t")[31].compareToIgnoreCase(l) == 0 ){
    hit++;
    }
    }
    if (hit == 0 ){
    printIt = false;
    }
    }
    if (CONSERVED) {
    Set con = new HashSet();
    if (tem.blast != null) {
    String[] bls = tem.blast.split(",");
    for (String d : bls) {
    System.out.println( ha.get(d.split("\\|")[0]) ) ;
    con.add(ha.get(d.split("\\|")[0]));
    }
    System.out.println( con.size() + "\t" +  psorts.size() ) ;
    if (con.size() < (psorts.size()-1) ){
    printIt = false;
    }
    } else {
    printIt = false;
    }
    }
    if (printIt) {
    out.write(tem.acc + "\t");
    out.write(tem.locus + "\t");
    if (tem.signalp) {
    out.write("YES\t");
    } else {
    out.write("NO\t");
    }
    if (tem.desc.split("\t")[33].contains("multiple localization sites")) {
    out.write(tem.desc.split("\t")[31] + "*\t");
    } else {
    out.write(tem.desc.split("\t")[31] + "\t");
    }
    out.write(tem.desc.split("\t")[1] + "\t");
    out.write(tem.blast + "\t");
    out.newLine();
    }
    }
    out.close();
    }*/

    public static String[] expandArray(String[] a) {
        String[] newArray = new String[a.length + (a.length / 2)];
        System.arraycopy(a, 0, newArray, 0, a.length);
        return newArray;
    }

    public static String fetchAnno(String proteinString, String database) throws MalformedURLException, IOException {
        String content = "";
        URLConnection connection = null;
        if (proteinString.contains(".")) {
            proteinString = proteinString.substring(0, proteinString.indexOf("."));
        }

        String ooo = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=" + database + "&id=" + proteinString + "&style=raw";
        threadMessage(ooo);
        connection = new URL(ooo).openConnection();
        Scanner scanner = new Scanner(connection.getInputStream());
        scanner.useDelimiter("\\z");
        while (scanner.hasNext()) {
            content += scanner.next();
        }
        return content;

    }

    public static File GetGBK(String acc, File wd) throws MalformedURLException, IOException, BioException {
        File gbk = new File(wd + File.separator + acc + ".gbk");
        if (!gbk.exists()) {
            String ooo = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=" + "refseqn" + "&id=" + acc + "&style=raw";
            threadMessage(ooo);
            URLConnection connection = new URL(ooo).openConnection();
            BufferedReader cont = new BufferedReader(new InputStreamReader(connection.getInputStream()));
            OutputStream gbkout = new FileOutputStream(gbk);
            RichSequenceIterator ss = RichSequence.IOTools.readGenbankDNA(cont, null);
            while (ss.hasNext()) {
                RichSequence zs = ss.nextRichSequence();
                RichSequence.IOTools.writeGenbank(gbkout, zs, null);
            }
        }
        return gbk;
    }

    public static String GenbankSubseq(String acc, int start, int stop, File wd) throws MalformedURLException, IOException, BioException {
        File gbk = new File(wd + File.separator + acc + ".gbk");
        String seq = "";
        if (!gbk.exists()) {
            String ooo = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=" + "refseqn" + "&id=" + acc + "&style=raw";
            threadMessage(ooo);
            URLConnection connection = new URL(ooo).openConnection();
            BufferedReader cont = new BufferedReader(new InputStreamReader(connection.getInputStream()));
            OutputStream gbkout = new FileOutputStream(gbk);
            RichSequenceIterator ss = RichSequence.IOTools.readGenbankDNA(cont, null);
            while (ss.hasNext()) {
                RichSequence zs = ss.nextRichSequence();
                RichSequence.IOTools.writeGenbank(gbkout, zs, null);
            }
        }
        BufferedReader br = new BufferedReader(new FileReader(gbk));
        RichSequenceIterator s = RichSequence.IOTools.readGenbankDNA(br, null);
        while (s.hasNext()) {
            RichSequence zz = s.nextRichSequence();
            seq = zz.subStr(start, stop);
        }
        return seq;
    }
}
// Create database of proteins for each table (if it doesnt exist)
 /*       threadMessage("Fetching protein sequences.");
List<File> fastaFiles = new ArrayList<File>();
int filec = 0;
for (String[] curFile : psort) {
// Check if file exists
File curz = new File(INPUTFILES.get(filec).toString() + ".faa");
fastaFiles.add(curz);
filec++;
if (curz.exists()) {
threadMessage("File " + curz.getName() + " already exists, is it complete?");
} else {
try {
threadMessage("Fetching sequences for " + curz.getName());
BufferedWriter proFast = new BufferedWriter(new FileWriter(curz));
// Fetch proteins
int proCount = 0;
String proNames = "";
for (String curLine : curFile) {
if (curLine != null) {
String[] lineArray = curLine.split("\t");
if (!lineArray[0].contains("SeqID")) {
try {
proNames += lineArray[0].split("\\|")[3] + ",";
// Bundle 50 at a time and fetch from intewebs
proCount++;
if (proCount % 50 == 0) {
String fas = fetchFasta(proNames, "refseqp");
if(fas.contains("ERROR")){
threadMessage(fas);
}
proFast.write(fas);
proFast.flush();
proNames = "";
threadMessage(".", true);
}
} catch (Exception e) {
e.printStackTrace();
}
}
}
}
if (proNames.length() > 0) {
proFast.write(fetchFasta(proNames, "refseqp"));
proFast.flush();
}
proFast.close();
} catch (Exception e) {
threadMessage("ERROR: " + e.getMessage());
}
}
}
// BLAST all-against-all (multi-thread jobs for each comparison)
threadMessage("Running BLASTp\n");
// Check if output already exists
List<Thread> blastThreads = new ArrayList<Thread>();
for (int i = 0; i < fastaFiles.size(); i++) {
File signalout = new File( fastaFiles.get(i) +".tab");
if (!signalout.isFile()) {
File[] dbs = new File[fastaFiles.size() ];
for (int j = 0; j < fastaFiles.size(); j++) {
if (fastaFiles.get(j).getName().compareTo(fastaFiles.get(i).getName()) != 0) {
dbs[j] = fastaFiles.get(j);
}
}
Blast si = new Blast(fastaFiles.get(i), dbs, OUT);
Thread t = new Thread(si, ("BLAST-" + i));
t.start();
blastThreads.add(t);
} else {
threadMessage("File " + signalout.getName() + " already exists, is it complete?\n");
}
}
for (Thread t : blastThreads) {
try {
t.join();
} catch (InterruptedException ex) {
threadMessage("Thread interupt: " + ex.getLocalizedMessage());
}
}
// run protein against signal p,
threadMessage("Running signalP\n");
List<Thread> signalThreads = new ArrayList<Thread>();
for (int i = 0; i < fastaFiles.size(); i++) {
File signalout = new File("signalp-" + fastaFiles.get(i));
if (!signalout.isFile()) {
signalP si = new signalP(fastaFiles.get(i), OUT);
Thread t = new Thread(si, ("SIGNALP-" + i));
t.start();
signalThreads.add(t);
} else {
threadMessage("File " + signalout.getName() + " already exists, is it complete?");
}
}
for (Thread t : signalThreads) {
try {
t.join();
} catch (InterruptedException ex) {
threadMessage("Thread interupt: " + ex.getLocalizedMessage());
}
}
// Add results to row.

// Apply filters, write to new file.
threadMessage("Applying filters");
threadMessage(".", true);
threadMessage("OK!\n", true);*/

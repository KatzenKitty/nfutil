/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package revar;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.io.FileProxyDNASequenceCreator;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;



/**
 *
 * @author nabil
 */
class revarSimple implements Runnable {

    File snpDir;
    File consenDir;
    File refSeq;
    File samDir;
    File wd;
    File sesh;
    List<Gene> genes = new ArrayList<Gene>();
    List<Gene> Badgenes = new ArrayList<Gene>();
    long runTime = System.currentTimeMillis();
    double percentNcutoff = 0.20;
    double duplicatePer = 0.60;
    String prefix = null;
    String crop = null ;

    public revarSimple(String session, String inprefix, String incrop) {
        try {
            prefix = inprefix;
            crop = incrop;
            time();
            sesh = new File(session);
            BufferedReader in = new BufferedReader(new FileReader(session));
            String line = "";
            while ((line = in.readLine()) != null) {
                if (!line.startsWith("#")) {
                    if (line.startsWith("ref")) {
                        refSeq = new File(line.split("=")[1]);
                    } else if (line.startsWith("consenseSeq")) {
                        consenDir = new File(line.split("=")[1]);
                    } else if (line.startsWith("snp")) {
                        snpDir = new File(line.split("=")[1]);
                    } else if (line.startsWith("wd")) {
                        wd = new File(line.split("=")[1]);
                    } else if (line.startsWith("sam")) {
                        samDir = new File(line.split("=")[1]);
                    } else if (line.startsWith("@")) {
                        try {
                            String[] lineArray = line.substring(1).split("\\s+");
                            Gene temp = new Gene(lineArray[0], Integer.parseInt(lineArray[1]), Integer.parseInt(lineArray[2]), null, null);
                            genes.add(temp);
                        } catch (Exception e) {
                            threadMessage("Could not read " + line + ": " + e.getMessage());
                        }
                    }
                }
            }
            in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void time() {
        threadMessage("Run time since call: " + Long.toString(System.currentTimeMillis() - runTime));
        runTime = System.currentTimeMillis();
    }

    public void run() {
        threadMessage("Begun REVAR-Simple");
        threadMessage("prefix: " +prefix);
        threadMessage("crop: " +crop);
       SortedMap snps = new TreeMap<Integer,String[]>();
        try {
            // Find List of genes to form consensus tree.
            // NB within sessionFile
            // For list of genes....
            // Load FASTA files in consensus folder
            // Touch gene files
            File fnaDir = new File(wd.getAbsolutePath() + File.separator + "fna");
            if (!fnaDir.exists()) {
                fnaDir.mkdir();
            }
            FastaReader<DNASequence, NucleotideCompound> fastaProxyReader =
                    new FastaReader<DNASequence, NucleotideCompound>(refSeq,
                    new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
                    new FileProxyDNASequenceCreator(refSeq, AmbiguityDNACompoundSet.getDNACompoundSet()));
            LinkedHashMap<String, DNASequence> a = fastaProxyReader.process();
            for (String h : a.keySet()) {
                for (Gene cur : genes) {
                    cur.seq = a.get(h).toString().substring(cur.start, cur.stop);
                }
            }
            ArrayList<String> conHead = new ArrayList<String>();
            String refname = refSeq.getName().split("\\.")[0].substring(0, 10);
            File[] snpFiles = snpDir.listFiles();
            int totalFiles = 1;
                        for (File snp : snpFiles) {
                if (snp.getName().endsWith(".vcf")) {
                    totalFiles++;
                }
            }
            HashMap snpLoc = new HashMap();
            int cure= 1;
            for (File snp : snpFiles) {
                if (snp.getName().endsWith(".vcf")) {
                    threadMessage("Adding snps from  " + snp.getName() + " Column " + (cure+1)) ;
                    BufferedReader in = new BufferedReader(new FileReader(snp.getAbsolutePath()));
                    String line = "";
                    while ((line = in.readLine()) != null) {
                        if (!line.startsWith("#")) {
                            String[] linArray = line.split("\\s+");
                            int pos = Integer.parseInt(linArray[1]);
                            double qual = Double.parseDouble(linArray[5]);
                            if (qual > 10) {
                                if (snps.containsKey(pos)) {
                                    String[] te = (String[]) snps.get(pos);
                                    te[cure] = linArray[4];
                                    if (linArray[4].contains(",")){
                                        te[cure] = linArray[4].split(",")[0];
                                    }
                                    snps.put(pos, te);
                                } else {
                                    String[] te = new String[totalFiles];
                                    for(int i =0; i<te.length; i++ ){
                                        te[i] = linArray[3];
                                    }
                                    te[cure] = linArray[4];
                                    if (linArray[4].contains(",")){
                                        te[cure] = linArray[4].split(",")[0];
                                    }                                    
                                    snps.put(pos, te);
                                }
                                snpLoc.put(pos, "SNP");
                            }
                        }
                    }
                    cure++;
                }
                
            }

            for (Object s :  snps.keySet() ) {
                System.out.print(s.toString() + "\t");
                for(String e : (String[]) snps.get(s) ){
                    System.out.print(e + "\t");
                }
                System.out.print("\n");
            }
            threadMessage("total number of snps: " + snps.keySet().size());

            time();
            String outfile = wd.getAbsolutePath() + File.separator + sesh.getName() + ".fna";
            BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
            out.write(">" + refname);
            out.newLine();
            int ch = 1;
            for (Object s : snps.keySet()) {
                out.write(((String[]) snps.get(s))[0]);
                if (ch % 80 == 0) {
                    out.newLine();
                }
                ch++;
            }
            out.newLine();
            cure = 1;
            for (File snp : snpFiles) {
                if (snp.getName().endsWith(".vcf")) {
                    String conname = snp.getName().split("\\.")[0];
                    if (crop != null) {
                        conname = conname.replaceAll(crop, "");
                    }
                    if (prefix != null) {
                        conname = prefix + conname;
                    }
                    out.write(">" + conname);
                    out.newLine();
                     ch = 1;
                    for (Object s : snps.keySet()) {
                        out.write(((String[]) snps.get(s))[cure]);
                        if (ch % 80 == 0) {
                            out.newLine();
                        }
                        ch++;
                    }
                    out.newLine();
                    cure++;
                }
            }
               out.close();

   /*             threadMessage("Loaded " + snpLoc.size()  + " snps");
            for (Gene cur : genes) {
                threadMessage("Preparing " + cur.geneName);
                time();
                String outfile = wd.getAbsolutePath() + File.separator + "fna" + File.separator + cur.geneName + ".fna";
                cur.fasta = new File(outfile);
                BufferedWriter out = new BufferedWriter(new FileWriter(cur.fasta));
                // Fetch reference

                out.write(">" + refname);
                if (firsta) {
                    conHead.add(refname);
                }
                out.newLine();
                char[] faz = cur.seq.toCharArray() ;
                int l = 0;
                int ch = 1;


                for (char f : faz) {
                    if (snpLoc.get(l ) != null) {
                        out.write(f);
                        if (ch % 80 == 0) {
                            out.newLine();
                        }
                        ch++;
                    }
                    l++;
                }
                out.newLine();
                // For each gene determine snps for each sample.
                for (File snp : snpFiles) {
                    if (snp.getName().endsWith(".vcf")) {
                        threadMessage("Reading " + snp.getName());
                        BufferedReader in = new BufferedReader(new FileReader(snp.getAbsolutePath()));
                        String line = "";
                        char[] fas = cur.seq.toCharArray();
                        while ((line = in.readLine()) != null) {
                            if (!line.startsWith("#")) {
                                String[] linArray = line.split("\\s+");
                                int pos = Integer.parseInt(linArray[1]);
                                double qual = Double.parseDouble(linArray[5]);
                                if (pos > cur.stop) {
                                    break;
                                }
                          //      System.out.println(  cur.start +"\t" + cur.stop + "\t" + pos + "\t" + qual);
                                if ((qual > 20) && pos >= cur.start && pos <= cur.stop && (linArray[4].length() == 1)) {
                                //    System.out.print( linArray[4].charAt(0));
                                    fas[pos ] = linArray[4].charAt(0);
                                }
                            }
                        }

                        String conname = snp.getName().split("\\.")[0];
                        if (crop != null ) {
                            conname = conname.replaceAll(crop, "");
                        }
                        if (prefix != null ){
                        conname = prefix + conname;
                        }
                        if (firsta) {
                            conHead.add(conname);
                        }
                        out.write(">" + conname);
                        out.newLine();
                         l = 0;
                         ch = 1;
                        for (char f : fas) {
                            if (snpLoc.get(l ) != null  ) {
                                out.write(f );
                                if (ch % 80 == 0 ) {
                                    out.newLine();
                                }
                                ch++;
                            }
                            l++;
                        }
                        out.newLine();
                        in.close();
                    }
                }
                firsta = false;
                out.close();
            }


            // --- You now have a pool of sequences MULTI-FASTA from variety of sources  --- //
            threadMessage("You now have a MULTI-FASTA pool of sequences  from variety of sources in " + wd.getAbsolutePath() + File.separator + "fna");
            // Merge all genes into single file - PHYLIP format

            String phyFile = wd.getAbsolutePath() + File.separator + sesh.getName() + ".phy";
            threadMessage("Formatting Phylip file to:  " + phyFile);
            // Merge Muscle results.
            // For each gene (check if bad gene) grab each sequence.

            if (conHead.size() > 0) {
                int first = 0;
                BufferedWriter out = new BufferedWriter(new FileWriter(phyFile));
                String myFile = wd.getAbsolutePath() + File.separator + sesh.getName() + ".fna";
                BufferedWriter out2 = new BufferedWriter(new FileWriter(myFile));

                for (Object head : conHead) {
                    threadMessage("Fetching " + head.toString().split(":")[0]);
                    time();
                    String hap = "";
                    for (Gene cur : genes) {
                        if (cur.fasta.exists()) {
                        FileInputStream inStream = new FileInputStream( refSeq );
                            LinkedHashMap<String, ProteinSequence> af = FastaReaderHelper.readFastaProteinSequence(cur.fasta);
                            // Check if sequence is just duplicated.
                            // Create gene Tree too!
                            // Run FastTreeMP
                            // Find correct header sequence.
                            Object[] tempHead = af.keySet().toArray();
                            for (Object temp : tempHead) {
                                if (temp.toString().compareTo(head.toString()) == 0) {
                                    hap += af.get( temp.toString() );
                                }
                            }
                        }
                    }
                    if (first == 0) {
                        out.write("   " + conHead.size() + "    " + hap.length());
                        out.newLine();
                        out.flush();
                    }

                    time();
                    String headerr = head.toString().split(":")[0];
                    if (headerr.length() < 10) {
                        threadMessage("Writing " + headerr.substring(0, headerr.length()).replaceAll(" ", "_") + "\t" + hap.length());
                        out2.write(">" + headerr.substring(0, headerr.length()).replaceAll(" ", "_"));
                        out.write(headerr.substring(0, headerr.length()).replaceAll(" ", "_"));
                        for (int i = headerr.length(); i < 10; i++) {
                            out2.write("_");
                            out.write("_");
                        }
                    } else {
                        threadMessage("Writing " + headerr.substring(0, 10).replaceAll(" ", "_") + "\t" + hap.length());
                        out2.write(">" + headerr.substring(0, 10).replaceAll(" ", "_"));
                        out.write(headerr.substring(0, 10).replaceAll(" ", "_"));
                    }
                    out2.newLine();
                    WordWrapFile(hap, 80, out2);
                    out2.newLine();
                    out.newLine();
                    out.flush();
                    FastaToPhyFile(hap, out);
                    out.newLine();
                    out.flush();
                    first++;
                }
                out.close();
                out2.close();
                threadMessage("Done (Now create the tree from the phy file)");
            }*/


        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    static void threadMessage(String message) {
        String threadName = Thread.currentThread().getName();
        System.out.format("%s: %s%n", threadName, message);
    }
    public static void FastaToPhyFile(String fasta, BufferedWriter out) throws IOException {
        char[] fastachar = fasta.toCharArray();
        for (int i = 0; i < fastachar.length; i++) {
            out.write(fastachar[i]);
            if (i % 60 == 0 && i != 0) {
                out.newLine();
                out.flush();
            } else if (i % 10 == 0 && i != 0) {
                out.write(" ");
            }
        }
    }

    public static void WordWrapFile(String unformat, int wrap, BufferedWriter out) throws IOException {
        char[] ar = unformat.toCharArray();
        for (int i = 0; i < ar.length; i++) {
            if (i % wrap == 0 && i != 0) {
                out.newLine();
                out.flush();
            }
            out.write(ar[i]);
        }

    }

    public String WordWrap(String unformat, int wrap) {
        String out = "";
        char[] ar = unformat.toCharArray();
        for (int i = 0; i < ar.length; i++) {
            if (i % wrap == 0 && i != 0) {
                out += "\n";
            }
            out += ar[i];
        }
        return out;
    }
}

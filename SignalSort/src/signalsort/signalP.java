/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package signalsort;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;

/**
 *
 * @author nabil
 */
public class signalP implements Runnable {

    static PrintStream OUT;
    File FILE;

    public signalP(File file, PrintStream inOut) {
        FILE = file;
        OUT = inOut;
    }

    public void run() {
        try {
            BufferedReader in = new BufferedReader(new FileReader(FILE));
            String line = "";
            int lineCount = 1;
            int headCount = 0;
            File tempfile = new File("temp-" + Thread.currentThread().getName());
            BufferedWriter proFast = new BufferedWriter(new FileWriter(tempfile));
            BufferedWriter proOut = new BufferedWriter(new FileWriter("signalp-" + FILE.getName()));
            threadMessage("Running SignalP on: " + FILE.getName() + "\n");
            while ((line = in.readLine()) != null) {
                if (line.contains(">") & lineCount != 1) {
                    proFast.close();
                    String exec = "signalp-3.0/signalp -t gram- -f summary temp-" + Thread.currentThread().getName();
                    Process q = Runtime.getRuntime().exec(exec);
                    InputStream istrm = q.getInputStream();
                    InputStreamReader istrmrdr = new InputStreamReader(istrm);
                    String data;
                    BufferedReader buffrdr = new BufferedReader(istrmrdr);
                    boolean NN = false;
                    int HMMcoun = 0;
                    String header = "" ;
                    while ((data = buffrdr.readLine()) != null) {
                        if (data != null) {
                       
                            if(data.startsWith(">")){
                                header = data.split("\\s+")[0].split("\\|")[1];
                            }
                            if (data.contains("Prediction: Signal peptide")) {
                                NN = true;
                            }
                            String[] linArray = data.split("\\s+");
                            if (!line.contains("#")) {
                                if (linArray.length == 6 && linArray[5].compareTo("YES") == 0 ) {
                                    HMMcoun++;
                                }
                                if (linArray.length == 7 && linArray[6].compareTo("YES") == 0 ) {
                                    HMMcoun++;
                                }
                            }
                        }
                    }
                    q.destroy();
                    proOut.write(header);
                    if (NN && (HMMcoun == 5)) {
                        proOut.write("\tYES");
                    }else{
                        proOut.write("\tNO");
                    }
                    proOut.newLine();
                    proOut.flush();
                    headCount++;
                    proFast = new BufferedWriter(new FileWriter("temp-" + Thread.currentThread().getName()));
                    if (headCount % 100 == 0) {
                        threadMessage(headCount + " complete.\n");
                    }
                }
                proFast.write(line);
                proFast.newLine();
                proFast.flush();
                lineCount++;
            }
            proOut.close();
            in.close();
            tempfile.delete();
        } catch (Exception e) {
            e.printStackTrace();
            threadMessage("ERROR: " + e.getLocalizedMessage());
        }
    }

    static void threadMessage(String message) {
        String threadName = Thread.currentThread().getName();
        OUT.print("[" + threadName + "] " + message);
    }
}

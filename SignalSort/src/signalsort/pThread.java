/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package signalsort;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.concurrent.Callable;

/**
 *
 * @author nabil
 */
public class pThread implements Callable<pStruct> {
    static PrintStream OUT;
    File FILE;
    File[] dbs;
    public pThread(File file, File[] db, PrintStream inOut) {
        FILE = file;
        OUT = inOut;
        dbs = db; 
    }


 /*   public void run() {
        try{
        File bls = new File( "bls/");
        bls.mkdir();
        String d = "";
        int dbcount = 0 ;
        for(File db : dbs){
            if( db != null ){
               d += db.getCanonicalPath()  + " ";
               dbcount++;
            }
        }
        if(dbcount > 1){
            d += "\""+d +"\"";
        }
        String exec = "makeblastdb -in "+d  +" -out " + bls.getAbsolutePath() +File.separator+ FILE.getName() + "-db";
        threadMessage("Formatting database of "+ d+"\n");
        execThis(exec, OUT);
        exec = "blastp -outfmt 6 -query "+ FILE.getCanonicalPath() +  " -db " + bls.getAbsolutePath() +File.separator+ FILE.getName() + "-db -out " + FILE.getName() + ".tab";
        threadMessage("Running BLAST: " + exec +"\n");
        execThis(exec, OUT);
        }catch(Exception e ){
            e.printStackTrace();
        }
    }*/
    
    public static void execThis(String exec, PrintStream str) {
        try {
            Process q = Runtime.getRuntime().exec(exec);
            InputStream istrm = q.getErrorStream();
            InputStreamReader istrmrdr = new InputStreamReader(istrm);
            String data;
            BufferedReader buffrdr = new BufferedReader(istrmrdr);
            while ((data = buffrdr.readLine()) != null) {
                threadMessage(data + "\n", str);
            }
            q.destroy();
        } catch (Exception e) {
            e.printStackTrace();
            threadMessage(e.getMessage() + "\n", str);
        }
    }

    static void threadMessage(String message) {
        String threadName = Thread.currentThread().getName();
        OUT.print("\n[" + threadName + "] " + message);
    }

    static void threadMessage(String message, boolean tick) {
        if (tick) {
            OUT.print(message);
        }
    }
    static void threadMessage(String message, PrintStream o ) {
        String threadName = Thread.currentThread().getName();
        o.print("\n[" + threadName + "] " + message);
    }
    public pStruct call() throws Exception {
        throw new UnsupportedOperationException("Not supported yet.");
            // BLAST against database, modify blast result for each entry

            // Load signalp datafile as a hash table.
            // check if signalp result in signalp hash.


            // if not run signalp.

            // add result to data structure
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package revar;

import java.io.File;

/**
 *
 * @author nabil
 */
public class Gene {

    public String geneName;
    public int start;
    public int stop;
    public File fasta;
    public String seq;
    boolean core;
    int counter;
    String genome;

    public Gene(String g, int s, int st, File f, String lseq) {
        geneName = g;
        start = s;
        stop = st;
        fasta = f;
        seq = lseq;
        core = true;
        counter = 0;
        genome = "" ;
    }
}
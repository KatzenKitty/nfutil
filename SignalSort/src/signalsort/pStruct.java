/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package signalsort;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author nabil
 */
public class pStruct {
    File location;
    LinkedHashMap<String, pEntry> entries;
    File fnafasta; 
    File faafasta;

    public pStruct(File pSortloc) throws FileNotFoundException, IOException {
        location = pSortloc;
        entries = new LinkedHashMap<String, pEntry>();
        buildEntries();
    }

    private void buildEntries() throws FileNotFoundException, IOException {
        BufferedReader in = new BufferedReader(new FileReader(location));
        String line = "";
        while ((line = in.readLine()) != null) {
            String[] linArray = line.split("\t");
            if (linArray[0].compareTo("SeqID") != 0) {
                pEntry temp = new pEntry();
                String get = linArray[0].split("\\|")[3];
                if (get.contains(".")) {
                    get = get.substring(0, get.indexOf("."));
                }
                temp.acc = get;
                temp.gi = linArray[0].split("\\|")[1];
                temp.desc = line ; 
                entries.put(temp.acc, temp);
            }
        }
            in.close();
    }


    
}


 class pEntry {
    String acc;
    String gi; 
    String locus;
    String desc = "" ;
    String fasta = "" ;
     int length = 0;
     String signalpHMMResult = "";
     String signalpHMMDesc = "";
     String signalpNNResult = "";
     String signalpNNDesc = "";
     boolean conserved = false;
    LinkedHashMap blastMatch = new LinkedHashMap<String,String>();


     public void addBlastMatch(String file, String inMatch) {
         if (!blastMatch.containsValue(file)) {
             blastMatch.put(file, inMatch);
         } else {
             String oldMatch = blastMatch.get(file).toString();
             if (Integer.parseInt(oldMatch.split("\t")[11]) < Integer.parseInt(inMatch.split("\t")[11])) {
                 blastMatch.put(file, inMatch);
             }
         }
     }
   
    
}



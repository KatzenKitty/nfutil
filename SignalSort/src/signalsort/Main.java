package signalsort;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 *
 * @author Nabil-Fareed Alikhan 2011.
 */
public class Main {
    final static String COMMAND = "java -jar signalSort.jar";
    final static String USAGE = "SignalSort version 0.5.\n Nabil-Fareed Alikhan 2011.\n"
            + "INPUT: Table of pSORTB cell localisation results & directories of genomes and databases (FASTA) to compare to\n"
            + "OUTPUT: Tab-delimited table of ACCESSION, LOCUS_TAG, SIGNAL-PEPTIDE, LOCALISATION, DESCRIPTION, BLAST MATCHES &"
            + " a seperate FASTA File for sequences.";
    final static String FOOTER = "";
    static boolean VERBOSE = false;
    static PrintStream OUT = (PrintStream) System.out;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        CommandLineParser parser = new GnuParser();
        final Options getOpt = constructOptions();
        CommandLine commandLine;
        try {
            commandLine = parser.parse(getOpt, args);
            if (commandLine.hasOption("h")) {
                printHelp(getOpt, 80, USAGE, FOOTER, 5, 3, true, OUT);
            } else {
                List inputFiles = commandLine.getArgList();
                if (inputFiles.size() != 1) {
                    threadMessage("BAD INPUT: You must specify only one pSORTB file\n");
                    printUsage(COMMAND, getOpt, OUT);
                } else {
                    if (commandLine.hasOption("v")) {
                        VERBOSE = true;
                    }

                    threadMessage("Checking Input");
                    String valid = validateInput(inputFiles);

                    if (valid.length() == 0) {
                        threadMessage("OK!\n", true);
                        boolean inSig = false;
                        String[] inFolders = null; 
                        if (commandLine.hasOption("s")) {
                            inSig = true;
                        }
                        String[] inlocal = null;
                        if (commandLine.hasOption("l")) {
                            inlocal = commandLine.getOptionValues("l");
                        }
                        if (commandLine.hasOption("f")) {
                             inFolders = commandLine.getOptionValues("f");
                        }
                        String inblast = "blastn";
                        if (commandLine.hasOption("b")) {
                            inblast = commandLine.getOptionValue("b");
                        }
                        Thread t = new Thread(new Signal(inputFiles, OUT, VERBOSE, inFolders, inSig, inlocal, inblast), "SIGNAL");
                        t.start();
                    } else {
                        threadMessage("\nBAD INPUT: " + valid +"\n");
                    }
                }
            }
        } catch (ParseException parseException) {
            threadMessage(
                    "BAD INPUT: "
                    + parseException.getMessage() +"\n");
             printUsage(COMMAND, getOpt, OUT);
        }

    }

    public static String validateInput(List inputFiles) {
        String out = "";
        for (Object input : inputFiles) {
            File inputFile = new File(input.toString());
            if (inputFile.isFile()) {
                try {
                    BufferedReader in = new BufferedReader(new FileReader(inputFile));
                    String line = "";
                    int lineCount = 1 ;
                    while ((line = in.readLine()) != null) {
                        String[] lineArray = line.split("\t");
                        if(lineArray.length != 37 ){
                            out += "Not valid a line: " + input.toString() + " [line: " +lineCount + "]\n";
                        }
                        lineCount++;
                    }
                } catch (FileNotFoundException ex) {
                    out += "Could not find file: " + input.toString() + " [" + ex.getMessage() + "]\n";
                } catch (IOException ex) {
                    out += "Error reading file: " + input.toString() + " [" + ex.getMessage() + "]\n";
                }
                threadMessage(".",true);
            } else {
                out += input.toString() + " is not a valid file\n";
            }
        }
        return out;
    }

    static void threadMessage(String message) {
        String threadName = Thread.currentThread().getName();
        OUT.print("[" + threadName + "] " + message);
    }

    static void threadMessage(String message, boolean tick) {
        if (tick) {
            OUT.print(message);
        }
    }
    public static Options constructOptions() {
        final Options gnuOptions = new Options();
        gnuOptions
                .addOption("h", "help", false, "print help message")
                .addOption("b", "blast", true, "Blast type (blastx, tblastn)")
                .addOption("s", "signal", false, "show only proteins with a signal peptide")
                .addOption("f", "folder", true, " compare psort file against gbk & fasta in this folder (you can specify multiple locations)")
                .addOption("l", "localisation", true, "show only proteins with the specified localisation (you can specify multiple locations)")
                .addOption("v", "verbose", false, "be extra verbose");
        return gnuOptions;
    }

    public static void printUsage(
            final String applicationName,
            final Options options,
            final OutputStream out) {
        final PrintWriter writer = new PrintWriter(out);
        final HelpFormatter usageFormatter = new HelpFormatter();
        usageFormatter.printUsage(writer, 80, applicationName, options);
        writer.close();
    }

    public static void printHelp(
            final Options options,
            final int printedRowWidth,
            final String header,
            final String footer,
            final int spacesBeforeOption,
            final int spacesBeforeOptionDescription,
            final boolean displayUsage,
            final OutputStream out) {
        final String commandLineSyntax = COMMAND;
        final PrintWriter writer = new PrintWriter(out);
        final HelpFormatter helpFormatter = new HelpFormatter();
        helpFormatter.printHelp(
                writer,
                printedRowWidth,
                commandLineSyntax,
                header,
                options,
                spacesBeforeOption,
                spacesBeforeOptionDescription,
                footer,
                displayUsage);
        writer.close();
    }
}

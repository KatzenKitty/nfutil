/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package revar;

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
 * @author nabil
 */
public class Main {

    final static String COMMAND = "java -jar revar.jar";
    final static String USAGE = "Revar version 0.4.\n Nabil-Fareed Alikhan 2011.\n";
    final static String FOOTER = "";
    static boolean VERBOSE = false;
    static PrintStream OUT = (PrintStream) System.out;

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
        gnuOptions.addOption("h", "help", false, "print help message").addOption("c", "crop", true, "crop this text from filenames").addOption("p", "prefix", true, "prefix filenames with this text").addOption("v", "verbose", false, "be extra verbose");
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

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws InterruptedException {
        CommandLineParser parser = new GnuParser();
        final Options getOpt = constructOptions();
        CommandLine commandLine;
        try {
            commandLine = parser.parse(getOpt, args);
            if (commandLine.hasOption("h")) {
                printHelp(getOpt, 80, USAGE, FOOTER, 5, 3, true, OUT);
            } else {
                List inputFiles = commandLine.getArgList();
                if (inputFiles.size() < 1) {
                    threadMessage("BAD INPUT: You must specify a config file\n");
                    printUsage(COMMAND, getOpt, OUT);
                } else {
                    if (commandLine.hasOption("v")) {
                        VERBOSE = true;
                    }
                    String crop = null;
                    String prefix = null;
                    //
                    if (commandLine.hasOption("c")) {
                        crop = commandLine.getOptionValue("c");
                    }
                    if (commandLine.hasOption("p")) {
                        prefix = commandLine.getOptionValue("p");
                    }
                    //          threadMessage("Checking Input");
                    String valid = validateInput(inputFiles);
                    if (valid.length() == 0) {
                        //     threadMessage("OK!\n", true);
                        Thread t = new Thread(new revarSimple(inputFiles.get(0).toString(), crop, prefix));
                        t.start();
                    } else {
                        threadMessage("BAD INPUT: " + valid + "\n");
                    }
                }
            }
        } catch (ParseException parseException) {
            threadMessage(
                    "BAD INPUT: "
                    + parseException.getMessage() + "\n");
            printUsage(COMMAND, getOpt, OUT);
        }

    }

    public static String validateInput(List inputFiles) {
        return "";
    }
}

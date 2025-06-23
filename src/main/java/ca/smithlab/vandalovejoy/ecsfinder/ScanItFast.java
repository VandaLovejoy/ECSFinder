package ca.smithlab.vandalovejoy.ecsfinder;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import static ca.smithlab.vandalovejoy.ecsfinder.ECSFinder.runExternalCommand;
import static java.lang.Math.log10;

public class ScanItFast implements Runnable {

    static boolean VERBOSE = false;
    private final String[] key;
    private final ArrayList<String[]> associativeList;

    private static String Path;
    private static String SSZBINARY;

    // Stats thresholds
    private double SSZR_THRESHOLD = -3.0;   // SISSIz threshold
    private int GAP_THRESHOLD = 50;         // Gap threshold

    // Some global measures we compute
    private double shannon       = 0.0;
    private double sampled_sd    = 0.0;
    private double sampled_MFE   = 0.0;
    private double zscore        = 0.0;
    private double sci           = 0.0;

    // Reverse complement measures
    private double sampled_sd_rc   = 0.0;
    private double sampled_MFE_rc  = 0.0;
    private double zscore_rc       = 0.0;
    private double sci_rc          = 0.0;

    /**
     * We track how many threads are currently writing to the CSV.
     * This is incremented before writing, and decremented after.
     */
    private static final AtomicInteger activeCsvWriters = new AtomicInteger(0);

    public ScanItFast(ArrayList<String[]> associativeList,
                      String[] key,
                      File outputDir,
                      String SSZBINARY,
                      boolean VERBOSE) {
        // We'll store Path with "/aln" subdir
        ScanItFast.Path      = outputDir + "/aln/";
        ScanItFast.SSZBINARY = SSZBINARY;
        ScanItFast.VERBOSE   = VERBOSE;
        this.associativeList = associativeList;
        this.key             = key;
    }

    @Override
    public void run() {
        if (VERBOSE) {
            System.out.println("- - -> Starting Scan");
        }


        Map<Character, Integer> letterMap = new HashMap<>();
        letterMap.put('A', 0);
        letterMap.put('T', 1);
        letterMap.put('C', 2);
        letterMap.put('G', 3);
        letterMap.put('N', 4);
        letterMap.put('-', 5);
        Map<Character, Integer> letterMapRC = new HashMap<>();
        letterMapRC.put('A', 1);
        letterMapRC.put('T', 0);
        letterMapRC.put('C', 3);
        letterMapRC.put('G', 2);
        letterMapRC.put('N', 4);
        letterMapRC.put('-', 5);
        // remove identical rows or those with too many gaps & N's
        Iterator<String[]> iter = associativeList.iterator();

        ArrayList<int[]> intTab = new ArrayList<>();
        ArrayList<int[]> intTabRC = new ArrayList<>();
        ArrayList<String> UniqueNames = new ArrayList<>();

// Check if the associative list is empty
        if (associativeList.isEmpty()) {
            System.err.println("Error: The associativeList is empty.");
            return;
        }

// Iterate over the associativeList
        while (iter.hasNext()) {
            String[] line = iter.next();
            String sequence = line[1].toUpperCase();
            int[] seqToInt = new int[sequence.length()];
            int[] seqToIntRC = new int[sequence.length()];

            // Convert the sequence characters to integers and their reverse complement
            for (int i = 0; i < sequence.length(); i++) {
                if (letterMap.get(sequence.charAt(i)) == null) {
                    System.err.println("Error: Invalid character in sequence: " + sequence.charAt(i));
                    continue;
                }
                // Convert char to int in the sequence
                seqToInt[i] = letterMap.get(sequence.charAt(i));
                seqToIntRC[sequence.length() - i - 1] = letterMapRC.get(sequence.charAt(i));
            }

            // Split species name and create unique name
            String[] species_part = line[0].split("_");
            String species_part_new = species_part[0] + "_" + species_part[1].split("\\.")[0];

            // If unique name not already present, add the name and corresponding sequences
            if (!UniqueNames.contains(species_part_new)) {
                UniqueNames.add(species_part_new);
                intTabRC.add(seqToIntRC);
                intTab.add(seqToInt);
            }
        }

// Remove rows where 50% or more of the elements are gaps (4) or Ns (5)
        int requiredCount = (int) (0.5 * intTab.get(0).length);

// Create new lists to store filtered rows
        ArrayList<int[]> newArray_gap = new ArrayList<>();
        ArrayList<int[]> newArrayRC_gap = new ArrayList<>();
        ArrayList<String> newUniqueNames = new ArrayList<>();

// Iterate through rows in intTab and intTabRC
        for (int i = 0; i < intTab.size(); i++) {
            int[] originalArray = intTab.get(i);
            int count4 = 0;
            int count5 = 0;

            // Count the number of gaps (4) and Ns (5) in the row
            for (int element : originalArray) {
                if (element == 4) {
                    count4++;
                }
                if (element == 5) {
                    count5++;
                }
            }

            // If fewer than 50% of the elements are gaps/Ns, keep the row and corresponding unique name
            if (count4 < requiredCount && count5 < requiredCount) {
                newArray_gap.add(originalArray);                // Keep row in intTab
                newArrayRC_gap.add(intTabRC.get(i));            // Keep corresponding row in intTabRC
                newUniqueNames.add(UniqueNames.get(i));         // Keep corresponding unique name
            }
        }

// Replace old lists with filtered ones
        intTab.clear();
        intTab.addAll(newArray_gap);

        intTabRC.clear();
        intTabRC.addAll(newArrayRC_gap);

        UniqueNames.clear();
        UniqueNames.addAll(newUniqueNames);

// Remove columns entirely made up of gaps (4) or Ns (5)
        ArrayList<Integer> columnsToRemove = new ArrayList<>();
        int numRows = intTab.size();
        if (numRows<= 3) {
            if (VERBOSE)
                System.out.println("Too many species with gappy sequences");
            return;
        }

        int numCols = intTab.get(0).length;

// Check each column
        for (int col = 0; col < numCols; col++) {
            boolean allGapsOrNs = true;

            for (int row = 0; row < numRows; row++) {
                if (intTab.get(row)[col] != 4 && intTab.get(row)[col] != 5) {
                    allGapsOrNs = false;
                    break;
                }
            }

            if (allGapsOrNs) {
                columnsToRemove.add(col);
            }
        }

// Create new ArrayLists without columns that are entirely made up of gaps/Ns
        ArrayList<int[]> newArray1 = new ArrayList<>();
        ArrayList<int[]> newArray2 = new ArrayList<>();

        for (int[] originalArray : intTab) {
            int[] modifiedArray = new int[numCols - columnsToRemove.size()];
            int colIndex = 0;

            for (int col = 0; col < numCols; col++) {
                if (!columnsToRemove.contains(col)) {
                    modifiedArray[colIndex] = originalArray[col];
                    colIndex++;
                }
            }

            newArray1.add(modifiedArray);
        }

        for (int[] originalArray : intTabRC) {
            int[] modifiedArray = new int[numCols - columnsToRemove.size()];
            int colIndex = 0;

            for (int col = 0; col < numCols; col++) {
                if (!columnsToRemove.contains(col)) {
                    modifiedArray[colIndex] = originalArray[col];
                    colIndex++;
                }
            }

            newArray2.add(modifiedArray);
        }

// Replace old arrays with the filtered ones
        intTab.clear();
        intTab.addAll(newArray1);

        intTabRC.clear();
        intTabRC.addAll(newArray2);

// Convert UniqueNames to an array if needed
        String[] nameTab = new String[UniqueNames.size()];
        nameTab = UniqueNames.toArray(nameTab);

// Check if enough sequences are left
        int goodSeqs = intTab.size();

        if (intTab.size() <= 3) {
            if (VERBOSE)
                System.out.println("-> Not Enough seqs ");
            return;
        }

        if (intTab.size() <= 3) {
            if(VERBOSE)
                System.out.println("-> Not Enough seqs in this window!");
            return;
        }

       // if (!(nameTab[0].contains("homo"))) {
      //      if(VERBOSE)
      //          System.out.println("-> No human in this alignment block");
       //     return;
     //   }



        /*********************************************************************
         calculate stats						*
         *********************************************************************/
        if (VERBOSE)
            System.out.println("- - -> calculating statistics");
        double uniqueSeqs = goodSeqs;
        int outCols = intTab.get(0).length; // Change last variable if CLUSTAL properties change
        double[] stats = new double[6];
        double[] totalChars;
        // boolean[] isNotUnique = new boolean[outCols]; // Update to match the number of columns
        double[] chars;
        double[] column = new double[outCols];

        // Calculate MPI
        for (int k = 0; k < outCols; k++) {
            double identicalNuc = 0.0;
            double totalNuc = 0.0;
            double[][] stats1 = new double[6][6];


            for (int i = 0; i < goodSeqs; i++) {
                for (int j = i + 1; j < goodSeqs; j++) {
                    int charI = intTab.get(i)[k];
                    int charJ = intTab.get(j)[k];

                    if (isValidNucleotide(charI) && isValidNucleotide(charJ)) {
                        stats1[charI][charJ] += 1.0;
                        totalNuc += 1.0;

                        if (charI == charJ) {
                            identicalNuc += 1.0;
                        }
                    }
                }
            }
            if (totalNuc > 0) {
                column[k] = identicalNuc / totalNuc;
            } else {
                column[k] = 0.0; // Avoid division by zero
            }

        }

        double sum = 0.0;
        for (double value : column) {
            sum += value;
        }

        double newMPI = sum / column.length;


        totalChars = new double[]{0.0, 0.0, 0.0, 0.0, 0.0};

        for (int i = 0; i < intTab.get(0).length; i++) {
            chars = new double[]{0.0, 0.0, 0.0, 0.0, 0.0};
            for (int j = 0; j < goodSeqs; j++) {
                if (intTab.get(j)[i] == 5) {
                    chars[4] += 1.0;
                    totalChars[4] += 1.0;
                } else {
                    chars[intTab.get(j)[i]] += 1.0;
                    totalChars[intTab.get(j)[i]] += 1.0;
                }
            }
            for (int z = 0; z != 5; z++) {
                double probz = chars[z] / uniqueSeqs;
                shannon = (chars[z] == 0) ? shannon + 0 : shannon + probz * (Math.log(probz) / Math.log(2));
            }

        }
        Map<Integer, Character> reverseMap = new HashMap<>();
        reverseMap.put(0, 'A');
        reverseMap.put(1, 'T');
        reverseMap.put(2, 'C');
        reverseMap.put(3, 'G');
        reverseMap.put(4, 'N');
        reverseMap.put(5, '-');

        // prepare clustalw file
        if (VERBOSE)
            System.out.println("- -> preparing Clustal format");
        String[] outAln = new String[goodSeqs];
        String[] outAlnRC = new String[goodSeqs];
        int iterate = 0;
        for (int seq = 0; seq < intTab.size(); seq++) { //removed x < goodseqs

            outAln[iterate] = nameTab[seq].substring(0, Math.min(nameTab[seq].length(), 20));
            outAlnRC[iterate] = nameTab[seq].substring(0, Math.min(nameTab[seq].length(), 20));
            for (int i = 0; i != 25 - Math.min(nameTab[seq].length(), 20); i++) {
                outAln[iterate] = outAln[iterate] + " ";
                outAlnRC[iterate] = outAlnRC[iterate] + " ";
            }
            for (int i = 0; i != intTab.get(0).length; i++) {

                outAln[iterate] = outAln[iterate] + reverseMap.get(intTab.get(seq)[i]);

                outAlnRC[iterate] = outAlnRC[iterate] + reverseMap.get(intTabRC.get(seq)[i]);

            }
            outAln[iterate] = outAln[iterate] + "\n";
            outAlnRC[iterate] = outAlnRC[iterate] + "\n";
            iterate++;


        }

        stats[0] = newMPI * 100;                                                                        // Mean Pairwise ID
        double var = 0.0;
        for (double v : column) {
            var = var + Math.pow(v - newMPI, 2);
        }
        double standard = Math.sqrt(var / (column.length - 1));
        stats[1] = standard;                        // Variance
        stats[2] = -1 * shannon / ((double) outCols);       // Normalized Shanon entropy
        stats[3] = 100 * (totalChars[2] + totalChars[3]) / (totalChars[0] + totalChars[1] + totalChars[2] + totalChars[3]);       // GC content
        stats[4] = 100 * totalChars[4] / (outCols * goodSeqs);// GAP content


        // save BED coords
        if (VERBOSE)
            System.out.println("- -> Calculating BED coords ");
        String BedFile = key[0] + "\t";
        BedFile = BedFile + key[1] + "\t" + key[2] + "\t";


        BedFile = BedFile + goodSeqs + "_" + ((double) (int) (10 * stats[0]) / 10) + "_"      // MPI
                + ((double) (int) (10 * stats[1]) / 10) + "_"            // STDEV
                + ((double) (int) (100 * stats[2]) / 100) + "_"                // SHANNON
                + ((double) (int) (10 * stats[3]) / 10) + "_"                  //      GC
                + ((double) (int) (10 * stats[4]) / 10);                     // GAPS


        if (VERBOSE)
            System.out.println("Pre SISSIz bed file: \n" + " " + BedFile);
        int random = (int) ((double) 10000 * Math.random());
        File Aln = new File(Path  + BedFile.replaceAll("\t", "_") + ".aln." + random),    //
                AlnRC = new File(Path  + BedFile.replaceAll("\t", "_") + "rc.aln." + random);  //
        // v v v v v v v v    INCLUSION STATS     v v v v v v v v v v v v v
        //MPI greater or equal than 50 and Gap content smaller or equal than 50
        if (stats[4] <= GAP_THRESHOLD && stats[0] >= 50) {
            // Write Sequences to ALN Format
            try (BufferedWriter WriteClustal = new BufferedWriter(new FileWriter(Aln));
                 BufferedWriter WriteClustalRC = new BufferedWriter(new FileWriter(AlnRC))) {
                WriteClustal.write("CLUSTAL format \n\n");
                WriteClustalRC.write("CLUSTAL format \n\n");

                for (int y = 0; y != goodSeqs; y++) {

                    WriteClustal.write(outAln[y]);
                    WriteClustalRC.write(outAlnRC[y]);
                }
                WriteClustal.close();
                WriteClustalRC.close();
            } catch (IOException Err) {
                if (VERBOSE) {
                    System.err.println("Couldn't write clustal file!");
                }
                Err.printStackTrace();
                Aln.delete();
                AlnRC.delete();
                return;
            }



        } else {
            if (VERBOSE) {
                System.out.println("---> rejected alignment");
                System.out.println("     outcols = " + outCols + "\tuniqueseqs = " + uniqueSeqs +
                        "\tGAPS = " + stats[4] + "\n    PID = " + stats[0]);
                if (stats[0] < 5)
                    System.out.println("-----> SUPER LOW PID");
            }
            Aln.delete();
            AlnRC.delete();
            return;
        }
        String FinalBedFile,
                FinalBedFileRC,
                Antisense = (key[3].equals("+")) ? "-" : "+";


        //***************** 	SISSIz scan & parse		******************
        String[] SissizOutTab = new String[12];


        try {

            SissizOutTab = ScanSSZ(Path, BedFile, random);

            if (SissizOutTab == null) { // timeout
                Aln.delete();
                return;
            }
        } catch (IOException | InterruptedException Err) {
            Err.printStackTrace();
            System.err.println("ScanSSZ failed with ");

            Aln.delete();
            return;
        }
        // delete empty alignments
        if (SissizOutTab == null || SissizOutTab[12] == null) {
            Aln.delete();
            return;
        } else {
            FinalBedFile = BedFile + "_" + (int) (Double.parseDouble(SissizOutTab[12]) * -100) + "_" + key[3];
            // delete low scoring alignments
            if (Double.parseDouble(SissizOutTab[12]) > SSZR_THRESHOLD) {

                Aln.delete();
                return;


            } else {
                sci= Double.parseDouble(SissizOutTab[7]);
                sampled_MFE = Double.parseDouble(SissizOutTab[10]);
                sampled_sd = Double.parseDouble(SissizOutTab[11]);
                zscore = Double.parseDouble(SissizOutTab[12]);
                String fileNameBed = FinalBedFile.replace("\t", "_");

                File NewFile = new File(Path+fileNameBed+".aln");
                int file_count = 0;
                while (NewFile.exists()) {
                    file_count++;
                    NewFile = new File(Path+fileNameBed + ".aln_" + file_count);
                }
                boolean result = Aln.renameTo( NewFile);
                // Run RNAalifold
                String clustalFilePath = NewFile.getAbsolutePath();

                runRNAalifold(clustalFilePath,fileNameBed);

                boolean rScapeResult = runRScape(Path+fileNameBed + ".stk");

// If R-scape returns false, skip alignment
                if (!rScapeResult) {
                    System.err.println("Skipping alignment because R-scape failed or encountered a fatal error.");

                    return;
                }
                String cleanFileNameBed = fileNameBed.replace('.', '_');




                // Default E-value = large number (treat as "not significant" if missing)
                double eval = 9999.0;
                File helixCovFile = new File(Path + cleanFileNameBed + "_0001.helixcov");
                if (helixCovFile.exists()) {
                    eval = FilterOutput.processFile(helixCovFile.getAbsolutePath(), "E-value: ");
                }
                // Default covarying base-pairs = 0 (assume none found)
                double cov = 0.0;
                File powerFile = new File(Path + cleanFileNameBed + "_0001.power");
                if (powerFile.exists()) {
                    cov = FilterOutput.processFile(powerFile.getAbsolutePath(), "# BPAIRS observed to covary ");
                }
                double[] energies = FilterOutput.processTxtFile(Path+fileNameBed+ ".txt");
                int len_prediction =Integer.valueOf(key[2])-Integer.valueOf(key[1]);
                double[] array_variates= new double[]{
                        energies[0] / len_prediction,
                        energies[1] / len_prediction,
                        log10(eval),
                        cov,
                        ((double) (int) (10 * stats[0]) / 10),
                        sampled_MFE / len_prediction,
                        sampled_sd,
                        zscore,
                        sci,
                };
                
                String path_csv = ECSFinder.OUT_PATH+"/csv/"+fileNameBed+".csv";

                // Instead of constructing path_csv per alignment, we just do:
                try {
                    writeFeaturesToCSV(array_variates, fileNameBed + ".aln");
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }


                deleteFilesWithPrefix(Path, cleanFileNameBed );
                deleteFilesWithPrefix(Path, fileNameBed );

            }

            
        }
        // * * * * * *  now for the RC  * * * * * *
        try {
            SissizOutTab = ScanSSZ(Path, BedFile + "rc", random);
            if (SissizOutTab == null) {
                AlnRC.delete();
            }
        } catch (IOException | InterruptedException Err) {
            Err.printStackTrace();
            System.err.println("ScanSSZ failed in RC with ");

            AlnRC.delete();
        }
        if (SissizOutTab == null || SissizOutTab[12] == null) {
            AlnRC.delete();
        } else {
            FinalBedFileRC = BedFile + "_" + (int) (Double.parseDouble(SissizOutTab[12]) * -100) + "_" + Antisense;
            // delete low scoring alignments
            if (Double.parseDouble(SissizOutTab[12]) > SSZR_THRESHOLD) {
                AlnRC.delete();

                //    System.out.println(FinalBedFileRC.replaceAll("_", "\t"));
            } else {
                sci_rc= Double.parseDouble(SissizOutTab[7]);
                sampled_MFE_rc = Double.parseDouble(SissizOutTab[10]);
                sampled_sd_rc = Double.parseDouble(SissizOutTab[11]);
                zscore_rc = Double.parseDouble(SissizOutTab[12]);
                String fileNameBedRc = FinalBedFileRC.replace("\t", "_");


                File NewFile = new File( Path+fileNameBedRc + ".aln");
                int file_count = 0;
                while (NewFile.exists()) {
                    file_count++;
                    NewFile = new File(Path+fileNameBedRc +  ".aln_" + file_count);
                }
                boolean result = AlnRC.renameTo(NewFile);
                // Run RNAalifold
                String clustalFilePath = NewFile.getAbsolutePath();



                runRNAalifold(clustalFilePath,fileNameBedRc);

                // Run R-scape

                boolean rScapeResult = runRScape(Path+fileNameBedRc + ".stk");

// If R-scape returns false, skip alignment
                if (!rScapeResult) {
                    System.err.println("Skipping alignment because R-scape failed or encountered a fatal error.");

                    return;
                }

                String cleanFileNameBedRC = fileNameBedRc.replace('.', '_');
                // Default E-value = large number (treat as "not significant" if missing)
                double eval = 9999.0;
                File helixCovFile = new File(Path + cleanFileNameBedRC + "_0001.helixcov");
                if (helixCovFile.exists()) {
                    eval = FilterOutput.processFile(helixCovFile.getAbsolutePath(), "E-value: ");
                }

                // Default covarying base-pairs = 0 (assume none found)
                double cov = 0.0;
                File powerFile = new File(Path + cleanFileNameBedRC + "_0001.power");
                if (powerFile.exists()) {
                    cov = FilterOutput.processFile(powerFile.getAbsolutePath(), "# BPAIRS observed to covary ");
                }
                double[] energies = FilterOutput.processTxtFile(Path+fileNameBedRc+ ".txt");
                int len_prediction =Integer.valueOf(key[2])-Integer.valueOf(key[1]);
                double[] array_variates= new double[]{

                        energies[0] / len_prediction,
                        energies[1] / len_prediction,
                        log10(eval),
                        cov,
                        ((double) (int) (10 * stats[0]) / 10),
                        sampled_MFE_rc / len_prediction,
                        sampled_sd_rc,
                        zscore_rc,
                        sci_rc

                };

               
                String path_csv = ECSFinder.OUT_PATH+"/csv/"+fileNameBedRc+".csv";
                // Instead of constructing path_csv per alignment, we just do:
                try {
                    writeFeaturesToCSV(array_variates, fileNameBedRc + ".aln");
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }

                deleteFilesWithPrefix(Path, cleanFileNameBedRC );
                deleteFilesWithPrefix(Path, fileNameBedRc );
            }
            return;
        }

    }

    /*********************************************************************
     SISSIz scan & parse						*
     //*********************************************************************/
    // sissiz-di       cluster.109999_step.aln  8       150     0.8759  0.8542  0.0094  -13.88  -8.20   3.48    -1.63
    /** SISSIz call and parse */
    protected static String[] ScanSSZ(String path, String bedFile, int id)
            throws IOException, InterruptedException {
        String name = path + bedFile.replaceAll("\t", "_") + ".aln." + id;

        List<String> commandList = new ArrayList<>();
        commandList.add(SSZBINARY);
        commandList.add("-j");
        commandList.add("-t");
        commandList.add("--sci");
        commandList.add(name);

        long timeoutMs = 300_000;  // 5 minutes
        List<String> outputLines = runExternalCommand(commandList, new File(Path), timeoutMs, VERBOSE);

        // We expect something that starts with "TREE..."
        String[] sissizOutTab = new String[12];
        for (String line : outputLines) {
            if (line != null && line.startsWith("TREE")) {
                String[] parts = line.split(";");
                if (parts.length > 1 && parts[1].trim().startsWith("sissiz")) {
                    if (VERBOSE) {
                        System.out.println("SISSIz output: " + parts[1]);
                    }
                    sissizOutTab = parts[1].trim().split("\\s+");
                    break;
                }
            }
        }
        return sissizOutTab;
    }

    /** Set SISSIz threshold */
    public void setSszR(double newValue) {
        SSZR_THRESHOLD = newValue;
    }

    /** Set gap threshold */
    public void setGap(int newGap) {
        GAP_THRESHOLD = newGap;
    }

    /** Simple check for valid nucleotides */
    private boolean isValidNucleotide(int c) {
        // Valid nucleotides are 0,1,2,3,5 (not 4 which is 'N')
        return (c >= 0 && c <= 5) && (c != 4);
    }

    /**
     * Run RNAalifold
     */
    private void runRNAalifold(String clustalFile, String noExt) {
        List<String> cmd = Arrays.asList(
                ECSFinder.RNAALIFOLD,
                "--noLP",
                "-r",
                "--noPS",
                "--aln-stk=" + noExt,
                "--id-prefix=" + noExt,
                clustalFile
        );

        long timeoutMs = 300_000;  // 5 minutes
        try {
            List<String> outputLines = runExternalCommand(cmd, new File(Path), timeoutMs, VERBOSE);
            File outputTextFile = new File(Path + noExt + ".txt");

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputTextFile))) {
                for (String line : outputLines) {
                    writer.write(line);
                    writer.newLine();
                }
            }
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }

    /**
     * Run R-scape
     * @return true if no fatal errors
     */
    private boolean runRScape(String stkFile) {
        List<String> cmd = Arrays.asList("R-scape", "--lancaster", "--nofigures", "-s", stkFile);
        long timeoutMs = 300_000;

        try {
            List<String> outputLines = runExternalCommand(cmd, new File(Path), timeoutMs, VERBOSE);

            boolean fatalErrorFound = false;
            for (String line : outputLines) {
                if (line.contains("Fatal exception") || line.contains("IncompleteGamma")) {
                    System.err.println("R-scape encountered a numerical issue, skipping alignment.");
                    fatalErrorFound = true;
                }
            }
            return !fatalErrorFound;
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * Deletes files that start with a prefix
     */
    public static void deleteFilesWithPrefix(String directory, String prefix) {
        if (prefix == null || prefix.isEmpty()) {
            System.err.println("Prefix is null or empty. Doing nothing to avoid accidental bulk delete.");
            return;
        }
        try (Stream<Path> paths = Files.list(Paths.get(directory))) {
            paths.filter(Files::isRegularFile)
                    .filter(path -> {
                        String fileName = path.getFileName().toString();
                        // Must start with the prefix
                        return fileName.startsWith(prefix) && fileName.startsWith("aln");
                    })
                    .forEach(path -> {
                        try {
                            Files.deleteIfExists(path);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    });
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Synchronized method to write data to the single CSV file
     */
    private static void writeFeaturesToCSV(double[] data, String nameFile) throws IOException {

        activeCsvWriters.incrementAndGet();

        String csvPath = ECSFinder.SINGLE_CSV_PATH; // Single CSV path

        synchronized (ECSFinder.CSV_LOCK) {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(csvPath, true))) {

                // e.g. strip any path portion from "nameFile"
                String fileName = nameFile.substring(nameFile.lastIndexOf('/') + 1);

                // Write: fileName, then the data columns
                writer.write(fileName + ",");
                for (int i = 0; i < data.length; i++) {
                    writer.write(String.valueOf(data[i]));
                    if (i < data.length - 1) {
                        writer.write(",");
                    }
                }
                writer.newLine();
            }
        }
        // Indicate we're done writing
        activeCsvWriters.decrementAndGet();
    }

    /**
     * If you want to call this in ECSFinder or after all tasks are done,
     * it will block until activeCsvWriters == 0, meaning no threads are writing.
     */
    public static void waitForCsvWritersToFinish() {
        while (activeCsvWriters.get() > 0) {
            try {
                Thread.sleep(100); // Sleep 100ms, then check again
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                break;
            }
        }
    }

}

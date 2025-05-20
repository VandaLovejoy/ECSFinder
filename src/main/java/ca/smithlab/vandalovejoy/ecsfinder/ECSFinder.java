package ca.smithlab.vandalovejoy.ecsfinder;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ECSFinder {

    // Use any object as a lock
    public static final Object CSV_LOCK = new Object();
    private static String RSCRIPT = "";

    public static String SINGLE_CSV_PATH = null;

    static int GAPS = 50, NTHREDS = 4;
    static boolean VERBOSE = false, MAFFT = false;
    static String FILENAME = "", OUT_PATH = "", dirProgram = "";
    static String SSZBINARY = "", ALIFOLDBINARY = "",
            RNAALIFOLD = "", R = "", RSCAPE = "", MAFFTBINARY = "";
    static double SSZR = -3.0;
    private static File path;  // Not used heavily, but left here if needed

    // 1) New global map to store the FASTA chunk coordinates for each block_part
    private static Map<String, CoordinateInfo> coordinateMap = new HashMap<>();

    public static void main(String[] args) throws IOException, InterruptedException {

        // usage info
        if (args.length == 0) {
            System.out.println("\n\t\t\t  version 1.0 \n" +
                    " ________    ______   ______    ________  _                 __                \n" +
                    "|_   __  | .' ___  |.' ____ \\  |_   __  |(_)               |  ]               \n" +
                    "  | |_ \\_|/ .'   \\_|| (___ \\_|   | |_ \\_|__   _ .--.   .--.| | .---.  _ .--.  \n" +
                    "  |  _| _ | |        _.____.    |  _|  [  | [ .-. |/ /'\\' |/ /__\\[ /'\\] \n" +
                    " _| |__/ |\\ .___.'\\| \\____) |  _| |_    | |  | | | || \\__/  || \\__., | |     \n" +
                    "|________| .____ .' \\______.' |_____|  [___][___||__]'.__.;__]'.__.'[___]    \n" +
                    "\t SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES\n\n" +
                    "Reads a set of maf files, calculates stats, scans with SISSIz and R-scape, outputs bed coordinates of high-confidence predictions\n\n" +
                    "Usage: java ECSFinder [options] -o output/directory -i input.maf (last parameter must be -i)\n\n" +
                    "Options:\n" +
                    " -c int number of CPUs for calculations (default 4)\n" +
                    " -g int max gap percentage of sequences for 2D prediction (default 50)\n" +
                    " -sszr double report SISSIz+RIBOSUM hits below this Z-score (default -3.0)\n" +
                    " -mafft realign aln using MAFFT (default FALSE)\n" +
                    " -v verbose (messy but detailed) output\n");
            System.exit(0);
        }

        // parse arguments
        parseArguments(args);

        // Now that OUT_PATH is known, define the single CSV path
        SINGLE_CSV_PATH = OUT_PATH + "/final.csv";

        // (Optional) Write header exactly once at start
        writeHeader(SINGLE_CSV_PATH);

        // get binary paths
        setBinaryPaths();
        // preprocess MAF file using MergeNFilter
        preprocessMafFiles();
        // run RNALalifold and process results
        runRNALalifoldAndProcessResults();
        // After feature CSV is complete, run R predictions
        String predictionsCsv = OUT_PATH + "/predictions.csv";
        List<Double> probabilities = callRScript(
                SINGLE_CSV_PATH,
                predictionsCsv
        );

        System.out.println("Predictions written to: " + predictionsCsv);
    }

    private static void writeHeader(String csvPath) {
        File csvFile = new File(csvPath);
        if (!csvFile.exists()) {
            synchronized (CSV_LOCK) {
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(csvFile, true))) {
                    writer.write("name_file,min_energy,pseudo_energy,log_min_evalue,covarying_bp,"
                            + "MPI,average_MFE_sample,sd_sample,zscore,sci");
                    writer.newLine();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private static void setBinaryPaths() throws IOException, InterruptedException {
        ALIFOLDBINARY = which("RNALalifold");
        SSZBINARY    = which("SISSIz");
        RNAALIFOLD   = which("RNAalifold");
        RSCRIPT      = which("Rscript");
        if (MAFFT) MAFFTBINARY = which("mafft-ginsi");
    }

    private static String which(String cmd) throws IOException, InterruptedException {
        ProcessBuilder pb = new ProcessBuilder("which", cmd);
        Process p = pb.start();
        try (BufferedReader r = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
            String path = r.readLine();
            p.waitFor();
            if (path == null) throw new IOException("Cannot find " + cmd);
            return path.trim();
        }
    }

    private static String getBinaryPath(String binaryName) throws IOException {
        List<String> command = Arrays.asList("which", binaryName);
        try {
            List<String> lines = runExternalCommand(command, null,10_000, VERBOSE);
            if (lines.isEmpty()) {
                throw new IOException("Cannot find " + binaryName + " in PATH");
            }
            return lines.get(0).trim();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted while locating " + binaryName, e);
        }
    }

    private static void parseArguments(String[] args) {
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-c":
                    NTHREDS = Integer.parseInt(args[++i]);
                    break;
                case "-g":
                    GAPS = Integer.parseInt(args[++i]);
                    break;
                case "-o":
                    OUT_PATH = System.getProperty("user.dir") + "/" + args[++i];
                    createDirectory(OUT_PATH);
                    dirProgram = System.getProperty("user.dir");
                    break;
                case "-v":
                    VERBOSE = true;
                    break;
                case "-mafft":
                    MAFFT = true;
                    break;
                case "-sszr":
                    SSZR = Double.parseDouble(args[++i]);
                    break;
                case "-i":
                    FILENAME = args[++i];
                    break;
                default:
                    System.err.println("Invalid argument: " + args[i]);
                    printUsageAndExit();
            }
        }

        // Additional argument validation
        if (FILENAME.isEmpty()) {
            System.err.println("Input MAF file (-i) is required.");
            printUsageAndExit();
        }
        if (OUT_PATH.isEmpty()) {
            System.err.println("Output directory (-o) is required.");
            printUsageAndExit();
        }
    }

    private static void printUsageAndExit() {
        System.out.println("Usage: java ECSFinder [options] -o output/directory -i input.maf");
        System.exit(1);
    }

    private static void createDirectory(String path) {
        File dir = new File(path);
        if (!dir.isDirectory()) {
            dir.mkdirs();
        }
    }

    private static void preprocessMafFiles() {
        File inputDir = new File(FILENAME);
        if (!inputDir.isDirectory()) {
            System.err.println("Error: Provided path is not a directory or does not exist.");
            System.exit(1);
        }

        File[] mafFiles = inputDir.listFiles((dir, name) -> name.endsWith(".maf"));
        if (mafFiles == null || mafFiles.length == 0) {
            System.err.println("Error: No .maf files found in the directory: " + FILENAME);
            System.exit(1);
        }

        String[] mafFilePaths = Arrays.stream(mafFiles).map(File::getAbsolutePath).toArray(String[]::new);

        try {
            MergeNFilter mergeNFilter = new MergeNFilter();
            mergeNFilter.process(mafFilePaths, OUT_PATH);
        } catch (IOException e) {
            System.err.println("Error processing MAF files: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static void runRNALalifoldAndProcessResults() throws IOException, InterruptedException {
        // Extract file path
        String inputFile = OUT_PATH + "/output.maf";
        String stockholmFolderPath = OUT_PATH + "/stockholm";
        File stockholmFolder = createStockholmFolder(stockholmFolderPath);

        // Handle MAFFT or direct MAF file processing
        if (!MAFFT) {
            runRNALalifold(inputFile);
        } else {
            processWithMafft(inputFile);  // MAFFT handling
        }

        // Process MAF file for alignment blocks
        try {
            processAlignmentBlocks(inputFile);
        } catch (ExecutionException e) {
            System.err.println("Error in executing task: " + e.getMessage());
            e.printStackTrace();
        }

        // Clean up the Stockholm folder
        deleteDirectoryRecursively(stockholmFolder);

        ScanItFast.waitForCsvWritersToFinish();
        System.out.println("All CSV writing operations have completed!");
    }

    private static File createStockholmFolder(String stockholmFolderPath) {
        File stockholmFolder = new File(stockholmFolderPath);
        if (!stockholmFolder.exists()) {
            stockholmFolder.mkdir();
        }
        return stockholmFolder;
    }

    public static void processWithMafft(String inputFile) throws IOException, InterruptedException {
        // Convert MAF to FASTA and realign with MAFFT
        convertMafToSeparateFastas(inputFile, OUT_PATH);  // Converts to smaller blocks with overlap

        // Realign each block individually to avoid memory issues
        File outputDir = new File(OUT_PATH + "/outputFastaDir");
        File[] fastaFiles = outputDir.listFiles((dir, name) -> name.endsWith(".fasta"));

        if (fastaFiles != null) {
            // Iterate through each FASTA block file
            for (File fastaFile : fastaFiles) {
                if(VERBOSE) {
                    System.out.println("Realigning file: " + fastaFile.getName());
                }
                // Realign each small block and save the result
                File realignedOutput = new File(fastaFile.getAbsolutePath().replace(".fasta", "_realigned.fasta"));
                realignSequences(fastaFile, outputDir);
                runRNALalifold(realignedOutput.getAbsolutePath());
            }
        }
    }

    private static void processAlignmentBlocks(String inputFile)
            throws IOException, InterruptedException, ExecutionException {
        BufferedReader readFile = new BufferedReader(new FileReader(inputFile));
        int blockAln = 0;
        StringBuilder temp = new StringBuilder();
        // Create the executor service
        ExecutorService multiThreads = Executors.newFixedThreadPool(NTHREDS);
        List<Future<?>> futures = new ArrayList<>();

        String path_aln = OUT_PATH + "/aln/";
        createDirectory(path_aln);

        String line;
        // Process the alignment blocks from the file
        while ((line = readFile.readLine()) != null) {
            if (line.length() >= 1 && line.charAt(0) != '#') {
                if (line.charAt(0) == 'a') {
                    blockAln++;
                } else if (line.charAt(0) == 's') {
                    temp.append(line).append("@");
                }
            } else if (!temp.toString().isEmpty()
                    && (temp.toString().split("@").length >= 3)
                    && line.equals("")) {
                // Submit the task for processing a block and store the Future
                iterateStockholm(temp.toString(), blockAln, futures, multiThreads);
                temp = new StringBuilder();
            } else {
                temp = new StringBuilder();
            }
        }
        readFile.close();

        for (Future<?> f : futures) {
            try {
                f.get();
            } catch (ExecutionException e) {
                // Log or handle the actual cause: e.getCause()
                e.printStackTrace();
            }
        }

        // Shut down the executor
        multiThreads.shutdown();
        if (!multiThreads.awaitTermination(10, TimeUnit.MINUTES)) {
            System.err.println("Executor did not terminate in time, forcing shutdown.");
            multiThreads.shutdownNow();
        }
    }

    private static void iterateStockholm(String block, int blockAln,
                                         List<Future<?>> futures,
                                         ExecutorService multiThreads) throws IOException {
        ArrayList<String[]> associativeList = new ArrayList<>();
        String[] mafTabTemp = block.split("@")[0].split("\\s+");

        // Create a prefix: "alifold_{blockAln}_"
        File dir = new File(OUT_PATH + "/stockholm");

        String prefix = "alifold_" + blockAln + "_";
        File[] stkFiles = dir.listFiles((dir1, name) ->
                name.startsWith(prefix) && name.endsWith(".stk")
        );

        if (stkFiles == null || stkFiles.length == 0) {
            // Instead of exiting, just warn and return
            System.out.println("No .stk files found for alignment block " + blockAln
                    + ". Skipping this block.\n");
            return;  // let the code proceed to the next blockAln
        }

        // Sort them by sub-block number (so e.g. alifold_1_0001.stk < alifold_1_0002.stk)
        Arrays.sort(stkFiles, (f1, f2) -> {
            int sub1 = parseSubBlockNumber(f1.getName(), blockAln);
            int sub2 = parseSubBlockNumber(f2.getName(), blockAln);
            return Integer.compare(sub1, sub2);
        });

        // For each .stk found, process it
        for (File stkFile : stkFiles) {
            if (VERBOSE) {
                System.out.println("Processing .stk for block #" + blockAln
                        + ": " + stkFile.getName());
            }
            processStockholmFile(stkFile, mafTabTemp, associativeList, futures, multiThreads);
        }
    }

    /**
     * Parse the sub-block number from a filename like "alifold_1_0002.stk".
     *
     * @param stkName   Something like "alifold_1_0002.stk"
     * @param blockAln  The integer block number, e.g. 1
     * @return          The integer sub-block, e.g. 2
     */
    private static int parseSubBlockNumber(String stkName, int blockAln) {
        // e.g. pattern "alifold_1_0002.stk"
        String patternStr = "^alifold_" + blockAln + "_(\\d+)\\.stk$";
        Pattern p = Pattern.compile(patternStr);
        Matcher m = p.matcher(stkName);
        if (m.matches()) {
            return Integer.parseInt(m.group(1));
        }
        // If we can't parse, return something large or handle differently
        return Integer.MAX_VALUE;
    }

    private static void processStockholmFile(File stockholmFile, String[] mafTabTemp,
                                             ArrayList<String[]> associativeList,
                                             List<Future<?>> futures,
                                             ExecutorService multiThreads) {
        System.out.println("[DEBUG] Processing Stockholm: " + stockholmFile.getName());
        String result = "";
        if(MAFFT) {
            String regex = "alifold_(\\d+)_(\\d+)";
            Pattern pattern = Pattern.compile(regex);
            Matcher matcher = pattern.matcher(stockholmFile.getName());

            if (matcher.find()) {
                String firstPart = matcher.group(1);
                // secondPart with leading zeros removed
                String secondPart = matcher.group(2).replaceFirst("^0+", "");
                result = firstPart + "_" + secondPart;
            } else {
                System.out.println("Error: Filename does not match expected pattern: " + stockholmFile.getName());
                return;
            }
        } else {
            String regex = "alifold_(\\d+)";
            Pattern pattern = Pattern.compile(regex);
            Matcher matcher = pattern.matcher(stockholmFile.getName());

            if (matcher.find()) {
                String secondPart = matcher.group(1).replaceFirst("^0+", "");
                result = secondPart;
            }
        }

        try (BufferedReader reader = new BufferedReader(new FileReader(stockholmFile))) {
            String currentLine;
            String[] arrayName = new String[5];
            String gcReference = "", gcSScons = "";

            while ((currentLine = reader.readLine()) != null) {
                if (currentLine.startsWith("#=GF ID ")) {
                    arrayName = currentLine.split("[_.]");
                    associativeList = new ArrayList<>();
                } else if (currentLine.startsWith("#=GC RF")) {
                    gcReference = extractValue(currentLine);
                } else if (currentLine.startsWith("#=GC SS_cons")) {
                    gcSScons = extractValue(currentLine);
                } else if (!currentLine.startsWith("#")
                        && !currentLine.equals("")
                        && !currentLine.startsWith("//")) {
                    associativeList.add(processSpeciesLine(currentLine));
                }

                if (!associativeList.isEmpty() && currentLine.startsWith("//")) {
                    processMotif(mafTabTemp, arrayName, associativeList,
                            gcReference, gcSScons, futures, multiThreads, result);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void processMotif(String[] mafTabTemp, String[] arrayName,
                                     ArrayList<String[]> associativeList,
                                     String gcReference, String gcSScons,
                                     List<Future<?>> futures,
                                     ExecutorService multiThreads,
                                     String result) {
        int[] cordMotif;
        if (MAFFT) {
            cordMotif = getRealCoordinates(Integer.parseInt(arrayName[4]),
                    mafTabTemp,
                    associativeList.get(0)[1],
                    result);

        } else {
            cordMotif = getRealCoordinates(Integer.parseInt(arrayName[3]),
                    mafTabTemp,
                    associativeList.get(0)[1],
                    result);
        }

        String loci = Arrays.toString(cordMotif);
        String chrom = mafTabTemp[1].substring((mafTabTemp[1].lastIndexOf(".") + 1));
        String lociChrm;
        if (MAFFT) {
            lociChrm = chrom + ", " + loci.substring(1, loci.length() - 1)
                    + ", " + mafTabTemp[4] + ", "
                    + arrayName[4] + ", " + arrayName[5]
                    + ", " + gcReference + ", " + gcSScons;
        } else {
            lociChrm = chrom + ", " + loci.substring(1, loci.length() - 1)
                    + ", " + mafTabTemp[4] + ", "
                    + arrayName[3] + ", " + arrayName[4]
                    + ", " + gcReference + ", " + gcSScons;
        }

        String[] arrayLociChrm = lociChrm.split(", ");
        // Skip short alignments
        if (Integer.parseInt(arrayLociChrm[2]) - Integer.parseInt(arrayLociChrm[1]) < 50) {
            return;
        }
        System.out.println("[DEBUG] Creating ScanItFast job for block = " + Arrays.toString(arrayName)
                + " on file result key = " + result);
        // Create the ScanItFast.java task
        ScanItFast aln = new ScanItFast(associativeList, arrayLociChrm,
                new File(OUT_PATH ), SSZBINARY, VERBOSE);
        aln.setSszR(SSZR);
        aln.setGap(GAPS);
        Future<?> future = multiThreads.submit(aln);
        futures.add(future);
    }

    private static String getBlockName(int blockAln) {
        String lastDigit = String.valueOf(blockAln);
        int numbZeros = 4 - lastDigit.length();
        StringBuilder finalName = new StringBuilder();
        for (int i = 0; i < numbZeros; i++) {
            finalName.append('0');
        }
        finalName.append(lastDigit);
        return finalName.toString();
    }

    private static String extractValue(String line) {
        String[] lineReference = line.split(" ");
        return lineReference[lineReference.length - 1];
    }

    private static String[] processSpeciesLine(String line) {
        String[] species = line.split(" ", 2);
        species[1] = species[1].trim();
        if (line == null || line.trim().isEmpty()) {
            return null;  // Handle invalid lines
        }
        return species;
    }

    // Method to delete a directory and all its contents, regardless of whether it contains files
    private static void deleteDirectoryRecursively(File directory) {
        if (directory.exists()) {
            File[] files = directory.listFiles();
            if (files != null) {
                for (File file : files) {
                    if (file.isDirectory()) {
                        deleteDirectoryRecursively(file);
                    } else {
                        file.delete();
                    }
                }
            }
            directory.delete();
        }
    }

    private static List<Double> callRScript(String inputCsv, String outputCsv)
            throws IOException, InterruptedException {
        String scriptDir = getRScriptPath();
        String rScriptFile = Paths.get(scriptDir, "predictions_ECSFinder.R").toString();
        List<String> cmd = Arrays.asList(
                RSCRIPT, rScriptFile, inputCsv, outputCsv, scriptDir
        );
        ProcessBuilder pb = new ProcessBuilder(cmd).directory(new File(scriptDir));
        Process p = pb.start();
        if (VERBOSE) inheritIO(p.getInputStream(), System.out);
        if (VERBOSE) inheritIO(p.getErrorStream(), System.err);
        if (p.waitFor() != 0) throw new IOException("Rscript failed");
        File out = new File(outputCsv);
        if (!out.exists()) throw new FileNotFoundException(outputCsv);
        List<Double> preds = new ArrayList<>();
        try (BufferedReader r = new BufferedReader(new FileReader(out))) {
            r.readLine(); String line;
            while ((line = r.readLine()) != null) preds.add(Double.parseDouble(line.split(",")[1]));
        }
        return preds;
    }

    // Utility to forward a stream
    private static void inheritIO(final InputStream src, final PrintStream dest) {
        new Thread(() -> {
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(src))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    dest.println(line);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }).start();
    }
    private static String getRScriptPath() throws IOException {
        File jarFile = new File(ECSFinder.class
                .getProtectionDomain()
                .getCodeSource()
                .getLocation()
                .getPath());
        String jarDir = jarFile.getParent();
        if (jarFile.isFile()) {
            return new File(jarDir, "RF/").getAbsolutePath();
        } else {
            return new File("RF/").getAbsolutePath();
        }
    }

    // Keep track of block offsets
    private static Map<String, Integer> blockStartMap = new HashMap<>();

    // Map block_part -> Realigned Homo sapiens sequence
    private static Map<String, String> homoSapiensSequences = new HashMap<>();

    public static void convertMafToSeparateFastas(String mafFilePath, String outputDirPath) {
        String fastaOutput = "outputFastaDir";
        int blockCount = 0;
        int overlapLength = 299;
        int maxBlockSize = 5000;
        int startAln = 0;
        int lengthAln = 0;
        int endAln = 0;
        String orientation = "+";

        try (BufferedReader reader = new BufferedReader(new FileReader(mafFilePath))) {
            Files.createDirectories(Paths.get(outputDirPath + "/" + fastaOutput));
            Map<String, StringBuilder> speciesSequences = new LinkedHashMap<>();

            String line;
            int currentBlockLength = 0;

            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("#") || line.isEmpty()) {
                    continue;
                }
                if (line.startsWith("a")) {
                    // Starting a new alignment block
                    if (!speciesSequences.isEmpty()) {
                        blockCount++;
                        splitAndWriteBlocks(outputDirPath, blockCount, speciesSequences,
                                maxBlockSize, overlapLength,
                                startAln, orientation, lengthAln, endAln);
                        speciesSequences.clear();
                    }
                    currentBlockLength = 0;
                } else if (line.startsWith("s")) {
                    String[] tokens = line.split("\\s+");
                    String sequenceId = tokens[1];
                    String sequence = tokens[tokens.length - 1];
                    if (line.contains("homo")) {
                        startAln = Integer.parseInt(tokens[2]);
                        lengthAln = Integer.parseInt(tokens[3]);
                        orientation = tokens[4];
                        endAln = Integer.parseInt(tokens[5]);
                    }
                    speciesSequences
                            .computeIfAbsent(sequenceId, k -> new StringBuilder())
                            .append(sequence);
                    currentBlockLength = Math.max(currentBlockLength,
                            speciesSequences.get(sequenceId).length());
                }
            }
            // Write the last block
            if (!speciesSequences.isEmpty()) {
                blockCount++;
                splitAndWriteBlocks(outputDirPath, blockCount, speciesSequences,
                        maxBlockSize, overlapLength,
                        startAln, orientation, lengthAln, endAln);
            }
        } catch (IOException e) {
            System.err.println("Error processing the MAF file: " + e.getMessage());
        }
    }

    private static void splitAndWriteBlocks(String outputDirPath, int blockCount,
                                            Map<String, StringBuilder> speciesSequences,
                                            int maxBlockSize, int overlapLength,
                                            int startAln, String orientation,
                                            int lengthAln, int chromosomeLength) throws IOException {

        int blockPart = 1;
        boolean isSplit = false;

        int fullLength = speciesSequences.values().iterator().next().length();

        String homoSequence = null;
        String homoSpeciesId = null;
        String chromosomeNumber = null;
        for (String key : speciesSequences.keySet()) {
            if (key.toLowerCase().contains("homo")) {
                homoSequence = speciesSequences.get(key).toString();
                homoSpeciesId = key;
                String[] parts = homoSpeciesId.split("\\.");
                if (parts.length > 1) {
                    chromosomeNumber = parts[1];
                } else {
                    throw new IllegalArgumentException("Chromosome number not found in Homo sapiens ID: " + homoSpeciesId);
                }
                break;
            }
        }
        if (homoSequence == null) {
            throw new IllegalArgumentException("Homo sapiens sequence not found in speciesSequences.");
        }

        for (int i = 0; i < fullLength; i += (maxBlockSize - overlapLength)) {
            int blockStartIndex = i;
            int blockEndIndex = Math.min(i + maxBlockSize, fullLength);

            String homoSubSequence = homoSequence.substring(blockStartIndex, blockEndIndex);
            int nucleotidesUpToStart = homoSequence
                    .substring(0, blockStartIndex)
                    .replaceAll("-", "")
                    .length();
            int blockNucleotideLength = homoSubSequence
                    .replaceAll("-", "")
                    .length();

            int genomicStart, genomicEnd;
            if (orientation.equals("+")) {
                genomicStart = startAln + nucleotidesUpToStart + 1;
                genomicEnd = genomicStart + blockNucleotideLength - 1;
            } else if (orientation.equals("-")) {
                genomicEnd = chromosomeLength - startAln - nucleotidesUpToStart;
                genomicStart = genomicEnd - blockNucleotideLength + 1;
            } else {
                throw new IllegalArgumentException("Invalid orientation: " + orientation);
            }

            int adjustedGenomicStart = genomicStart;
            int adjustedGenomicEnd = genomicEnd;
            if (adjustedGenomicStart > adjustedGenomicEnd) {
                int temp = adjustedGenomicStart;
                adjustedGenomicStart = adjustedGenomicEnd;
                adjustedGenomicEnd = temp;
            }

            // Construct the FASTA filename including coords
            String blockFileName = outputDirPath + "/outputFastaDir/block_"
                    + blockCount + "_part_" + blockPart + "_chr"
                    + chromosomeNumber + "_" + adjustedGenomicStart + "_"
                    + adjustedGenomicEnd + "_" + orientation + ".fasta";

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(blockFileName))) {
                for (Map.Entry<String, StringBuilder> entry : speciesSequences.entrySet()) {
                    String speciesId = entry.getKey();
                    String sequence = entry.getValue().toString();
                    String subSequence = sequence.substring(blockStartIndex, blockEndIndex);
                    writer.write(">" + speciesId + "\n");
                    writer.write(subSequence + "\n");
                }
            }

            if (VERBOSE) {
                System.out.println("Wrote block part: " + blockFileName
                        + " starting at position: " + blockStartIndex);
            }

            blockPart++;
            if (i + maxBlockSize < fullLength) {
                isSplit = true;
            }
        }

        if (isSplit && VERBOSE) {
            System.out.println("Block " + blockCount + " was split into "
                    + (blockPart - 1) + " parts.");
        }
    }

    /**
     * realignSequences(...) reads the FASTA block filename, captures the
     * chromosome, start, end, and orientation, then puts them in coordinateMap.
     */
    public static void realignSequences(File inputFilePath, File outputFilePath)
            throws IOException, InterruptedException {
        String realignedFilePath = inputFilePath.getAbsolutePath()
                .replace(".fasta", "_realigned.fasta");

        // Regex to parse e.g. block_1_part_1_chr1_1000_2000_+.fasta
        String regex = "block_(\\d+)_part_(\\d+)_chr(\\w+)_(\\d+)_(\\d+)_([+-])(?:_\\w+)?\\.fasta";
        Pattern pattern = Pattern.compile(regex);
        Matcher matcher = pattern.matcher(inputFilePath.getAbsolutePath());

        String result = "";
        if (matcher.find()) {
            // block_(1) _part_(1) _chr(1) _(1000) _(2000) _(+)
            String blockNumber = matcher.group(1);
            String partNumber  = matcher.group(2);
            String chromName   = matcher.group(3);
            int    fastaStart  = Integer.parseInt(matcher.group(4));
            int    fastaEnd    = Integer.parseInt(matcher.group(5));
            String strand      = matcher.group(6);

            // This becomes the key used in getRealCoordinates
            result = blockNumber + "_" + partNumber;

            // Save these coords in coordinateMap
            CoordinateInfo coordInfo = new CoordinateInfo(chromName, fastaStart, fastaEnd, strand);
            coordinateMap.put(result, coordInfo);

            if (VERBOSE) {
                System.out.println("Parsed FASTA coords: block_part = " + result
                        + ", chr=" + chromName
                        + ", start=" + fastaStart
                        + ", end=" + fastaEnd
                        + ", strand=" + strand);
            }
        } else {
            System.out.println("No match found in realignSequences for file: "
                    + inputFilePath.getName());
        }

        // Next, run MAFFT
        List<String> command = Arrays.asList(
                MAFFTBINARY,
                "--quiet",
                "--thread", String.valueOf(NTHREDS),
                inputFilePath.getAbsolutePath()
        );
        ProcessBuilder pb = new ProcessBuilder(command);
        Process mafftProcess = pb.start();

        try (BufferedReader reAligned = new BufferedReader(new InputStreamReader(mafftProcess.getInputStream()));
             BufferedWriter writer = new BufferedWriter(new FileWriter(realignedFilePath))) {

            StringBuilder sequence = new StringBuilder();
            String speciesReal = "";
            String line;
            while ((line = reAligned.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (!speciesReal.isEmpty()) {
                        writer.write(">" + speciesReal + "\n");
                        writer.write(sequence.toString().toUpperCase() + "\n");

                        // If this is Homo sapiens, store the realigned seq
                        if (speciesReal.toLowerCase().contains("homo")) {
                            homoSapiensSequences.put(result, sequence.toString());
                        }
                        sequence = new StringBuilder();
                    }
                    speciesReal = line.substring(1);
                } else {
                    sequence.append(line);
                }
            }
            // Last sequence
            if (!speciesReal.isEmpty()) {
                writer.write(">" + speciesReal + "\n");
                writer.write(sequence.toString().toUpperCase() + "\n");

                if (speciesReal.toLowerCase().contains("homo")) {
                    homoSapiensSequences.put(result, sequence.toString());
                }
            }
        }

        int exitCode = mafftProcess.waitFor();
        if (exitCode != 0) {
            throw new IOException("MAFFT exited with error code: " + exitCode);
        }
        try (BufferedReader errorReader = new BufferedReader(
                new InputStreamReader(mafftProcess.getErrorStream()))) {
            String errorLine;
            while ((errorLine = errorReader.readLine()) != null) {
                if(VERBOSE) {
                    System.err.println("MAFFT-ERR: " + errorLine);
                }
            }
        }

        // Clean up
        File tempInputFile = new File(inputFilePath.getParent(), "temp_gap_stripped.fasta");
        if (tempInputFile.exists()) {
            tempInputFile.delete();
        }
    }

    private static void runRNALalifold(String inputFilePath)
            throws IOException, InterruptedException {
        System.out.println("[DEBUG] About to run RNALalifold on FASTA: " + inputFilePath);
        String outputFilePath = OUT_PATH + "/stockholm";
        List<String> command;

        if (inputFilePath.endsWith(".fasta")) {
            // We parse blockNumber & partNumber for the alifold name,
            // but the real coords were already saved in coordinateMap by realignSequences
            String regex = ".*/block_(\\d+)_part_(\\d+)_chr(\\w+)_(\\d+)_(\\d+)_([+-])(?:_realigned)?\\.fasta$";
            Pattern pattern = Pattern.compile(regex);
            Matcher matcher = pattern.matcher(inputFilePath);

            int blockNumber = 0;
            int partNumber  = 0;
            if (matcher.find()) {
                blockNumber = Integer.parseInt(matcher.group(1));
                partNumber  = Integer.parseInt(matcher.group(2));
            } else {
                System.err.println("Block/Part not found in the input file path: "
                        + inputFilePath + ". Skipping file.");
                return;
            }
            command = Arrays.asList(
                    ALIFOLDBINARY,
                    "--id-prefix=alifold_" + blockNumber,
                    "--id-start=" + partNumber,
                    "--noLP",
                    "--maxBPspan=300",
                    "--ribosum_scoring",
                    "--aln-stk",
                    inputFilePath
            );
        } else {
            // If not a FASTA, just do a simpler approach
            command = Arrays.asList(
                    ALIFOLDBINARY,
                    "--id-prefix=alifold",
                    "--noLP",
                    "--maxBPspan=300",
                    "--ribosum_scoring",
                    "--aln-stk",
                    inputFilePath
            );
        }

        ProcessBuilder pb = new ProcessBuilder(command);
        pb.directory(new File(outputFilePath));
        Process process = pb.start();

        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        BufferedReader errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));

        Thread outputThread = new Thread(() -> reader.lines().forEach(line -> {}));
        Thread errorThread = new Thread(() -> errorReader.lines().forEach(line -> {}));
        outputThread.start();
        errorThread.start();
        int exitCode = process.waitFor();
        outputThread.join();
        errorThread.join();

        if (exitCode != 0) {
            System.err.println("RNALalifold exited with error code: " + exitCode
                    + " for file: " + inputFilePath + ". Skipping.");
            return;
        }
        if (VERBOSE && inputFilePath.endsWith(".fasta")) {
            System.out.println("RNALalifold output saved to: " + outputFilePath);
        }
    }

    /**
     * The final method that calculates absolute coordinates.
     * We combine the offset of the motif in the realigned Homo sapiens sequence
     * with the chunk's actual genomic start/end from coordinateMap.
     */
    private static int[] getRealCoordinates(int start, String[] mafCord,
                                            String motifHuman,
                                            String blockPartKey) {

        CoordinateInfo info = coordinateMap.get(blockPartKey);
        if (info == null) {
            return oldMafBasedCoordinates(start, mafCord, motifHuman);
        }

        // 2) Count how many real bases appear up to 'start' in the realigned Homo sapiens row
        String homoAligned = homoSapiensSequences.get(blockPartKey);
        if (homoAligned == null) {
            return oldMafBasedCoordinates(start, mafCord, motifHuman);
        }
        String upToStart = homoAligned.substring(0, start-1);
        int offsetInChunk = upToStart.replaceAll("-", "").length();
        int motifLen = motifHuman.replaceAll("-", "").length();


        // 4) Combine offset with the chunk's known start or end (both are 1-based in your code)
        int finalStart, finalEnd;
        if (info.strand.equals("+")) {
            // 1-based + offsetInChunk
            finalStart = info.start + offsetInChunk ;  // no extra +1 needed
            finalEnd   = finalStart + motifLen - 1;   // also 1-based
        } else {
            // Negative strand block, still in the forward coordinate system
            finalEnd   = info.end - offsetInChunk;  // remove the +1
            finalStart = finalEnd - motifLen + 1;
        }


        if (finalStart > finalEnd) {
            int temp = finalStart;
            finalStart = finalEnd;
            finalEnd   = temp;
        }

        // 6) Return [start, end], still 1-based inclusive.
        return new int[]{ finalStart, finalEnd };
    }

    /**
     * If we didn't parse coords from FASTA, we do the old MAF-based approach.
     */
    /**
     * If we didn't parse coords from FASTA, we do the older MAF-based approach.
     */
    private static int[] oldMafBasedCoordinates(int startIndexInAlignment,
                                                String[] mafCord,
                                                String motifHuman) {

        // mafCord[2] = start
        // mafCord[3] = size
        // mafCord[4] = strand
        // mafCord[5] = srcSize
        int mafStart  = Integer.parseInt(mafCord[2]);
        int mafSize   = Integer.parseInt(mafCord[3]);
        int mafSrcLen = Integer.parseInt(mafCord[5]);
        String strand = mafCord[4];

        // The offset from the left side of the alignment block:
        // Count the aligned (non-gap) nucleotides in the "homo" row up to startIndexInAlignment.
        // Or however you have 'nuc' / 'nuclStockholm' computed.
        // For simplicity, I'll call it offset:
        int offset = countRealBases(mafCord[6].substring(0, startIndexInAlignment));

        // Convert block to forward coords
        int forwardStart, forwardEnd;
        if (strand.equals("+")) {
            forwardStart = mafStart;
            forwardEnd   = mafStart + mafSize - 1;
        } else {
            forwardStart = mafSrcLen - (mafStart + mafSize);
            forwardEnd   = mafSrcLen - mafStart - 1;
        }

        // Now place the motif offset within that block
        int motifLen = motifHuman.replaceAll("-", "").length();
        int finalStart = forwardStart + offset;
        int finalEnd   = finalStart + motifLen - 1;

        // Return [start, end] in a 0-based system. If you want 1-based, add +1
        return new int[]{ finalStart, finalEnd };
    }

    private static int countRealBases(String seq) {
        // Count A/C/G/T (non-gap) in 'seq'
        return seq.replaceAll("-", "").length();
    }


    /**
     * Helper class for storing the chunk's genomic coords from the FASTA filename.
     */
    public static class CoordinateInfo {
        public final String chrom;
        public final int start;
        public final int end;
        public final String strand;

        public CoordinateInfo(String chrom, int start, int end, String strand) {
            this.chrom = chrom;
            this.start = start;
            this.end   = end;
            this.strand = strand;
        }
    }

    public static void mergeLogFiles(String logDirPath, String finalCsvPath) {
        File logDir = new File(logDirPath);
        File[] logFiles = logDir.listFiles((dir, name) -> name.endsWith(".csv"));
        if (logFiles != null && logFiles.length > 0) {
            try (BufferedWriter finalWriter = new BufferedWriter(new FileWriter(finalCsvPath))) {
                finalWriter.write("name_file,min_energy,pseudo_energy,log_min_evalue,covarying_bp,"
                        + "MPI,average_MFE_sample,sd_sample,zscore,sci\n");
                for (File logFile : logFiles) {
                    try (BufferedReader logReader = new BufferedReader(new FileReader(logFile))) {
                        String line;
                        while ((line = logReader.readLine()) != null) {
                            finalWriter.write(line);
                            finalWriter.newLine();
                        }
                    } catch (IOException e) {
                        System.err.println("Error reading log file: " + logFile.getName());
                        e.printStackTrace();
                    }
                    logFile.delete();
                }
            } catch (IOException e) {
                System.err.println("Error writing to final CSV file.");
                e.printStackTrace();
            }
        }
    }

    public static List<String> runExternalCommand(
            List<String> command,
            File workingDir,
            long timeoutMs,
            boolean verbose
    ) throws IOException, InterruptedException {

        List<String> outputLines = new ArrayList<>();
        ProcessBuilder pb = new ProcessBuilder(command);
        if (workingDir != null && workingDir.isDirectory()) {
            pb.directory(workingDir);
            if (verbose) {
                System.out.println("Setting working directory to " + workingDir.getAbsolutePath());
            }
        }

        Process process = null;
        try {
            process = pb.start();
            process.getOutputStream().close();

            try (BufferedReader stdOut = new BufferedReader(new InputStreamReader(process.getInputStream()));
                 BufferedReader stdErr = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {

                long startTime = System.currentTimeMillis();
                boolean finished = false;

                while (true) {
                    while (stdOut.ready()) {
                        String line = stdOut.readLine();
                        if (line == null) break;
                        if (verbose) {
                            System.out.println("CMD-OUT: " + line);
                        }
                        outputLines.add(line);
                    }
                    while (stdErr.ready()) {
                        String line = stdErr.readLine();
                        if (line == null) break;
                        if (verbose) {
                            System.err.println("CMD-ERR: " + line);
                        }
                    }
                    try {
                        int exitVal = process.exitValue();
                        finished = true;
                        if (exitVal != 0) {
                            throw new IOException("Process exited with code: " + exitVal);
                        }
                        break;
                    } catch (IllegalThreadStateException e) {
                        // still running
                    }
                    if (System.currentTimeMillis() - startTime > timeoutMs) {
                        if (verbose) {
                            System.err.println("Process timed out, destroying...");
                        }
                        process.destroyForcibly();
                        throw new IOException("Process timed out after " + timeoutMs + " ms: " + command);
                    }
                    Thread.sleep(100);
                }

                // Drain any remaining
                while (stdOut.ready()) {
                    String line = stdOut.readLine();
                    if (line != null) {
                        outputLines.add(line);
                        if (verbose) {
                            System.out.println("CMD-OUT: " + line);
                        }
                    }
                }
                while (stdErr.ready()) {
                    String line = stdErr.readLine();
                    if (verbose) {
                        System.err.println("CMD-ERR: " + line);
                    }
                }
                if (!finished) {
                    throw new IOException("Process ended unexpectedly without a proper finish.");
                }
            }
        } finally {
            if (process != null && process.isAlive()) {
                process.destroyForcibly();
            }
        }

        return outputLines;
    }
}

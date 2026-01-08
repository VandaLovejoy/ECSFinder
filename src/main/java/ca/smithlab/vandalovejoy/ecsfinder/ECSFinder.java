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

    // Reference species: raw argument + compiled regex pattern
    public static String  REF_SPECIES_RAW     = "homo";  // default: Homo sapiens (substring)
    public static Pattern REF_SPECIES_PATTERN = Pattern.compile("homo", Pattern.CASE_INSENSITIVE);

    public static String SINGLE_CSV_PATH = null;

    static int GAPS = 50, NTHREDS = 4, MIN_MPI = 50;
    static boolean VERBOSE = false, MAFFT = false;
    static String FILENAME = "", OUT_PATH = "", dirProgram = "";
    static String SSZBINARY = "", ALIFOLDBINARY = "",
            RNAALIFOLD = "", R = "", RSCAPE = "", MAFFTBINARY = "";
    static double SSZR = -3.0;

    // Map block_part -> realigned reference-species sequence (MAFFT pipeline)
    private static Map<String, String> refSpeciesSequences = new HashMap<>();

    // Map block_part -> FASTA chunk coordinates (MAFFT pipeline)
    private static Map<String, CoordinateInfo> coordinateMap = new HashMap<>();

    static final int MIN_SEQS_PER_BLOCK = 5;

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
                    "Usage: java ECSFinder [options] -o output/directory -i input.maf_or_dir\n\n" +
                    "Options:\n" +
                    " -c int     number of CPUs for calculations (default 4)\n" +
                    " -g int     max gap percentage of sequences for 2D prediction (default 50)\n" +
                    " -sszr dbl  report SISSIz+RIBOSUM hits below this Z-score (default -3.0)\n" +
                    " -mpi int   minimum MPI (percent identity) cutoff (default 50)\n" +
                    " -mafft     realign aln using MAFFT (default FALSE)\n" +
                    " -ref str   reference species (substring or /regex/) (default \"homo\")\n" +
                    " -v         verbose (messy but detailed) output\n"+
            "Requirements / input expectations:\n" +
                    " - Alignments should contain at least " + MIN_SEQS_PER_BLOCK + " sequences to provide enough\n" +
                    "   evolutionary signal for structural conservation/covariation. Blocks with fewer sequences are skipped" +
                    " or alignments with fewer than 50 nucleotides per species.\n" );
            System.exit(0);
        }

        // parse arguments
        parseArguments(args);

        // Now that OUT_PATH is known, define the single CSV path
        SINGLE_CSV_PATH = OUT_PATH + "/final.csv";

        // Write header exactly once at start
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

        keepTPAlignments(predictionsCsv, OUT_PATH + "/aln");

        System.out.println("Predictions written to: " + predictionsCsv);
    }

    /* ==============================
       Reference species helpers
       ============================== */

    /**
     * Configure the reference species from the user argument.
     * If arg is of the form /.../, treat the inside as a regex.
     * Otherwise, treat arg as a literal substring (case-insensitive).
     */
    public static void setRefSpecies(String arg) {
        REF_SPECIES_RAW = arg;

        String patternBody;
        if (arg.startsWith("/") && arg.endsWith("/") && arg.length() > 2) {
            // Full regex, strip surrounding slashes
            patternBody = arg.substring(1, arg.length() - 1);
        } else {
            // Literal substring: quote it and use .find()
            patternBody = Pattern.quote(arg);
        }

        REF_SPECIES_PATTERN = Pattern.compile(patternBody, Pattern.CASE_INSENSITIVE);

        if (VERBOSE) {
            System.out.println("Reference species set to: " + REF_SPECIES_RAW +
                    " (regex: " + REF_SPECIES_PATTERN.pattern() + ")");
        }
    }

    /* ==============================
       CSV header
       ============================== */

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

    /* ==============================
       Binaries
       ============================== */

    private static void setBinaryPaths() throws IOException, InterruptedException {
        ALIFOLDBINARY = which("RNALalifold");
        SSZBINARY     = which("SISSIz");
        RNAALIFOLD    = which("RNAalifold");
        RSCRIPT       = which("Rscript");
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
            List<String> lines = runExternalCommand(command, null, 10_000, VERBOSE);
            if (lines.isEmpty()) {
                throw new IOException("Cannot find " + binaryName + " in PATH");
            }
            return lines.get(0).trim();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted while locating " + binaryName, e);
        }
    }

    /* ==============================
       Argument parsing
       ============================== */

    private static void parseArguments(String[] args) {
        // Default reference pattern
        setRefSpecies("homo");

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-c":
                    NTHREDS = Integer.parseInt(args[++i]);
                    break;
                case "-g":
                    GAPS = Integer.parseInt(args[++i]);
                    break;
                case "-o": {
                    String outArg = args[++i];
                    OUT_PATH = normalizePath(outArg);
                    createDirectory(OUT_PATH);
                    dirProgram = System.getProperty("user.dir");
                    break;
                }
                case "-v":
                    VERBOSE = true;
                    break;
                case "-mafft":
                    MAFFT = true;
                    break;
                case "-sszr":
                    SSZR = Double.parseDouble(args[++i]);
                    break;
                case "-mpi":
                    MIN_MPI = Integer.parseInt(args[++i]);
                    break;
                case "-i":
                    FILENAME = normalizePath(args[++i]);
                    break;
                case "-ref":
                    setRefSpecies(args[++i]);
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
        System.out.println("Usage: java ECSFinder [options] -o output/directory -i input.maf_or_dir");
        System.exit(1);
    }

    private static void createDirectory(String path) {
        File dir = new File(path);
        if (!dir.isDirectory()) {
            dir.mkdirs();
        }
    }

    /**
     * Normalize a user-provided path:
     *  - expand "~" and "~/..." to the user's home directory
     *  - resolve relative paths against the current working directory
     *  - return an absolute path string
     */
    private static String normalizePath(String arg) {
        if (arg == null || arg.isEmpty()) {
            return arg;
        }

        String expanded = arg;

        // Handle "~" and "~/..."
        if (arg.equals("~")) {
            expanded = System.getProperty("user.home");
        } else if (arg.startsWith("~" + File.separator) || arg.startsWith("~/")) {
            expanded = System.getProperty("user.home") + arg.substring(1);
        }

        File f = new File(expanded);
        if (!f.isAbsolute()) {
            f = new File(System.getProperty("user.dir"), expanded);
        }
        return f.getAbsolutePath();
    }

    /* ==============================
       MAF preprocessing
       ============================== */

    private static void preprocessMafFiles() {
        // FILENAME is already normalized to an absolute path (and ~ expanded)
        File input = new File(FILENAME);

        if (input.isDirectory()) {
            // Case 1: directory containing .maf files
            File[] mafFiles = input.listFiles((dir, name) -> name.endsWith(".maf"));
            if (mafFiles == null || mafFiles.length == 0) {
                System.err.println("Error: No .maf files found in the directory: " + input.getAbsolutePath());
                System.exit(1);
            }

            String[] mafFilePaths = Arrays.stream(mafFiles)
                    .map(File::getAbsolutePath)
                    .toArray(String[]::new);

            try {
                MergeNFilter mergeNFilter = new MergeNFilter();
                mergeNFilter.process(mafFilePaths, OUT_PATH);
            } catch (IOException e) {
                System.err.println("Error processing MAF files: " + e.getMessage());
                e.printStackTrace();
                System.exit(1);
            }

        } else if (input.isFile()) {
            // Case 2: single .maf file
            if (!input.getName().endsWith(".maf")) {
                System.err.println("Error: Input file does not have .maf extension: " + input.getAbsolutePath());
                System.exit(1);
            }

            String[] mafFilePaths = new String[]{ input.getAbsolutePath() };

            try {
                MergeNFilter mergeNFilter = new MergeNFilter();
                mergeNFilter.process(mafFilePaths, OUT_PATH);
            } catch (IOException e) {
                System.err.println("Error processing MAF file: " + e.getMessage());
                e.printStackTrace();
                System.exit(1);
            }

        } else {
            System.err.println("Error: Input path does not exist: " + input.getAbsolutePath());
            System.exit(1);
        }
    }

    /* ==============================
       RNALalifold + block processing
       ============================== */

    private static void runRNALalifoldAndProcessResults() throws IOException, InterruptedException {
        String inputFile = OUT_PATH + "/output.maf";
        String stockholmFolderPath = OUT_PATH + "/stockholm";
        File stockholmFolder = createStockholmFolder(stockholmFolderPath);

        // MAFFT pipeline uses RNALalifold on realigned FASTAs
        if (MAFFT) {
            processWithMafft(inputFile);
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
        convertMafToSeparateFastas(inputFile, OUT_PATH);

        // Realign each block individually to avoid memory issues
        File outputDir = new File(OUT_PATH + "/outputFastaDir");
        File[] fastaFiles = outputDir.listFiles((dir, name) -> name.endsWith(".fasta"));

        if (fastaFiles != null) {
            for (File fastaFile : fastaFiles) {
                if (VERBOSE) {
                    System.out.println("Realigning file: " + fastaFile.getName());
                }
                File realignedOutput = new File(
                        fastaFile.getAbsolutePath().replace(".fasta", "_realigned.fasta")
                );
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

        ExecutorService multiThreads = Executors.newFixedThreadPool(NTHREDS);
        List<Future<?>> futures = new ArrayList<>();

        String path_aln = OUT_PATH + "/aln/";
        createDirectory(path_aln);
        createDirectory(OUT_PATH + "/stockholm");

        String line;
        while ((line = readFile.readLine()) != null) {
            if (line.length() >= 1 && line.charAt(0) != '#') {
                if (line.charAt(0) == 'a') {
                    blockAln++;
                } else if (line.charAt(0) == 's') {
                    temp.append(line).append("@");
                }
            } else if (line.equals("") && !temp.toString().isEmpty()) {
                long nSeq = Arrays.stream(temp.toString().split("@"))
                        .map(String::trim)
                        .filter(s -> !s.isEmpty())
                        .count();

                if (nSeq < MIN_SEQS_PER_BLOCK) {
                    if (VERBOSE) {
                        System.out.println("[INFO] Skipping block " + blockAln +
                                " (only " + nSeq + " sequences; need >= " + MIN_SEQS_PER_BLOCK + ")");
                    }
                    temp = new StringBuilder();
                    continue;
                }


                // For non-MAFFT: run RNALalifold per block on a FASTA version of this block
                if (!MAFFT) {
                    try {
                        runRNALalifoldOnBlock(temp.toString(), blockAln);
                    } catch (IOException | InterruptedException e) {
                        System.err.println("RNALalifold failed on block " + blockAln + ": " + e.getMessage());
                        e.printStackTrace();
                    }
                }

                iterateStockholm(temp.toString(), blockAln, futures, multiThreads);
                temp = new StringBuilder();
            } else {
                temp = new StringBuilder();
            }
        }
        readFile.close();

// handle last block if file doesn't end with blank line
        if (!temp.toString().isEmpty()) {
            long nSeq = Arrays.stream(temp.toString().split("@"))
                    .map(String::trim)
                    .filter(s -> !s.isEmpty())
                    .count();

            if (nSeq >= MIN_SEQS_PER_BLOCK) {
                if (!MAFFT) {
                    try {
                        runRNALalifoldOnBlock(temp.toString(), blockAln);
                    } catch (IOException | InterruptedException e) {
                        System.err.println("RNALalifold failed on last block " + blockAln + ": " + e.getMessage());
                        e.printStackTrace();
                    }
                }
                iterateStockholm(temp.toString(), blockAln, futures, multiThreads);
            } else if (VERBOSE) {
                System.out.println("[INFO] Skipping last block " + blockAln +
                        " (only " + nSeq + " sequences; need >= " + MIN_SEQS_PER_BLOCK + ")");
            }
        }


        for (Future<?> f : futures) {
            try {
                f.get();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }

        multiThreads.shutdown();
        if (!multiThreads.awaitTermination(10, TimeUnit.MINUTES)) {
            System.err.println("Executor did not terminate in time, forcing shutdown.");
            multiThreads.shutdownNow();
        }
    }

    /**
     * Run RNALalifold on a single MAF block (non-MAFFT pipeline).
     * Convert the block (s-lines) to a temporary FASTA and call RNALalifold
     * with a block-specific id-prefix: alifold_<blockAln>.
     */
    /**
     * Run RNALalifold on a single MAF block (non-MAFFT pipeline).
     * Converts the block (s-lines) to a temporary FASTA and calls RNALalifold.
     * Skips blocks with too few sequences (MIN_SEQS_PER_BLOCK).
     */
    private static void runRNALalifoldOnBlock(String block, int blockAln)
            throws IOException, InterruptedException {

        String stockholmDir = OUT_PATH + "/stockholm";
        createDirectory(stockholmDir);

        File fastaFile = new File(stockholmDir, "block_" + blockAln + ".fasta");

        // Split once
        String[] sLines = block.split("@");

        // Count sequences (MAF 's' lines)
        int nS = 0;
        for (String sLine : sLines) {
            if (sLine == null) continue;
            sLine = sLine.trim();
            if (sLine.startsWith("s")) { // handles "s ..." or "s\t..."
                nS++;
            }
        }

        // Enforce minimum sequences
        if (nS < MIN_SEQS_PER_BLOCK) {
            if (VERBOSE) {
                System.out.println("[INFO] Not running RNALalifold on block " + blockAln +
                        " (only " + nS + " sequences; need >= " + MIN_SEQS_PER_BLOCK + ")");
            }
            return;
        }

        // Write the MAF block to FASTA
        try (BufferedWriter w = new BufferedWriter(new FileWriter(fastaFile))) {
            for (String sLine : sLines) {
                if (sLine == null) continue;
                sLine = sLine.trim();
                if (sLine.isEmpty()) continue;
                if (!sLine.startsWith("s")) continue;

                String[] tok = sLine.split("\\s+");
                if (tok.length < 7) continue;

                String id  = tok[1];
                String seq = tok[tok.length - 1];

                w.write(">" + id);
                w.newLine();
                w.write(seq);
                w.newLine();
            }
        }

        List<String> cmd = Arrays.asList(
                ALIFOLDBINARY,
                "--id-prefix=alifold_" + blockAln,
                "--noLP",
                "--maxBPspan=500",
                "--ribosum_scoring",
                "--aln-stk",
                fastaFile.getAbsolutePath()
        );

        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.directory(new File(stockholmDir));
        Process process = pb.start();

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
             BufferedReader err    = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {

            // Drain stdout
            while (reader.readLine() != null) { /* ignore */ }

            // Drain stderr (optionally print)
            String errLine;
            while ((errLine = err.readLine()) != null) {
                if (VERBOSE) {
                    System.err.println("RNALalifold(block " + blockAln + ") ERR: " + errLine);
                }
            }
        }

        int exitCode = process.waitFor();

        // Always try to delete temp FASTA
        if (!fastaFile.delete() && VERBOSE) {
            System.err.println("Warning: could not delete temporary FASTA " + fastaFile.getAbsolutePath());
        }

        if (exitCode != 0) {
            throw new IOException("RNALalifold exited with " + exitCode + " on block " + blockAln);
        }
    }


    private static void iterateStockholm(String block, int blockAln,
                                         List<Future<?>> futures,
                                         ExecutorService multiThreads) throws IOException {
        ArrayList<String[]> associativeList = new ArrayList<>();
        String[] mafTabTemp = block.split("@")[0].split("\\s+");

        File dir = new File(OUT_PATH + "/stockholm");

        // For both MAFFT and non-MAFFT pipelines, .stk files are now named alifold_<blockAln>_NNNN.stk
        String prefix = "alifold_" + blockAln + "_";

        File[] stkFiles = dir.listFiles((dir1, name) ->
                name.startsWith(prefix) && name.endsWith(".stk")
        );

        if (stkFiles == null || stkFiles.length == 0) {
            System.out.println("No .stk files found for alignment block " + blockAln
                    + ". Skipping this block.\n");
            return;
        }

        Arrays.sort(stkFiles, (f1, f2) -> {
            int sub1 = parseSubBlockNumber(f1.getName(), blockAln);
            int sub2 = parseSubBlockNumber(f2.getName(), blockAln);
            return Integer.compare(sub1, sub2);
        });

        for (File stkFile : stkFiles) {
            if (VERBOSE) {
                System.out.println("Processing .stk for block #" + blockAln
                        + ": " + stkFile.getName());
            }
            processStockholmFile(stkFile, mafTabTemp, associativeList, futures, multiThreads);
        }
    }

    private static int parseSubBlockNumber(String stkName, int blockAln) {
        String patternStr = "^alifold_" + blockAln + "_(\\d+)\\.stk$";
        Pattern p = Pattern.compile(patternStr);
        Matcher m = p.matcher(stkName);
        if (m.matches()) {
            return Integer.parseInt(m.group(1));
        }
        return Integer.MAX_VALUE;
    }

    private static void processStockholmFile(File stockholmFile, String[] mafTabTemp,
                                             ArrayList<String[]> associativeList,
                                             List<Future<?>> futures,
                                             ExecutorService multiThreads) {
        if(VERBOSE) {
            System.out.println("[DEBUG] Processing Stockholm: " + stockholmFile.getName());
        }
        String result = "";

        if (MAFFT) {
            String regex = "alifold_(\\d+)_(\\d+)";
            Pattern pattern = Pattern.compile(regex);
            Matcher matcher = pattern.matcher(stockholmFile.getName());

            if (matcher.find()) {
                String firstPart = matcher.group(1);
                String secondPart = matcher.group(2).replaceFirst("^0+", "");
                result = firstPart + "_" + secondPart;
            } else {
                System.out.println("Error: Filename does not match expected pattern: " + stockholmFile.getName());
                return;
            }
        } else {
            // Non-MAFFT: coordinateMap/refSpeciesSequences are not used,
            // result doesn't need to be a valid key.
            result = "0";
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

    /**
     * Parse motif start index from GF ID tokens robustly.
     * We collect all integers from arrayName and assume the last two
     * are (start, end). Returns start, or null if nothing parseable.
     */
    private static Integer parseMotifStartFromGFID(String[] arrayName) {
        List<Integer> nums = new ArrayList<>();
        for (String token : arrayName) {
            try {
                nums.add(Integer.parseInt(token));
            } catch (NumberFormatException ignored) {
            }
        }
        if (nums.isEmpty()) {
            return null;
        }
        if (nums.size() == 1) {
            return nums.get(0);
        }
        // typical: alifold_<block>_aln_<start>_<end>
        return nums.get(nums.size() - 2);
    }

    private static void processMotif(String[] mafTabTemp, String[] arrayName,
                                     ArrayList<String[]> associativeList,
                                     String gcReference, String gcSScons,
                                     List<Future<?>> futures,
                                     ExecutorService multiThreads,
                                     String result) {

        // 1) robustly parse motif start index in alignment coords
        Integer startIndexInAlignment = parseMotifStartFromGFID(arrayName);
        if (startIndexInAlignment == null) {
            if (VERBOSE) {
                System.err.println("[WARN] Could not parse start index from GF ID: "
                        + String.join(" ", arrayName) + " — skipping motif.");
            }
            return;
        }

        int[] cordMotif = getRealCoordinates(
                startIndexInAlignment,
                mafTabTemp,
                associativeList.get(0)[1],
                result
        );

        if (cordMotif[0] < 0 || cordMotif[1] < 0) {
            if (VERBOSE) {
                System.err.println("[WARN] Invalid coordinates for motif from GF ID: "
                        + String.join(" ", arrayName) + " — skipping.");
            }
            return;
        }

        String loci = Arrays.toString(cordMotif);
        String chrom = mafContig(mafTabTemp[1]);;
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
        if (Integer.parseInt(arrayLociChrm[2]) - Integer.parseInt(arrayLociChrm[1]) < 50) {
            return;
        }
        if(VERBOSE) {
            System.out.println("[DEBUG] Creating ScanItFast job for block = " + Arrays.toString(arrayName)
                    + " on file result key = " + result);
        }

        ScanItFast aln = new ScanItFast(associativeList, arrayLociChrm,
                new File(OUT_PATH), SSZBINARY, VERBOSE);
        aln.setSszR(SSZR);
        aln.setGap(GAPS);
        aln.setMinMPI(MIN_MPI);
        Future<?> future = multiThreads.submit(aln);
        futures.add(future);
    }

    /* ==============================
       Utility helpers
       ============================== */

    private static String extractValue(String line) {
        String[] lineReference = line.split(" ");
        return lineReference[lineReference.length - 1];
    }

    private static String[] processSpeciesLine(String line) {
        String[] species = line.split(" ", 2);
        species[1] = species[1].trim();
        if (line == null || line.trim().isEmpty()) {
            return null;
        }
        return species;
    }

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

    /* ==============================
       Rscript call + TP-only aln cleanup
       ============================== */

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
            r.readLine();
            String line;
            while ((line = r.readLine()) != null) preds.add(Double.parseDouble(line.split(",")[1]));
        }
        return preds;
    }

    /**
     * Delete any .aln file whose predicted probability is below threshold.
     * predictionsCsv is expected to have header, then "name_file,probability".
     */
    private static void keepTPAlignments(String predictionsCsv, String alnDirPath)
            throws IOException {

        File predFile = new File(predictionsCsv);
        if (!predFile.exists()) {
            if (VERBOSE) {
                System.err.println("[WARN] predictions file not found: " + predictionsCsv +
                        " — skipping TP-based cleanup.");
            }
            return;
        }

        // 1) Parse header and find name_file + Predicted_Class columns
        String[] headerCols;
        try (BufferedReader br = new BufferedReader(new FileReader(predFile))) {
            String header = br.readLine();
            if (header == null) return;
            headerCols = header.split(",");
        }

        int nameIdx = -1;
        int labelIdx = -1;
        for (int i = 0; i < headerCols.length; i++) {
            String col = headerCols[i].trim();
            if (col.equals("name_file") || col.equals("name")) {
                nameIdx = i;
            }
            // Important: use Predicted_Class from your R script
            if (col.equals("Predicted_Class")) {
                labelIdx = i;
            }
        }

        if (nameIdx == -1 || labelIdx == -1) {
            // No explicit TP label → nothing to do, preserve original behavior
            if (VERBOSE) {
                System.out.println("[INFO] Could not find name_file and Predicted_Class "
                        + "in " + predictionsCsv + ". Not deleting any .aln files.");
            }
            return;
        }

        // 2) Collect base names (e.g. "5_...rc.aln") that are TP
        Set<String> tpBases = new HashSet<>();
        try (BufferedReader br = new BufferedReader(new FileReader(predFile))) {
            br.readLine(); // skip header
            String line;
            while ((line = br.readLine()) != null) {
                String[] parts = line.split(",");
                if (parts.length <= Math.max(nameIdx, labelIdx)) continue;

                String baseName = parts[nameIdx].trim();   // e.g. "5_RF00017_...rc.aln"
                String labelVal = parts[labelIdx].trim();  // "TP" or "FP"

                if (labelVal.equals("TP")) {
                    tpBases.add(baseName);
                }
            }
        }

        File alnDir = new File(alnDirPath);
        if (!alnDir.isDirectory()) return;

        // 3) Visit all files containing ".aln" (final + temp, rc, etc.)
        File[] alnFiles = alnDir.listFiles((d, name) -> name.contains(".aln"));
        if (alnFiles == null) return;

        for (File f : alnFiles) {
            String fname = f.getName();
            int idx = fname.indexOf(".aln");
            if (idx == -1) continue;

            // Canonical base: everything up to and including ".aln"
            String base = fname.substring(0, idx + 4);

            if (!tpBases.contains(base)) {
                if (VERBOSE) {
                    System.out.println("Deleting non-TP alignment: " + fname +
                            " (base=" + base + ")");
                }
                if (!f.delete() && VERBOSE) {
                    System.err.println("Warning: could not delete " + f.getAbsolutePath());
                }
            }
        }
    }

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

        // jarDir is the directory containing the JAR, or the directory itself if not a file
        File jarDir = jarFile.isFile() ? jarFile.getParentFile() : jarFile;

        // Candidate locations for the RF directory
        List<File> candidates = new ArrayList<>();
        // 1) RF next to the JAR: target/RF
        candidates.add(new File(jarDir, "RF"));
        // 2) RF one level above the JAR: project-root/RF (your current layout)
        if (jarDir.getParentFile() != null) {
            candidates.add(new File(jarDir.getParentFile(), "RF"));
        }
        // 3) RF in the current working directory
        candidates.add(new File("RF"));

        for (File c : candidates) {
            if (c.isDirectory()) {
                if (VERBOSE) {
                    System.out.println("Using RF directory: " + c.getAbsolutePath());
                }
                return c.getAbsolutePath();
            }
        }

        StringBuilder sb = new StringBuilder("Could not find RF directory with predictions_ECSFinder.R and final_rf_model.rds. Checked:\n");
        for (File c : candidates) {
            sb.append("  - ").append(c.getAbsolutePath()).append("\n");
        }
        throw new IOException(sb.toString());
    }


    /* ==============================
       MAF -> FASTA splitting (MAFFT pipeline)
       ============================== */

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
                    String sequence   = tokens[tokens.length - 1];

                    // Use the reference species line to get start/len/orientation/srcSize
                    if (REF_SPECIES_PATTERN.matcher(sequenceId).find()) {
                        startAln    = Integer.parseInt(tokens[2]);
                        lengthAln   = Integer.parseInt(tokens[3]);
                        orientation = tokens[4];
                        endAln      = Integer.parseInt(tokens[5]);
                    }

                    speciesSequences
                            .computeIfAbsent(sequenceId, k -> new StringBuilder())
                            .append(sequence);
                    currentBlockLength = Math.max(currentBlockLength,
                            speciesSequences.get(sequenceId).length());
                }
            }
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

        String refSequence = null;
        String refSpeciesId = null;
        String chromosomeNumber = null;

        for (String key : speciesSequences.keySet()) {
            if (REF_SPECIES_PATTERN.matcher(key).find()) {
                refSequence  = speciesSequences.get(key).toString();
                refSpeciesId = key;
                chromosomeNumber = mafContig(refSpeciesId);
                break;
            }
        }
        if (refSequence == null) {
            throw new IllegalArgumentException(
                    "Reference species sequence not found in speciesSequences (pattern: "
                            + REF_SPECIES_RAW + ").");
        }

        for (int i = 0; i < fullLength; i += (maxBlockSize - overlapLength)) {
            int blockStartIndex = i;
            int blockEndIndex = Math.min(i + maxBlockSize, fullLength);

            String refSubSequence = refSequence.substring(blockStartIndex, blockEndIndex);
            int nucleotidesUpToStart = refSequence
                    .substring(0, blockStartIndex)
                    .replaceAll("-", "")
                    .length();
            int blockNucleotideLength = refSubSequence
                    .replaceAll("-", "")
                    .length();

            int genomicStart, genomicEnd;
            if (orientation.equals("+")) {
                genomicStart = startAln + nucleotidesUpToStart + 1;
                genomicEnd   = genomicStart + blockNucleotideLength - 1;
            } else if (orientation.equals("-")) {
                genomicEnd   = chromosomeLength - startAln - nucleotidesUpToStart;
                genomicStart = genomicEnd - blockNucleotideLength + 1;
            } else {
                throw new IllegalArgumentException("Invalid orientation: " + orientation);
            }

            int adjustedGenomicStart = genomicStart;
            int adjustedGenomicEnd = genomicEnd;
            if (adjustedGenomicStart > adjustedGenomicEnd) {
                int temp = adjustedGenomicStart;
                adjustedGenomicStart = adjustedGenomicEnd;
                adjustedGenomicEnd   = temp;
            }

            String blockFileName = outputDirPath + "/outputFastaDir/block_"
                    + blockCount + "_part_" + blockPart + "_chr"
                    + chromosomeNumber + "_" + adjustedGenomicStart + "_"
                    + adjustedGenomicEnd + "_" + orientation + ".fasta";

            try (BufferedWriter writer = new BufferedWriter(new FileWriter(blockFileName))) {
                for (Map.Entry<String, StringBuilder> entry : speciesSequences.entrySet()) {
                    String speciesId = entry.getKey();
                    String sequence  = entry.getValue().toString();
                    String subSequence = sequence
                            .substring(blockStartIndex, blockEndIndex)
                            .replace(".", "");
                    writer.write(">" + speciesId);
                    writer.newLine();
                    writer.write(subSequence);
                    writer.newLine();
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

    /* ==============================
       MAFFT realignment
       ============================== */

    public static void realignSequences(File inputFilePath, File outputFilePath)
            throws IOException, InterruptedException {
        String realignedFilePath = inputFilePath.getAbsolutePath()
                .replace(".fasta", "_realigned.fasta");

        String regex = "block_(\\d+)_part_(\\d+)_chr(\\w+)_(\\d+)_(\\d+)_([+-])(?:_\\w+)?\\.fasta";
        Pattern pattern = Pattern.compile(regex);
        Matcher matcher = pattern.matcher(inputFilePath.getAbsolutePath());

        String result = "";
        if (matcher.find()) {
            String blockNumber = matcher.group(1);
            String partNumber  = matcher.group(2);
            String chromName   = matcher.group(3);
            int    fastaStart  = Integer.parseInt(matcher.group(4));
            int    fastaEnd    = Integer.parseInt(matcher.group(5));
            String strand      = matcher.group(6);

            result = blockNumber + "_" + partNumber;

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

                        if (REF_SPECIES_PATTERN.matcher(speciesReal).find()) {
                            refSpeciesSequences.put(result, sequence.toString());
                        }

                        sequence = new StringBuilder();
                    }
                    speciesReal = line.substring(1);
                } else {
                    sequence.append(line);
                }
            }
            if (!speciesReal.isEmpty()) {
                writer.write(">" + speciesReal + "\n");
                writer.write(sequence.toString().toUpperCase() + "\n");

                if (REF_SPECIES_PATTERN.matcher(speciesReal).find()) {
                    refSpeciesSequences.put(result, sequence.toString());
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
                if (VERBOSE) {
                    System.err.println("MAFFT-ERR: " + errorLine);
                }
            }
        }

        File tempInputFile = new File(inputFilePath.getParent(), "temp_gap_stripped.fasta");
        if (tempInputFile.exists()) {
            tempInputFile.delete();
        }
    }

    /* ==============================
       RNALalifold call (MAFFT pipeline)
       ============================== */

    private static void runRNALalifold(String inputFilePath)
            throws IOException, InterruptedException {
        if(VERBOSE) {
            System.out.println("[DEBUG] About to run RNALalifold on FASTA: " + inputFilePath);
        }
        String outputFilePath = OUT_PATH + "/stockholm";
        List<String> command;

        if (inputFilePath.endsWith(".fasta")) {
            String fileName = new File(inputFilePath).getName();
            String regex = "block_(\\d+)_part_(\\d+)_chr(\\w+)_(\\d+)_(\\d+)_([+-])(?:_realigned)?\\.fasta$";
            Pattern pattern = Pattern.compile(regex);
            Matcher matcher = pattern.matcher(fileName);

            int blockNumber = 0;
            int partNumber  = 0;
            if (matcher.find()) {
                blockNumber = Integer.parseInt(matcher.group(1));
                partNumber  = Integer.parseInt(matcher.group(2));
            } else {
                System.err.println("Block/Part not found in the input FASTA name: "
                        + fileName + ". Skipping file.");
                return;
            }
            command = Arrays.asList(
                    ALIFOLDBINARY,
                    "--id-prefix=alifold_" + blockNumber,
                    "--id-start=" + partNumber,
                    "--noLP",
                    "--maxBPspan=500",
                    "--ribosum_scoring",
                    "--aln-stk",
                    inputFilePath
            );
        } else {
            // Not used in new non-MAFFT pipeline
            command = Arrays.asList(
                    ALIFOLDBINARY,
                    "--id-prefix=alifold",
                    "--noLP",
                    "--maxBPspan=500",
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

    /* ==============================
       Coordinate calculation
       ============================== */

    private static int[] getRealCoordinates(int start, String[] mafCord,
                                            String motifHuman,
                                            String blockPartKey) {

        CoordinateInfo info = coordinateMap.get(blockPartKey);
        if (info == null) {
            return oldMafBasedCoordinates(start, mafCord, motifHuman);
        }

        String refAligned = refSpeciesSequences.get(blockPartKey);
        if (refAligned == null) {
            return oldMafBasedCoordinates(start, mafCord, motifHuman);
        }
        String upToStart = refAligned.substring(0, start - 1);
        int offsetInChunk = upToStart.replaceAll("-", "").length();
        int motifLen = motifHuman.replaceAll("-", "").length();

        int finalStart, finalEnd;
        if (info.strand.equals("+")) {
            finalStart = info.start + offsetInChunk;
            finalEnd   = finalStart + motifLen - 1;
        } else {
            finalEnd   = info.end - offsetInChunk;
            finalStart = finalEnd - motifLen + 1;
        }

        if (finalStart > finalEnd) {
            int temp = finalStart;
            finalStart = finalEnd;
            finalEnd   = temp;
        }

        return new int[]{ finalStart, finalEnd };
    }

    private static int[] oldMafBasedCoordinates(int startIndexInAlignment,
                                                String[] mafCord,
                                                String motifHuman) {

        int mafStart  = Integer.parseInt(mafCord[2]);
        int mafSize   = Integer.parseInt(mafCord[3]);
        int mafSrcLen = Integer.parseInt(mafCord[5]);
        String strand = mafCord[4];

        String alnSeq = mafCord[6];
        int alnLen = alnSeq.length();

        if (startIndexInAlignment < 0 || startIndexInAlignment > alnLen) {
            if (VERBOSE) {
                System.err.println(
                        "[WARN] oldMafBasedCoordinates: startIndexInAlignment=" +
                                startIndexInAlignment + " outside alignment length=" + alnLen +
                                " for " + mafCord[1] + "; skipping this motif.");
            }
            return new int[]{-1, -1};
        }

        int offset = countRealBases(alnSeq.substring(0, startIndexInAlignment));

        int forwardStart, forwardEnd;
        if (strand.equals("+")) {
            forwardStart = mafStart;
            forwardEnd   = mafStart + mafSize - 1;
        } else {
            forwardStart = mafSrcLen - (mafStart + mafSize);
            forwardEnd   = mafSrcLen - mafStart - 1;
        }

        int motifLen = motifHuman.replaceAll("-", "").length();
        int finalStart = forwardStart + offset;
        int finalEnd   = finalStart + motifLen - 1;

        return new int[]{ finalStart, finalEnd };
    }

    private static int countRealBases(String seq) {
        return seq.replaceAll("-", "").length();
    }

    /* ==============================
       CoordinateInfo
       ============================== */

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

    /* ==============================
       Log merging & external commands
       ============================== */

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
    /** Return "species" part from a MAF src like species.contig. If no dot, return full src. */
    private static String mafSpecies(String src) {
        if (src == null) return "";
        int dot = src.indexOf('.');
        return (dot > 0) ? src.substring(0, dot) : src;
    }

    /** Return "contig/chrom" part from a MAF src like species.contig. If no dot, return full src. */
    private static String mafContig(String src) {
        if (src == null) return "";
        int dot = src.indexOf('.');
        return (dot >= 0 && dot < src.length() - 1) ? src.substring(dot + 1) : src;
    }

}

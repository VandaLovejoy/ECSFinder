package ca.smithlab.vandalovejoy.ecsfinder;

import java.io.*;
import java.util.*;

/**
 * This class filters multiple sequence alignments and removes duplicate species within alignment blocks.
 */
public class MergeNFilter {

    public void process(String[] args, String outputPath) throws IOException {
        BufferedWriter segDups = null;
        BufferedWriter out = new BufferedWriter(new FileWriter(outputPath + "/output.maf"));

        out.write("##maf version=2 \n" +
                "# original dump date: 2020-12-10 13:44:00\n# ensembl release: 103\n" +
                "# emf comment: Alignments: 46 eutherian mammals EPO\n# emf comment: Region:" +
                " Homo sapiens chromosome:GRCh38\n\n\n");

        boolean hasRemovedLines = false;

        for (String file : args) {
            BufferedReader entry = new BufferedReader(new FileReader(file));
            String line;
            boolean firstAlignment = true;
            while ((line = entry.readLine()) != null) {
                if (line.length() != 0 && line.charAt(0) == 's') {
                    LinkedHashMap<String, String[]> speciesSequences = new LinkedHashMap<>();
                    ArrayList<String[]> duplicateSequences = new ArrayList<>();
                    while (line != null && line.length() != 0 && line.charAt(0) == 's') {
                        String[] arraySequenceInfo = line.split("\\s+");
                        String nameSpeciesWithChro = arraySequenceInfo[1];
                        String nameSpeciesOnly = nameSpeciesWithChro.substring(0, nameSpeciesWithChro.indexOf("."));
                        if (!nameSpeciesOnly.equals("ancestral_sequences")) {
                            if (!speciesSequences.containsKey(nameSpeciesOnly)) {
                                speciesSequences.put(nameSpeciesOnly, arraySequenceInfo);
                            } else {
                                duplicateSequences.add(arraySequenceInfo);
                            }
                        }
                        line = entry.readLine();
                    }
                    if (speciesSequences.size() == 1) {
                        duplicateSequences.add(speciesSequences.get("homo_sapiens"));
                    }
                    if (!duplicateSequences.isEmpty()) {
                        if (segDups == null) {
                            segDups = new BufferedWriter(new FileWriter(outputPath + "/removedLines.txt"));
                        }
                        if (duplicateSequences.indexOf(duplicateSequences.get(0)) == 0) {
                            segDups.write("\na score=0\n");
                        }
                        for (String[] arraySpecies : duplicateSequences) {
                            segDups.write(String.join("\t", arraySpecies) + "\n");
                        }
                        hasRemovedLines = true;
                    }
                    if (speciesSequences.size() > 1) {
                        if (!firstAlignment) {
                            out.write("\n\n");
                        }
                        out.write("\n\na\n");
                        for (String key : speciesSequences.keySet()) {
                            String[] value = speciesSequences.get(key);
                            out.write(String.join("\t", value) + "\n");
                        }
                        firstAlignment = false;
                    }
                }
            }
            entry.close();
        }
        out.write("\n\n");
        out.close();

        if (segDups != null) {
            segDups.close();
            if (!hasRemovedLines) {
                File removedLinesFile = new File(outputPath + "/removedLines.txt");
                removedLinesFile.delete();
            }
        }
    }
}

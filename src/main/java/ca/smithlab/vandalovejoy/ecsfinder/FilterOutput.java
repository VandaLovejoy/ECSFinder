package ca.smithlab.vandalovejoy.ecsfinder;


import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.regex.*;
import java.util.stream.*;

public class FilterOutput {
    double min_eval;
    public static double processFile(String filePath, String searchPattern) {
        double minEval = Double.MAX_VALUE;
        try {
            List<Double> eValues = new ArrayList<>();
            BufferedReader reader = new BufferedReader(new FileReader(filePath));
            String line;

            // Extract values based on the search pattern
            while ((line = reader.readLine()) != null) {
                Matcher matcher = Pattern.compile(searchPattern + "([^ ]*)").matcher(line);
                if (matcher.find()) {
                    eValues.add(Double.valueOf(matcher.group(1)));
                }
            }
            reader.close();
            if (!eValues.isEmpty()) {
                minEval = Collections.min(eValues);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return minEval;
    }


    public static double[] processTxtFile(String filePath) {
        double[] energies = new double[2]; // Default to 0.0 for both elements
        try {
            BufferedReader br = new BufferedReader(new FileReader(filePath));
            String line;
            Pattern pattern = Pattern.compile("\\(([^=]+)=([^+]+)\\+([^\\)]+)\\)");

            while ((line = br.readLine()) != null) {
                Matcher matcher = pattern.matcher(line);
                if (matcher.find()) {
                    String part2 = matcher.group(2).trim().replaceAll("\\)", "");
                    String part3 = matcher.group(3).trim().replaceAll("\\)", "");
                    energies = new double[]{Double.parseDouble(part2), Double.parseDouble(part3)};
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return energies;
    }
}

//Jack Newman
//Date: 2024-12-4: Assigment 3. 
import java.util.*;
import java.util.stream.*;

public class newFunctionsLearned {
    static int size = 6000;
    static int colLen = 5;

    public static void main(String[] args) {
        // Generate 6000 lines of input with String[] format
        List<String[]> lines = new ArrayList<>(size);
        for (int i = 0; i < size; i++) {
            lines.add(new String[]{"1.0", "2.0", "3.0", "4.0", "5.0"});
        }

        // Declare variables for each time measurement
        int iter = 1000;
        long loopTime = 0, loopWithSetAllTime = 0, streamTime = 0, streamWithSetAllTime = 0, paraStreamWithSetAllTime = 0;

        // Warm-up phase to ensure JIT optimizations are active
        for (int i = 0; i < 1; i++) {
//     	   loopTime += loop(lines);
//   	   loopWithSetAllTime += loopSetAll(lines);
//       streamTime += stream(lines);
//       streamWithSetAllTime += streamSetAll(lines);
//       paraStreamWithSetAllTime += paraStreamSetAll(lines);
        }

        // Run each method and accumulate the total time for each
        for (int i = 0; i < iter; i++) {
//      	   loopTime += loop(lines);
//      	   loopWithSetAllTime += loopSetAll(lines);
//          streamTime += stream(lines);
      	  streamWithSetAllTime += streamSetAll(lines);
//          paraStreamWithSetAllTime += paraStreamSetAll(lines);
            doExtraneousWork(iter);
            doExtraneousWork(iter);
        }

        // Output average times
        if (loopTime != 0) System.out.printf("%-18s: %.8f ms%n", "Avg Loop", loopTime / 1e6 / iter);
        if (loopWithSetAllTime != 0) System.out.printf("%-18s: %.8f ms%n", "Avg Loop/SetAll", loopWithSetAllTime / 1e6 / iter);
        if (streamTime != 0) System.out.printf("%-18s: %.8f ms%n", "Avg Stream", streamTime / 1e6 / iter);
        if (streamWithSetAllTime != 0) System.out.printf("%-18s: %.8f ms%n", "Avg Stream/SetAll", streamWithSetAllTime / 1e6 / iter);
        if (paraStreamWithSetAllTime != 0) System.out.printf("%-18s: %.8f ms%n", "Avg Parallel Stream/SetAll", paraStreamWithSetAllTime / 1e6 / iter);
    }

    // Measure loop-based approach
    private static long loop(List<String[]> lines) {
        long start = System.nanoTime();
        double[][] arr = new double[size][colLen];
        for (int i = 0; i < size; i++) {
            String[] parts = lines.get(i);
            double[] row = arr[i];
            for (int j = 0; j < colLen; j++) {
                row[j] = Double.parseDouble(parts[j]);
            }
        }
        return System.nanoTime() - start;
    }

    // Measure loop with Arrays.setAll()
    private static long loopSetAll(List<String[]> lines) {
        long start = System.nanoTime();
        double[][] arr = new double[size][colLen];
        for (int i = 0; i < size; i++) {
            String[] parts = lines.get(i);
            double[] tmp = arr[i];
            Arrays.setAll(tmp, j -> Double.parseDouble(parts[j]));
        }
        return System.nanoTime() - start;
    }

    // Measure stream-based approach
    private static long stream(List<String[]> lines) {
        long start = System.nanoTime();
        double[][] arr = new double[size][colLen];
        IntStream.range(0, size).forEach(i -> {
            String[] parts = lines.get(i);
            arr[i] = Arrays.stream(parts)
                           .mapToDouble(Double::parseDouble)
                           .toArray();
        });
        return System.nanoTime() - start;
    }

    // Measure stream with Arrays.setAll()
    private static long streamSetAll(List<String[]> lines) {
        long start = System.nanoTime();
        double[][] arr = new double[size][colLen];
        IntStream.range(0, size).forEach(i -> {
            String[] parts = lines.get(i);
            double[] tmp = arr[i];
            Arrays.setAll(tmp, j -> Double.parseDouble(parts[j]));
        });
        return System.nanoTime() - start;
    }

    // Measure parallel stream with Arrays.setAll()
    private static long paraStreamSetAll(List<String[]> lines) {
        long start = System.nanoTime();
        double[][] arr = new double[size][colLen];
        IntStream.range(0, size).parallel().forEach(i -> {
            String[] parts = lines.get(i);
            Arrays.setAll(arr[i], j -> Double.parseDouble(parts[j]));
        });
        return System.nanoTime() - start;
    }

    private static void doExtraneousWork(int iter) {
        Random rand = new Random();
        List<Integer> intList = new ArrayList<>();
        for (int i = 0; i < iter * 25; i++) {
            intList.add(rand.nextInt(1000));
        }
        Map<Integer, List<String>> map = new HashMap<>();
        for (int i = 0; i < iter * 50; i++) {
            map.put(i, Arrays.asList("a", "b", "c", "d"));
        }
        Map<Integer, Integer> lookupMap = new HashMap<>();
        for (int i = 0; i < iter * 50; i++) {
            lookupMap.put(i, rand.nextInt(1000));
        }
        helperMethod(intList, map, lookupMap);
    }

    private static void helperMethod(List<Integer> intList, Map<Integer, List<String>> map, Map<Integer, Integer> lookupMap) {
        int sum = intList.stream().mapToInt(Integer::intValue).sum();
        int mapSize = map.size();
        int lookupResult = lookupMap.getOrDefault(mapSize, 0);
        double average = sum / (double) intList.size();
        moreComplexMethod(average, mapSize, lookupResult);
    }

    private static double moreComplexMethod(double average, int mapSize, int lookupResult) {
        double result = (average * mapSize + Math.pow(average * mapSize, 2)) * lookupResult;
        if (lookupResult % 3 == 0) {
            result = Math.pow(result, 2);
        }
        if (result % 2 == 0) {
            result = Math.pow(result, 2);
        }
        return result;
    }
}

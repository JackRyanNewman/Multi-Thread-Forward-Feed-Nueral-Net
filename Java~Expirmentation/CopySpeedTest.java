//Jack Newman
//Date: 2024-12-4: Assigment 3. 
import java.util.Random;

public class CopySpeedTest {
   public static void main(String[] args) {
      final int ROWS = 10000; // Number of rows
      final int COLS = 1000;  // Number of columns
      final int ITERATIONS = 100000; // Number of iterations for the test

      // Generate a sample 2D array
      int[][] array = new int[ROWS][COLS];
      Random rand = new Random();
      for (int i = 0; i < ROWS; i++) {
         for (int j = 0; j < COLS; j++) {
            array[i][j] = rand.nextInt(100);
         }
      }

      // Test System.arraycopy (deep copy)
      long startArrayCopy = System.nanoTime();
      for (int it = 0; it < ITERATIONS; it++) {
         int[] copiedRow = new int[COLS];
         System.arraycopy(array[it % ROWS], 0, copiedRow, 0, COLS);
      }
      long endArrayCopy = System.nanoTime();
      long arrayCopyTime = endArrayCopy - startArrayCopy;

      // Test reference copy (shallow copy)
      long startRefCopy = System.nanoTime();
      for (int it = 0; it < ITERATIONS; it++) {
         int[] refRow = array[it % ROWS];
      }
      long endRefCopy = System.nanoTime();
      long refCopyTime = endRefCopy - startRefCopy;

      // Print results
      System.out.println("Time taken by System.arraycopy (deep copy): " + arrayCopyTime + " ns");
      System.out.println("Time taken by reference copy (shallow copy): " + refCopyTime + " ns");
   }
}

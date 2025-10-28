//Jack Newman
//Date: 2024-12-4: Assigment 3. 
import java.util.Random;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

public class VectorMultiTest {

    // Initialize random 2D arrays and an array of vectors
    static int Rows = 6000, Cols = 10;
    static double[][] inputs = new double[Rows][Cols];
    static double[][] weights = new double[Rows][Cols];
    static DMatrixRMaj matrix = new DMatrixRMaj(Rows, Cols); // Matrix initialization
    static DMatrixRMaj scalarsMatrix = new DMatrixRMaj(1, Cols); // 1xCols matrix initialization
    
    public static void main(String[] args) {
        Random random = new Random();
        for (int i = 0; i < Rows; i++) {
            double[] inputRow = new double[Cols];
            double[] weightRow = new double[Cols];
            for (int j = 0; j < Cols; j++) {
                inputRow[j] = random.nextDouble();
                weightRow[j] = random.nextDouble();
            }
        }
        
        // Fill the matrix with 1.0 values
        CommonOps_DDRM.fill(matrix, 1.0);
        
        // Fill the scalarsMatrix with values
        for (int i = 0; i < Cols; i++) {
            scalarsMatrix.set(0, i, i + 2);  // Setting values for row 0 (since it's 1xCols)
        }
        
        // Perform time-based tests
        long start = System.nanoTime();
        multiplyWeightsAndInputsLoop(); // Perform loop-based multiplication
        long loopTime = System.nanoTime() - start;

        start = System.nanoTime();
        test();
        long vectorTime = System.nanoTime() - start;

        System.out.println("Time (Loop-Based): " + loopTime / 1_000_000 + " ms");
        System.out.println("Time (Vectorized): " + vectorTime / 1_000_000 + " ms");
    }

    public static void test() {
        // Element-wise multiplication
   	 CommonOps_DDRM.elementMult(matrix, scalarsMatrix);
//       DMatrixRMaj rowSums = new DMatrixRMaj(Rows, 1);
//       CommonOps_DDRM.sumRows(matrix, rowSums);
    }
    
    // Loop-based multiplication for comparison
    public static void multiplyWeightsAndInputsLoop() {
        for (int row = 0; row < Rows; row++) {
            for (int col = 0; col < Cols; col++) {
                // Multiply weights with inputs
                weights[row][col] = weights[row][col] * inputs[row][col];
            }
        }
    }
}

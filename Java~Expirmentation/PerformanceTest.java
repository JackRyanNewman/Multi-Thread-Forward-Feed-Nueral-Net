//Jack Newman
//Date: 2024-12-4: Assigment 3. 
public class PerformanceTest {

	static final int iterations = 1_000_000;
   static final double div = 1_000_000_000.0; // Set divisor for per iteration time, adjust as needed 
   static final int type = 0; 
   
	public static void main(String[] args) {
        long totalInlineTime = 0;
        long totalMethodTime = 0;
        long totalLambdaTime = 0;
        double[] resultInline = {0};
        final java.util.function.BiFunction<Double, Double, Double> lambdaFunction = (a, b) -> a + b * b + a * a; 
        
        

        // Measure the total time for each computation type
        for (int iter = 0; iter < iterations; iter++) {
            long start = System.nanoTime();
//            switch(type) {
//            	case 0 -> resultInline[0] = 2 + 3 * 3 + 2 * 2;
//            	case 1 -> resultInline[0] = 2 + 3 * 3 + 2 * 2;
//            	case 2 -> resultInline[0] = 2 + 3 * 3 + 2 * 2;
//            	case 3 -> resultInline[0] = 2 + 3 * 3 + 2 * 2;
//            }
            resultInline[0] = 2 + 3 * 3 + 2 * 2;
            totalInlineTime += System.nanoTime() - start;

            start = System.nanoTime();
            resultInline[0] = computeMethod(2); // Method call
            totalMethodTime += System.nanoTime() - start;

            start = System.nanoTime();
            resultInline[0] = lambdaFunction.apply(2.0, 3.0); // Lambda call
            totalLambdaTime += System.nanoTime() - start;
        }

        // Output the total times and average times per iteration
        System.out.printf("T=Total time, Iters = %d, %n%n", iterations);
        System.out.printf("Inline computation T: %.8f  T/I: %.8f%n", totalInlineTime / (double) div, totalInlineTime / (double) iterations / div); 
        System.out.printf("Method call        T: %.8f  T/I: %.8f%n", totalMethodTime / (double) div, totalMethodTime / (double) iterations / div);
        System.out.printf("Lambda function    T: %.8f  T/I: %.8f%n", totalLambdaTime / (double) div, totalLambdaTime / (double) iterations / div);
    }

    // Method for computation (used in method call example)
    public static double computeMethod(double a) {
   	 switch(type) {
	    	case 0 -> a = 2 + 3 * 3 + 2 * 2;
	    	case 1 -> a = 2 + 3 * 3 + 2 * 2;
	    	case 2 -> a = 2 + 3 * 3 + 2 * 2;
	    	case 3 -> a = 2 + 3 * 3 + 2 * 2;
   	 }
   	 return a;
    }
}

//Jack Newman
//Date: 2024-12-4: Assigment 3. 
import java.util.*;
import java.util.stream.IntStream;

public class MinMaxSpeedTest {
	
    static final int iter = 1000, IBal = 100, InputLen = 10, setSize = 10000, ClassLen = InputLen;
    static double[][] BseSet = new double[setSize][InputLen], MinMax = new double[InputLen][2];
    static double div = 1_000_000_000.0;
    
    public static void main(String[] args) {
        for (int i = 0; i < setSize; i++) {
            for (int j = 0; j < InputLen; j++) { 
            	BseSet[i][j] = Math.random() * iter;
            }
        }
        long totalSeqTime = 0, totalParTime = 0, totalR = 0;
        for (int i = 0; i < iter; i++) {
            
//            long start = System.nanoTime();
//            sequentialMinMax(); totalSeqTime += System.nanoTime() - start;
            
            long start = System.nanoTime();
            rabbit2(); totalR = System.nanoTime()-start;
            
            Arrays.fill(MinMax, new double[]{0, 0}); 
            doExtraneousWork(100);
            
            start = System.nanoTime();
            parallelMinMax();
            
            totalParTime += System.nanoTime() - start;
            
            Arrays.fill(MinMax, new double[]{0, 0});
        }
        
        
        System.out.printf("T=Total time, Iters = %d, %n%n",iter);
        System.out.printf("Randomize  T: %.8f  T/I: %.8f%n", totalR/div, totalR  / iter / div); 
        System.out.printf("Parallel   T: %.8f  T/I: %.8f%n", totalParTime / div, totalParTime / iter / div);
        System.out.printf("Sequential T: %.8f  T/I: %.8f%n",totalSeqTime / div, totalSeqTime / iter / div);
             
      
    }

    public static void sequentialMinMax() {
        for (int col = 0; col < InputLen; col++) {
      	  double currMin = Double.MAX_VALUE;
           double currMax =  Double.MIN_VALUE;
            for (int row = 0; row < setSize; row++) {
                double tmp = BseSet[row][col];
                if (tmp < currMin) {
                    currMin = tmp;
                    MinMax[col][0] = tmp;
                } else if (tmp > currMax) {
                    currMax = tmp;
                    MinMax[col][1] = tmp;
                }
            }
        }
    }

    public static void parallelMinMax() {
        IntStream.range(0, InputLen).parallel().forEach(col -> {
      	  double currMin = Double.MAX_VALUE;
           double currMax =  Double.MIN_VALUE;
            for (int row = 0; row < setSize; row++) {
                double tmp = BseSet[row][col];
                if (tmp < currMin) {
                    currMin = tmp;
                    MinMax[col][0] = tmp;
                } else if (tmp > currMax) {
                    currMax = tmp;
                    MinMax[col][1] = tmp;
                }
            }
        });
    }
    
    public static void rabbit() {
   	 Random rand = new Random();
		 int rowCount = BseSet.length; 
   	 IntStream.range(0, rowCount).parallel().forEach(i -> {
			// Synchronize on BseSet[i] and BseSet[j] to avoid race conditions 
			  int j = rand.nextInt(rowCount);
	         synchronized (BseSet[i]) {
	             synchronized (BseSet[j]) {
	                 double[] tmp = BseSet[i];
	                 BseSet[i] = BseSet[j];
	                 BseSet[j] = tmp;
	             }
	         }
	     });
    }
    
    
    public static void rabbit2() {
   	 Random rand = new Random();
		 int rowCount = BseSet.length; 
   	 IntStream.range(0, rowCount).parallel().forEach(i -> {
			// Synchronize on BseSet[i] and BseSet[j] to avoid race conditions 
			  int j = rand.nextInt(rowCount);
	         synchronized (BseSet[i]) {
	             synchronized (BseSet[j]) {
	            	 double[] tmp = new double[BseSet[i].length];
	                System.arraycopy(BseSet[i], 0, tmp, 0, tmp.length);
	                System.arraycopy(BseSet[j], 0, BseSet[i], 0, BseSet[i].length);
	                System.arraycopy(tmp, 0, BseSet[j], 0, tmp.length);
	             }
	         }
	     });
    }
    
    
    public static void randomize() {
		 //1.5-2.5x faster than collections shuffle
		 Random rand = new Random();
		 int rowCount = BseSet.length;
		 for (int i = 0; i < rowCount; i++) {
		    int j = rand.nextInt(rowCount); 
		    double[] tmp = BseSet[i];
		    BseSet[i] = BseSet[j];
		    BseSet[j] = tmp;
		}
 }
    
    public static void randomize2() {
       Random rand = new Random();
       int rowCount = BseSet.length;
       for (int i = 0; i < rowCount; i++) {
           int j = rand.nextInt(rowCount);
           double[] tmp = new double[BseSet[i].length];
           System.arraycopy(BseSet[i], 0, tmp, 0, tmp.length);
           System.arraycopy(BseSet[j], 0, BseSet[i], 0, BseSet[i].length);
           System.arraycopy(tmp, 0, BseSet[j], 0, tmp.length);
       }
   }  
    
    
//    IntStream.range(0, InputLen).parallel().forEach(col -> {
//       double currMin = MinMax[col][0];
//       double currMax = MinMax[col][1];
//       for (int row = 0; row < setSize; row++) {
//           double tmp = BseSet[row][col];
//           if (tmp < currMin) {
//               currMin = tmp;
//               MinMax[col][0] = tmp;
//           } else if (tmp > currMax) {
//               currMax = tmp;
//               MinMax[col][1] = tmp;
//           }
//       }
//		});
//		IntStream.range(0, InputLen).parallel().forEach(col -> {
//			double min = MinMax[col][0];
//       double div = min-MinMax[col][1]; 
//       if(div==0) {
//       	for (int i= 0; i < bseSize; i++) {
//       		RegSet[i][col]=-1;
//       	}
//       }
//       else {
//       	for (int row = 0; row < bseSize; row++) {
//       		RegSet[row][col]-= 1 + 2*((BseSet[row][col]-min)/div);
//       	}
//       }
//		});
    
    
    
    
    
 	private void setMinMax() {
 		 for(int x = 0; x < ClassLen; x++) {
   			MinMax[x][0] = Double.MAX_VALUE;
   			MinMax[x][1] = Double.MIN_VALUE;
   		}
	}
    

    private static void doExtraneousWork(int iter) {
       Random rand = new Random();
       List<Integer> intList = new ArrayList<>();
       Map<Integer, List<String>> map = new HashMap<>();
       Map<Integer, Integer> lookupMap = new HashMap<>();
       for (int i = 0; i < iter * 25; i++) intList.add(rand.nextInt(iter));
       for (int i = 0; i < iter * 50; i++) map.put(i, Arrays.asList("a", "b", "c", "d"));
       for (int i = 0; i < iter * 50; i++) lookupMap.put(i, rand.nextInt(iter));
       helperMethod(intList, map, lookupMap);
   }

   private static void helperMethod(List<Integer> intList, Map<Integer, List<String>> map, Map<Integer, Integer> lookupMap) {
       int sum = intList.stream().mapToInt(Integer::intValue).sum();
       int mapSize = map.size(), lookupResult = lookupMap.getOrDefault(mapSize, 0);
       double average = sum / (double) intList.size();
       moreComplexMethod(average, mapSize, lookupResult);
   }

   private static double moreComplexMethod(double average, int mapSize, int lookupResult) {
       double result = (average * mapSize + Math.pow(average * mapSize, 2)) * lookupResult;
       if (lookupResult % 3 == 0) result = Math.pow(result, 2);
       if (result % 2 == 0) result = Math.pow(result, 2);
       return result;
   }
}

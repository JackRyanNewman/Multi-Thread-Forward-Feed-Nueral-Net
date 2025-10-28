//Jack Newman
//Date: 2024-12-4: Assigment 3. 
import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.ArrayList;
import java.util.Arrays;


public class UnRoll {
	static int rowCount = 6000, colCount = 10; 
	static Random rand = new Random();

	
	
	 public static void main(String[] args) {
       long shuffle = 0, unroll = 0, rabbit = 0;
       long startTime;
       int iter = 40;
       double[][] array1 = new double[rowCount][colCount];
       double[][] array2 = new double[rowCount][colCount];
       double[][] array3 = new double[rowCount][colCount];
       for(int i=0; i < iter; i++) {
      	 
      	 
//        startTime = System.nanoTime();
//        Collections.shuffle(Arrays.asList(array2));
//        shuffle += System.nanoTime() - startTime;
//        

        doExtraneousWork(100);
        startTime = System.nanoTime();
        unroll2(array1);
        unroll += System.nanoTime() - startTime;
    	  //doExtraneousWork(100);
    	 
    	 
        startTime = System.nanoTime();
        para2(array3);
        rabbit += System.nanoTime() - startTime;
       }
       System.out.printf("Shuffled duratiom:  %,.8f%n", shuffle/1e6/iter);
       System.out.printf("Unrolled duration:  %,.8f%n", unroll/1e6/iter);
       System.out.printf("Parallel duratiom:  %,.8f%n", rabbit/1e6/iter);

   }
	
	
	 public static void para(double[][] arr) {
   	 Random rand = new Random();
		 int rowCount = arr.length; 
   	 IntStream.range(0, rowCount).parallel().forEach(i -> {
			// Synchronize on arr[i] and arr[j] to avoid race conditions 
			  int j = rand.nextInt(rowCount);
	         synchronized (arr[i]) {
	             synchronized (arr[j]) {
	                 double[] tmp = arr[i];
	                 arr[i] = arr[j];
	                 arr[j] = tmp;
	             }
	         }
	     });
    }
    
    
    public static void para2(double[][] arr) {
   	 Random rand = new Random();
		 int rowCount = arr.length; 
   	 IntStream.range(0, rowCount).parallel().forEach(i -> {
			// Synchronize on arr[i] and arr[j] to avoid race conditions 
			  int j = rand.nextInt(rowCount);
	         synchronized (arr[i]) {
	             synchronized (arr[j]) {
	            	 double[] tmp = new double[arr[i].length];
	                System.arraycopy(arr[i], 0, tmp, 0, tmp.length);
	                System.arraycopy(arr[j], 0, arr[i], 0, arr[i].length);
	                System.arraycopy(tmp, 0, arr[j], 0, tmp.length);
	             }
	         }
	     });
    }
    
    
    public static void unroll(double[][] arr) {
		 //1.5-2.5x faster than collections shuffle
		 Random rand = new Random();
		 int rowCount = arr.length;
		 for (int i = 0; i < rowCount; i++) {
		    int j = rand.nextInt(rowCount); 
		    double[] tmp = arr[i];
		    arr[i] = arr[j];
		    arr[j] = tmp;
		}
  }
    
    public static void unroll2(double[][] arr) {
       Random rand = new Random();
       int rowCount = arr.length;
       for (int i = 0; i < rowCount; i++) {
           int j = rand.nextInt(rowCount);
           double[] tmp = new double[arr[i].length];
           System.arraycopy(arr[i], 0, tmp, 0, tmp.length);
           System.arraycopy(arr[j], 0, arr[i], 0, arr[i].length);
           System.arraycopy(tmp, 0, arr[j], 0, tmp.length);
       }
   }  
    
	
	
	
	  
  
	  
	 
	  
	  
	  
//	  public static void rabbit(double[][] array) {
//		  ForkJoinPool forkJoinPool = new ForkJoinPool();
//        forkJoinPool.submit(() -> {
//            for (int i = 0; i < rowCount; i++) {
//                int j = rand.nextInt(i + 1);
//                // Swap rows i and j
//                double[] tmp = array[i];
//                array[i] = array[j];
//                array[j] = tmp;
//            }
//        }).join();
//	  }
//	  
	  
	  
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

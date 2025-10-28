//Jack Newman
//Date: 2024-12-4: Assigment 3. 
//==================================================================

//NeuralNetWorkTester: 
//Description: This class is used to run expriements by either calling the nueralNetController, and running many expirments 
//Or by calling specific tester, which grabs a test case from a plethora of options. 

public class NeuralNetworkTester {
	
    public static void main(String[] args) throws Exception {
    	//specificTester(0, 8, true);
    	NeuralNetController.main();
    }

    //This takes in 3 arugment, type is what type of methodd you want to call, opt is the switch case number you want to select, and par is if you want it to run in parallel. 
    public static void specificTester(int type, int opt, boolean par) throws Exception {
   	 String tmp = switch(type) {
	   	 case 0 -> quickTest(opt);
	   	 case 1 -> documentCheck(opt);
	   	 default -> "";
   	 };
   	 if(par) tmp+=" -p"; 
   	 System.out.println(tmp);  
   	 double start = System.nanoTime();
   	 NeuralNetController.main(tmp.split(" "));
   	 System.out.printf("\n•› Total runtime in %.3fms%n", (System.nanoTime() - start) / 1_000_000);
    }
      
    
    //This was used as a quick place to put in any sort of test. 
    public static String quickTest(int option) {
       return switch (option) {
           case 0 -> "-f ../457-ML-03-JN/a03-data/iris-shuf.dat -v 2";
           case 1 -> "-f ../457-ML-03-JN/a03-data/image-05.dat -v 2";
           case 2 -> "-f ../457-ML-03-JN/a03-data/image-10.dat -v 2";
           case 3 -> "-f ../457-ML-03-JN/a03-data/image-15.dat -v 2";
           case 4 -> "-f ../457-ML-03-JN/a03-data/image-20.dat -v 0 -h 5 10 7 20 10 5 -% 20";
           case 5 -> "-f ../457-ML-03-JN/a03-data/mnist.dat -v 0 -h 5 10 7 20 10 5 -% 20";
           case 6 -> "-f ../457-ML-03-JN/a03-data/image-05.dat -v 3 -w 0.01 -a 0.005 -e 500 -l 0.0 -% 10 -g softplus -h 3 6 6 6";
           case 7 -> "-f ../457-ML-03-JN/a03-data/mnist.dat -v 1 -h 2 70 175 -e 100 -w 0.1 -l 0.0001 -a 0.04";
           case 8 -> "-f ../457-ML-03-JN/a03-data/image-10.dat -e 1000 -a 0.01 -l 0.005 -w 0.05 -% 0 -h 5 15 20 15 15 20";
           default -> ""; // Default case for invalid option
       };
    }
      
    //This was to test all my output agiest your examples to confirm if they were right. 
    public static String documentCheck(int option) {
   	 return switch (option) {
       case 0 -> "-f ../457-ML-03-JN/a03-data/iris-shuf.dat";
       case 1 -> "-f ../457-ML-03-JN/a03-data/iris-shuf.dat -v 2";
       case 2 -> "-f ../457-ML-03-JN/a03-data/iris-shuf.dat -v 3 -h 1 2";
       case 3 -> "-f ../457-ML-03-JN/a03-data/iris-shuf.dat -v 3 -e 5000 -h 2 3 2 -l 0.0001";
       case 4 -> "-f ../457-ML-03-JN/a03-data/iris-small.dat -v 4 -h 1 2 -e 1 -w 0 -a 0.5";
       default -> ""; // Default case for invalid option
   };
   }
 

    
    
}

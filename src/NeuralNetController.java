//Jack Newman
//Date: 2024-12-4: Assigment 3. 
//Project Overview: This program trains and evaluates a multi-class classification neural network. It allows users to define the structure of a multi-layer feed-forward network, then trains it using backpropagation
//with mini-batch gradient descent. The program supports parallel calculations for mini-batch processing, loss, and accuracy. It also includes a dynamic hyperparameter generator for testing 
//different configurations quickly, running them in parallel, logging results, and saving the top-performing models.
//==================================================================
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ConcurrentLinkedQueue;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;


//NeuralNetController: 
//Descripition: This file is multi purpose, one is it can take in command line arguments to do any configuration, or you tweak varibles found in HP.D and run this classes main without any arguments. 
//Which then runs a bunch of different networks. The primary purpose of this class however, is to parse the dataset, minimax scale it, and create the minibatches, and creating the neuralnetwork set, and then it

public class NeuralNetController {
	static double[][] BseSet;     //Holds the orginal Data set. 
	//RegSet is innfiecent, but it was made for the case, in which, you wanted to start over, and rescale the set. It adds another pointer
	//And double the data. But I was late in development, so I kept it. 
	static double[][] RegSet;     //regularized data of bse set, if not randomize then this points to bseSet, and it is regualrized
	static short IOLen, OPLen;    //IOlen = amount of inputs per row, O 
	static AtomicInteger Trial;   //Thread safe counter for my main method
	static boolean Store = false; //Controlls if you want to write to a file 
	private static ConcurrentLinkedQueue<ConfigEntry> topEntries = new ConcurrentLinkedQueue<>(); 
	//This will save the top 20 configs, but you have to currently uncomment on where they exist. 
  
	
//==============================================================================================
//Hyperparamter Expirmentation Methods. Class constructor 
	
	
	//Main: Generates many perumtations and then runs mutplie network in parallel in linear fashion. 
	public static void main() throws Exception {		
	//My notes and intinal observsations. 
	//Note: x-y (a-b) means uniformly x through y features all have the min of a, and max b
		
		final byte V = 0;
		final HP HP = new HP();
		final boolean pSrch= true;  
		final boolean R = true, P = false; 		
		Trial = new AtomicInteger(0);  // count of correct predictions
	 
		final String[] fp = {  //Files to hyperparamterized. 								
			 "../457-ML-03-JN/a03-data/iris-shuf.dat", //0: 4 IO, 3OP. 
			 "../457-ML-03-JN/a03-data/image-05.dat",  //1: 65 IO,  14 OP, 0-40 (0-1), 41-64 (0-255) : sz = 2516
			 "../457-ML-03-JN/a03-data/image-10.dat",  //2: 140 IO, 14 OP, continue 41-140 (0-255)	 : sz = 2516
			 "../457-ML-03-JN/a03-data/image-15.dat",  //3: 265 IO, 14 OP, continue 41-255 (0-255)	 : sz = 2516
			 "../457-ML-03-JN/a03-data/image-20.dat",  //4  440 IO, 14 OP, continue 41-440 (0-255)	 : sz = 2516
			 "../457-ML-03-JN/a03-data/mnist.dat"      //5  784 IO, 10 OP, Major ranges 506 (0-255), 130(0-0) 32 0(0-254), 26(0-253): sz 5000
		}; 
		if(Store) HP.openDualStream();
		String fileName = fp[4];
		createBseSet(fileName, R, V);
		ArrayList<String[]> trials = HP.createTrials(0, true);
		System.out.println("Search For Best Parameters for file:"+fileName);
		System.out.println("› Searching over "+ trials.size()+" trials\n");
		double start = System.nanoTime();
		
		if(!pSrch) for(String[] config: trials) trainNeuralNet(config, R, V, P);
		trials.parallelStream().forEach(config -> trainNeuralNet(config, R, V, P));
		double end = (System.nanoTime() - start)/1_000_000; 	
		System.out.printf("•› %d Completed in %.3fms, (%.1f per trial)%n%n", trials.size(), end, end/trials.size());
		//printConfigs();
		if(Store) HP.closeDualStream();
	} 
		
	//Modfied version of my other trainNueralNet method, and this one however takes in arguments, and constructs data as needed, instead of passsing the data. 
	//This is because this trainNueralNet is called iteraively or in parallel, so it needs all the information in one spot. I would suggest following the other main
	//method before as its more organized 
	private static void trainNeuralNet(String[] args, final boolean r, final byte v, boolean P) {
		double initWE = 0.1;   				// -w <DOUBLE>   	Weight initialization parameter; default 0.1 
		double alpha = 0.01;     	   	// -a <DOUBLE>		Learning rate for gradient descent; default 0.01
		double lambda = 0;         		// -l <DOUBLE>		Regularization parameter; default 0.0
		int btchSz = 1;           		// -m <INTEGER>	Batch Sz; 0 for full batch; default 1
		int eLim = 1000;       				// -e <INTEGER>	Epoch limit for gradient descent; default 1000: 
	  int[] hiddenLayers = {0}; 		// -h <NH> <S1>	Number and Szs of each hidden layers; default 0 layers, max is 10 layers, 500 nuerons per layer.  
	  byte AF = 0; 	   							// -g <byte> 		activation Function, value will come as its name. 
	  byte bDiv = -1; 							// -% <byte> 	  !CUSTOM FLAG!, btchSize will now be determined by, (bseSetLen*.8) * bScalar. Default is -1 which uses default btchSz
	  double[][] regSet = r? new double[RegSet.length][RegSet[0].length]: RegSet;
	  if(r) for (int i = 0; i < regSet.length; i++) regSet[i] = Arrays.copyOf(RegSet[i], RegSet[i].length); 	
	 	
		for (int i = 0; i < args.length; i++) {
		   switch (args[i]) {
		      case "-w" -> initWE = Double.parseDouble(args[++i]);
		      case "-a" -> alpha = Double.parseDouble(args[++i]);
		      case "-e" -> eLim = Integer.parseInt(args[++i]);
		      case "-l" -> lambda = Double.parseDouble(args[++i]);
		      case "-m" -> btchSz = Integer.parseInt(args[++i]);
		      case "-g" -> AF = setActivationFunction(args[++i]);
		      case "-%" -> bDiv = Byte.parseByte(args[++i]);
		      case "-h" -> {
		      	int nh = Integer.parseInt((args[++i]));
		         if(nh!=0) { //If some reason, he species the default 
		        	 hiddenLayers = new int[nh];
		        	 for (int j = 0; j < nh; j++) {
		               hiddenLayers[j] = Short.parseShort(args[++i]);
		            }
		         }
		      }
		   }  
		}
		
		
	
		if(bDiv>0) btchSz = ((int) (BseSet.length * .8)) / bDiv;
		else if(btchSz==0 || bDiv == 0 || btchSz > regSet.length) btchSz = (int) (BseSet.length * .8);
		
		int[][] ret = initMB(btchSz, regSet.length, P);
		int[] TrainSet = ret[0], FastCalc = ret[1];
		NeuralNet neuralNet = new NeuralNet(hiddenLayers, initWE, alpha, lambda, ret[2][0], IOLen, OPLen, v, AF, P);
		
		String printVal = "";
		int eCount = 1, iterCount = 0; 
		boolean stop = false, print = false;
		int pInterval = eLim <= 10? 1: eLim/10; 
		int iters = TrainSet.length-1; //how many iteratations occur per minibatch. 
		double start = System.nanoTime();
		
		for(; eCount <= eLim; eCount++){
			if(r) randomize(regSet); //Caught an error I put in RegSet?
			if(P) neuralNet.paralleltTrain(regSet, TrainSet, btchSz);
			else  neuralNet.linearTrain(regSet, TrainSet);
			iterCount+=iters;
			
			if(!P && v>1 && eCount%pInterval==0) {
				printVal = String.format("    After  %6d epochs (%7d iter.)", eCount, iterCount);
				print = true;
			}
			if(P) stop =  neuralNet.pCalc(regSet, TrainSet, FastCalc, print,printVal);
			else stop = neuralNet.lCalc(regSet, TrainSet, print,printVal);
			print = false;
			if(stop) break; 
		}
		
		float[] results = new float[2]; 
		double end =(System.nanoTime() - start)/1_000_000; 		
		
		int trialNum = Trial.incrementAndGet();
		String header = "•›TRIAL: " + trialNum; 
		String massPrint = "\n ›Flags: " + String.join(" ", args);
		massPrint+=String.format("\n ›Time %.0f ms, %d epochs, %d iters:(%.4f ms/iter)", end, eCount-1, iterCount, end/iterCount);
		String stopCon = stop? "Abouslute Error at every output neuron was <= .01": "Epoch Limit";
		massPrint+=": Stopped at: "+stopCon;
		if(P) massPrint+=neuralNet.pAcc(regSet, TrainSet, FastCalc, results);
		else  massPrint+=neuralNet.lAcc(regSet,TrainSet, results);
		if(results[1]!=0) {
			System.out.printf("%s •› Completed in %.3fms%s%n",header, end, massPrint);
			//topEntries.add(new ConfigEntry(args, results));
		}
		else System.out.printf("%s » Completed in %.3fms%n%n",header, end);
	} 
  
	
	//Optional method, if I wanted to print the best configrations at the end. 
	private static void printConfigs(){
	 List<ConfigEntry> sortedEntries = new ArrayList<>(topEntries);
	 sortedEntries.sort((e1, e2) -> Float.compare(e2.score[1], e1.score[1]));
	 sortedEntries = sortedEntries.stream().limit(20).collect(Collectors.toList());
	 sortedEntries.forEach(entry -> 
	    System.out.println("Params: " + Arrays.toString(entry.params) + ", Score: " + Arrays.toString(entry.score))
	  );
	}
	
//==============================================================================================
//Main call either through arguments or through my NeuralNetWork tester class. 
	
	public static void main(String[] args) throws Exception{
		String fileName = "";     // -f <FILENAME>	File to read data from
		double initWE = 0.1;   	  // -w <DOUBLE>   	Weight initialization parameter; default 0.1 ///////////////////////////////// TWO PLACES he say two different things
		double alpha = 0.01;     	// -a <DOUBLE>		Learning rate for gradient descent; default 0.01
		double lambda = 0;        // -l <DOUBLE>		Regularization parameter; default 0.0
		int btchSz = 1;           // -m <INTEGER>	Batch Sz; 0 for full batch; default 1
		int epochLimit = 1000;    // -e <INTEGER>	Epoch limit for gradient descent; default 1000: 
	  int[] hiddenLayers = {0}; // -h <NH> <S1>	Number and Szs of each hidden layers; default 0 layers, max is 10 layers, 500 nuerons per layer.  
	  byte v = 1; 		  				// -v <INTEGER>  !MODFIED!		,  vrb level; default 1 | added 0 for silent mode, for hyperparamter tuning. 
	  byte AF = 0; 	   					// -g <byte> 		activation Function, value will come as its name. 
	  boolean r = false;  			// -r <Boolean>   rnds data if specified; default is false
	  boolean P = false; 			  // -p <boolean>  !CUSTOM FLAG!, items in mb are done in parallel. WARNING: Depending on NueralNetwork set up this may crash the program. 
	  byte bDiv = -1; 					// -% <byte> 	  !CUSTOM FLAG!, btchSize will now be determined by, (bseSetLen*.8) * bScalar. Default is -1 which uses default btchSz
   	
   	
		for (int i = 0; i < args.length; i++) {
		   switch (args[i]) {
		   		case "-f" -> fileName = args[++i];	
		      case "-w" -> initWE = Double.parseDouble(args[++i]);
		      case "-a" -> alpha = Double.parseDouble(args[++i]);
		      case "-e" -> epochLimit = Integer.parseInt(args[++i]);
		      case "-l" -> lambda = Double.parseDouble(args[++i]);
		      case "-m" -> btchSz = Integer.parseInt(args[++i]);
		      case "-g" -> AF = setActivationFunction(args[++i]);
		      case "-%" -> bDiv = Byte.parseByte(args[++i]);
		      case "-v" -> v = Byte.parseByte(args[++i]);
		      case "-r" -> r = true;  
		      case "-p" -> P = true;  
		      
		      case "-h" -> {
		      	int nh = Integer.parseInt((args[++i]));
		         if(nh!=0) { //If some reason, he species the default 
		        	 hiddenLayers = new int[nh];
		        	 for (int j = 0; j < nh; j++) {
		               hiddenLayers[j] = Short.parseShort(args[++i]);
		            }
		         }
		      }
		   }  
		}		
		
		if(!fileName.equals("")) {
			if(P) System.out.println("Parallel mode: large networks may cause memory overload and crash.");
			createBseSet(fileName, r, v); //inits data, and minimax scales it.
			if(bDiv>0) btchSz = ((int) (BseSet.length * .8)) / bDiv;
			else if(btchSz==0 || bDiv == 0 || btchSz > RegSet.length) btchSz = (int) (BseSet.length * .8); //adjusting btchSize basedd off requirements. 
			int[][] ret = initMB(btchSz, RegSet.length, P); //See method for further explnation of items 
			int[] TrainSet = ret[0]; //The minibatch Set indices
			int[] FastCalc = ret[1]; //FastCalc, if parallel is active, it is holds indices of a minibatch, if the batchsize were .2 of the data set, or it matches the trainset if the batch size is larger. 
			//By setting the fastCaalc By ensuring the minibatch size is at least 0.2% of the total dataset, we can efficiently compute the validation set in parallel
			
			NeuralNet neuralNet = new NeuralNet(hiddenLayers, initWE, alpha, lambda, ret[2][0], IOLen, OPLen, v, AF, P); //Creates my nueral net objects nessarcy for training. 
			trainNeuralNet(neuralNet, TrainSet, FastCalc, alpha, lambda, epochLimit, btchSz, r, v, P);							
		}
		else System.out.print("The file name was incorrect");
	}	
	
	//This handles training a nueralNet for the main 
	private static void trainNeuralNet(NeuralNet neuralNet,int[] TrainSet, int[] FastCalc, double alpha, double lambda, int eLim, int btchSz, final boolean r, final byte v, final boolean P) {
		if(v!=0) {
			System.out.println("* Training network (using " + (int) (BseSet.length * .8) + " examples)");
			if(v>1) {
				int minSz = neuralNet.NetWork[0].Aj.length; 
				System.out.println("  * Beginning mini-batch gradient descent");
				if(P) System.out.printf("    (batchSize=%d, epochLimit=%d, learningRate=%.4f, lambda=%.4f, quick-calc-size=%d)%n" ,btchSz, eLim, alpha, lambda, minSz);
				else  System.out.printf("    (batchSize=%d, epochLimit=%d, learningRate=%.4f, lambda=%.4f)%n" ,btchSz, eLim, alpha, lambda);
			}
			if(v>2) {
				if(P) neuralNet.pCalc(RegSet, TrainSet, FastCalc, true, "    Initial model with random weights  :");
				else neuralNet.lCalc(RegSet, TrainSet, true, "    Initial model with random weights  :");
			}
		}
		
		int eCount = 1; 	//epoch count 
		int iterCount = 0; 
		boolean stop = false; 
		boolean print = false;
		String printVal = "";
		int pInterval = eLim <= 10? 1: eLim/10; //Determines when we should print during epochs 
		int iters = TrainSet.length-1; //how many iteratations occur per minibatch. 
		double start = System.nanoTime();
		
		
		for(; eCount <= eLim; eCount++){
			if(r) randomize(RegSet); //Re-orders orginal dataset, to be split up again. 
			if(P) neuralNet.paralleltTrain(RegSet, TrainSet, btchSz);
			else  neuralNet.linearTrain(RegSet, TrainSet);
			iterCount+=iters;
			
			if(v>1 && eCount%pInterval==0) {
				print = true;
				printVal = String.format("    After  %6d epochs (%7d iter.)", eCount, iterCount);
			}
			if(P) stop =  neuralNet.pCalc(RegSet, TrainSet, FastCalc, print, printVal);
			else stop = neuralNet.lCalc(RegSet, TrainSet, print,printVal);
			print = false;
			if(stop) break; 
		}
		double end =(System.nanoTime() - start)/1_000_000; 		
		if(v>1) {
			if(stop) {
				printVal = String.format("    After  %6d epochs (%7d iter.)", eCount, iterCount);
				if(P) neuralNet.pCalc(RegSet, TrainSet, FastCalc, print, printVal);
				else neuralNet.lCalc(RegSet, TrainSet, print,printVal);
			}
			String stopCon = stop? "Abouslute Error at every output neuron was <= .01": "Epoch Limit";
			System.out.println("  * Done with fitting!");
			System.out.printf ("    Training took %.0f ms, %d epochs, %d iterations ( %.4f ms/ iteration)\n", end, eCount-1, iterCount, end/iterCount);
			System.out.println("    GD Stop condition: "+stopCon);
		}
		if(v>0) {
			System.out.println("* Evaluating accuracy");
			if(P) neuralNet.pAcc(RegSet,TrainSet, FastCalc);
			else neuralNet.lAcc(RegSet,TrainSet);
		}
		if(v==0) {
			System.out.printf ("    Training took %.0f ms, %d epochs, %d iterations ( %.4f ms/ iteration)\n", end, eCount-1, iterCount, end/iterCount);
			if(P) neuralNet.pAcc(RegSet,TrainSet, FastCalc);
			else neuralNet.lAcc(RegSet,TrainSet);
			System.out.printf("•› Completed in %.3fms%n%n", end);
		}
	}
	
//==============================================================================================
//Data set configurations.  
	
	private static byte setActivationFunction(String name){
		byte tmp = 0;
		switch(name.toLowerCase()) {
			case "logistic" -> tmp=0;
			case "softplus" -> tmp=1;
			case "relu" -> tmp=2;
			case "tanh" -> tmp=3;
		}
		return tmp;
	}
	
	private static void createBseSet(String fileName, boolean r, byte v) throws Exception {	
		if(v!=0) System.out.println("* Reading"+fileName);
		String line;  																		//tmp string, to traverse non data. 
   	ArrayList<double[]> parser = new ArrayList<>();  					//Keeps, track of every row of data.   
   	BufferedReader br = new BufferedReader(new FileReader(fileName)); //Buffered reader to read 
   	do line = br.readLine(); while(line.charAt(0) != '('); 				//I read until the first entry of data. 
   	
   	int colLen = line.replaceAll("[()]", "").split(" ").length;	 //Gets us how many columns are needed. I first remove all ( or ), and then blank spaces, then get Sz of resulting array. 
   	OPLen = (short) ((line.split("\\) \\(")[1]).split(" ")).length;// ) ( is indeintifer, to find y postion, [1] would be the postion, and length would be the amount of outputs. 
   	//Lesson learned, its expecting a regular expression, and () need //
   	IOLen = (short) (colLen - OPLen); 							//The amount of inputs, and above is amount of classes. 
   	int trueColLen = IOLen + 1; 
   	
   	
   	do { //I insert all the datum, and remove all whitespace and parens from them.  
   		String[] tmp = line.replaceAll("[()]", "").split(" "); //I choose to convert them to doubles right away
   		double[] row = new double[trueColLen]; 							 //Considering the storage cost, and double accesses 
   		int j = 0; 
   		for(; j < IOLen; j++) {							 //Required. 
           row[j] = Double.parseDouble(tmp[j]);
   		}	
   		for(int loc=0; j < colLen; j++, loc++) {	//I reduce output column to be only its valid output 
            double tmp2 = Double.parseDouble(tmp[j]);
   			if(tmp2!=0) {
   				row[trueColLen-1] = loc; //You save  where the 
   				break; //You save the postion of the the element! cause that will corespond to the correct nueron! 
   			}
    		}
   		parser.add(row);
   	} while((line = br.readLine()) != null); 
   	if(r) Collections.shuffle(parser);
   	br.close();
   	
   	int rowLen = parser.size();
   	BseSet = new double[rowLen][trueColLen];
   	RegSet = r? new double[rowLen][trueColLen]: BseSet; 
   	
   	for(int i=0; i < rowLen; i++) { 
   		System.arraycopy(parser.get(i), 0, BseSet[i], 0 , trueColLen);
   	  if (r) System.arraycopy(parser.get(i), 0, RegSet[i], 0 , trueColLen);
   	}	
   	minMaxScale(v);
   	
   }
		
//Apply min-max normalization to scale feature values to the range [-1, 1]. For each feature j, 
//use its min (Lj) and max (Uj) in the training set to transform the values: x'_ij = -1 + 2 * ((x_ij - Lj) / (Uj - Lj))
	private static void minMaxScale(byte v){
		int bseSz = BseSet.length;
		int setSz = (int) (bseSz * .8);
		final double[][] MinMax = new double[IOLen][2]; 	
		for (int col = 0; col < IOLen; col++) {
          double currMin = Double.MAX_VALUE;
          double currMax =  Double.MIN_VALUE;
          for (int row = 0; row < setSz; row++) {
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
	
		if(v!=0) {
			System.out.println("* Doing train/validation split");
			System.out.println("* Scaling features");
		}
		if(v>1) {
			System.out.println("  * min/max values on training set");
			for(int i=0; i < IOLen; i++) {
				System.out.printf("      Feature %d: %.4f, %.4f%n", i + 1, MinMax[i][0], MinMax[i][1]);
			}
		}
		
		
		for(int c = 0; c < IOLen; c++) {
			double min = MinMax[c][0];
         double div = MinMax[c][1]-min; 
         if(min==MinMax[c][1]) {
         	for (int i= 0; i < bseSz; i++) {
         		RegSet[i][c]=-1;
         	}
         }
         else {
         	for (int row = 0; row < bseSz; row++) {
         		RegSet[row][c]= -1 + 2*((BseSet[row][c]-min)/div);
         	}
         }
		}
		
	}
 	
	
//Randomizes passed in set. 	
	private static void randomize(double[][] RegSet) {
	//1.5-2.5x faster than collections shuffle
	 Random rand = new Random();
  	 int bseSz = (int) (BseSet.length * .8);
  	 for (int i = 0; i < bseSz; i++) {
		    int j = rand.nextInt(bseSz); 
		    double[] tmp = RegSet[i];
		    RegSet[i] = RegSet[j];
		    RegSet[j] = tmp;
		}
	}
	
//Intilizes minibatch set and quickCalc set. 	
	private static int[][] initMB(int btchSz, int bseSz, boolean P) {
		int ret = 1; 
		int[] TrainSet, FastCalc = null; 
		int trnSz = (int) (bseSz * .8);
		if(btchSz==trnSz) TrainSet = new int[] {0, trnSz};
		else { //else its minibatch grd descent or stoachastic. 
			int arrSz = trnSz/btchSz; //How many minibatches we will have. 
			TrainSet = trnSz%btchSz!=0? new int[arrSz+2]: new int[arrSz+1];
			for(int i=1, loc=btchSz; i< TrainSet.length-1; i++, loc+=btchSz) {
				TrainSet[i]=loc; 
			}
			TrainSet[TrainSet.length-1] = trnSz; 
		}
		if(P) { //If parallel exists then make another set of indices so we can 
			int mbSz = bseSz-trnSz;
			int arrSz = trnSz/mbSz; //How many minibatches we will have. 
			FastCalc = trnSz%mbSz!=0? new int[arrSz+2]: new int[arrSz+1];
			for(int i=1, l=mbSz; i< FastCalc.length-1; i++, l+=mbSz) FastCalc[i]=l; 
			FastCalc[FastCalc.length-1] = trnSz;
			ret = mbSz>btchSz? mbSz: btchSz; 
		//By setting the fastCaalc By ensuring the minibatch size is at least 0.2% of the total dataset, we can efficiently compute the validation set in parallel
		}
		return new int[][]{TrainSet, FastCalc, {ret}};
	}
	
//This class, stores the results of one nueral net train. 	
  public static class ConfigEntry {
    final String[] params; //The string paramaters
    final float[] score;   //The scores 

    ConfigEntry(String[] params, float[] score) {
        this.params = params;
        this.score = score;
    }
}
	
	
}

	

//Jack Newman
//Date: 2024-12-4: Assigment 3. 
import java.util.Random;
import java.util.stream.IntStream;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

//NeuralNet:
//Description: This class creates and trains the nueral net. Note this class has two seperate sets of methods that related to each other. 
//if a method begins with an l its means the nueral network calculations are happening linearly, while if a methodd starts with a p, which stands for parallel, 
//Then it means the neural net work calculations are being done in parallel. 

public class NeuralNet {
	final double Alpha, Lambda; //Just globals varibles I dont want to exesively pass through. 
	final Layer[] NetWork; 			//The object the represent my nueral net. 
	final byte AF, V;						//AF stands for activation function, v stands for verbosity.
	final boolean P; 						//P stands for parallel. 
	final int NwLen;						//NWlen stands for netWork length. 
	
  //==============================================================================================	
	
	//NueralNet: Intilizes the class and sets all global varibles, and does some printing. 
	public NeuralNet(int[] HL, double initWE, double alpha, double lambda, int minSz, short IOLen, short OPLen, byte v, byte af, boolean p) {
		NwLen = HL[0]==0? 1: HL.length+1;
		V=v; AF=af; P=p;
		Alpha = alpha; 
		Lambda = lambda; 
		NetWork = new Layer[NwLen]; 		
	
		int tmp = 1; 
		if(V!=0) System.out.println("* Building network");
		if(V>1) {
			if(V>1)  	System.out.println("  * Layer sizes (excluding bias neuron(s)):");
			if(V>1)  	System.out.printf("     Layer %-2d (input) : %d\n", tmp++, IOLen); 
			if(NwLen>1) System.out.printf("     Layer %-2d (hidden): %d\n", tmp++, HL[0]);
		}
		if(NwLen==1) NetWork[0] = new Layer(IOLen, OPLen, initWE, minSz); //If only one layer, then 
		else NetWork[0] = new Layer(IOLen, HL[0], initWE, minSz);  
		for(int n=1; n < NwLen-1; n++) { //add hidden layers. 
			int prvLayerSz = HL[n-1];  //Incoming weights per neuron, or amount of input. 
			int nxtLayer = HL[n];      //Amount of Nuerons for next Layer
			NetWork[n] = new Layer(prvLayerSz, nxtLayer, minSz, initWE); //calls hidden layer constructor. 
			if(V>1)System.out.printf("     Layer %-2d (hidden): %d\n", tmp++, nxtLayer);
		}	
		if(V>1)	 System.out.printf("     Layer %-2d (output): %d\n", NwLen+1, OPLen);
		if(NwLen>1) { //NwLen != 0, and if it not 1, then i have to make output class. 
			NetWork[NwLen-1] = new Layer(HL[NwLen-2], OPLen, minSz, initWE); 
		}
		
	}
	
	//Layer: Represents each layer in my nueral net. 
	public class Layer{
		final double[] BiasWeights;//Each bias weight, corresponds to one nueron. 
		final double[][] Neurons; 	//Neuron Storage to storage currTotal
		final double[][] Aj; 		    //incoming IO*wieghts to get Inj.  
		final double[][] Inj; 		  //Result of Aj*weights, or weight update value.  
		final double[][] NUpdates;  //In order, 0=w14=^4a1, ^4a2, ^4a3, 3=w15=^5a1, ^5a2, ^5a3, Assume for that example, input layer is 3, and output is 2.  
		final double[][] BUpdates; 
		
		
		//First layer, could be treated as input, or output. 
		public Layer(int prvLayerSz, int newLayerSz, double initWE, int maxBhSz) { 
			NUpdates = new double[maxBhSz][newLayerSz*prvLayerSz];
			Neurons = new double[newLayerSz][prvLayerSz];
			BUpdates = new double[maxBhSz][newLayerSz]; 
			Aj = new double[maxBhSz][];
			Inj = new double[maxBhSz][newLayerSz];
			BiasWeights =  new double[newLayerSz]; 
			intilizeWeights(initWE, prvLayerSz);
		}
		
		//Hidden Layer Nueron. 
		public Layer(int prvLayerSz, int newLayerSz, int maxBhSz, double initWE){
			NUpdates = new double[maxBhSz][newLayerSz*prvLayerSz];
			Neurons = new double[newLayerSz][prvLayerSz];
			BUpdates = new double[maxBhSz][newLayerSz]; 
			Aj = new double[maxBhSz][prvLayerSz];
			Inj = new double[maxBhSz][newLayerSz];
			BiasWeights =  new double[newLayerSz]; 
			intilizeWeights(initWE, prvLayerSz);
		}
		
		
		public void intilizeWeights(double initWE, int weightLen) {
			Random rand = new Random();		
			for(int n=0; n < Neurons.length; n++) {
				double[] neuron = Neurons[n];
				BiasWeights[n] = (rand.nextFloat() * 2 * initWE) - initWE;
				for(int w=0; w < weightLen; w++) neuron[w] = (rand.nextFloat() * 2 * initWE) - initWE;		
			}
		}
	}

	//A printer object to handle verbsoity 4 printing, more spefically it was created to handle printing verbsoity 4 in parallel. 
	private class Printer{
		String MassPrint; //The main object that holds all my print info while im doing a parallel train run. I save it then print it at the end so nothing coming out of order
		String Inj; //The INJ,AJ,Y,G are are all saved due that fact how i caclucate and store these values that are done in way not good for printing. 
		String Aj;
		String Y;  //Y value
		String G;  //Gradient value 
	
		//Sets the intinial things needed for a print. 
		public Printer(double[] row, int ex){
			Inj = String.format("        Layer %-2d (output):    in_j:", NwLen+1);
			Aj = "\n                               a_j:";
			Y = "                 examples actaul y:";
			G = "      * Backward Propagation on example "+ex+"\n"+
			    "        Layer "+NetWork.length+" (output): Delta_j:";
			
			String s = "      * Foward Propgation on example "+ex+"\n"+
			  		     "        Layer 1  (input) :     a_j:  1.000";
			for(int aj=0; aj<row.length-1; aj++) {
				s=s+String.format("  %.3f", row[aj]);
			}
			MassPrint = P? s: "";
			if(!P) System.out.print(s); 
		}
				
		private void printHidden(int layerNum, double[] inj, double[] aj){
			  String s = String.format("\n        Layer %-2d (hidden):    in_j:", layerNum);
			  for(int in=0; in<inj.length; in++) {
				  s=s+String.format("  %.3f", inj[in]);
			  }
			  s+="\n                               a_j:";
			  for(int a=0; a< aj.length; a++) {
				  s=s+String.format("  %.3f", aj[a]);
			  }
			  if(P) MassPrint+=s+"\n";
			  else System.out.println(s);
		 }
		 
		private void addOutput(double in , double aj, double y, double gradient) {
			Inj+=String.format("  %.3f", in); 
			Aj+=String.format("  %.3f", aj); 
			Y+=String.format("  %.3f", y); 
			G+=String.format("  %.3f", gradient); 
		}
		
		private void addGradients(double[] inj, int layer) {
			G+= String.format( "\n        Layer %d (output): Delta_j:", layer);
			for(int in=0; in < inj.length; in++) {
				G+=String.format("  %.3f", inj[in]); 
			}
		}
		
		private void printAll() {
			System.out.print(MassPrint+Inj+Aj+"\n"+Y+"\n"+G+"\n\n");
		}
		
	}

 //==============================================================================================
 //Training controllers. 
	 public void linearTrain(double[][] RegSet, int[] TrainSet) {
		Printer printer = null; 
		int yI = RegSet[0].length-1;
		
		for(int mb=0; mb< TrainSet.length-1; mb++) { //Mini Batch Controller. 
			int strt = TrainSet[mb], end = TrainSet[mb+1], btchSz = end-strt;
			for(int row=strt; row<end; row++) {
				if(V==4) printer = new Printer(RegSet[row], row+1);
				NetWork[0].Aj[0] = RegSet[row];
				fwdProp(0, printer);
				bwdProp(0, (int) RegSet[row][yI], printer); 	
			}
			linearUpdateWeights(Alpha, Lambda, NwLen, btchSz);
		}
	}
	
	 public void paralleltTrain(double[][] RegSet, int[] TrainSet, int btchSz) {
			int yI = RegSet[0].length-1;
			for(int mb=0; mb< TrainSet.length-1; mb++) { //Mini Batch Controller. 
				final int safeMb = mb;
				int strt = TrainSet[mb], end = TrainSet[mb+1];
				IntStream.range(strt, end).parallel().forEach(row -> {	
					Printer printer = V==4? new Printer(RegSet[row], row+1): null;
					int nsr = row >= btchSz? row-(btchSz*safeMb): row; 	//NSR=nueron storage row, loc for thread safety.
					NetWork[0].Aj[nsr] = RegSet[row]; //First incoming Input. 
					fwdProp(nsr,printer);
					bwdProp(nsr, (int) RegSet[row][yI],printer); 
				});
				pUpdateWeights(Alpha, Lambda, NwLen, btchSz);
			}
		}
 //===========================================================================================
 //Propatation methods. 
	 
	 private void fwdProp(int nsr, Printer printer) {
		 	for(int l=0; l < NwLen-1; l++) {	//This does foward propatagion, however not for the last layer. 
				Layer layer = NetWork[l];				//Layer hold set of nuerons.
				double[][] neurons = layer.Neurons; 			//Nuerons holds set of weights. 
				double[] IO = layer.Aj[nsr];						  //output of prev layer, each output maps to each i in wij.
				double[] Inj = layer.Inj[nsr];					  //Result of w*IO = inj, each results maps to each nueron in next layer.
				double[] nextAJ = NetWork[l+1].Aj[nsr];		//activaton result of inj. 
				double[] biasWeights = layer.BiasWeights; //bias for each nueron. 
				for(int n=0; n < neurons.length ; n++) { //For each nueron, we multiply its set of weights.. 
					double op = biasWeights[n];
					double[] WE = neurons[n];								//Nueron Weights
					for(int w=0; w<WE.length;w++) op+=WE[w]*IO[w]; 	//input*Nuerons
					Inj[n] = op;		 	//Store neurons Inj
					nextAJ[n] = AF(op);	//Store nuerons next AJ. 
				}
				if(V==4) printer.printHidden(l+2, Inj, nextAJ);
			}
		}
	
	 private void bwdProp(int nsr, int actaulY, Printer printer) {
		 //This section does last bit of foward prop, gradient calculation, and weight update storage. I do this to avoid extra storage. 
			Layer OPL = NetWork[NwLen-1]; 				//Output layer. 
			double[] biasWeights = OPL.BiasWeights;	//bias for each nueron. 
			double[] NUpdates = OPL.NUpdates[nsr]; 	//Each wieght update for whole layer.  
			double[] bUpdates = OPL.BUpdates[nsr];    //Tmp storage for BUpdates
			double[][] neurons = OPL.Neurons; 			//Each nuerons set of weights. 
			double[] AJ = OPL.Aj[nsr];				 		//Output of prev layer/or data. 
			double[] inj = OPL.Inj[nsr];			 		//Storage for sum of weights*input, and then storage for gradient. 
			int wLen = neurons[0].length; 			   //The amount of weights corresponds to the amount of Aj's. 
			for(int n=0, u=0; n < neurons.length ; n++) { //I grabb all neurons, and their incoming weights
				double y = n==actaulY? 1: 0; 						 //The y value holds the index of the nueron that should output 1. 
				double in = biasWeights[n]; 						 //Output of nueron
				double[] WE = neurons[n];							 //Nueron Weights
				for(int w=0;w<WE.length;w++) in+=WE[w]*AJ[w]; //Mupltiy input*Nuerons
				double aj = 1 / (1 + Math.exp(-in));			 //Output, value. 
				double gradient = -2*(y-aj)*AFP(in);	 		 //cacluate grd 
				inj[n] = gradient; 	 //Inj maps to one neuron, thus I use it to map to one gradient. 
				bUpdates[n]+=gradient; 
				for(int i=0; i < wLen; i++, u++) NUpdates[u]+=gradient*AJ[i]; //For each output. wij = ^j*ai
				if(V==4) printer.addOutput(in, aj, y, gradient);
			}
			
		//BackProp part 2: This edits last hidden layer to inputlayer, and accesses ouput layer to inpuut layer. My implementation, each nueron is a set of incoming weights. NOT outgoing weightss. 
	   //I need access to 2 layerrs. The layer i holds gradients( now held in its inj), that infleunce layer i-1, and layer i holds the "outgoing" weights that layer i-1 should be connected too. 
			for(int l=NwLen-1; l > 1; l--) { 
				double[][] prvNeurons = NetWork[l].Neurons; //Layer i holds layer i-1 set of outgoing weights. This needed for gradient calculations. 
				double[] prvGrd = NetWork[l].Inj[nsr];  	  //Layer i holds the previous gradients
				wLen = prvNeurons.length;
				Layer layer = NetWork[l-1];  	//The Layer we aim to caculate gradients for. 
				AJ = layer.Aj[nsr]; 			  	//Values needed for calculaute this layers weight updates. I can do this because this layer represents a set of incoming weightts. 
				inj = layer.Inj[nsr];		  	//the input values that were used to create layer i's, aj, which are needed to caclulaute layers i-1 gradient. 
				NUpdates = layer.NUpdates[nsr];
				bUpdates = layer.BUpdates[nsr];
				for(int g=0, u=0; g< inj.length; g++) { //inj becomes new gradient storage, each inj also maps tto postion in wieghts. 
					double grd = 0;  					  //Represents tmp sum for gradient. g'(inj) * sigma wjj`^j`	
					for(int n =0; n< wLen;n++) grd+= prvNeurons[n][g]*prvGrd[n]; 
					grd*=AFP(inj[g]); //^i = g`(in(layer i-1)) * layer(i)(prvneuron[0][g]*prvGrd[g]+prvnueron[1][g]*grd[g]				
					inj[g] = grd;    //then we save the gradient. 
					bUpdates[g]+=grd; 
					for(int i=0; i < AJ.length; i++, u++) NUpdates[u] += grd*AJ[i]; //For each output. wij = ^j*ai
				}
				if(V==4) printer.addGradients(inj, l+1);
			}
			if(NwLen>1) {
				double[][] prvNeurons = NetWork[1].Neurons; //Layer i holds layer i-1 set of outgoing weights. This needed for gradient calculations. 
				double[] prvGrd = NetWork[1].Inj[nsr];  	  //Layer i holds the previous gradients
				wLen = prvNeurons.length;
				Layer layer = NetWork[0];  	//The Layer we aim to caculate gradients for. 
				AJ = layer.Aj[nsr]; 			  	//Values needed for calculaute this layers weight updates. I can do this because this layer represents a set of incoming weightts. 
				inj = layer.Inj[nsr];		  	//the input values that were used to create layer i's, aj, which are needed to caclulaute layers i-1 gradient. 
				NUpdates = layer.NUpdates[nsr];
				bUpdates = layer.BUpdates[nsr];
				for(int g=0, u=0; g< inj.length; g++) { //inj becomes new gradient storage, each inj also maps tto postion in wieghts. 
					double grd = 0;  					  //Represents tmp sum for gradient. g'(inj) * sigma wjj`^j`	
					for(int n =0; n< wLen;n++) grd+= prvNeurons[n][g]*prvGrd[n]; 
					grd*=AFP(inj[g]); //^i = g`(in(layer i-1)) * layer(i)(prvneuron[0][g]*prvGrd[g]+prvnueron[1][g]*grd[g]				
					inj[g] = grd;    //then we save the gradient. 
					bUpdates[g]+=grd;
					for(int i=0; i < AJ.length-1; i++, u++) NUpdates[u] += grd*AJ[i]; //For each output. wij = ^j*ai
				}
				if(V==4) printer.addGradients(inj, 2);
			}
			if(V==4) printer.printAll();
		}
	 
 //===========================================================================================
 //Weight Update mehtods
	 private void linearUpdateWeights(double alpha, double lambda, int NwLen, int btchSz){
		for(int l=0; l < NwLen; l++) {
			Layer layer = NetWork[l];
			double[][] nuerons = layer.Neurons; //grabs the nueron for each layer
			double[] sums = layer.NUpdates[0]; //I use this postion, to sum update all tmp wij=ai*gradientj
			double[] biasWeights = layer.BiasWeights;		//bias for each nueron. 
			double[] grds = layer.BUpdates[0];
			for(int n=0, uI=0; n < nuerons.length; n++) { //Here I do the actaul weight updates. 
				double[] weights = nuerons[n];				 //I access each nuerons weights. 
				for(int w=0; w < weights.length; w++, uI++) { //Then grab their corresponding weight update and do the math. 					
					weights[w] = weights[w] -((alpha*(sums[uI]/ btchSz))  - 2*alpha*lambda*weights[w]);
					sums[uI] = 0; 
				}
				biasWeights[n] = biasWeights[n] - ((alpha*( grds[n]/ btchSz)) - 2*alpha*lambda*biasWeights[n]);
				grds[n] = 0;
			}
		
		}
	 }
	 
	 private void pUpdateWeights(double alpha, double lambda, int NwLen, int btchSz){
		for(int l=0; l < NwLen; l++) {
			Layer layer = NetWork[l];
			double[][] nuerons = layer.Neurons; 		//grabs the nueron for each layer
			double[] biasWeights = layer.BiasWeights;	//bias for each nueron.
			
			double[][] updates = layer.NUpdates; //if 2 weights, N[1][1]=updates[1], N[2][1] = updates[3]
			double[] sums = layer.NUpdates[0]; //I use this postion, to sum update all tmp wij=ai*gradientj
			double updatesLen = sums.length;
			
			double[][] grdPointer = layer.Inj; 
			double[] grds = grdPointer[0];     //Use this as the same postion to sum all gradients. 
			double grdLen = grds.length; 
		
			for(int nsr=1; nsr < btchSz; nsr++) { //Could mabye optimization here by different access patterns. 				
				double[] tmpU = updates[nsr]; 	//nsr represents the access positon. for a certain chain of inputs
				double[] tmpG = grdPointer[nsr];	
				for(int u=0; u < updatesLen; u++) sums[u]+=tmpU[u]; //WE START AT ONE, because we use sum as our base. Each weight coresponds to a postion wihtin the updates row. 
				for(int g=0; g<grdLen; g++) grds[g] += tmpG[g];
			}
			
			for(int n=0, uI=0; n < nuerons.length; n++) { //Here I do the actaul weight updates. 
				double[] weights = nuerons[n];				 //I access each nuerons weights. 
				for(int w=0; w < weights.length; w++, uI++) { //Then grab their corresponding weight update and do the math. 
					weights[w] = weights[w] -((alpha*(sums[uI]/ btchSz))  - 2*alpha*lambda*weights[w]);
					sums[uI] = 0; 
				}
				biasWeights[n] = biasWeights[n] - (alpha*( grds[n]/ btchSz)) - 2*alpha*lambda*biasWeights[n];
				grds[n] = 0; 
			}	
		}
	 }
		
 //===========================================================================================
 //Calculation Methods 
 //My caculation methods are specially designed to do early stop if i dont have to expclity log 
 //Furthermore they are different than foward prop methtodds as I dont need to store extra values.
 //Thus making these methods expoetinally less expesnive. 
 	
  //=================================================
	//These calculate the loss, cost, and accuarcy during gradient descent.  
	public boolean lCalc(double[][] RegSet, int[] TrainSet, boolean log, String print){	
		 int sz = TrainSet[TrainSet.length-1];	
		 float[] loss = lCalcLoss(RegSet, 0, sz,  0, RegSet[0].length-1, log); //I pass in log, as a way to evaulute if we can consider stopping doing cacluations
		 boolean stop = loss[1] <= .01? true: false;  //This tells us if need to stop doing epochs. 
		 if(V>2 && log) {
			 float acc =   loss[2]/sz;
			 double cost = (loss[0]/(double)sz) + Lambda*regularize();
			 System.out.printf("%s Cost = %.6f; Loss = %.6f; Acc = %.4f%n", print, cost, loss[0]/sz, acc);		
		 }
		 return stop; 
	 }
	
	public boolean pCalc(double[][] RegSet, int[] TrainSet, int[] FastCalc, final boolean log, String print){
		int sz = TrainSet[TrainSet.length-1];	
		int valEnd = RegSet.length; 
		float[] loss = pCalcLoss(RegSet,  FastCalc, RegSet[0].length-1, valEnd-sz, log);
		boolean stop = loss[1] <= .01? true: false;  //This tells us if need to stop doing epochs. 
		if(V>2 && log) {
			float acc =   loss[2]/sz;
			double cost = (loss[0]/(double)sz) + Lambda*regularize();
			System.out.printf("%s Cost = %.6f; Loss = %.6f; Acc = %.4f%n", print, cost, loss[0]/sz, acc);		
		}
		return stop; 
	}
	
	//=================================================
	//These calculate the accuarcy after the stopping condition is met.   	
	public void lAcc(double[][] RegSet,int[] TrainSet){
		int sz = TrainSet[TrainSet.length-1];	
		int valEnd = RegSet.length; 
		float train = lCalcLoss(RegSet, 0, sz,  0, RegSet[0].length-1, true)[2]/sz;
		float val = lCalcLoss(RegSet, sz, valEnd,  0, RegSet[0].length-1, true)[2]/(valEnd-sz);
		System.out.printf("  TrainAcc: %.6f\n  ValidAcc: %.6f\n", train, val);
	}

	public void pAcc(double[][] RegSet, int[] TrainSet, int[] FastCalc){
		int sz = TrainSet[TrainSet.length-1];	
		int valEnd = RegSet.length; 
		float train = pCalcLoss(RegSet,  FastCalc, RegSet[0].length-1, valEnd-sz, false)[2]/sz;
		float val = pCalcLossVal(RegSet, sz, valEnd, RegSet[0].length-1, false)/(valEnd-sz);
		System.out.printf("  TrainAcc: %.6f\n  ValidAcc: %.6f\n", train, val);
		
	}

 //=================================================
 //These are the methods lCalc,pCalc, lAcc, and pAcc use to calculate.  			
	private float[] lCalcLoss(double[][] RegSet, int strt, int end, int nsr, int yI, final boolean log){
		float loss = 0, stop = 0, acc=0; //stop tells us the largest different of nueron output. 
		boolean earlyStop = false;    
		for(int row=strt; row < end; row++) {
			NetWork[0].Aj[nsr] = RegSet[row];  //First incoming Input. 
			int actaulY = (int) RegSet[row][yI]; 
			//This simply does foward prop till last layer. 
			for(int l=0; l < NwLen-1; l++) {			//Foward prop till last layer.  
				Layer layer = NetWork[l];				//Layer hold set of nuerons.
				double[][] neurons = layer.Neurons; //Nuerons holds set of weights. 
				double[] Aj = layer.Aj[nsr];		   	 //OP of prev layer, each output maps to each i in wij. 
				double[] biasWeights = layer.BiasWeights;		//bias for each nueron. 
				double[] nextAJ = NetWork[l+1].Aj[nsr]; //Next Layer storage. 
				for(int n=0; n < neurons.length ; n++) {//For each nueron, we multiply its set of weights.. 
					double op = biasWeights[n];
					double[] WE = neurons[n];								//Nueron Weights
					for(int w=0; w<WE.length;w++) op+=WE[w]*Aj[w]; 	//input*Nuerons
					nextAJ[n] = AF(op);	//Store nuerons next AJ. 										
				}
			}	
			int index = -1; //Location of largest nueron val, and indicator if result is correct. 
			double largest = Double.MIN_VALUE;
			Layer OPL = NetWork[NwLen-1];
			double[] Aj = OPL.Aj[nsr];
			double[][] neurons = OPL.Neurons;
			double[] biasWeights = OPL.BiasWeights;		
			
			outer:
			for(int n=0; n < neurons.length; n++) {
				int y = n==actaulY? 1: 0;
				double[] WE = neurons[n];								
				double inj = biasWeights[n]; 
				for(int w=0; w<WE.length;w++) inj+=WE[w]*Aj[w];
				inj = 1 / (1 + Math.exp(-inj)); //inj now becomes aj. 
				if (inj > largest) {
					index = n;
					largest = inj;
				}
				float abs =(float) Math.abs(y-inj);
				if(abs>stop) { 
					stop = abs;
					if(!log && abs > .01) { 
						earlyStop=true; 
						break outer; 
					}
				} 
				loss+=Math.pow(y-inj, 2);
			}
			if(index == actaulY) acc++; 
			if(earlyStop) loss = .11f;
		}	
		return new float[] {loss, stop, acc};
	}
			
	private float[] pCalcLoss(double[][] RegSet, int[] FastCalc, int yI, int btchSz, final boolean log){		
		AtomicReference<Double> sync0 = new AtomicReference<>(0.0); // loss accumulator 
		AtomicReference<Double> sync1 = new AtomicReference<>(0.0); // max stop value
		AtomicInteger sync2 = new AtomicInteger(0);  // count of correct predictions
		AtomicBoolean earlyStop = new AtomicBoolean(false);
		
		 for (int fc = 0; fc < FastCalc.length - 1; fc++) {
		     final int safeFc = fc;
		     int strt = FastCalc[fc], end = FastCalc[fc + 1];
		     IntStream.range(strt, end).parallel().forEach(row -> {
		   	  if(!earlyStop.get()) {  
		   	 	int nsr = row >= btchSz ? row - (btchSz * safeFc) : row;
		   	 	if(nsr==632) {
		   	 		int a = 0;
		   	 	}
		   	 	
		         NetWork[0].Aj[nsr] = RegSet[row]; // First incoming Input.
		         int actaulY = (int) RegSet[row][yI]; 
		
		         // Forward propagation through the layers
		         for (int l = 0; l < NwLen - 1; l++) {
		             Layer layer = NetWork[l];
		             double[][] neurons = layer.Neurons;
		             double[] Aj = layer.Aj[nsr];
		             double[] biasWeights = layer.BiasWeights;
		             double[] nextAJ = NetWork[l + 1].Aj[nsr];
		             for (int n = 0; n < neurons.length; n++) {
		                 double op = biasWeights[n];
		                 double[] WE = neurons[n];
		                 for (int w = 0; w < WE.length; w++) op += WE[w] * Aj[w];  // input * neurons
		                 nextAJ[n] = AF(op);  // Store neurons' next AJ.
		             }
		         }
		         double loss = 0;
		         int index = -1; // Location of largest neuron value, and indicator if result is correct.
		         double stop = 0;
		         double largest = Double.MIN_VALUE;
		         Layer OPL = NetWork[NwLen - 1];
		         double[] Aj = OPL.Aj[nsr];
		         double[][] neurons = OPL.Neurons;
		         double[] biasWeights = OPL.BiasWeights;
		         outer: 
		         for (int n = 0; n < neurons.length; n++) {
		             double y = n == actaulY ? 1 : 0;
		             double[] WE = neurons[n];
		             double inj = biasWeights[n];
		             for (int w = 0; w < WE.length; w++) inj += WE[w] * Aj[w];
		             inj = 1 / (1 + Math.exp(-inj)); // inj now becomes aj.
		             if (inj > largest) {
		                 index = n;
		                 largest = inj;
		             }
		             float abs = (float) Math.abs(y - inj);
		             if(abs>stop) stop = abs;
		             if(abs>stop) { 
		 					stop = abs;
		 					if(!log && abs > .01) { 
		 						earlyStop.set(true);
		 						break outer; 
		 					}
		 				} 
		            loss += Math.pow(y - inj, 2); // Add loss
		         }
		         final double l = loss, s = stop;  
		         sync0.updateAndGet(prev -> prev + l); // Add loss atomically
		         sync1.updateAndGet(prev -> Math.max(prev, s)); // Update max stop value atomically
		         if (index == actaulY) sync2.incrementAndGet(); // Increment count of correct predictions atomicall
		   	  }
	        });
		    }
		    return new float[] {sync0.get().floatValue(), sync1.get().floatValue(), sync2.get()};
	}
	
	//Special calc for loss.
	private float pCalcLossVal(double[][] RegSet, int strt, int end, int yI, final boolean log){		
		int btchSz = end-strt;
		int tmp = end-btchSz; 
		AtomicInteger sync2 = new AtomicInteger(0);  // count of correct predictions
		     IntStream.range(strt, end).parallel().forEach(row -> {
		   	 	int nsr = row >= btchSz ? row - tmp : row;
		         NetWork[0].Aj[nsr] = RegSet[row]; // First incoming Input.
		         int actaulY = (int) RegSet[row][yI]; 
		
		         // Forward propagation through the layers
		         for (int l = 0; l < NwLen - 1; l++) {
		             Layer layer = NetWork[l];
		             double[][] neurons = layer.Neurons;
		             double[] Aj = layer.Aj[nsr];
		             double[] biasWeights = layer.BiasWeights;
		             double[] nextAJ = NetWork[l + 1].Aj[nsr];
		             for (int n = 0; n < neurons.length; n++) {
		                 double op = biasWeights[n];
		                 double[] WE = neurons[n];
		                 for (int w = 0; w < WE.length; w++) op += WE[w] * Aj[w];  // input * neurons
		                 nextAJ[n] = AF(op);  // Store neurons' next AJ.
		             }
		         }
		         double loss = 0;
		         int index = -1; // Location of largest neuron value, and indicator if result is correct.
		         double stop = 0;
		         double largest = Double.MIN_VALUE;
		         Layer OPL = NetWork[NwLen - 1];
		         double[] Aj = OPL.Aj[nsr];
		         double[][] neurons = OPL.Neurons;
		         double[] biasWeights = OPL.BiasWeights;
		         for (int n = 0; n < neurons.length; n++) {
		             double y = n == actaulY ? 1 : 0;
		             double[] WE = neurons[n];
		             double inj = biasWeights[n];
		             for (int w = 0; w < WE.length; w++) inj += WE[w] * Aj[w];
		             inj = 1 / (1 + Math.exp(-inj)); // inj now becomes aj.
		             if (inj > largest) {
		                 index = n;
		                 largest = inj;
		             }
		         }
		         if (index == actaulY) sync2.incrementAndGet(); // Increment count of correct predictions atomicall
	        });
		    return  sync2.get();
	}

	private double regularize() {
		double sum = 0; 
		for(Layer layer: NetWork) {
			for(double[] weights: layer.Neurons) {
				for(double weight: weights) {
					sum+=Math.pow(weight, 2);
				}
			}
		}
		return sum;
	}
 
 //================================================================================================
 //Activation Functions
	
	private double AF(double t){
		return switch(AF) {
			case 0 -> 1 / (1 + Math.exp(-t)); 		//logistic
			case 1 -> Math.log(1 + Math.exp(-t));	//softplus
			case 2 -> Math.max(0, t);			 		//Relu
			default -> {									//tanh
				double eat = Math.exp(t); 	//postive t
				double est = Math.exp(-t); //negative t
				yield (eat-est)/(eat+est);
			}
		};
	}

	private double AFP(double t) {
		return switch(AF) {
			case 0 -> { //logistic	
				double funct = 1 / (1 + Math.exp(-t));
				yield funct*(1.0-funct);
			}
			case 1 -> { //softplus
				yield 1 / (1 + Math.exp(-t));
			}
			case 2 -> { //Relu
				if(t>0) yield 1;
				else yield 0;	 		
			}
			default -> { //tanh								
				double eat = Math.exp(t); 	//postive t
				double est = Math.exp(-t); //negative t
				yield 1 - Math.pow((eat-est)/(eat+est), 2) ;
			}
		};
	}
	
 //================================================================================================
 //Special non print calc statements for HPT 
	public String lAcc(double[][] RegSet,int[] TrainSet, float[] results){
		int sz = TrainSet[TrainSet.length-1];	
		int valEnd = RegSet.length; 
		float train = lCalcLoss(RegSet, 0, sz,  0, RegSet[0].length-1, true)[2]/sz;
		float val = lCalcLoss(RegSet, sz, valEnd,  0, RegSet[0].length-1, true)[2]/(valEnd-sz);
		results[0] = train; results[1]=val;
		return String.format("\n  TrainAcc: %.6f  ValidAcc: %.6f\n", train, val);
	}

	public String pAcc(double[][] RegSet, int[] TrainSet, int[] FastCalc, float[] results){
		int sz = TrainSet[TrainSet.length-1];	
		int valEnd = RegSet.length; 
		float train = pCalcLoss(RegSet,  FastCalc, RegSet[0].length-1, valEnd-sz, false)[2]/sz;
		float val = pCalcLossVal(RegSet, sz, valEnd, RegSet[0].length-1, false)/(valEnd-sz);
		results[0] = train; results[1]=val;
		return String.format("\n  TrainAcc: %.6f  ValidAcc: %.6f\n", train, val);
	}

	
	
	
	
}

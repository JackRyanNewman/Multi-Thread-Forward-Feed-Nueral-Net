//Jack Newman
//Date: 2024-12-4: Assigment 3. 
import java.util.Random;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.IntStream;

import NeuralNet.Layer;






public class atomitcs {
	Layer[] NetWork; 
	int NwLen;
	byte AF;
	
	private float[] pCalcLoss(double[][] RegSet, int[] FstClc, int yI, final boolean log) {
		 AtomicReference<Double> sync0 = new AtomicReference<>(0.0); // loss accumulator
		 AtomicReference<Double> sync1 = new AtomicReference<>(0.0); // max stop value
		 AtomicInteger sync2 = new AtomicInteger(0);  // count of correct predictions
		 AtomicBoolean earlyStop = new AtomicBoolean(false);
		 for (int fc = 0; fc < FstClc.length - 1; fc++) {
		     final int safeFc = fc;
		     int strt = FstClc[fc], end = FstClc[fc + 1], btchSz = end - strt;
		     IntStream.range(strt, end).parallel().forEach(row -> {
		   	  if (!earlyStop.get()) {  
		   	  		int nsr = row >= btchSz ? row - (btchSz * safeFc) : row;
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
	
	
	private float[] pCalcLoss(double[][] RegSet, int[] FstClc, int yI) {
		 AtomicReference<Double> sync0 = new AtomicReference<>(0.0); // loss accumulator
		 AtomicReference<Double> sync1 = new AtomicReference<>(0.0); // max stop value
		 AtomicInteger sync2 = new AtomicInteger(0);  // count of correct predictions
	
		 for (int fc = 0; fc < FstClc.length - 1; fc++) {
		     final int safeFc = fc;
		     int strt = FstClc[fc], end = FstClc[fc + 1], btchSz = end - strt;
		     IntStream.range(strt, end).parallel().forEach(row -> {
	   	  		int nsr = row >= btchSz ? row - (btchSz * safeFc) : row;
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
		             float abs = (float) Math.abs(y - inj);
		             if(abs>stop) stop = abs;
		 					
		 			
		             loss += Math.pow(y - inj, 2); // Add loss
		         }
		         final double l = loss, s = stop;  
		         sync0.updateAndGet(prev -> prev + l); // Add loss atomically
		         sync1.updateAndGet(prev -> Math.max(prev, s)); // Update max stop value atomically
		         if (index == actaulY) sync2.incrementAndGet(); // Increment count of correct predictions atomicall
	        });
		    }
		    return new float[] {sync0.get().floatValue(), sync1.get().floatValue(), sync2.get()};
		}
	
	private float[] PCalcLoss(double[][] RegSet, int[] FstClc, int yI, final boolean log){		
		float[] sync = {0,0,0};  
		
		
		for(int fc=0; fc< FstClc.length-1; fc++) {
			final int safeFc = fc; 
			int strt = FstClc[fc], end = FstClc[fc+1], btchSz = end-strt;
			IntStream.range(strt, end).parallel().forEach(row -> {
				int nsr = row >= btchSz? row-(btchSz*safeFc): row;
				NetWork[0].Aj[nsr] = RegSet[row];  //First incoming Input. 
				int actaulY = (int) RegSet[row][yI]; 
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
				double loss = 0; 
				int index = -1; //Location of largest nueron val, and indicator if result is correct. 
				double stop = 0; 
				double largest = Double.MIN_VALUE;
				Layer OPL = NetWork[NwLen-1];
				double[] Aj = OPL.Aj[nsr];
				double[][] neurons = OPL.Neurons;
				double[] biasWeights = OPL.BiasWeights;		
				for(int n=0; n < neurons.length; n++) {
					double y = n==actaulY? 1: 0;
					double[] WE = neurons[n];								
					double inj = biasWeights[n]; 
					for(int w=0; w<WE.length;w++) inj+=WE[w]*Aj[w];
					inj = 1 / (1 + Math.exp(-inj)); //inj now becomes aj. 
					if (inj > largest) {
						index = n;
						largest = inj;
					}
					float abs =(float) Math.abs(y-inj);
					if(abs>stop) stop = abs;
					loss+=Math.pow(y-inj, 2);
				}
				synchronized(sync){
					sync[0]+=loss; 
					if(stop>sync[1]) sync[1] = (float) stop; 
					if(index == actaulY) sync[2]++; 
				}
			});
			
			
		}
		return sync;
	}
	
	
	
	
	private class Layer{
		final double[] BiasWeights;//Each bias weight, corresponds to one nueron. 
		final double[][] Neurons; 	//Neuron Storage to storage currTotal
		final double[][] Aj; 		//incoming IO*wieghts to get Inj.  
		final double[][] Inj; 		//Result of Aj*weights, or weight update value.  
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
	
	
	
	
}

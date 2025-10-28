//Jack Newman
//Date: 2024-12-4: Assigment 3. 
//==================================================================
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;

//HP:
//Description: This class was designed to generate different dymanic permuntations. The whole generation is based off the stattic varibles 
//you set in defaults. Essestinally it allows for rapid testing, and its called by my main method in NueralNetController. 

public class HP { //This class'es contents are NWCTyye meant to be collasped. 
	final defaults D;
	private FileOutputStream fileOutputStream;
	
	public HP() throws Exception {
		fileOutputStream = new FileOutputStream("../457-ML-03-JN/a03-data/output.txt", false);
		D = new defaults();
	}
	
	//This returns arraylist<String[]> which is arraylist that holds all arguments createdd. 
	public ArrayList<String[]> createTrials(int type, boolean all){
		ArrayList<String> netWorks = createNetWorks(type, all);
		ArrayList<String[]> configs = new ArrayList<>();
		for (String netWork : netWorks) {
			 for (int epoch : D.EpochLimit) {
			  for (double alpha : D.Alpha) {
			   for (double lambda : D.Lambda) {
			    for (double InitWEight : D.InitWE) {
			     for (int batch : D.BtchSz) {
			      String tmp = "-e " + epoch + " -a " + alpha + " -l " + lambda + " -w " + InitWEight + " -% " + batch + " -h" + " " + netWork;
			      configs.add(tmp.split(" "));
			     }
			    }
			   }
			  }
			 }
			}
	  return configs; 
	}
	
	public ArrayList<String> createNetWorks(int type, boolean all) {
		ArrayList<String> netWorks = new ArrayList<>(); 
		for (int t = 0; t < D.NWCType.length; t++) {
		  if (all || type == t) {
		    for (Deque<String> set : D.NWCType[t]) {
		      for (int size : D.nwLens) {
		        String net = "" + size;
		        for (int h = 0; h < size; h++) {
		          String val = set.removeFirst();
		          net += " " + val;
		          set.addLast(val);
		        }
		        netWorks.add(net);
		      }
		    }
		  }
		}
		String TNK = ""+D.tSz; //TNWCTyye and Skinny ttmp string. 
		 for(int t=0; t < D.tSz; t++) TNK+=" "+D.TNWCTyye;
		 netWorks.add(TNK);
		 TNK = ""+D.sSz;
		 for(int t=0; t < D.sSz; t++) TNK+=" "+D.Skinny;
		 netWorks.add(TNK);
		 return netWorks;
	}
	

//==================================================================
//Helper classes for storage. 	
	
	public class defaults{
		final int[] EpochLimit = {500,1000};  	 // -e <INTEGER> Epoch limit for gradient descent; default 1000: 
		final double[] Alpha =   {.005, 0.01};   // -a <DOUBLE>	 Learning rate for gradient descent; default 0.01
		final double[] Lambda =  {0,.005, 0.01, .1}; // -l <DOUBLE>		Regularization parameter; default 0.0
		final double[] InitWE =  {0.05,.1,};	 // -w <DOUBLE>  Weight initialization parameter; default 0.1   
		final int[] BtchSz =     {0, 50};    	 	 // -% <byte> 	 !CUSTOM FLAG!, btchSize will now be determined by, (bseSetLen*.8) * bScalar
	  final int TNWCTyye = 6, tSz = 2; 
	  final int Skinny = 4, sSz = 3; 
	  final int[] nwLens = new int[] {3, 5, 8};
	  @SuppressWarnings({ "unchecked", "serial" })
		final ArrayList<Deque<String>>[] NWCType = new ArrayList[] {
//        new ArrayList<Deque<String>>() {{ //Small network configurations
//            add(new ArrayDeque<String>() {{ add("2"); add("3"); }}); //even
//            add(new ArrayDeque<String>() {{ add("3"); add("3"); }}); //Odd
//        }},
//        new ArrayList<Deque<String>>() {{ //Medium network configurations
//          add(new ArrayDeque<String>() {{ add("4"); add("5"); }});
//          add(new ArrayDeque<String>() {{ add("5"); add("5"); }});
//      }},
      new ArrayList<Deque<String>>() {{ //Large network configurations.
          add(new ArrayDeque<String>() {{ add("5"); add("15"); }}); 
          add(new ArrayDeque<String>() {{ add("50"); add("25");  add("20"); }});
      }},
      new ArrayList<Deque<String>>() {{
        add(new ArrayDeque<String>() {{ add("100"); add("50"); add("25"); add("10"); add("5"); }});
        add(new ArrayDeque<String>() {{ add("20"); add("15"); add("10"); }});
      }},
    };
	  final String[] ActivationFuncs = { 		 //-g <String>	 activation functions     
				"logistic", 
				//"softplus", 
				//"relu", 
				//"tanh" 
	  }; 	  
		defaults(){}
	}
	
 //=================================================================
 //Way to observe and track Configs. 
	
    public void openDualStream() throws Exception {
        PrintStream originalOut = System.out;
        OutputStream dualStream = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                fileOutputStream.write(b);       // Write to file
                originalOut.write(b);            // Write to console
            }
        };
        System.setOut(new PrintStream(dualStream));
    }
    
    public void closeDualStream() throws IOException { 
       fileOutputStream.close();  // Close the file output stream
    }
}

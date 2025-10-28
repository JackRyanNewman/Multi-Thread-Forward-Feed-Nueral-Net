//Jack Newman
//Date: 2024-12-4: Assigment 3. 
import java.util.function.Function;

public class Lambdas {
	final Function<Double, Double> ACTF = Lambdas::logistic; 
	final Function<Double, Double> ACTPF = Lambdas::logisticPrime;
	
	 //================================================================================================
	  //Activation Functions
		 public static double logistic(double t) {
			 return 1 / (1 + Math.exp(-t));
		 }
		 public static double logisticPrime(double t) {
			 double funct = 1 / (1 + Math.exp(-t));
			 return funct*(1-funct);
		 }
	
	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

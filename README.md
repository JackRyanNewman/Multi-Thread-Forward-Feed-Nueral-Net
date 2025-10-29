# Multi-Class Classification Neural Network with Threaded Hyperparameter Tester (Java)

**Link:** [Project Repository]  

This project implements a fully customizable **multi-layer feed-forward neural network** for **multi-class classification**, written entirely in Java, with a focus on performance, flexibility, and real-world experimentation. It features a **threaded hyperparameter tester** that allows parallel exploration of parameter combinations, accelerating model development and tuning.  

The network supports **configurable architecture**, including multiple hidden layers, diverse activation functions (logistic, ReLU, softplus, tanh), and mini-batch gradient descent with optional **threaded mini-batch processing** for multi-core acceleration. Designed with **data-oriented principles** and **JVM profiling research**, this project emphasizes both performance and reproducibility.  

## Key Features

### Core Neural Network
- Fully custom multi-layer perceptron with user-defined hidden layers.
- Mini-batch gradient descent with stochastic, mini-batch, and full-batch options.
- Forward and backward propagation with detailed logging capabilities.
- L2 regularization to reduce overfitting.
- Feature scaling via min-max normalization learned from training data.  

### Advanced Optimization
- Threaded hyperparameter tester for parallel experimentation.
- Optional threaded mini-batch processing.
- Deterministic weight initialization for reproducible debugging.
- Flexible stopping conditions (epoch limits, error thresholds).
- Performance profiling with iteration timing metrics.  

### Experimental Framework
- Automated train/validation splitting (80/20) with optional randomization.
- Adjustable verbosity levels from silent to detailed propagation traces.
- Designed for experimentation on complex datasets such as MNIST and image classification.
- Systematic tracking of results and best-performing configurations.  

## Custom Flags

The program includes a variety of configurable flags for precise control over network behavior and hyperparameter testing:  

- `-f <FILENAME>`: Reads input data from the specified file.  
- `-h <NH> <S1> <S2> ...`: Defines the number of hidden layers `<NH>` and their respective sizes `<S1>, <S2>, ...` (0 ≤ NH ≤ 10, 1 ≤ Si ≤ 500).  
- `-a <DOUBLE>`: Learning rate α for mini-batch gradient descent (default: 0.01).  
- `-e <INTEGER>`: Epoch limit (default: 1000).  
- `-m <INTEGER>`: Batch size (default: 1 for stochastic gradient descent; 0 indicates full-batch).  
- `-l <DOUBLE>`: L2 regularization hyperparameter λ (default: 0.0).  
- `-r`: Enables randomization of data for train/validation split and batch construction. Default is off.  
- `-w <DOUBLE>`: Weight initialization value ϵ (default: 0.1).  
- `-v <INTEGER>`: Verbosity level; 0 = near silent, 1 = default, higher values for detailed logging.  
- `-p <BOOLEAN>`: Enables optional threaded
- `% <byte>`: batch size divisor. batchSize= (dataSet*.8)/ batchSzDiv. -1 is default, which is off, 0 is full batch.


##READ 00-References\0-assignment-03.pdf For more info
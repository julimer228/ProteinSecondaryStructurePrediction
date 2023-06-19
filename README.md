# ProteinSecondaryStructurePrediction
The protein secondary structure prediction is an essential problem in bioinformatics. The structure mostly depends on the primary amino acid sequence of the protein. Secondary structure prediction belongs to the group of pattern recognition and classification problems. The secondary structure of a given instance is predicted based on its sequence features. One of the known solutions is using the Support Vector Machine (SVM) to predict the secondary structure, which has been described in [1] and [2]. The aim of this work was the implementation of the protein secondary structure predictor based on a logistic regression model. To do this we implemented the algorithm described in the mentioned articles. The project was implemented in R programming language.
<h1>Methods</h1>

<h2>Dataset</h2>
The dataset consists of three text files: the training dataset, the testing dataset and the validation dataset. Each file has the following structure: the first line provides information about the sequence identification code. The second line contains the sequence of amino acids. In the third line, the secondary protein structure is written. The proteins are separated with an empty line. There are no missing values in this dataset.
<h3>Working with the large dataset</h3>
Due to the large size of the dataset, we decided to use R libraries that allowed us to perform calculations on multiple cores: parallel and doParallel. In addition, we saved the trained binary classifiers in RDS files. This operation allowed us to remove models from the workspace and free the memory, which was very important for performing further calculations. We used google drive to store the models, which also allowed us to transfer data between two computers.
<h2>Measures</h2>
To evaluate achieved results we used two commonly used measures Q3 and SOV.
<h3>Q3</h3>
The secondary structure prediction is usually evaluated by Q3 accuracy, which measures the percent of residues for which a 3-state secondary structure is correctly predicted. 
<h3>SOV</h3>
The segment overlap score (SOV) is used to evaluate the predicted protein secondary structures, a sequence composed of helix (H), strand (E), and coil (C), by comparing it with the native or reference secondary structures. The main advantage of SOV is that it can consider the size of continuous overlapping segments and assign extra allowance to longer continuous overlapping segments instead of only judging from the percentage of overlapping individual positions as Q3 score does. 

<h2>Logistic Regression</h2>
Logistic regression analyzes the relationship between multiple independent variables and a categorical dependent variable and estimates the probability of occurrence of an event by fitting data to a logistic curve. Binary logistic regression is commonly used when the outcome variable is binary and the predictor variables are either continuous or categorical.

<h2>Algorithm Description</h2>
After investigating the dataset, we created the input groups for the logistic regression classifier. We followed the instructions described in the articles. First, we had to implement the sliding window scheme. This method allows to preserve the information about the local interactions among neighbouring residues. In the beginning, we had to choose the size of the window. To predict the structure of the amino acid in the middle we need to use the sequence whose size is equal to the size of the window. As can be seen, the problem of missing amino acids at the ends of the sequence had to be solved. We complemented the missing values with the empty character "-" to keep the correct window size. The described algorithm for the window size 5 is presented in the image.
<h3>Sliding window coding scheme</h3>
After investigating the dataset, we created the input groups for the logistic regression classifier. We followed the instructions described in the articles. First, we had to implement the sliding window scheme. This method allows to preserve the information about the local interactions among neighbouring residues. In the beginning, we had to choose the size of the window. To predict the structure of the amino acid in the middle we need to use the sequence whose size is equal to the size of the window. As can be seen, the problem of missing amino acids at the ends of the sequence had to be solved. We complemented the missing values with the empty character "-" to keep the correct window size. The described algorithm for the window size 5 is presented in the image.

<p align="center">
  <img  width="500" src="https://github.com/julimer228/ProteinSecondaryStructurePrediction/assets/56163818/740640b6-8b88-43d9-8e27-4bbdda3d4125" alt="Sublime's custom image"/>
</p>
<h3>Orthogonal Input profile</h3>

The next step was to use orthogonal encoding to assign a unique binary vector to each residue. The weights of all the residues in the window have the value 1 and the rest have the value 0. As a result, we obtain the input matrix with 21 columns (there are 20 different amino acids and one value assigned to the empty character) and with a number of rows equal to the size of the window. Next, we reshaped the orthogonal input to get the one-dimensional input vector. The size of the vector is equal to the result of the multiplication of the number of rows of the matrix by the number of columns. The following rows were written to the vector one by one. Encoding for the input window of the size 5 is presented in the image.

<p align="center">
  <img  width="500" src="https://github.com/julimer228/ProteinSecondaryStructurePrediction/assets/56163818/dca757d4-a0e2-4ac0-af7c-f9662ca1f3f2" alt="Sublime's custom image"/>
</p>

<h3>Constructing the binary classifiers</h3>
We constructed six binary classifiers: three one-versus-one classifiers (H/E, C/E, C/H) and three one-versus-rest classifiers (C/~C, E/~E, H/~H). For each classifier, we trained the logistic regression model. 



<h3>Constructing tertiary classifier</h3>

The binary classifiers were used to create the different tertiary classifiers. We created three tree classifiers described in the articles. (C/~C & H/E,  E/~E & C/H, H/~H & C/E). For example for the second classifier when the first binary classifier classifies the sample as C its predicted value is C, otherwise the class is predicted by the second one-versus-one classifier H/E. Their structures are presented in figures. We also tested the classifier based on three one-versus-one classifiers (C/~C & E/~E & H/~H). The sample is assigned the class with the highest probability.

<p align="center">
  <img width="300" src="https://github.com/julimer228/ProteinSecondaryStructurePrediction/assets/56163818/b143e4ec-d107-4af3-9462-f5cd6766fb01" alt="Sublime's custom image"/>
</p>

<p align="center">
  <img width="300" src="https://github.com/julimer228/ProteinSecondaryStructurePrediction/assets/56163818/c006c433-37e7-4325-b1c7-bbd4fb4c53fc" alt="Sublime's custom image"/>
</p>
<p align="center">
  <img  width="300" src="https://github.com/julimer228/ProteinSecondaryStructurePrediction/assets/56163818/a7a29a10-aef0-4650-bb12-890bc878ec7b" alt="Sublime's custom image"/>
</p>

<h1>Results</h1>
We tested different window sizes (from 5 to 13 amino acids, only odd sizes). The table 3.1 presents the accuracy scores obtained for each binary classifier on the test dataset. The best results were obtained for the window of size 13.

<p align="center">
  <img  width="600" src="https://github.com/julimer228/ProteinSecondaryStructurePrediction/assets/56163818/a200ac25-8c38-44a8-bce5-abbf81bb5370" alt="Sublime's custom image"/>
</p>

Then we compared the Q3 and SOV results we got for each of the tertiary classifiers. To do that, we saved the predicted structure in a FASTA format. For each classifier, we used the window that provided the best accuracy for the binary classification. Results are presented in the table 3.2.

<h1>Bibliography</h1>
[1] Mayuri Patel and Hitesh Shah. ‘Protein Secondary Structure Prediction Using Support Vector Machines (SVMs)’. In: 2013 International Conference on Machine Intelligence and Research Advancement. 2013, pp. 594–598. doi: 10.1109/ICMIRA.2013.124.
[2] Hua Sujun, and Sun Zhirong. ‘A novel method of protein secondary structure prediction with high segment overlap measure: support vector machine approach’. In:Journal of molecular biology. 2001, 397––407.

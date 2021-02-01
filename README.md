# DTLDephos- Deep Transfer Learning based approach to predict dephosphorylation sites

DTLDephos is a transfer learning based approach employing deep learning to predict dephosphorylation sites of S,T and Y.

DTLDephos is a deep learning based method for Dephosphorylation sites of S,T and Y. It utilizes phosphorylation data through transfer learning method to address prediction to present scarce dephosphorylation data.It is implemented using Keras (version 2.2.4) and Tensorflow (version 1.15) backend and has been tested on both in Windows and Linux OS.
 


# Pre-requisites
  Python 3.6<br/>
  Tensorflow (>version 1.15) - (Optional: Tensorflow GPU for GPU based machine)<br/>
  Keras (>version 2.2.4) - (Optional: Keras GPU for GPU based machine)<br/>
  Numpy <br/>
  Biopython <br/>
  Sklearn <br/>
  Imblearn <br/>
  
 # Running on CPU or GPU
 To run in CPU, installation of Tensorflow and Keras will suffice. However, to run in GPU, further Tensorflow-gpu and keras-gpu must be installed. 
 Tensorflow GPU and Keras GPU version utilizes cuda cores in our GPU for faster training time. However, running in GPU is not mandatory.
 
 # Dataset
 Dataset is in fasta format. Both training and testing datasets are provided which are independent (one does not include others). Training dataset for positive and negative are X.fasta and train_X.fasta respectively. Testing dataset for positive and negative are test_X.fasta and test_X.fasta respectively. Training dataset is made available so that future models can be trained for the comparison purpose.
 
 # Model
 The model learned through transfer learning from the phosphorylation data, for residue ST and Y are provided. The ComDephos_ST.h5 and ComDephos_Y are the optimized model for ST and Y respectively.

# Code
Independent Testing code is provided. The model provided can be used to predict input dephosphorylation of S,T and Y sites for the given window sequences. The code will take input of window size 31.

# Prediction for given test dataset (Procedure)
  - Download test datasets, test_Pos_ST.fasta and test_Neg_ST.fasta(from dataset folder), and python code test_model.py.
    Keep them in the same folder as model files.
  - Run test_model.py and you will get output mentioned in our research.
  - In linux code will be, $python3 test_model.py
  
  # Prediction for your dataset
  If you would like to use DTLDephos to predict dephosphorylation sites in the protein of your interest, you should prepare your dataset in the same format as the test dataset which is basically a FASTA format. This model works for window size 31 only, meaning for the residue of your interest you should provide 25 resiudes downstream and 25 residues upstream. e.g. if you want to predict whether the target residue(S/T/Y) in Position 735 in protein Q4KWH8 is dephosphorylated or not, the input file should contain 16 residues upstream of the target residue (position 735 in protein Q4KWH8) and 16 residues downstream of the target residue.

The general format for your dataset should be:

sp|Q4KWH8|PLCH1_HUMAN%730%755
PKKQLILKVISGQQLPKPPDSMFGDSGEIIDPFVEVEIIGLPVDCCKDQTR

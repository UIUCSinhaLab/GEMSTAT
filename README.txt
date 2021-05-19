Thermodynamics-based models of transcriptional regulation by enhancers: the roles of synergistic activation, cooperative binding and short-range repression
Author: Xin He <xinhe2@illinois.edu>
--------------------------------------------------------------------------
INSTALLATION

The program needs GSL (GNU Scientific Library). After installing GSL, go to the main distribution directory, and type:
./configure
make

The main executable will be generated: 
src/seq2expr: the main program to fit a sequence-to-expression model

Type the program name without parameters will print usage information. 
--------------------------------------------------------------------------
PROGRAM USAGE

The program takes as input: sequences, the expression profiles of these sequences, the PWMs of the relevant TFs and the expression profiles of these TFs, and computes the paramaters of the underlying sequence-to-expression model as well as the predicted expression patterns.

The data/ directory contain the example files. A simple command of running the program: 
$INSTALL_DIR/src/seq2expr -s seqs.fa -e expr.tab -m factors.wtmx -f factor_expr.tab -fo obs_pre.txt -i factor_info.txt -o Direct

For more examples of running the program with various options, see the run.sh script (need to modify the installation directory of the script). 

Explanation of parameters: 
-s <seq_file>: required, the sequence file in FASTA format. See data/seqs.fa. 

-e <expr_file>: required, the expression data of the sequences. The first line specifies the name of the expression conditions. See data/expr.tab. 

-m <motif_file>: required, the PWM (motif) of the relevant TFs. See data/factors.wtmx. 

-f <factor_expr_file>: required, the expression data of the TFs. Must match the format of expr_file, and the order of TFs in motif_file. See data/factor_expr.tab. 

-fo <output_file>: required, the output file, the predicted expression patterns of all sequences as well as the observed expression patterns (alternating, the first row is the observed and the second the predicted). 

-o <model_option>: the sequence-to-expression model. Options are: Logistic, Direct (DirectInt model), ChrMod_Unlimited (SRR model with N_MA = inf), ChrMod_Limited (SRR model with finite N_MA). 

-c <coop_file>: the list of cooperative interactions. One line per cooperative pair. If not specified, then no cooperative interaction is allowed. See data/coop_file. 

-i <factor_info_file>: the role of TFs (activators or repressors). This would be required for SRR models. See data/factor_info.txt. The second column indicates whether the TF is an activator and the third whether repressor (in theory, an activator could have two roles, thus we have two columns). 

-oo <obj_option>: the option of objective function. Options are: SSE - sum of squared error (default), Corr - average correlation coefficient. 

-mc <max_contact>: the N_MA parameter of the SRR model. 

-p <par_file>: the parameter file. When this option is specified, the values of the parameters in par_file will be used as initial parameter values of the optimizer. In particular, if no parameter estimation is performed (with the option: -na 0, see below), the parameter values will be used for predicting expression patterns. 

-na <nAlternations>: a parameter of the optimizer. The number of alternations between two optimization methods (GSL simplex method and GSL BFGS method). If it is 0, no parameter estimation. Typically 3 to 5. 

-a <annotation_file>: the file of sequences represented by a set of TFBSs (site annotation). With this option, only the specified sites will be used in the sequence-to-expression model. See data/seqs_p002.ann. Note that in the file, the first column is the start position of a site (from 1), the second is the strand of the site, and the last the factor name. 

-ct <coopDistThr>: the distance threshold of cooperative interactions. Default = 50 bp. 

-rt <repressionDistThr>: the distance threshold of short range repression. Default = 150 bp. 

--------------------------------------------------------------------------
SOURCE CODE

The main classes and functions of the program: 

Tools.*
The utility classes and functions, e.g. Matrix class, mathematical functions, I/O classes. 

class Sequence (SeqAnnotator.h): 
The DNA sequences, represented as a vector of A,C,G,T. The functions that read and write sequences in FASTA format are also included. 

class Motif (SeqAnnotator.h): 
The PWM representation of binding profiles. Defined the methods for computing the LLR score and mismatch energy of any sequence element. Also included read/write functions. 

class Site (SeqAnnotator.h): 
The TF binding sites. The data members are: position, which TF the site is associted with (index), and its mismatch energy and binding affinity (:= exp(-energy)). 

class SeqAnnotator (SeqAnnotator.h): 
Create the site representation of a sequence, i.e. a Sequence object can now be represented as a vector of Site objects. All sites that exceed certain energy thresholds will be extracted. 

class FactorIntFunc (ExprPredictor.h): 
The function that computes the interaction between two occupied TFBSs, according to the TF pairs, the distance and orientation of the sites. 

class ExprPar (ExprPredictor.h): 
The parameters of a sequence-to-expression model, including the binding parameters (K of the strongest site multiplied by TF_max, see the paper), the activation parameters (alpha), the repression parameters (beta), the basal transcription, and the TF-TF interaction matrix. Because some models only have a subset of these parameters (e.g. beta is only used in the SRR model),  methods are defined to: 1) create an ExprPar object by a vector of free parameters (constructor); 2) extract all free parameters from the ExprPar object: getFreePars() function. Exactly how these methods are implemented depend on the current model option (static ModelType modelOption). 

class ExprFunc (ExprPredictor.h): 
The class that predicts expression level of a sequence (represented by a vector of Site objects)according to TF expression levels and the model parameters. The main function is predictExpr(). For thermodynamic models, this main function is implemented through compPartFuncOff() and compPartFuncOn(). These two functions implemented the dynamic programming algorithms in the paper. Also note that all model parameters are supposed to known. 

class ExprPredictor (ExprPredictor.h): 
This is the main class of the program. 
* It contains data members that store the input data of the program: sequences and their expression profiles, the PWMs and expression profiles of TF. It computes the best model parameters through the method train() (saved in the data member ExprPar par_model), and can be applied to a new sequence through the method predict(). 
* The main optimization methods are: simplex_minimize() and gradient_minimize(). These two methods use the GSL functions, see GSL manual for the style of using the library. Also the two functions print information when doing optimization. This could be turned off by commenting out the lines in the main do-while loop. 
* Also note that several functions are defined to deal with parameters: testPar() for testing if the parameter ranges are valid and randomSamplePar() sample randomly parameter values (for random starts of the optimizer). 

seq2expr.cpp
The main() function. The default values of many parameters are defined here. Also control the output. 


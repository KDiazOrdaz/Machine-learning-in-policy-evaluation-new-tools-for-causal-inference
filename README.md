# Machine-learning-in-policy-evaluation-new-tools-for-causal-inference
This is the repository holding the code used to perform the analysis used in the manuscript "Machine learning in policy evaluation: new tools for causal inference"


It contains R code to perform parallel super-learning of both propensity score models and outcome regression models, used in the mauscript to obtain Average Treatment Effects (ATE) estimates using Inverse probability of Treatment Weighting (IPTW) ATE, as well as Augmented IPTW and TMLE. 

Other files in this repository perform Double Machine Learning with cross-fitting (Chernozhukov 2018). A separate R file allows the user to perform CTMLE analysis.

Finally, a file implementing a double-lasso for variable selection, including all 2-way interactions is also available. 

Note that the data cannot have missing values. 
You can deal with this in a pre processing step by either using a single imputation (using for example Predictive Mean Matching) *and* adding the missingness indicator columns to the set of potential factors included in all models. This gives valid inferences under the assumption that the counfounders are only confounders when measured. 
Other alternatives are to use Multiple Imputation or Inverse probability of Censoring Weights.

All factor variables have to be "expanded" into dummies.


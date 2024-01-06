# MODS analysis
Gene Expression Signatures Associated with Multiple Organ Dysfunction Syndrome (MODS) Trajectories.

## Table of Contents

- [Description](#description)
- [Overall Workflow](#overall-workflow)
- [Data](#data)
- [Code Description](#code-Description)
- [Dependencies](#dependencies)
- [Citation](#citation)

## Description
Multiple organ dysfunction syndrome (MODS) disproportionately drives sepsis morbidity and mortality among critically ill patients. The biology of this heterogeneous syndrome is complex, dynamic, and incompletely understood. Gene-expression signatures associated with MODS trajectories may improve the prediction of at-risk patients and inform their underlying biology. In the current study, we leveraged publicly available datasets to identify the gene signatures associated with a persistent MODS trajectory among critically ill patients and unraveled biological mechanisms at play. We implemented supervised machine learning (ML) approaches to identify a
parsimonious set of genes predictive of the outcome of interest, trained and validated a model to identify those at high risk of MODS reliably, and demonstrated the reproducibility of our approach across test datasets irrespective of the cause of organ dysfunctions.
Finally, we tested whether our limited set of genes we identified improved upon previously published gene sets that have been demonstrated to predict sepsis mortality in identifying patients at risk of MODS.

## Overall Workflow
![workflow](https://github.com/banerjeeshayantan/CC_MODS_codes/blob/main/Figure%201-1.png)

## Data
![workflow](https://github.com/banerjeeshayantan/CC_MODS_codes/blob/main/datasets-1.png)

**Training data:** Sweeney TE, Shidham A, Wong HR, Khatri P. A comprehensive time-course-based multicohort analysis of sepsis and sterile inflammation reveals a robust diagnostic gene set. Sci Transl Med 2015 May 13;7(287):287ra71. PMID: 25972003

**Validation data:** Snyder A, Jedreski K, Fitch J, Wijeratne S, Wetzel A, Hensley J, Flowers M, Bline K, Hall MW, Muszynski JA. Transcriptomic Profiles in Children With Septic Shock With or Without Immunoparalysis. Front Immunol. 2021 Oct 1;12:733834. doi: 10.3389/fimmu.2021.733834. PMID: 34659221; PMCID: PMC8517409.

**Test data:**  
* Claudia P. Cabrera, Joanna Manson, Joanna M. Shepherd, Hew D. Torrance, David Watson, M. Paula Longhi, Mimoza Hoti, Minal B. Patel, Michael O'Dwyer, Sussan Nourshargh, Daniel J. Pennington, Michael R. Barnes, Karim Brohi. Signatures of inflammation and impending multiple organ dysfunction in the hyperacute phase of trauma: A prospective cohort study
* Shankar R, Leimanis ML, Newbury PA, Liu K et al. Gene expression signatures identify pediatric patients with multiple organ dysfunction who require advanced life support in the intensive care unit. EBioMedicine 2020 Dec;62:103122. PMID: 33248372

## Code Description  
* **MODS_DEG_analysis.R** contains codes required to reproduce the differential expression analysis using the training dataset (GSE66099). It generates a list of differentially expressed genes, a volcano plot, a heatmap, and plots for gene ontology analysis.
* **Prepare_datasets_for_training_testing.R** contains codes to prepare the training, validation, and test sets for machine learning analysis. It performs normalization (if not done already), probes to gene mapping, and averages the expression of duplicate genes. It repeats the same process for both MODS and mortality labels. The mortality labels are needed to compare the machine learning performances between our list of MODS and already published mortality features available from https://www.synapse.org/#!Synapse:syn5612563
* **Stable_features.py** contains code to run the repeated cross-validation experiments to generate a list of 111 stable features for further analysis.
* **validation_and_test.py** contains code to perform validation experiments to determine the best parameters, i.e., top n stable features, sampling-classifier combination, classification threshold, and scaling technique using the dataset E-MTAB-10938. The same parameters are tested on independent validation sets, E-MTAB-5882 and GSE144406. A module to plot the combined ROC curve is also included. 


## Dependencies
scikit-learn - 1.1.1  
pandas - 1.4.1  
numpy - 1.18.5  
imblearn - 0.9.1  
ggplot2 - 3.3.2  
reshape2 - 1.4.4   
stringr - 1.4.0  
tidyr - 1.1.2  
readr - 1.4.0  
caret - 6.0.86

## Citation
Atreya, M. R.^, **Banerjee, S.^**, Lautz, A. J., Alder, M. N., Varisco, B. M., Wong, H. R., Muszynski, J. A., Hall, M. W., Sanchez-Pinto, L. N., Kamaleswaran, R., & Genomics of Pediatric Septic Shock Investigators (2023). Machine learning-driven identification of the gene-expression signature associated with a persistent multiple organ dysfunction trajectory in critical illness. EBioMedicine, 99, 104938. Advance online publication. https://doi.org/10.1016/j.ebiom.2023.104938. [***^ denotes equal contribution***]

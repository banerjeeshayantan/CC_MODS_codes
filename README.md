# MODS analysis
Gene Expression Signatures Associated with Multiple Organ Dysfunction Syndrome (MODS) Trajectories.

## Table of Contents

- [Description](#description)
- [Overall Workflow](#overall-workflow)
- [Data](#data)
- [Dependencies](#dependencies)
- [Preprint Link](#links)
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
Atreya MR* , Banerjee S* , Lautz AJ, et al. Machine learning-driven identification of gene-expression signatures correlated with multiple organ dysfunction trajectories and complex sub-endotypes of pediatric septic shock.
Research Square; 2022. DOI: 10.21203/rs.3.rs-2093663/v1 [* denotes equal contribution]

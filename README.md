# MODS analysis
Gene Expression Signatures Associated with Multiple Organ Dysfunction Syndrome (MODS) Trajectories.

## Table of Contents

- [Description](#description)
- [Overall Workflow](#overall-workflow)
- [Data](#data)
- [Dependencies](#dependencies)
- [Preprint Link](#links)
- [Acknowledgements](#acknowledgements)

## Description
Multiple organ dysfunction syndrome (MODS) disproportionately drives sepsis morbidity and mortality among critically ill patients. The biology of this heterogeneous syndrome is complex, dynamic, and incompletely understood. Gene-expression signatures associated with MODS trajectories may improve prediction of at-
risk patients and inform their underlying biology. In the current study, we leveraged publicly available datasets to identify the gene signatures associated with a persistent MODS trajectory among critically ill patients and unraveled biological mechanisms at play. We implemented supervised machine learning (ML) approaches to identify a
parsimonious set of genes predictive of the outcome of interest, trained and validated a model to identify those at high risk of MODS reliably, and demonstrated the reproducibility of our approach across test datasets irrespective of the cause of organ dysfunctions.
Finally, we tested whether our limited set of genes we identified improved upon previously published gene sets that have been demonstrated to predict sepsis mortality in identifying patients at risk of MODS.

## Overall Workflow
![workflow](https://github.com/banerjeeshayantan/CC_MODS_codes/blob/main/Figure%201-1.png)

## Data
Training data was derived from a study by Brown et al., where they published mutation data from experimental assays labelled as drivers/passengers.
><cite>Brown AL, Li M, Goncearenco A, Panchenko AR (2019) Finding driver mutations in cancer: Elucidating the role of background mutational processes. PLOS Computational Biology 15(4): e1006981. https://doi.org/10.1371/journal.pcbi.1006981</cite>  

Independent test dataset from a benchamrking study by Martelotto et al. consisted of 989 labelled driver and passenger mutations. 
><cite>Martelotto, L.G., Ng, C.K., De Filippo, M.R. et al. Benchmarking mutation effect prediction algorithms using functionally validated cancer-related missense mutations. Genome Biol 15, 484 (2014). https://doi.org/10.1186/s13059-014-0484-1</cite>  

## Dependencies
scikit-learn - 0.22.1  
pandas - 0.25.3  
numpy - 1.18.5  
imblearn - 0.5.0  
ggplot2 - 3.3.2  
reshape2 - 1.4.4   
stringr - 1.4.0  
tidyr - 1.1.2  
readr - 1.4.0  
caret - 6.0.86

## Citation
Banerjee, S.; Raman, K.; Ravindran, B. Sequence Neighborhoods Enable Reliable Prediction of Pathogenic Mutations in Cancer Genomes. Cancers 2021, 13, 2366. https://doi.org/10.3390/cancers1310236

### Acknowledgements
* [Initiative for Biological Systems Engineering](https://ibse.iitm.ac.in/)
* [Robert Bosch Centre for Data Science and Artificial Intelligence (RBCDSAI)](https://rbcdsai.iitm.ac.in/)

<img title="IBSE logo" src="https://github.com/RBC-DSAI-IITM/rbc-dsai-iitm.github.io/blob/master/images/IBSE_logo.png" height="200" width="200"><img title="RBC-DSAI logo" src="https://github.com/RBC-DSAI-IITM/rbc-dsai-iitm.github.io/blob/master/images/logo.jpg" height="200" width="351">

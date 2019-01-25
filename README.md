# Tianyi Liu's Senior Honors Thesis
Repository for the R project for my senior honors thesis. I will constantly add to or edit the files. Questions and concerns regarding the codes should be directed to me [through this link](mailto:tianyi96@live.unc.edu). 

## Research Synopsis
Multiple studies have been published extracting gene expression information from cancer patients via RNA sequencing. In tandem with known patient subtypes, machine learning tools can build on such data to assist in predicting cancer subtype in new patients. It has been previously demonstrated that parametric approaches such as penalized logistic regression can perform well in predicting the tumor subtypes from trial-generated RNA-seq data.  However, non-parametric models such as random forests and support vector machines may offer more robustness to issues such as non-linearity in variable effects and potential interactions between genes. In addition, cross-study variability in the effect of each gene may further impact the accuracy of prediction.  I am currently using existing machine learn approaches, such as random forests and support vector machines, to examine the accuracy of several strategies for training such models to account for cross-study heterogeneity in predicting cancer subtype, including rank-based transformation of expression datasets and horizontal data integration prior to model training. I will also apply recently developed multi-task learning techniques to better account for between-study heterogeneity jointly across multiple datasets and also demonstrate the improvement in prediction accuracy over the prior approaches.

***
## Notes

* Datasets used in this project were uploaded to the 'dataset_used' folder.
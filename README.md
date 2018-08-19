# RJClust
Maintained by: Shahina Rahman.

Authors of the main manuscript: Shahina Rahman and Valen E. Johnson

To use the RJ clustering, the following steps are needed. 

1. All the .R files and .csv file should be downloaded in the same folder RJ. 
2. Set the current folder as the working directory.  Example: setwd("~/Documents/Clustering/RJ")
3. Run the main_RJ.R file. 



The validated labels of the Dataset Dyrskot has two sources. 

Source1:  De Souto et al. 2008, BMC Bioinformatics. Authors maintain the repository of all the 32 datasets used in this article. The label obtained from their processed and filtered dataset "Dyrskjot-2003.txt" is available in https://schlieplab.org/Static/Supplements/CompCancer/datasets.htm. We call their labels as group1.

Source2:  Dyrskjot et al. 2003, Nature Genetics. This is the paper when the data was first published and the data respository is available in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE089 . Along with T2+ and T1 Bladder Carcinoma tumors, authors have categorized Ta tumor cells into three further groups, Ta-grade2, Ta-grade3 and Ta-grade3-CIS. We denote the labels (obtained directly from Table 3 of Dyrskjot et al. 2003, Nature Genetics) as group2.

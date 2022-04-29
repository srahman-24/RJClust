# GMac 
Maintained by: Shahina Rahman.

Authors of the main manuscript: Shahina Rahman, Valen E. Johnson and Suhasini Subba Rao

To use the GMac clustering, the following steps are needed. 

1. All the .R  and .cpp files and .csv file should be downloaded in the same folder RJ. 
2. Set the current folder as the working directory.  Example: setwd("~/Documents/Clustering/GMac")
3. Install 2 external packages "mclust" and "infotheo" by the commands: 
           
           install.packages("mclust")
           install.packages("infotheo")
           
4. Run the main_GMac.R file. 



The validated labels of the dataset "dyrskjot-2003_database.csv" used in main_RJ.R has two sources. 

Source1:  De Souto et al. 2008, BMC Bioinformatics. Authors maintain the repository of all the 32 datasets used in this article. The label obtained from their processed and filtered dataset "dyrskjot-2003_database.txt" is available in https://schlieplab.org/Static/Supplements/CompCancer/Affymetrix/dyrskjot-2003/. We call their labels as group1.

Source2:  Dyrskjot et al. 2003, Nature Genetics. This is the paper when the data was first published and the data respository is available in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE089 . Along with T2+ and T1 Bladder Carcinoma tumors, authors have categorized Ta tumor cells into three further groups, Ta-grade2, Ta-grade3 and Ta-grade3-CIS. We denote the labels (obtained directly from Table 3 of Dyrskjot et al. 2003, Nature Genetics) as group2.

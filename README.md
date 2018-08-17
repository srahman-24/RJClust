# RJClust
Authors: Shahina Rahman and Valen E. Johnson. 

The validated labels has two source. 
Source1:  De Souto et al 2008, BMC Bioinformatics. Authors maintains the repository of all the 32 datasets used in this article. The label is obtained from their processed and filtered dataset "dyrskjot-2003_database.csv". We call it group1.

Source2:  Dyrskot et al 2003, Nature Genetics. This is the original paper when the data was first published and the data respository is available in ? . Along with T2 and T1 Bladder Carcinoma tumors, authors have categorized T2+ tumor cells into three further groups, T2+, T2-grade2 and T2-grade3. We denote the labels obtained directly from Page ? of original article by group2. 

The following Rcode should be copy pasted in the main.R file along with the data to get the original labels denoted by "group"

group1         = c(1:nrow(D));
group1[c(1,2,5,6,8,9,11,22,34)] = 1;
group1[c(3,4,7,10,12,13,14,16,17,18,19,23,26:33)] = 2;
group1[c(15,20,21,24,25,35:40)] = 3;
group2 = group1;
group2[c(3,13,14,30)] = 4;
group2[c(4,7,12,16,18,19,23,31,33)]= 5;

group = group1; 
group = group2

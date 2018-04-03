Table of contents
=================
   * [Model](#model)
   * [Parameters](#parameters)
   * [Summary statistics](#summary-statistics)
   * [Usage](#usage)
      * [Testing for the directionality of introgression](#testing-for-the-directionality-of-introgression)
      * [Variable importance](#variable-importance)


# Model  
**ABC_4pop** is made to investigate various models of speciation between four populations/species/gene-pools.
The following topology of the species tree is the only assumed to date: ( (A, B), (C, D) ). Alternative models only differ by their temporal + geographical patterns of gene flow.  
Thus, gene flow can be unidirectional ( A |--> C ) or bidirectional ( A <--> C ).  


![model](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/model.png)
# Parameters  
_**Tsplit_AB**_ : time of split between population **A** and population **B**. Units of times are in 4.N generations where N is the effective population size of the reference population (better explained in the original document describing _ms_ [found here](https://snoweye.github.io/phyclust/document/msdoc.pdf))    
_**Tsplit_CD**_ : time of split between population **C** and population **D**.  
_**Tsplit**_ : time of split between clade **AB** and clade **CD**.  
  
_**T_SC_AC**_ : time of secondary contact between population **A** and population **C**. Gene flow occured at rates M_AC and M_CA between the present time and _**T_SC_AC**_. Those two migration rates are equal to zero for periods older than _**T_SC_AC**_ (backward in time). Units of migration rates are 4.N.m (better explained in the original document describing _ms_ [found here](https://snoweye.github.io/phyclust/document/msdoc.pdf))  
_**T_SC_BD**_ : time of secondary contact between population **B** and population **D**.  
    
_**N_popA**_, _**N_popB**_, _**N_popC**_ and _**N_popD**_ are the effective population sizes of the current populations **A**, **B**, **C** and **D**.  
_**Na_AB**_ is the effective population size of the ancestral population of **A** and **B**.  
_**Na**_ is the effective population size of the ancestral population of **A**, **B**, **C** and **D**.  

Prior distributions will feed the coalescent simulator [msnsam](https://github.com/rossibarra/msnsam), and are generated by [priorgen](https://github.com/popgenomics/ABC_4pop/blob/master/priorgen_4pop.py). Priorgen does not have to be executed by users, but can be easily modified by them in order to change the prior boundaries, some probability distributions (uniform, Beta, normal, etc ... ).   
Values for all parameters are found in the file **priorfile.txt**, one line per multilocus simulations.  
    
# Summary statistics  
Summary statistics are directly computed from the [msnsam](https://github.com/rossibarra/msnsam)'s output. For each locus, [mscalc](https://github.com/popgenomics/ABC_4pop/blob/master/mscalc_4pop.py) will compute:
![ABCstatistics](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/statistics.png)  
An array of statistics corresponding to the average statistics computed over loci and their standard deviation will be returned every multilocus simulation and written in the file **ABCstat.txt**.  
  
  
# Usage  
**To run the simulations, simply use the following command:**  
```
ABC_4pop.py [model] [migration] [number of multilocus simulations]  
  
Ex: ABC_4pop.py SC_2M_2N AC 10  
```
  
This example will run 10 multilocus simulations, with a secondary contact between population **A** and **C**.  
Migration is heterogeneous over the genome (2M) as well as the effective population size (2N).  
All of the genomic heterogeneities are modeled by a rescaled Beta distribution.  
  
Model in [SI, SC_1M_1N, SC_2M_1N, SC_1M_2N, SC_2M_2N]  
Migration in [none, A, B, C, D, AC, BD, ACBD]  
  
   
**This package requires**:  
 [pypy](https://pypy.org) (has to be linked to the user's bin)  
 [numpy](http://www.numpy.org/)  
 [mscalc_4pop.py](https://github.com/popgenomics/ABC_4pop/blob/master/mscalc_4pop.py) (has to be linked to the user's bin)  
 [priorgen_4pop.py](https://github.com/popgenomics/ABC_4pop/blob/master/priorgen_4pop.py) (has to be linked to the user's bin)  
 [ABC_4pop.py](https://github.com/popgenomics/ABC_4pop/blob/master/ABC_4pop.py) (has to be linked to the user's bin)  
 [msnsam](https://github.com/rossibarra/msnsam)  
  
  
Statistical comparison between "observation" and "simulations" can be made using various R libraries ([_abc_](https://cran.r-project.org/web/packages/abc/abc.pdf), [_abcrf_](https://cran.r-project.org/web/packages/abcrf/abcrf.pdf)).  
  
In the following example, I measured the classification error among four models:  
![confusion matrix](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/confusion_matrix.png)  
Probabilities of each model was inferred for 5,000 of simulated datasets under each model. Errors are lying between 1.5% and 2.76%.  
  
## Testing for the directionality of introgression  
For a given pair of gene pools A and C, one can easily test whether the introgression occurred unidirectionally (A→C or A←C) *versus* bidirectionally (A↔C).
  
## Confusion matrix: all statistics
|                    | SC A←C | SC A↔C | SC A↔C and B↔D | SC B←D | SC B↔D | SC A→C | SC B→D | SI   | class.error |
|:-------------------|:--------|:---------|:-------------------|:---------|:---------|:---------|:---------|:-----|:------------|
| __SC A←C__            | 4693    | 108      | 28                 | 3        | 4        | 65       | 1        | 98   | 0.0614      |
| __SC A↔C__           | 307     | 4298     | 16                 | 3        | 3        | 338      | 2        | 33   | 0.1404      |
| __SC A↔C and B↔ D__ | 25      | 30       | 4844               | 17       | 32       | 13       | 22       | 17   | 0.0312      |
| __SC B←D__           | 6       | 6        | 20                 | 4734     | 105      | 2        | 50       | 77   | 0.0532      |
| __SC B↔D__           | 8       | 6        | 19                 | 327      | 4282     | 1        | 318      | 39   | 0.1436      |
| __SC A→C__           | 60      | 101      | 24                 | 2        | 3        | 4734     | 2        | 74   | 0.0532      |
| __SC B→D__           | 13      | 7        | 14                 | 48       | 105      | 9        | 4725     | 79   | 0.055       |
| __SI__                 | 20      | 13       | 19                 | 13       | 12       | 12       | 9        | 4902 | 0.0196      |

Adding the unidirectional models increased the classification error for bidirectional models compared to the previous analysis. Cases where the difference between the two migration rates is large will be interpreted as unidirectional.  
The following figures show the join distributions for M(A←C) and M(A→C) of the datasets simulated under the model "A↔C" but supported as being:  
"A↔C"  
"A←C"  
"A→C"  
![space_param](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/space_param.png)  

The ABC analysis tends to consider an "A↔C" model as being "A→C" when the migration rates M(A→C) is much greater than M(A←C).  
Even if the model is labelled as bidirectional, randomly chosen values for migration rates produce some simulated datasets for which migration is biologically unidirectional.  

   
Removing some statistics will increase the classification error. 
  
## confusion matrix: effects of discarding ABBA-BABA statistics (__D__ and __fd__)
|                | SC A←C | SC A↔C | SC A↔C and B↔D | SC B←D | SC B↔D | SC A→C | SC B→D | SI   | class.error |
|:---------------|:-------|:-------|:---------------|:-------|:-------|:-------|:-------|:-----|:------------|
| __SC A←C__         | 4606   | 177    | 30             | 7      | 3      | 88     | 3      | 86   | 0.0788      |
| __SC A↔C__         | 351    | 4194   | 15             | 4      | 3      | 399    | 5      | 29   | 0.1612      |
| __SC A↔C and B↔D__ | 32     | 38     | 4833           | 18     | 34     | 16     | 15     | 14   | 0.0334      |
| __SC B ← D__       | 7      | 5      | 20             | 4632   | 172    | 4      | 80     | 80   | 0.0736      |
| __SC B ↔ D__       | 9      | 3      | 19             | 352    | 4196   | 2      | 377    | 42   | 0.1608      |
| __SC A → C__       | 81     | 187    | 25             | 6      | 2      | 4624   | 2      | 73   | 0.0752      |
| __SC B → D__       | 10     | 9      | 18             | 82     | 163    | 7      | 4630   | 81   | 0.074       |
| __SI__             | 26     | 13     | 15             | 21     | 14     | 11     | 17     | 4883 | 0.0234      |
  
  
## confusion matrix: effects of discarding Gmin, minDiv, Gmax and maxDiv statistics 
|                | SC A←C | SC A↔C | SC A↔C and B↔D | SC B←D | SC B↔D | SC A→C | SC B→D | SI   | class.error |
|:---------------|:-------|:-------|:---------------|:-------|:-------|:-------|:-------|:-----|:------------|
| __SC A←C__         | 4680   | 117    | 16             | 3      | 5      | 61     | 1      | 117  | 0.064       |
| __SC A↔C__         | 361    | 4207   | 18             | 4      | 3      | 361    | 2      | 44   | 0.1586      |
| __SC A↔C and B↔D__ | 39     | 48     | 4763           | 28     | 50     | 23     | 25     | 24   | 0.0474      |
| __SC B ← D__       | 3      | 7      | 16             | 4700   | 119    | 0      | 43     | 112  | 0.06        |
| __SC B ↔ D__       | 6      | 4      | 16             | 361    | 4200   | 2      | 352    | 59   | 0.16        |
| __SC A → C__       | 51     | 138    | 13             | 5      | 2      | 4673   | 3      | 115  | 0.0654      |
| __SC B → D__       | 16     | 5      | 18             | 46     | 106    | 5      | 4680   | 124  | 0.064       |
| __SI__             | 22     | 12     | 16             | 19     | 15     | 16     | 6      | 4894 | 0.0212      |
  
  
## confusion matrix: without ABBA-BABA and Gmin statistics
|                | SC A←C | SC A↔C | SC A↔C and B↔D | SC B←D | SC B↔D | SC A→C | SC B→D | SI   | class.error |
|:---------------|:-------|:-------|:---------------|:-------|:-------|:-------|:-------|:-----|:------------|
| __SC A←C__         | 4529   | 214    | 23             | 11     | 7      | 86     | 7      | 123  | 0.0942      |
| __SC A↔C__         | 390    | 4116   | 18             | 8      | 5      | 416    | 10     | 37   | 0.1768      |
| __SC A↔C and B↔D__ | 29     | 43     | 4779           | 35     | 47     | 21     | 25     | 21   | 0.0442      |
| __SC B← D__       | 13     | 5      | 20             | 4531   | 225    | 6      | 85     | 115  | 0.0938      |
| __SC B↔ D__       | 9      | 6      | 25             | 375    | 4143   | 7      | 390    | 45   | 0.1714      |
| __SC A→ C__       | 74     | 232    | 20             | 17     | 6      | 4517   | 8      | 126  | 0.0966      |
| __SC B→ D__       | 11     | 11     | 23             | 81     | 218    | 10     | 4522   | 124  | 0.0956      |
| __SI__             | 39     | 15     | 19             | 44     | 19     | 35     | 27     | 4802 | 0.0396      |


  
## Variable importance
The importance of variables in model classifications can be represented as follows (top 30 over 186):  
![variable_importance](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/variable_importance.png)    
__ABBA-BABA__ _**fd**_ and _**D**_ statistics as well as _**Gmin**_ are the statistics with the greatest power to classify models, _i.e_, they reduce the dispersion of models along branches of each decision tree (measured by the Gini index).  
  
However, it is important to note that these statistics are either sensitive to an available outgroup (_**ABBA-BABA**_), or are dependent on good quality inferences of haplotype phases (_**Gmin**_).  
Although they are the most informative on simulated data, biases in their measurements on real data can induce biases in inferences. Excluding these statistics does not greatly reduce inferential power. The remaining 111 statistics do not depend as much on external groups or third-party inferences, and they contain enough combined information to make model selection.  
The choice of whether or not to include them in the inferences should belong only to the experimenter.



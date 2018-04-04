ABC 4 pop
=================
   * [Model](#model)
   * [Parameters](#parameters)
   * [Summary statistics](#summary-statistics)
   * [Usage](#usage)
   * [Results](#results)
      * [Testing for 'ongoing migration' _versus_ 'current isolation'](#testing-for-ongoing-migration-versus-current-isolation)
      * [Testing for the directionality of introgression](#testing-for-the-directionality-of-introgression)
      * [Effects of statistics on model confusion](#effects-of-statistics-on-model-confusion)
      * [Estimating the parameters](#estimating-the-parameters)

# Model  
**ABC_4pop** is made to investigate various models of speciation between four populations/species/gene-pools.
The following topology of the species tree is the only assumed to date: ( (A, B), (C, D) ). Alternative models only differ by their temporal + geographical patterns of gene flow.  
Thus, gene flow can be unidirectional ( A→C ) or bidirectional ( A↔C ).  


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
Values for all parameters are written in the file **priorfile.txt**, one line per multilocus simulations.  
    
# Summary statistics  
Summary statistics are directly computed from the [msnsam](https://github.com/rossibarra/msnsam)'s output. For each locus, [mscalc](https://github.com/popgenomics/ABC_4pop/blob/master/mscalc_4pop.py) will compute:

| Statistics         | Description                                                                 |
|:-------------------|:----------------------------------------------------------------------------|
| __bialsites__          |  number of SNPs in the alignment                                            |
| __sf XY__               |  number of fixed differences between species X and Y / locus length         |
| __sx X__                |  number of exclusively polymorphic positions in species X / locus length    |
| __ss XY__               |  number of shared biallelic positions between species X and Y/ locus length |
| __pi X__                |  Tajima’s Theta within species X                                            |
| __theta X__             |  Watterson’s Theta witin species X                                          |
| __pearson_r_pi XY__    |  correlation’s coefficient for pi over orthologs between X and Y            |
| __pearson_r_theta XY__ |  correlation’s coefficient for theta over orthologs between X and Y         |
| __Dtaj X__              |  Tajima’s D for species X                                                   |
| __div XY__              |  raw divergence Dxy measured between X and Y                                |
| __netdiv XY__           |  net divergence Da measured between X and Y                                 |
| __minDiv XY__           |  smallest divergence measured between one individual from X and one from Y  |
| __maxDiv XY__           |  highest divergence measured between one individual from X and one from Y   |
| __Gmin XY__             |  minimum divergence between one sequence from X and one from Y __minDivXY__ divided by the average __divXY__                                                             |
| __Gmax XY__             |  __maxDivXY/divXY__                                                             |
| __FST XY__             |  FST between X and Y compute as 1-(pi_X + pi_Z) / (2*pi_XY)                                                        |
| __D XY_Z__             |  ABBA-BABA’s D statistics where pop_1 = X, pop_2 = Y and pop_3 = Z          |
| __fd XY_Z__            |  ABBA-BABA’s fd statistics where pop_1 = X, pop_2 = Y and pop_3 = Z         |

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
  
ABC_4pop.py will simply execute the pipeline  
```
priorgen | msnsam | mscalc
```
The output files are __ABCstat.txt__ (containing the computed summary statistics) and __priorfile.txt__ (containing the parameter values used to simulate the data from which the summary statistics were calculated).  

  
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
 
 
# Results  
## Testing for 'ongoing migration' _versus_ 'current isolation'  
Four models were compared: secondary contacts (_i)_ between A and C, _ii)_ between A/C and B/D, _iii)_ between B and D) and strict isolation (__SI__).  
The classification error of the ABC approach when it comes to specifically comparing these 4 models are shown in the following table:  
![confusion matrix](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/confusion_matrix.png)  
Classification error corresponds to the percentage over 5,000 of simulated datasets under the model _Mi_ that were not statistically supported as being _Mi_. Those errors are lying between 1.5% and 2.76%.  
  
  
## Testing for the directionality of introgression  
For a given pair of gene pools A and C, one can easily test whether the introgression occurred unidirectionally (A→C or A←C) *versus* bidirectionally (A↔C).
  
## Confusion matrix (using all statistics)
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

Adding the unidirectional models increased the classification error for bidirectional models compared to the previous analysis.   
The following figures show the join distributions for M(A←C) and M(A→C) of the datasets simulated under the model "A↔C" but supported as being:  
"A↔C"  
"A←C"  
"A→C"  
![space_param](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/space_param.png)  

__Figure A:__ Combinations of the two migration rates (M(A←C) and M(A→C)) producing datasets that were correctly supported as being the A↔C model.  
__Figure B:__ Combinations of the two migration rates (M(A←C) and M(A→C)) producing datasets that were uncorrectly supported as being the A←C model.  
__Figure C:__ Combinations of the two migration rates (M(A←C) and M(A→C)) producing datasets that were uncorrectly supported as being the A→C model.  
Cases where the difference between the two migration rates is large will be interpreted as unidirectional. Thus, the ABC analysis tends to consider an "A↔C" model as being "A→C" when the migration rate M(A→C) is much greater than M(A←C).  
Even if the model is labelled as bidirectional, randomly chosen values for migration rates produce some simulated datasets for which migration is biologically unidirectional.  

   
## Effects of statistics on model confusion  
The same model comparison was also realized by removing statistics depending on an available sequenced outgroup (_no ABBA-BABA stats_ analysis), by removing statistics depending on the identification of gametic phases (_no Gmin stats_ analysis) or by only keeping statistics describing a folded site-frequency-spectrum (_no ABBA-BABA and Gmin stats_ analysis).  
For each analysis, the error rate in model classification was also measured for 5,000 of randomly simulated datasets under each of the 8 alternative models. This error rate simply corresponds to the rate of misclassification.  
  
| Targets            | All stats | No ABBA-BABA stats | No Gmin stats | No ABBA-BABA and Gmin stats |
|:-------------------|:----------|:-------------------|:--------------|:----------------------------|
| SC A← C            | 0.0614    | 0.0788             | 0.064         | 0.0942                      |
| SC A ↔ C           | 0.1404    | 0.1612             | 0.1586        | 0.1768                      |
| SC A ↔ C and B ↔ D | 0.0312    | 0.0334             | 0.0474        | 0.0442                      |
| SC B ← D           | 0.0532    | 0.0736             | 0.06          | 0.0938                      |
| SC B ↔ D           | 0.1436    | 0.1608             | 0.16          | 0.1714                      |
| SC A → C           | 0.0532    | 0.0752             | 0.0654        | 0.0966                      |
| SC B → D           | 0.055     | 0.074              | 0.064         | 0.0956                      |
| SI                 | 0.0196    | 0.0234             | 0.0212        | 0.0396                      |
  
Removing ABBA-BABA or Gmin statistics will individually increase the classification error by a low factor, however, this reduction begins to be worrisome when both categories of statistics are neglected. It should always be borne in mind that model errors concern border line cases between two models.  
  
__ABBA-BABA__ _**fd**_ and _**D**_ statistics as well as _**Gmin**_ have an important power to classify models, _i.e_, they reduce the dispersion of models along branches of each decision tree (measured by the Gini index).  
![variable_importance](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/variable_importance.png)    
  
However, it is important to note that these statistics are either sensitive to an available outgroup (_**ABBA-BABA**_), or are dependent on good quality inferences of haplotype phases (_**Gmin**_).  
Although they are the most informative on simulated data, biases in their measurements on real data can induce biases in inferences. Excluding these statistics does not greatly reduce inferential power. The remaining 111 statistics do not depend as much on external groups or third-party inferences, and they contain enough combined information to make model selection.  
The choice of whether or not to include them in the inferences should belong only to the experimenter.

## Estimating the parameters  
### SI model  
Parameters of the SI model were estimated for 500 simulated datasets.  
![parameters_SI](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/parameters_SI.png)  
Each point represent a simulated dataset.  
The x-axis shows the real value to estimate.  
The y-axis shows the estimated value by using a random-forest approach.  
The red line is a line of slope 1 passing through the origin.  
  
Results for current populations are shown only for population A, and shows similar relationships between _real values_ and _estimated values_ for B, C and D.  
In the same way for Na_AB and Na_CD, and for Tsplit_AB and Tsplit_CD.  

### Effect of gene flow on paramters estimates  
![parameters_tsplit_5models](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/tsplit_5models.png)  


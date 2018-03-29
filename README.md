# Model  
**ABC_4pop** is made to investigate various models of speciation between four populations/species/gene-pools.
The following topology of the species tree is the only assumed to date: ( (A, B), (C, D) ).

![model](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/model.png)
### Parameters  
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
    
### Summary statistics  
Summary statistics are directly computed from the [msnsam](https://github.com/rossibarra/msnsam)'s output. For each locus, [mscalc](https://github.com/popgenomics/ABC_4pop/blob/master/mscalc_4pop.py) will compute:
![ABCstatistics](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/statistics.png)  
An array of statistics corresponding to the average statistics computed over loci and their standard deviation will be returned every multilocus simulation and written in the file **ABCstat.txt**.  
  
  
# Usage  
This package requires:  
[pypy](https://pypy.org) (has to be linked to the user's bin)  
[numpy](http://www.numpy.org/)  
[mscalc_4pop.py](https://github.com/popgenomics/ABC_4pop/blob/master/mscalc_4pop.py) (has to be linked to the user's bin)  
[priorgen_4pop.py](https://github.com/popgenomics/ABC_4pop/blob/master/priorgen_4pop.py) (has to be linked to the user's bin)  
[ABC_4pop.py](https://github.com/popgenomics/ABC_4pop/blob/master/ABC_4pop.py) (has to be linked to the user's bin)  
  
**To run the simulations, simply use the following command:  **
ABC_4pop.py [model] [migration] [number of multilocus simulations]  
  
**Ex:** ABC_4pop.py SC_2M_2N AC 10  
  
This example will run 10 multilocus simulations, with a secondary contact between population **A** and **C**.  
Model in [SI, SC_1M_1N, SC_2M_1N, SC_1M_2N, SC_2M_2N]  
Migration in [none, AC, BD, ACBD]  
  
Statistical comparison between "observation" and "simulations" can be made using various R libraries ([_abc_](https://cran.r-project.org/web/packages/abc/abc.pdf), [_abcrf_](https://cran.r-project.org/web/packages/abcrf/abcrf.pdf)).  
  
In the following example, I measured the classification error among four models:  
SC_1M_1N with migration between **A** and **C**.  
SC_1M_1N with migration between **A** and **C**, and migration between **B** and **D**.  
SC_1M_1N with migration between **B** and **D**.  
SI (no migration)  

![confusion matrix](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/confusion_matrix.png)  
  
In this example, I measured over 5,000 of simulated data under each model the error rate in classification. Errors are lying from 1.5% to 2.76%.  
  

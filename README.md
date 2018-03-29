**ABC_4pop** is made to investigate various models of speciation between four populations/species/gene-pools.
The following topology of the species tree is the only assumed to date: ( (A, B), (C, D) ).

![model](https://github.com/popgenomics/ABC_4pop/blob/master/pictures/model.png)

_**Tsplit_AB**_ : time of split between population **A** and population **B**. Units of times are in 4.N generations where N is the effective population size of the reference population (better described in the original document describing _ms_ [found here](https://snoweye.github.io/phyclust/document/msdoc.pdf))    
_**Tsplit_CD**_ : time of split between population **C** and population **D**.  
_**Tsplit_ABCD**_ : time of split between clade **AB** and clade **CD**.  
_**T_SC_AC**_ : time of secondary contact between population **A** and population **C**. Gene flow occured at rates M_AC and M_CA between the present time and _**T_SC_AC**_. Those two migration rates are equal to zero for periods older than _**T_SC_AC**_ (backward in time). Units of migration rates are 4.N.m (better described in the original document describing _ms_ [found here](https://snoweye.github.io/phyclust/document/msdoc.pdf))  


# PyChange - Change detection with Python

## Quickstart

### Example1: Differences in two sequences.  

```
python PyChange.py --filename=random.csv --cell=B --values=A --time=T --method=MaChaMP --preprocessing=diff
```

Computes the changes in trend in two sequences stored in `random.csv` with the `MaChaMP` method. One of the sequences has a significant trend change, the other does not.   

### Example2: Comparsion of methods.   

```
python Example.py
```

Creates [`Example.pdf`](./Example.pdf), where all available methods are applied to the same step function.    

## Use as a script   

Store your time series in a .csv file and run PyChange.  

`--filename`: Name of `.csv` File (e.g. from standard QtFy output).  
`--cell`: Name of unique time series identification column.     
`--values`: Name of column of time-series values    
`--time`: Name of column of timepoints.    

`--method`: Type of change detection method.   
- `PELT`: **default** Linear method based on dynamic programming from the popular r-package [changepoint](https://cran.r-project.org/web/packages/changepoint/index.html). 
- `SMUCE`: For changes at different scales using the r-package [stepR](https://cran.r-project.org/web/packages/stepR/index.html).    
- [`WBS`](https://cran.r-project.org/web/packages/wbs/index.html): Stochastic method. 
- `E-Divise`: Non-parametric offline from the r-package [ecp](https://cran.r-project.org/web/packages/ecp/index.html).    
- `Lepage`: Non-parametric online method using [cpm](https://cran.r-project.org/web/packages/cpm/index.html).   
- [`Segmentor3IsBack`](https://cran.r-project.org/web/packages/Segmentor3IsBack/index.html): An r-package for Seq-data.  
- [`Fpop`](https://r-forge.r-project.org/R/?group_id=1851): OP using functional pruining for large data.   
- [`BCP`](https://cran.r-project.org/web/packages/bcp/index.html): Bayesian changepoint analysis.   

`--preprocessing`: Transformation of sequence before applying change detection. 
- `none`: **default**  Raw sequence, for change in mean.     
- `diff`: Change detection on gradient sequence, to remove linear dependencies   
- `logdiff`: To remove multiplicative dependencies in `diff`.  
- `percdiff`: Percentile differences for changes in intensities.    
- `logpercdiff`: To remove multiplicative dependiencies of `percdiff`. 

If `--filename=Name.csv`, the output is a `ChangesName.csv` file. The columns are: CellID, Changepoint location, and timepoint of change.    

## Use as a module  

```
import nump as np
from PyChange import PyChange

r = np.random.RandomState(42)
seq = list(r.randn(200)) + list(r.randn(200) + 1.) + list(r.randn(200)) + list(r.randn(200) + 1.)

cp = PyChange(seq)
print cp
```

The function `PyChange` has three input variables:   
- `seq`: **required** The time seires as a list.  
- `transform`: the preprocessing as above.    
- `method`: the method as above.   

## Remarks  

Be carful:    
- `r` and the desired package needs to be installed.      
- `E-Divise` is slow.    
- The Bayesian methods return posterior probabilities for each point being a changepoint.  
- All methdos other are set to the settings suggested in repective the documentation.   

## Any questions?  

Come to my office BSB 2.01.   
Or email: thomas.gumbsch@bsse.ethz.ch 

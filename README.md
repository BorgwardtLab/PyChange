# PyChange - Change detection with Python

## Quickstart

### Example1: Differences in two sequences.  

```
python3 PyChange.py --filename=random.csv --cell=B --values=A --time=T --method=CUSUM --preprocessing=diff
```

Computes the changes in trend in two sequences stored in `random.csv` with the `MaChaMP` method. One of the sequences has a significant trend change, the other does not.   

### Example2: Comparsion of methods.   

```
python3 Example.py
```

Creates [`Example.pdf`](./Example.pdf), where all available methods are applied to the same step function.    

## Use as a script   

Store your time series in a .csv file and run PyChange.  

`--filename`: Name of `.csv` File (e.g. from standard QtFy output).  
`--cell`: Name of unique time series identification column.     
`--values`: Name of column of time-series values    
`--time`: Name of column of timepoints.    

`--method`: Type of change detection method.   
- `CUSUM`: Pages Cummulative sum control chart. (h=8.01, k=0.25).     
- `EWMA`: Expoentially weighted moving average control chart. (burn-in period 30, lambda = 0.05, L=3.7).  
- `QChart`: For start up processes and short runs (3o3, maxlike=2.).    

`--preprocessing`: Transformation of sequence before applying change detection. 
- `none`: **default**  Raw sequence, for change in mean.     
- `diff`: Change detection on gradient sequence, to remove linear dependencies   
- `logdiff`: To remove multiplicative dependencies in `diff`.  
- `percdiff`: Percentile differences for changes in intensities.    
- `logpercdiff`: To remove multiplicative dependiencies of `percdiff`. 

If `--filename=Name.csv`, the output is a `ChangesName.csv` file. The columns are: CellID, Changepoint location, and timepoint of change.    

## Use as a module  

```
import numpy as np
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

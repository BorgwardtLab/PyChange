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

Creates the plot `Example.pdf` where each method is applied to the same step function.   


## Use as a script   

Store your time series in a .csv file and run PyChange.  

`--filename`: Name of `.csv` File (e.g. from standard QtFy output).  
`--cell`: Name of unique time series identification column.     
`--values`: Name of column of time-series values    
`--time`: Name of column of timepoints.    

`--method`: Type of change detection method. 
- `MaChaMP`: **default** reports the biggest configuration of changes using Welch's t-Test and dimensionalty analysis.     
- `PELT`: Linear method from the popular r-package [changepoint](https://cran.r-project.org/web/packages/changepoint/index.html). 
- `SMUCE`: For changes at different scales from the r-package [stepR](https://cran.r-project.org/web/packages/stepR/index.html).    
- `WBS`: Stochastic method implemented in [wbs](https://cran.r-project.org/web/packages/wbs/index.html)   
- `E-Divise`: Non-arametric method from the r-package [ecp](https://cran.r-project.org/web/packages/ecp/index.html).   

`--preprocessing`: Transformation of sequence before applying change detection. 
- `none`: **default**  Raw sequence, for change in mean.   
- `diff`: Change detection of gradient sequecne.   
- `logdiff`: To remove skewness in `diff`.  
- `percdiff`: Percentile differences. Use when time series incements are interesting wrt on abs value, i.e. financial derrivatives.    
- `logpercdiff`:To remove skewness of `percdiff`. 

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
- Changes in differences have a high FN rate.     
- `E-Divise` is slow.    

## Any questions?  

Come to my office BSB 2.01.   
Or email: thomas.gumbsch@bsse.ethz.ch 

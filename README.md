# PyChange - Multiple change detection with python

## Quickstart:

```
python PyChange.py --filename=random.csv --cell=B --values=A --time=T --method=MaChaMP --preprocessing=diff
```

computes the changes in trend in two sequences stored in `random.csv` with the `MaChaMP` method. One of the sequences has a significant trend change, the other does not.   

## Use as a script:   

Store your time series in a .csv file and run PyChange.  

`--filename`: Name of `.csv` File (e.g. from standard QtFy output).  
`--cell`: Name of unique time series identification column.     
`--values`: Name of column of time-series values    
`--time`: Name of column of timepoints.    

`--method`: Type of change detection method. 
- `MaChaMP`: reports the biggest combination of changes with Welch's t-Test and a window method.   **default**. 
- `PELT`: popular r-package based on dynamic programming.   
- `SMUCE`: for changes at different scales.  
- `WBS`: stochastic method, n^2.    
- `E-Divise`: Based on dp, n^3   **carful: this is slow!!!**. 

`--preprocessing`: Transformation of sequence before applying change detection. 
- `none`: raw sequence, for change in mean.  **default**. 
- `diff`: Change detection of gradient sequecne.   
- `logdiff`: To remove skewness in `diff`.  
- `percdiff`: Percentile differences. Use when time series incements are interesting wrt on abs value, i.e. financial derrivatives.    
- `logpercdiff`:To remove skewness of `percdiff`. 

If `--filename=Name.csv`, the output is a `ChangesName.csv` file. The columns are: CellID, Changepoint location, and timepoint of change.    

## Use as a Python module:  

```
from PyChange import PyChange
```




## Remarks  

Be carful:  
- Changes in differences have a high FN rate 

## Any questions?  

Come to my office BSB 2.01.   
Or email: thomas.gumbsch@bsse.ethz.ch 

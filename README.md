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

`--method`: Type of change detection method
- `MaChaMP`: reports the biggest combination of changes with Welch's t-Test and a window method.   **default**
- `PELT`: popular r-package based on dynamic programming.  
- `WBS`: stochastic method, n^2.  
- `SMUCE`: other method  
- `E-Divise`: dp.   **carful: this is slow!!!**


## Output

Is a `.csv` file. Columns are: CellID, Changepoint location, (Fishers) combined *p*-value, Timepoint of change.  
**Warning:** The *p*-values are not corrected for multiple hypothesis testing when looking for more than one change!  

## Remarks  

Be carful:  
- Either a change in trend or a change in mean - this is due to a two-sample t-Test behind the routines.  
- For many changes, the *p*-values are not corrected for multiple hypothesis testing. They can therefore not be interpreted as *p*-values. 
- General rule of thumb: Do not input sequences longer than 100!  

## Any questions?  

Come to my office BSB 2.01.   
Or email: tgumbsch@gmail.com. 

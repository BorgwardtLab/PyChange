# PyChange - Multiple change detection with python

## Quickstart:

```
python PyChange.py random.csv B A single_diff
```

computes the single biggest change in trend in two sequences stored in `random.csv`. One of the sequences has a significant trend change, the other does not.   

## Arguments explained. 

```
python PyChange.py Filename CellID Values Method
```

`Filename`: Name of `.csv` File (e.g. from standard QtFy output).   
`CellID`: Name of unique cell identification column.   
`Values`: Column of intensities.  

`Method`: Type of change detection method
- `single_diff`: *standard* reports single biggest change on gradient sequence.  
- `single`: single biggest change in mean on raw sequence.  
- `DP`: Dynamic programming method mines for multiple changes.  
- `exact`: Exact method mines for multiple changes **carful: this is slow!!!**.   
- `ANOVA`: Find all pairwise dissimilar subsequences.  **careful: This is not what you usually want**
- all methods can be combined with `_diff` for applying differences to the sequence.  
- Try also `exact_diff_3` to find the three biggest changes accoring to the exact method.  


## Output

Is a `.csv` file. Columns are: CellID, Changepoint location, (Fishers) combined *p*-value.  
**Warning:** The *p*-values are not corrected for multiple hypothesis testing when looking for more than one change!  

## Remarks  

Be carful:  
- Either change in trend or change in mean - this is due to a two-sample TTest behind all the routines.  
- For many changes, meaningful *p*-value can only be assigned to a combination of changes.   
- General rule of thumb: Do not input sequences longer than 100!  

## Any questions?  

Come to my office BSB 2.01.   
Or email: tgumbsch@gmail.com. 

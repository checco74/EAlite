# EAlite

## Quick start

Get GIT repository: 
```
git clone https://github.com/scimerc/EAlite
cd EAlite
```

Install EAtestslite R-package:
```
R CMD INSTALL EAtestslite
```

Run ealite script on example data: 
```
Rscript ealite.Rscript --annot example_data/annot-sorted.txt example_data/pvalues.txt --qq
```

This will generate enrichment plots of the pvalues in 
`example_data/pvalues.txt` for the categories in 
`example_data/annot-sorted.txt`.


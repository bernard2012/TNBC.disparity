# TNBC.disparity

## Running Linear Mixed Model to get race-specific cell-cell interactions

Prerequisite: R and lme4 package.

### Step 1:
```
cd groups
python3 extract.py
```

This will split the `full.matrix.txt` into a series of files, corresponding to interactions.
For example, the file "7 CD152 CD8a--7 CD152 CD8a.txt" contains:

```
7 CD152 CD8a--7 CD152 CD8a  AA  patient5    1.55560
7 CD152 CD8a--7 CD152 CD8a  AA  patient6    1.28396
7 CD152 CD8a--7 CD152 CD8a  AA  patient6    -0.11636
7 CD152 CD8a--7 CD152 CD8a  AA  patient6    2.22011
7 CD152 CD8a--7 CD152 CD8a  AA  patient4    1.26598
7 CD152 CD8a--7 CD152 CD8a  AA  patient4    1.02704
7 CD152 CD8a--7 CD152 CD8a  AA  patient7    -0.27023
7 CD152 CD8a--7 CD152 CD8a  AA  patient7    1.42868
7 CD152 CD8a--7 CD152 CD8a  AA  patient1    0.22347
7 CD152 CD8a--7 CD152 CD8a  AA  patient1    2.80255
7 CD152 CD8a--7 CD152 CD8a  AA  patient2    1.97753
7 CD152 CD8a--7 CD152 CD8a  AA  patient2    2.31938
7 CD152 CD8a--7 CD152 CD8a  AA  patient3    1.71408
7 CD152 CD8a--7 CD152 CD8a  AA  patient18   1.82975
7 CD152 CD8a--7 CD152 CD8a  AA  patient14   1.78164
7 CD152 CD8a--7 CD152 CD8a  AA  patient14   -0.00288
...
```
Where the first column is the interaction, followed by race, then patient ID, and finally interaction z-score (positive if enrichment or negative if depleted).


### Step 2:

Make sure `command1.sh` is in the directory. Then run:
```
./command1.sh
```

Content of `command1.sh`:
```
IFS=$'\n'; for i in `ls -1 *.txt`; do Rscript do.one.R "$i"; done
```

Content of `do.one.R`:
```
options(echo=T)
library(lme4)
args<-commandArgs(trailingOnly = TRUE)
print(args)
x<-read.table(args[1], sep="\t")
colnames(x)=c("interaction", "group", "patient", "score")
mixed=lmer(score ~ 1+ group + (1|patient), data=x)
reduced.mixed=lmer(score ~ 1 + (1|patient), data=x)
an=anova(mixed, reduced.mixed)
#an$"Pr(>Chisq)"
tt<-cbind(t(fixef(mixed)), t(an$"Pr(>Chisq)"))
write.table(tt, file=paste0("stats/", args), sep="\t", quot=F, row.names=F)
```

This R script `do.one.R` performs the linear mixed mode on an interaction file, e.g. 7 CD152 CD8a--7 CD152 CD8a.txt, and outputs the model coefficient and P-value significance for the statistical comparison: interaction in AA versus interaction in EA. Please see our BioRxiv paper for explanation and rationale for using linear mixed model.

The `command1.sh` iterates through all interactions there are in the directory.
The outputs are in `stats` directory.

### Step 3:

```
cd stats
source ../command.sh
```

Content of `command.sh`:
```
IFS=$'\n'; for i in `ls -1`; do echo $i; cat "$i"|sed "1d"; done|paste - -|sort -t" " -g -k2 -r|less
```

This will combine all the testing results from all interactions into a summary table, like the following:
```
9 GranzymeB CD152 HIF1a--9 GranzymeB CD152 HIF1a.txt    1.37754484199537        0.328769288177209       NA      0.150836842476321
8 PLK1 PanCK Ki67--9 GranzymeB CD152 HIF1a.txt  -0.634285221486421      0.252390878398147       NA      0.214136418469996
8 PLK1 PanCK Ki67--8 PLK1 PanCK Ki67.txt        1.05368343509042        0.0013638392215268      NA      0.990284913892305
7 CD152 CD8a--9 GranzymeB CD152 HIF1a.txt       -0.0166020416521173     0.292592851817632       NA      0.049433886078644
```

Where the first column is the interaction, second column is the coefficient (AA), third column is the coefficient (EA), and the last column is the Padj value.
For example, the interaction "7 CD152 CD8a--9 GranzymeB CD152 HIF1a.txt" is significant (Padj=0.049) and is higher in EA than AA.



# MSGene
R package for multistate modeling

Welcome to our package for estimating lifetime probabilities iterated over a multistate model!! 

To install:

```{r}
devtools::install_github("surbut/MSGene");
library(MSGene)
```

to run this code with minimal example:
```{r}
data("testdf")
```

# set a training set

Now, we select a portion of the data for training
```{r}
train=df2[1:(nrow(df2)*0.80),]
```

# Parameter choice

Here we choose the ages we will consider and the states we will move through.

```{r}
ages=c(40:80)
nstates=c("Health", "Ht","HyperLip","Dm","Cad","death","Ht&HyperLip","HyperLip&Dm","Ht&Dm","Ht&HyperLip&Dm")
```

## Fit our model

Now we're ready to fit our model! Sit tight, this step takes between 60-90 seconds  ... but it's worth it.

```{r}
modelfit = fitfunc2(
  df_frame=data.table::data.table(train),
  ages = ages,
  nstates = nstates,
  mode = "binomial",
  covariates = "cad.prs+f.31.0.0+smoke+antihtn_now+statin_now")
```

For more details, including state to state coefficient extraction, lifetime probabilities of disease under treated and untreated strategy, please see our test vignette at https://surbut.github.io/MSGene/vignette.html



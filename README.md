# MSGene
R package for multistate modeling

Welcome to our package for estimating lifetime probabilities iterated over a multistate model!! 

To install:

```{r}
devtools::install_github("surbut/MSGene");
library(MSGene)
```

To run this code with minimal example and compute the model fit, we create a sample data frame with the following columns. Note that cases are written as 2, controls as 1. If an individual does not have a condition, his censor age is the age of last follow up or death. Please see [ukbpheno](https://github.com/niekverw/ukbpheno/tree/master) for directions on deriving case control status from your data. Here statin and antihtn refer to the use of statins or antithn, and age (or NA) of prescription. This is to produce the training and computation listed in Urbut et al. 

# Load the dummy data

```{r}
df2 <- readRDS(system.file("data", "msgene_sampledf.rds", package = "MSGene"))
head(df2)
```


# set a training set

Now, we select a portion of the data for training, and convert to a data table
```{r}
train=data.table(df2[1:(nrow(df2)*0.80),])
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
  print(attributes(modelfit))
```

This produces the model fit. For more details, including state to state coefficient extraction, lifetime probabilities of disease under treated and untreated strategy, please see our test vignette at 
https://surbut.github.io/MSGene/vignette.html 


For more flexible states inputting an array of state positions and possible transitions with any names, please see:
https://surbut.github.io/MSGene/usingMarginal.html



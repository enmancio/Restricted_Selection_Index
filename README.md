# Restricted Selection Index


## Upload or download all packages and built-in functions:

```r
set.seed(1526)
source("function/JDS_function.R") # load all built-in function
pkg = c("AlphaSimR","MASS","ggplot2","pedigreemm","sparseinv","reshape2")
libraries(pkg)
```

### Use AlphaSimR to simulate the first population

In this example, we mimic a population with similar parameters to the one described in Mancin et al. 2022 (Economic weights for restriction of selection index as the optimal strategy for combining multiple traits ). Note that the genetic variance and phenotype mean have been changed for a matter of clarity and a better visual appraisal of the graphs. Indeed, here we assumed a genetic variance equal to 1 and a phenotype mean equal to 0.

```r 
pop_size = 10000

# we simulate traits similar to those present in the article 
# using the same heritability and genetic correlations 
# see Table 1 on the mauscript 
n_traits = 5
#  number of candidates males
n_candidate_male = as.integer(pop_size/50)
#  number of candidates females
n_candidate_female = as.integer(pop_size/4)
# we compute selection intensity using function of AlphaSimR
selInt((n_candidate_female+n_candidate_male)/pop_size)
# genetic correlations matrix
G_cor = matrix(c(1,0.758,0.845,0.069,-0.458,
                0.758,1,0.824,0.067,-0.413,
                0.845,0.824,1,0.088,-0.397,
                0.069,0.067,0.088,1,-0.156,
               -0.458,-0.413,-0.397,-0.156,1),ncol=5)
# heritability 
h2 = c(0.219,0.178,0.125,0.13,0.328)

var_G = rep(1,n_traits)
meanP = rep(0,n_traits)
# economic weights  (same values as in the manuscript)
econWt = c(0,0.3,0.4,0,0.3)
```

In this part, we simulate founder haplotype and general population parameters

```r
# create founder haplotype
founderPop = quickHaplo(nInd=pop_size, nChr=10, segSites=1000)
# set population parameters
SP = SimParam$new(founderPop)
# generate the additive traits with the previous parameters
SP$addTraitA(nQtlPerChr = 1000, mean=meanP, var=var_G, corA=G_cor)
SP$setVarE(h2=h2)
SP$setSexes("yes_sys")
SP$setTrackPed(TRUE)
 
```
In the following example **six** different scenarios have been simulated using the initial economic weights proposed in the manuscript (**econWt**) and also defined above:

1. Selection using **unrestricted** phenotypic selection index.
2. Selection using **restricted** phenotypic selection index.
3. Selection using **restricted** phenotypic selection index with sex-limited traits, **without accounting** for different accuracy across candidates.
4. Selection using **restricted** phenotypic selection index with sex-limited traits **accounting** for different accuracies across candidates.
5. Selection using **restricted genetic** selection index with sex-limited traits, **without accounting** for different accuracies across candidates for the estimation of restricted economic weight
6. Selection using **restricted genetic** selection index with sex-limited traits using the adaptation of Lin (1990) for **accounting** of different accuracies across candidates, as proposed in our manuscript.


### Scenario 1.

```r 
# create first generation
pop = newPop(founderPop, simParam=SP)
G = varG(pop)
P = varP(pop)

# calcolate unrestricted b
b = smithHazel(econWt, G, P)

genMean = data.frame(t(meanG(pop)))  # mean of the breeding values at present generation
names(genMean) = paste0("traits",1:n_traits)

# perform selection for 4 generations
for(generation in 1:4){
  # select the best females as candidate dams
  damsc = selectInd(pop,nInd=n_candidate_female, trait=selIndex,sex="F", b=b)
  # select the best males as candidate sires
  sirec = selectInd(pop,nInd=n_candidate_male, trait=selIndex,sex="M", b=b)
  # create population of candidates to selection
  candidate=mergePops(list(damsc,sirec)) 
  # create next generation by random crossing sires and dams
  pop=randCross(candidate,nCrosses=pop_size)
  cat("generation ",generation,"\n")
  # collect genetic mean for each trait  
  genMean = rbind(genMean,meanG(pop))
}

# plot genetic progress per generation
plot_g_prog(genMean)
```
[![name](/enmancio/Restricted_Selection_Index/)]

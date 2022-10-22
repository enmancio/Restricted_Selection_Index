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
![Scenario1](https://github.com/enmancio/Restricted_Selection_Index/files/9836803/selection.uno-1.pdf?raw=true "Optional Title")



### Scenario 2.

```r 
restricted_traits = c(4,5)
```

we use the **res_sel_I** built-in function to restrict the economic weights of traits 4 and 5:

```r 
# create first generation
pop = newPop(founderPop, simParam=SP)
# get genetic and phenotypic variance
G = varG(pop)
P = varP(pop)

# restricted selection index (built-in function created from Appendix 1)
# that permit to obtain the restricted economic weights and 
# restricted b vector
restriction=res_sel_I(econWt,G,P,tr=restricted_traits)
b_r=restriction$res_b
genMean = data.frame(t(meanG(pop)))  # mean of the breeding values at present generation
names(genMean) = paste0("traits",1:n_traits)


for(generation in 1:4){
  # select the best females as candidate dams
  damsc = selectInd(pop,nInd=n_candidate_female, trait=selIndex,sex="F", b=b_r)
  # select the best males as candidate sires
  sirec = selectInd(pop,nInd=n_candidate_male, trait=selIndex,sex="M", b=b_r)
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

### Scenario 3.

```r 
sex_lim_traits  = c(2,5)
```

here we use same selection differential (b) for both males and females

```r

pop = newPop(founderPop, simParam=SP)


pop=remove_phenotype_male(pop,sex_lim_traits) # phenotypes 2 and 5 as sex limited traits
head(pop@pheno) # shows the missing first 6 animals

P=varP(pop)
G=varG(pop)

# restricted selection index (built-in function created from code in Appendix 1)
# that function permit to obtain the restricted economic weights and thus
# restricted b vector
restriction=res_sel_I(econWt,G,P,tr=restricted_traits )
b_r=restriction$res_b

# mean of the breeding values at present generation
genMean = data.frame(t(meanG(pop)))  
names(genMean) = paste0("traits",1:n_traits)


for(generation in 1:4){
  cat("generation ",generation,"\n")
  # select the best females as candidate dams
  damsc = selectInd(pop,sex="F",nInd=n_candidate_female, trait=selIndex, b=b_r)
  # select the best males as candidate sire
  sirec = selectInd(pop,sex="M",nInd=n_candidate_male, trait=selIndex, b=b_r)
  # create population of candidates to selection
  candidate=mergePops(list(damsc,sirec))
  # create next generation by random crossing sires and dams
  pop=randCross(candidate,nCrosses=pop_size)
  # remove sex-limited traits
  pop=remove_phenotype_male(pop,sex_lim_traits)
  # collect genetic means for each trait 
  genMean = rbind(genMean,meanG(pop))
}

# plot genetic progress per generation
plot_g_prog(genMean)
```

### Scenario 4. 

```r
# create first generation
pop = newPop(founderPop, simParam=SP)

# restricted selection index (built-in function created from Appendix 1)
# that permit to obtain the restricted economic weight and 
# restricted b vector

pop=remove_phenotype_male(pop,sex_lim_traits)
P=varP(pop)
G=varG(pop)

G_m = varG(pop[pop@sex=="M"])
P_m = varP(pop[pop@sex=="M"])
G_f = varG(pop[pop@sex=="F"])
P_f = varP(pop[pop@sex=="F"])


genMean = data.frame(t(meanG(pop)))  # mean of the breeding values at present generation
names(genMean) = paste0("traits",1:n_traits)

#calculate separate economic weights for the females and  males
restriction_m=res_sel_I(econWt,G_m,P_m,tr=restricted_traits )
br_m=restriction_m$res_b

restriction_f=res_sel_I(econWt,G_f,P_f,tr=restricted_traits )
br_f=restriction_f$res_b


for(generation in 1:4){
  cat("generation ",generation,"\n")
  # select the best females as candidate dams
  damsc = selectInd(pop,sex="F",nInd=n_candidate_female, trait=selIndex, b=br_f)
  # select the best males as candidate sires
  sirec = selectInd(pop,sex="M",nInd=n_candidate_male, trait=selIndex, b=br_m)
  # create population of candidates to selection
  candidate=mergePops(list(damsc,sirec))
  # create next generation by random crossing sires and dams
 pop=randCross(candidate,nCrosses=pop_size)
  # remove sex-limited traits
  pop=remove_phenotype_male(pop,sex_lim_traits)
  # collect genetic means for each trait  
  genMean = rbind(genMean,meanG(pop))
}

plot_g_prog(genMean )
```

### Scenario 5
For computational reasons we have reduced the size of the population and due to the need to accurately 
estimate the genetic value of the animals (deep pedigree), before selection, we perform 3 generations of 
random mating.
The intensity of the selection was increased for better visual evaluation


```r 

pop_size = 1300
n_candidate_male = as.integer(pop_size/50)
n_candidate_female = as.integer(pop_size/4)

founderPop = quickHaplo(nInd=pop_size, nChr=10, segSites=1000)

# create founder haplotype
pop = newPop(founderPop, simParam=SP)

# generate pedigree for 3 generation of random mating
for (i in 1:3){
    candidate = selectInd(pop,nInd=pop@nInd*1/5, use = "rand",sex = "B")
    new_p =randCross(candidate,nCrosses=pop_size/4)
    pop = mergePops(list(pop,new_p))
    pop=remove_phenotype_male(pop,sex_lim_traits)
}

#  restricted economic weights
a=res_sel_I(econWt,varG(pop),varP(pop),c(4,5))[["res_a"]]

genMean=as.data.frame(t(meanG(pop)))

for(generation in 1:3){
  cat("generation ",generation,"\n")
  ebv = get_ebv(pop) # calculate breeding values with BLUP
  # select the best females as candidate dams by combing a with EBV
  damsc = my_sel_id(pop,SEX="M",n_candidate_female/2,a,ebv)
  # select the best males as candidate sire by combing a with EBV
  sirec = my_sel_id(pop,SEX="F",n_candidate_male/2,a,ebv)
  # create population of candidate to selection 
  candidate=mergePops(list(damsc,sirec))
  # generate next population by random mating
  new_p=randCross(candidate,nCrosses=pop_size)
  # combine the populations to trace-back pedigree infomation
  pop = mergePops(list(pop,new_p))
  pop=remove_phenotype_male(pop,sex_lim_traits)
  # create next generation by random crossing sire and dams
  genMean = rbind(genMean,meanG(pop))
}

plot_g_prog(genMean )


```

### Scenario 6

```r 
#  restricted economic weights
pop = newPop(founderPop, simParam=SP)


# generate pedigree by 3 generation of random mating
for (i in 1:4){
    candidate = selectInd(pop,nInd=pop@nInd*1/5, use = "rand",sex = "B")
    new_p =randCross(candidate,nCrosses=pop_size/4)
    pop = mergePops(list(pop,new_p))
    pop=remove_phenotype_male(pop,sex_lim_traits)
}

genMean=as.data.frame(t(meanG(pop)))

for(generation in 1:3){
  cat("generation ",generation,"\n")
# built-in function to solve 5 traits mixed model and obtained  
# the Caa values in the correct order
  MME = estimateg3(pop)
  G_hat = MME[["G_hat"]] # obtain Ghat = G - Caa
  EBV=MME[["SOL"]] # get the EBV
  G_hat = as.matrix(G_hat)
  # ranking animals based on the restricted economic weights
  af=res_weigth_lin(pop,econWt,G_hat,restricted_traits)
  # select the best females as candidate dams by combing a with EBV
  damsc = my_sel_id2(pop,SEX="F",n_candidate_female/2,af,EBV)
  # select the best males as candidate dams by combing a with EBV
  sirec = my_sel_id2(pop,SEX="M",n_candidate_male/2,af,EBV)
  # create population of candidates to selection 
  candidate=mergePops(list(damsc,sirec))
  # generate next population by random mating
  new_p=randCross(candidate,nCrosses=pop_size)
  #merge previous and new population to trace-back pedigree information  
  pop = mergePops(list(pop,new_p))
  pop=remove_phenotype_male(pop,sex_lim_traits)
  genMean = rbind(genMean,meanG(pop))
}


plot_g_prog(genMean )

```

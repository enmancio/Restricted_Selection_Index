# Restricted Selection Index


## Upload or download all packages and built-in functions:

```r
set.seed(1234)
source("function/JDS_function.R") # load all built-in function
pkg = c("AlphaSimR","MASS","ggplot2","pedigreemm","sparseinv","Rcpp")
libraries(pkg)
```

### Use AlphaSimR to simulate the first population

In this example, we mimic a population with similar parameters to the one described in [Mancin et al. (2022)](https://reader.elsevier.com/reader/sd/pii/S0022030222006130?token=D918083579920C8D0B9ADAC59C8D893F94465DE9A79C795EB679FFED6584E63371702C17D9EE5F09BC84C52D439C23B1&originRegion=eu-west-1&originCreation=20221026133930). 
Note that the genetic variance and phenotype mean have been changed for a matter of clarity and a better visual appraisal of the graphs. Indeed, here we assumed a genetic variance equal to 1 and a phenotype mean equal to 0.

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
![Scenario1](https://user-images.githubusercontent.com/93042853/197342581-375f4b10-43c5-468b-9c48-f5c6290f0780.png)


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

![Scenario2](https://user-images.githubusercontent.com/93042853/197343046-81449382-fb64-4c0b-b016-4a5fe6261c52.png)

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

![Scenario3](https://user-images.githubusercontent.com/93042853/197343040-25a9e9a8-b6e6-467a-8233-30b403a0373d.png)

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

![Scenario4](https://user-images.githubusercontent.com/93042853/197343043-bfb347fc-0ecf-430d-82cb-8f51c48e3ad4.png)

### Scenario 5
For computational reasons we have reduced the size of the population and due to the need to accurately 
estimate the genetic value of the animals (deep pedigree), before selection, we perform 3 generations of 
random mating.
The intensity of the selection was increased for better visual evaluation


```r 

pop_size = 1700
n_candidate_male = as.integer(pop_size/50)
n_candidate_female = as.integer(pop_size/4)

founderPop = quickHaplo(nInd=pop_size, nChr=10, segSites=1000)

# create founder haplotype
pop = newPop(founderPop, simParam=SP)

# generate pedigree for 3 generation of random mating
for (i in 1:3){
    candidate = selectInd(pop,nInd=pop@nInd*1/2, use = "rand",sex = "B")
    new_p =randCross(candidate,nCrosses=pop_size/2)
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
  damsc = my_sel_id(pop,SEX="M",n_candidate_female,a,ebv)
  # select the best males as candidate sire by combing a with EBV
  sirec = my_sel_id(pop,SEX="F",n_candidate_male,a,ebv)
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
![Scenario5](https://user-images.githubusercontent.com/93042853/197704265-04b589cb-ff6a-41ea-a21b-c2ed5e7e20b6.png)


### Scenario 6

```r 
#  restricted economic weights
pop = newPop(founderPop, simParam=SP)


# generate pedigree by 3 generation of random mating
for (i in 1:3){
    candidate = selectInd(pop,nInd=pop@nInd*1/2, use = "rand",sex = "B")
    new_p =randCross(candidate,nCrosses=pop_size/2)
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
  damsc = my_sel_id2(pop,SEX="F",n_candidate_female,af,EBV)
  # select the best males as candidate dams by combing a with EBV
  sirec = my_sel_id2(pop,SEX="M",n_candidate_male,af,EBV)
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
![Scenario6](https://user-images.githubusercontent.com/93042853/197704272-20360d7b-1781-4aa7-82b1-2a4be93c7af0.png)

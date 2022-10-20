
require("sparseinv")

# aggiungere locazione

corgen = function(size=3){
#  set.seed(1234) 
  A <- matrix(runif(size^2)*2-1, ncol=size) 
  return(t(A) %*% A)
}

get_ebv = function(pop){
    VARG=varG(pop)
    VARP=varP(pop)
    VARE=VARP-VARG
    n = pop@nInd
    p = pop@nTraits

    ped=data.frame("id"=pop@iid,"fid"=pop@father,"mid"=pop@mother,"u"=1,
                  "P"=pop@pheno[,],"sex"=pop@sex)

    #head(ped)
    #tail(ped)
    #A_noi=AGHmatrix::Amatrix(ped[,c(1:3)])
    library("pedigreemm")
    ped2 <- with(ped, pedigree(sire=fid, dam=mid, label=id))
    U <- relfactor(ped2)
    A_noi <- crossprod(U)
    A=solve(A_noi)
    e=solve(VARE)
    VARGi= solve(VARG)

    require(Matrix)
    X <- Matrix(model.matrix(~ u-1, data=ped))
    X=as(X, "dgTMatrix")
    VARGi=as(VARGi, "dgTMatrix")
    e=as(e, "dgTMatrix")

    Z1 = Matrix(diag(rep(1,pop@nInd)))
    #Z1=as(Z1, "dgTMatrix")
    Z3 = Z1
    Z4 = Z1
    Z2 = Matrix(model.matrix(~ sex-1, data=ped))
    Z2= diag(Z2[,1])
    Z2= as(Z2, "dgTMatrix")
    Z5=Z2

    Rcpp::sourceCpp("mme_simplify.cpp")
    TMP = fast_five_traitsMM(A,Z1,Z2,VARGi,e)
    #TMP[TMP<1e-2]=0
    
    #Takahashi equations
    if (isSymmetric(TMP)) {
    TMP=sparseinv::Takahashi_Davis(TMP)
    } else {
      break
    }
    TMP=as.matrix(TMP)
    pheno = as.matrix(ped[,5:9])
    sol = SOLVEMME(pheno,as.matrix(e),as.matrix(Z1),as.matrix(Z2),TMP)


    print("re-order the solutions traits and then animals")
    sol = data.frame("val"=c(as.matrix(sol)))
    l=nrow(sol)/p
    sol_fix = data.frame("id"=ped$id)
    std = 1:n
    dove=std
      for (i in 1:p){
            sol_fix[,paste0("tr_",i)]   = sol[dove,1]
            dove = std + max(dove)
      }

    ok = as.matrix(sol_fix)[,-1]
    ok = as.matrix(sol_fix)[,-1]
return( sol_fix)
}



my_sel_id = function(pop,SEX,n,a,ebv) {

  ebv=as.data.frame(ebv)
  ebv = ebv[,paste0("tr_",1:pop@nTraits)]
  rank = as.matrix(ebv)%*%(a) 
  pop@ebv = rank
  mask=pop@sex==SEX
  db=data.frame(id=pop@id,"aw"=pop@ebv)
  db = db[mask,]  
  top=tail(db[order(db$aw),],n)
  pop2=pop[pop@id %in% top$id]
  return(pop2)
}

my_sel_id2 = function(pop,SEX,n,a,EBV) {
  ebv=EBV
#  res_weigth_lin = function(a,G,restricted_traits)
  ebv=as.data.frame(ebv)
  ebv = ebv[,paste0("tr_",1:pop@nTraits)]
  tot_rank=NULL
  for (i in 1:pop@nInd){
  rank = as.matrix(ebv[i,])%*%(af[i,]) 
  tot_rank = c(tot_rank,rank)
}
  #print(ebv[,])
  pop@ebv = as.matrix(tot_rank)
  mask=pop@sex==SEX
  db=data.frame(id=pop@id,"aw"=pop@ebv)
  db = db[mask,]  
  top=tail(db[order(db$aw),],n)
  pop2=pop[pop@id %in% top$id]
  return(pop2)
}

res_weigth_lin = function(pop,a,G,restricted_traits) {
      i=0
      x=1
      af = matrix(NA,nrow=pop@nInd,ncol=pop@nTraits)
      while(x != pop@nInd+1) {
      A1 = G[i+(1:pop@nTraits),i+(1:pop@nTraits)]
      C = A1[restricted_traits,]
      C = t(C)
      A1i = MASS::ginv(A1)
      K=C%*%(t(C)%*% A1i %*%C)%*%t(C)
      I = diag(rep(1,ncol(A1)))
      a_S= (I - A1i%*%K)%*%a
      af[x,] <- as.matrix(a_S,ncol=pop@nTraits)
      #print(a_S)
      i=i+5
      x=x+1
      }
      return(af)
      }



res_weigth_lin_TBV = function(pop,a,G,restricted_traits) {
      A1 = G[i+(1:5),i+(1:5)]
      C = A1[restricted_traits,]
      C = t(C)
      K=C%*%(t(C)%*%solve(A1)%*%C)%*%t(C)
      I = diag(rep(1,ncol(A1)))
      a_S= (I - solve(A1)%*%K)%*%a
      return(a_S)
      }




remove_phenotype_male = function(pop,traits) {
  mask = pop@sex == "M"
  pop@pheno[mask,traits] <-0
  return(pop)
}



res_sel_I = function(econom_wt,G,P,tr) {
  a = econWt
  ds_G = diag(G)
  a_std = a/(ds_G)
  Pinv = MASS::ginv(P)
  b = Pinv %*% G %*% (a_std) 
  ds_I = sqrt(t(b) %*% P %*% b)
  gp = (1/ds_I) %*% t(b) %*% G
  C = G[tr,]
  C = t(C)
  I = diag(rep(1,ncol(G)))
  K = I - (Pinv %*% ((C) %*% MASS::ginv(t(C) %*% Pinv %*% (C)) %*% t(C)))
  br = K %*% Pinv %*% G %*% (a_std) # rturn qua 
  
  gpr = ((1/(ds_I))%*%t(br)%*%G)
  new_a = solve(G)%*%(P)%*%(br)
  new_a_std = c(new_a)/(sum(abs(new_a)))
  #cat(paste(paste0("traits"),"=","G.P.","G.P. res", sep="\t"),sep="\n")
  #cat(paste(paste0("tr", 1:nrow(G)),"=",round(gp,2),round(gpr,2), sep="\t"),sep="\n")
  out =list()
  out[["res_a"]] = new_a
  out[["res_b"]] = br
  out[["expected_genetic_progress"]] = gpr
  return(out)
}



select_sire = function(pop,n=100,b) {
  n_traits = ncol(pop@pheno)
  dat = data.frame("iid"=pop@iid,"id"=pop@id,
                   "sex" = pop@sex)
  dat[,paste0("pheno",1:n_traits)] =pop@pheno
  
  dat_m = dat[dat$sex == "M",]
  dat_m$totalmert= as.matrix(dat_m[,paste0("pheno",1:n_traits )]) %*% b
  top_sire = dat_m[order(dat_m$totalmert*-1),][1:n,]
  nrow(dat_m[dat_m$id %in% as.character(top_sire$id),])
  sire_pop = pop[pop@id %in% as.character(top_sire$id),]
  return(sire_pop)
}



libraries <- function(x){
  for( i in x ){
    
    if( ! require( i , character.only = TRUE ) ){
      
      install.packages( i , dependencies = TRUE )
      
      require( i , character.only = TRUE )
    }
  }
}

plot_g_prog = function(genMean) {
  genMean$generation = rownames(genMean)
  names(genMean) = c(paste0("trait",1:5),"generation")
  to_plot = reshape2::melt(genMean)
  plotg=ggplot(to_plot,aes(x=as.integer(generation),y=value,group=variable,color=variable))+
    geom_line(size=1.5) +
    geom_point(size=3,color="black")+
    geom_point(size=2)+
    theme_light()+
    xlab("generation")+
    ylab("genetic value")+
    
    scale_color_brewer(name = "Traits:", palette = "Set1")+
    theme(title =element_text(size=16),
          legend.position="bottom",
          #plot.margin = margin(1, 1, 1, 1, "cm"),
          legend.text = element_text(colour="black",size = 10,color='black'),
          axis.text.y = element_text(size=13, color='black'),
          axis.text.x = element_text(size=13, color='black'),
          axis.title.x = element_text(size=13,color='black'),
          axis.title.y = element_text(size=13,color='black'),
          axis.line.x = element_line(colour = "black", size = 1),
          axis.line.y = element_line(colour = "black", size = 1),
          axis.ticks = element_line(colour = "black", size = 1))
  return(plotg)}




estimateg3=function(pop,VARG,VARE){
        VARG=varG(pop)
        VARP=varP(pop)
        VARE=VARP-VARG

        ped=data.frame("id"=pop@iid,"fid"=pop@father,"mid"=pop@mother,"u"=1,
                      "P"=pop@pheno[,],"sex"=pop@sex)

        #A_noi=AGHmatrix::Amatrix(ped[,c(1:3)])
        library("pedigreemm")
        ped2 <- with(ped, pedigree(sire=fid, dam=mid, label=id))
        U <- relfactor(ped2)
        A_noi <- crossprod(U)
        A=solve(A_noi)

        e=solve(VARE)
        VARGi= solve(VARG)

        require(Matrix)
        X <- Matrix(model.matrix(~ u-1, data=ped))
        X=as(X, "dgTMatrix")
        VARGi=as(VARGi, "dgTMatrix")
        e=as(e, "dgTMatrix")

        Z1 = Matrix(diag(rep(1,pop@nInd)))
        Z3 = Z1
        Z4 = Z1
        Z2 = Matrix(model.matrix(~ sex-1, data=ped))
        Z2= diag(Z2[,1])
        Z2= as(Z2, "dgTMatrix")
        Z5=Z2

        Rcpp::sourceCpp("mme_simplify.cpp")
        TMP = fast_five_traitsMM(A,Z1,Z2,VARGi,e)
        print("inversion")
        TMP=solve(TMP)
        TMP=as.matrix(TMP)
        pheno = as.matrix(ped[,5:9])
        sol = SOLVEMME(pheno,as.matrix(e),as.matrix(Z1),as.matrix(Z2),TMP)

        n = pop@nInd
        p = pop@nTraits

        dim(TMP)
        A_LIST=re_orderCaa(TMP,p,n)

        print("order done")

        G= A_noi %x% varG(pop)

        PEV=G
        PEV[,] <- 0
        i=1
        std = 1:p
        dove=std
        while(i < n) {
          PEV[dove,dove] <- A_LIST[[i]] #*VARE
        # print(dove)
          dove = std + max(dove)
        # print(dove)
          i=i+1
        }

        G_hat = G-(PEV)

        #if(F){
        print("re-order the solutions")
        sol = data.frame("val"=c(as.matrix(sol)))
        l=nrow(sol)/p
        sol_fix = data.frame("id"=ped$id)
        std = 1:n
        dove=std
          for (i in 1:p){
                sol_fix[,paste0("tr_",i)]   = sol[dove,1]
                dove = std + max(dove)
          }
        #}

        ii=list(G_hat=G_hat,A_LIST=A_LIST,SOL=sol_fix)

        return(ii)
}






#Rcpp::sourceCpp("mme_simplify.cpp")
#TMP = fast_five_traitsMM(A,Z1,Z2,VARGi,e)

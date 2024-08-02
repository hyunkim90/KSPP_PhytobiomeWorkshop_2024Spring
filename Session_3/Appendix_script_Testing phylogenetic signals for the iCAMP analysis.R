###### Phylogenetic signal test when the env object is not available ######
#### 1. Install R packages ####
install.packages("ape")
install.packages("picante")
install.packages("phytools")

#### 2. Load R packages ####
library(ape)
library(picante)
library(phytools)

#### 3. Prepare data ####
set.seed(1234) ## for reproducibility
#### 3-1. Create continuous trait data (use this if you do not have any trait data) ####
trait_data <- ape::rTraitCont(tree, model="BM", sigma=sqrt(2)) ### The tree object: phylogenetic tree of ASVs or OTUs
names(trait_data) <- tree$tip.label
K_continuous <- phytools::phylosig(tree, trait_data, method="K")
print(K_continuous)

#### 3-2. When considering mean abundances of ASVs as a continuous trait ####
trait_data<-apply(comm,2,mean) ### The variable comm is the ASV or OTU abundance table data.
names(trait_data) <- tree$tip.label
K_continuous <- phytools::phylosig(tree, trait_data, method="K")
print(K_continuous)

#### 4. Test phylogenetic signal (Blomberg's K) for each bin size ####
# If you want to see how different bin sizes (number of bins) affect the results,
# you could loop through different numbers of bins and calculate Blomberg's K for each
K_binned_results <- sapply(2:40, function(nbins) {
  binned_data <- cut(trait_data, breaks=quantile(trait_data, probs=seq(0, 1, length.out=nbins+1)), labels=FALSE, include.lowest=TRUE)
  phytools::phylosig(tree, binned_data, method="K")
})

### In the sapply(2:40, function(x)...), the range from 2 to 40 (2:40) indicates the size of phylogenetic bins that will be tested.
### The range of the tested bin sizes could be adjusted depending on your data.

#### 5. Print the results for different bin sizes ####
print(K_binned_results)
test.result <- data.frame(bin.size = c(2:40), physig.value = K_binned_results)

### 6. Choose the bin size which K value is the highest and apply it to the iCAMP analysis ####
test.result$bin.size[which(test.result$physig.value == max(test.result$physig.value))]


############### Determine the minimum bin size using pNST

##### Install package NST if not yet
if(!("NST" %in% installed.packages()[,"Package"])){install.packages("NST")}
library(NST)

##### Set a grouping value
i=1
treat.use=treat[,i,drop=FALSE]

##### Calculate stochasticity using pNST
pnstout=NST::pNST(comm=comm, pd.desc=pd.big$pd.file, pd.wd=pd.big$pd.wd, 
                  pd.spname=pd.big$tip.label, group=treat.use, abundance.weighted=TRUE,
                  rand=rand.time, phylo.shuffle=TRUE, nworker=nworker,
                  output.rand = TRUE, SES=FALSE, RC=FALSE)

write.csv(pnstout$index.grp,file = paste0(prefix,".pNST.summary.",colnames(treat)[i],".csv"))
write.csv(pnstout$index.pair.grp,file = paste0(prefix,".pNST.pairwise.",colnames(treat)[i],".csv"))

pnst.bt=NST::nst.boot(nst.result=pnstout, group=treat.use,
                      rand=rand.time, nworker=nworker)
write.csv(pnst.bt$summary,file = paste0(prefix,".pNST.bootstr.",colnames(treat)[i],".csv"))
write.csv(pnst.bt$compare,file = paste0(prefix,".pNST.compare.",colnames(treat)[i],".csv"))


pnst.res<-pnst.bt$summary
pnst.res<-pnst.res[pnst.res$Index == "NST",]
print(mean(pnst.res$mean))


####### Calculate iCAMP-based stochasticity with different bin size value
bin.size.test <- data.frame(bin.size=c(12:48,"pNST"), stochasticity=0)
bin.size.test$stochasticity[bin.size.test$bin.size == "pNST"] <- mean(pnst.res$mean)
for (i in 12:48){
  bin.size.limit = i 
  sig.index="Confidence" 
  icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                         prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                         phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                         phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                         nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                         qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                         correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                         ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
  
  icbin <- icamp.bins(icamp.detail = icres$detail,treat = treat.use,
                      clas=clas,silent=FALSE, boot = TRUE,
                      rand.time = rand.time,between.group = TRUE)
  
  icres.bintest<-icbin$Pt
  icres.bintest <- icres.bintest[!(grepl("vs",icres.bintest$Group)),]
  j<-c(6:8)
  icres.bintest[, j] <- apply(icres.bintest[ , j], 2,function(x) as.numeric(as.character(x)))
  icres.bintest$stochasticity <- rowSums(icres.bintest[c(6:8)])
  
  bin.size.test$stochasticity[bin.size.test$bin.size == i] <- mean(icres.bintest$stochasticity)
  
  print(i)
  print(bin.size.test$stochasticity[bin.size.test$bin.size == i])
}


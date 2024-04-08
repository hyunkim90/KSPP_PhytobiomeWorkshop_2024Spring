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

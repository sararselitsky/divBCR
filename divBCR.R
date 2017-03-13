library(plyr)
args <- commandArgs(TRUE)
igh <- read.delim(args[1], sep="\t",header=TRUE)

# step 1, BCR Counts are normalized to total RNA-seq counts #
countsForMed<-unique(data.frame(igh$sample,igh$total_count))
igh$normFactor<-igh$total_count/median(countsForMed$igh.total_count)
igh$normed_counts<-igh$normFactor*igh$expected_counts

# step 2, Calculations of biodiversity and other metrics at the sequence level#
chararcterizeBCR=function(counts,sampleList,sequence,cdr3,vregion){
  char=matrix(NA,nrow=length(levels(sampleList)),ncol=12)
  i=0
  for(sample in levels(sampleList)){
    i=i+1
    filt<-which(sample==sampleList)
    # total counts
    char[i,1]= sum(counts[filt])
    # avg cdr3 length
    char[i,2] = mean(nchar(as.character(cdr3[filt])))
    # weighted avg cdr3 length
    char[i,3] = sum(nchar(as.character(cdr3[filt]))*counts[filt])/sum(counts[filt])
    # mean v-region identity
    char[i,4] = mean(vregion[filt])
    # weighted mean v-region identity
    char[i,5] = sum(vregion[filt]*counts[filt])/sum(counts[filt])
    # evenness by sequence
    char[i,6] = -sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt])) ) /log(length(sequence[filt]))
    # shannon entropy by sequence
    char[i,7] = -sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt])))
    # true diversity by sequence
    char[i,8] = exp(-sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt]))))
    # gini simpson by sequence
    char[i,9] = 1-(sum((counts[filt]/sum(counts[filt]))^2))
    # top clone prop by sequence
    char[i,10] = max(counts[filt])/sum(counts[filt])
    # second top clone by sequence
    if(length(sequence[filt])>1){
      char[i,11] = (sort(counts[filt],partial=length(counts[filt])-1)[length(counts[filt])-1])/sum(counts[filt])
    }
    # number of contigs by cluster
    char[i,12] = length(sequence[filt])
  }
  char<-as.data.frame(char)
  rownames(char)<-levels(sampleList)
  return(char)
}

# step 3, sums of clusters #
sumCountsCluster<-function(data)
{c(clusterSums = with(data, sum(normed_counts) ))}
sumCountsCluster <-ddply(igh, .variables=c("sample","cluster"), .fun=sumCountsCluster)

# step 4, calculations of biodiversity at the cluster level #
diversityBCR=function(counts,samplesList,clusters){
  div=matrix(NA,nrow=length(levels(samplesList)),ncol=7)
  i=0
  for(sample in levels(samplesList)){
    i=i+1
    filt<-which(sample==samplesList)
    # evenness
    div[i,1] = -sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt])) ) /log(length(clusters[filt]))
    # shannon entropy
    div[i,2] = -sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt])))
    # true diversity
    div[i,3] = exp(-sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt]))))
    # gini simpson
    div[i,4] = 1-(sum((counts[filt]/sum(counts[filt]))^2))
    # top clone prop
    div[i,5] = max(counts[filt])/sum(counts[filt])
    # second top clone
    if(length(clusters[filt])>1){
      div[i,6] = (sort(counts[filt],partial=length(counts[filt])-1)[length(counts[filt])-1])/sum(counts[filt])
    }
    # number of contigs
    div[i,7] = length(clusters[filt])
  }
  div<-as.data.frame(div)
  rownames(div)<-levels(samplesList)
  return(div)
}

# step 5, executing the functions, merging the two data frames and saving the output to a file. #
char<-chararcterizeBCR(counts = igh$normed_counts,sampleList = igh$sample,sequence= igh$sequence,cdr3 = igh$aa_cdr3,vregion = igh$vregion_identity)
div<-diversityBCR(counts = sumCountsCluster$clusterSums,clusters = sumCountsCluster$cluster,samplesList = sumCountsCluster$sample)
char<-merge(char,div,by.x = "row.names",by.y = "row.names")
colnames(char)<-c("sample","total_counts","mean_cdr3_aa_len","mean_cdr3_aa_len_w","mean_v_region_ident","mean_v_region_ident_w","evenness_by_seq","shannon_entropy_by_seq","true_div_by_seq","gini_simpson_by_seq","top_clone_prop_by_seq","second_top_clone_by_seq","num_unique_seq","evenness_by_clus","shannon_entropy_by_clus","true_div_by_clus","gini_simpson_by_clus","top_clone_prop_by_clus","second_top_clone_by_clus","num_unique_clus")
write.table(char,args[2],sep="\t",quote = FALSE,row.names = FALSE)


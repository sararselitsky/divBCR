library(plyr)
args <- commandArgs(TRUE)
igh <- read.delim(args[1], sep="\t",header=TRUE)

# step 1, BCR Counts are normalized to total RNA-seq counts #
countsForMed<-unique(data.frame(igh$sample,igh$total_count))
igh$normFactor<-igh$total_count/median(countsForMed$igh.total_count)
igh$normed_counts<-igh$normFactor*igh$expected_counts

# step 2, calculations of biodiversity and other metrics at the sequence level for each sample, for each isotype#
chararcterizeBCR=function(counts,sampleList,sequence,cdr3,vregion,isotype){
  char=matrix(NA,nrow=length(unique(paste(as.factor(sampleList),as.factor(isotype)))),ncol=15)
  i=0
  for(sample in levels(sampleList)){
    for(iso in levels(isotype)){
      filt_sam<-which(sample==sampleList)
      filt_iso<-which(iso==isotype)
      filt<-intersect(filt_sam,filt_iso)
      if(length(filt)>0){
        i=i+1
        sampleIso=paste(sample,iso)
        char[i,1] = sampleIso
        char[i,2] = sample
        char[i,3] = iso
        # total counts
        char[i,4]= sum(counts[filt])
        # avg cdr3 length
        char[i,5] = mean(nchar(as.character(cdr3[filt])))
        # weighted avg cdr3 length
        char[i,6] = sum(nchar(as.character(cdr3[filt]))*counts[filt])/sum(counts[filt])
        # mean v-region identity
        char[i,7] = mean(vregion[filt])
        # weighted mean v-region identity
        char[i,8] = sum(vregion[filt]*counts[filt])/sum(counts[filt])
        # evenness by sequence
        char[i,9] = -sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt])) ) /log(length(sequence[filt]))
        # shannon entropy by sequence
        char[i,10] = -sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt])))
        # true diversity by sequence
        char[i,11] = exp(-sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt]))))
        # gini simpson by sequence
        char[i,12] = 1-(sum((counts[filt]/sum(counts[filt]))^2))
        # top clone prop by sequence
        char[i,13] = max(counts[filt])/sum(counts[filt])
        # second top clone by sequence
        if(length(sequence[filt])>1){
          char[i,14] = (sort(counts[filt],partial=length(counts[filt])-1)[length(counts[filt])-1])/sum(counts[filt])
        }
        # number of contigs by cluster
        char[i,15] = length(sequence[filt])
      }
    }
  }
  char<-as.data.frame(char)
  rownames(char)<-unique(paste(as.factor(sampleList),as.factor(isotype)))
  colnames(char)<-c("sampleIso","sample","isotype","total_counts","mean_cdr3_aa_len","mean_cdr3_aa_len_w","mean_v_region_ident","mean_v_region_ident_w","evenness_by_seq","shannon_entropy_by_seq","true_div_by_seq","gini_simpson_by_seq","top_clone_prop_by_seq","second_top_clone_by_seq","num_unique_seq")
  return(char)
}

# step 3, sums of clusters #
sumCountsCluster<-function(data)
{c(clusterSums = with(data, sum(normed_counts) ))}
sumCountsCluster <-ddply(igh, .variables=c("sample","cluster","isotype"), .fun=sumCountsCluster)

# step 4, calculations of biodiversity at the cluster level,for each sample, for each isotype #
diversityBCR=function(counts,sampleList,clusters,isotype){
  div=matrix(NA,nrow=length(unique(paste(as.factor(sampleList),as.factor(isotype)))),ncol=8)
  i=0
  for(sample in levels(sampleList)){
    for(iso in levels(isotype)){
      filt_sam<-which(sample==sampleList)
      filt_iso<-which(iso==isotype)
      filt<-intersect(filt_sam,filt_iso)
      if(length(filt)>0){
        i=i+1
        sampleIso=paste(sample,iso)
        div[i,1] = sampleIso
        # evenness
        div[i,2] = -sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt])) ) /log(length(clusters[filt]))
        # shannon entropy
        div[i,3] = -sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt])))
        # true diversity
        div[i,4] = exp(-sum( (counts[filt]/sum(counts[filt])) * log(counts[filt]/sum(counts[filt]))))
        # gini simpson
        div[i,5] = 1-(sum((counts[filt]/sum(counts[filt]))^2))
        # top clone prop
        div[i,6] = max(counts[filt])/sum(counts[filt])
        # second top clone
        if(length(clusters[filt])>1){
          div[i,7] = (sort(counts[filt],partial=length(counts[filt])-1)[length(counts[filt])-1])/sum(counts[filt])
        }
        # number of contigs
        div[i,8] = length(clusters[filt])
      }
    }
  }
  div<-as.data.frame(div)
  rownames(div)<-unique(paste(as.factor(sampleList),as.factor(isotype)))
  colnames(div)<-c("sampleIso","evenness_by_clus","shannon_entropy_by_clus","true_div_by_clus","gini_simpson_by_clus","top_clone_prop_by_clus","second_top_clone_by_clus","num_unique_clus")
  return(div)
}

# step 5, executing the functions, merging the two data frames and saving the output to a file. #
char<-chararcterizeBCR(counts = igh$normed_counts,sampleList = igh$sample,sequence= igh$sequence,cdr3 = igh$aa_cdr3,vregion = igh$vregion_identity,isotype = igh$isotype)
div<-diversityBCR(counts = sumCountsCluster$clusterSums,clusters = sumCountsCluster$cluster,sampleList = sumCountsCluster$sample,isotype = sumCountsCluster$isotype)
char<-merge(char,div,by.x = "sampleIso",by.y = "sampleIso")
char$sampleIso<-NULL
write.table(char,args[2],sep="\t",quote = FALSE,row.names = FALSE)


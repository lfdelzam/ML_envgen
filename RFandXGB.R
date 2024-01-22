rm(list = ls())
options(warn=-1)




#1. Installing packages

if ( ! "argparser" %in% installed.packages()[,"Package"]) {
  installed.packages("argparser")
}

if ( ! "DiagrammeRsvg" %in% installed.packages()[,"Package"]) {
devtools::install_github('rich-iannone/DiagrammeRsvg')
}


bioconduc_packages=c("ape")
new.packages <- bioconduc_packages[!(bioconduc_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  for (s in new.packages) { suppressPackageStartupMessages(BiocManager::install(s))}
}

#libraries
list_of_packages <- c("tidyverse","anytime","viridis","reshape2","vegan","ggpubr","ape","lubridate",
                      "dplyr","plyr","svglite","gridExtra","geosphere","corrplot","ggthemes",
                      "cowplot", "data.table", "mlr", "tidymodels", "randomForest", "xgboost", "DiagrammeR",
                      "Ckmeans.1d.dp","DiagrammeRsvg", "rsvg")

new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  for (s in new.packages) { suppressPackageStartupMessages(install.packages(s))}
}



suppressPackageStartupMessages(library("argparser"))

# arguments
p <- arg_parser("clean_dataset.R")
p <- add_argument(p, "-b", help="16S_ASV sequences", default="Raw_data/ASVs_16S.txt")
p <- add_argument(p, "-e", help="18S_ASV sequences", default="Raw_data/ASVs_18S.txt")
p <- add_argument(p, "-sb", help="16S_count table", default="Raw_data/seqtab_16S.txt")
p <- add_argument(p, "-se", help="18S_count table", default="Raw_data/seqtab_18S.txt")
p <- add_argument(p, "-m", help="Metadata", default="Raw_data/metadata.txt")
p <- add_argument(p, "-a", help="Abiotic parameters", default="contextual_data/physical_chemical_data_2019-2020_2022-03-16.txt")
p <- add_argument(p, "-tb", help="taxa_16S table", default="Raw_data/taxa_16S.txt")
p <- add_argument(p, "-te", help="taxa_18S table", default="Raw_data/taxa_18S.txt")
p <- add_argument(p, "-w", help="working directory", default="/Users/luisdelgado/Documents/VAE")
p <- add_argument(p, "-d", help="Minimum depth", default=10)
p <- add_argument(p, "-o", help="output directory", default="Test_xgb_mlr")



argv <- parse_args(p)

for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}


setwd(argv$w)

#### Anders code
# load sequence tables with counts
seqtab_16S = as.matrix(read.delim(argv$sb))
seqtab_18S = as.matrix(read.delim(argv$se))
# load asv tables
asv_16S = as.matrix(read.delim(argv$b))
asv_18S = as.matrix(read.delim(argv$e))
# load tax tables
taxa_16S = as.matrix(read.table(argv$tb))
taxa_18S = as.matrix(read.table(argv$te))


#load metadata file
metadata = read.delim(paste(argv$m))
#remove row starting with "#"
metadata = metadata[2:nrow(metadata),]

#load physical chemical data
psdata = read.delim(argv$a)
#remove values of less than 10m depth 
psdata = psdata[(which(psdata$sample_depth_m <= argv$d)),]
#merge stations (are very close)
psdata$station_name[psdata$station_name == 'B3'] = 'B7'
#sum salinity
ps_mean = aggregate(Salinity_CTD_o_oo_psu ~ station_name_alt, data = psdata, FUN = mean)
ps_sd = aggregate(Salinity_CTD_o_oo_psu ~ station_name_alt, data = psdata, FUN = sd)
colnames(ps_mean)[2] <- 'mean_salinity'
colnames(ps_sd)[2] <- 'sd_salinity'
ps_sal = merge(ps_mean, ps_sd)



# For 16S
metadata_16S = metadata
meta_ix_16S = match(colnames(seqtab_16S), metadata_16S$sample_alias.library_name)
metadata_16S = metadata_16S[meta_ix_16S,]
#metadata_16S$station_name[metadata_16S$station_name == 'B3'] = 'B7' # Merge stations B3 and B7 (practically one location)
## Remove special samples
# remove blanks from the dataset
ix_16S = sapply('blank', grep, metadata_16S$test_type, invert = TRUE)
# remove volumetest
ix_16S = setdiff(ix_16S, which(metadata_16S$test_type == 'volume_test'))
# Remove undersequenced samples
ix_16S = intersect(which(colSums(seqtab_16S) > 40000), ix_16S) 


# For 18S
metadata_18S = metadata
meta_ix_18S = match(colnames(seqtab_18S), metadata_18S$sample_alias.library_name)
metadata_18S = metadata_18S[meta_ix_18S,]
#metadata_18S$station_name[metadata_18S$station_name == 'B3'] = 'B7' # Merge stations B3 and B7 (practically one location)
## Remove special samples
# remove blanks from the dataset
ix_18S = sapply('blank', grep, metadata_18S$test_type, invert = TRUE)
# remove volume test
ix_18S = setdiff(ix_18S, which(metadata_18S$test_type == 'volume_test'))
# Remove undersequenced samples
#ix_18S = intersect(which(colSums(seqtab_18S) > 20000), ix_18S) 

# for 18S:
# remove samples not included in 16S
ix_18S = ix_18S[which(!is.na(match(metadata_18S$sample_title[ix_18S], metadata_16S$sample_title[ix_16S])))]
# choose just one of replicates from the same time and station
dates_stations = paste(metadata_18S$collection.date[ix_18S], metadata_18S$station_name[ix_18S])
ix_18S = ix_18S[match(unique(dates_stations), dates_stations)]

# for 16S again:
# remove samples not included in 18S
ix_16S = ix_16S[which(!is.na(match(metadata_16S$sample_title[ix_16S], metadata_18S$sample_title[ix_18S])))]
# order samples as for 18S
ix_16S = ix_16S[match(metadata_18S$sample_title[ix_18S], metadata_16S$sample_title[ix_16S])]

# select only relevant data for further use, 18S
metadata_18S = metadata_18S[ix_18S,]
seqtab_18S = seqtab_18S[,ix_18S]
sample_id_18S = metadata_18S$sample_title
colnames(seqtab_18S) = sample_id_18S
# remove ASVs with sum zero for remaining samples
iz_18S =which(rowSums(seqtab_18S) > 0) #seqtab_18S_temp = seqtab_18S[which(rowSums(seqtab_18S) > 0),]
asv_18S = asv_18S[iz_18S]
#seqtab_18S = seqtab_18S[iz_18S,]
taxa_18S = taxa_18S[iz_18S,]

# select only relevant data for further use, 16S
metadata_16S = metadata_16S[ix_16S,]
seqtab_16S = seqtab_16S[,ix_16S]
sample_id_16S = metadata_16S$sample_title
colnames(seqtab_16S) = sample_id_16S
# remove ASVs with sum zero for remaining samples
iz_16S =which(rowSums(seqtab_16S) > 0) #seqtab_16S_temp = seqtab_16S[which(rowSums(seqtab_16S) > 0),]
asv_16S = asv_16S[iz_16S]
#seqtab_16S = seqtab_16S[iz_16S,]
taxa_16S = taxa_16S[iz_16S,]

## Now samples come in the same order for 16S and 18S, so we can make joint sample_id and metadata
identical(metadata_16S$sample_title, metadata_18S$sample_title) # should be TRUE
metadata = metadata_16S

#check for duplicate samples
duplicated(metadata_18S$sample_title)

sample_id = metadata$sample_title
sample_station = metadata$station_name
sample_date = anydate(metadata$collection.date)
sample_date_text = str_remove(metadata$sample_title, '_.+_.+$')
sample_replicate = metadata$replicate

#these are the spike DNA sequences
spike_16S = 'CGGGCAGCTCTCGATAACCGGCGGAAGGTGGTAGCCACGGACAGGATCAGAACAATTAGAAGTGCCGCAGGTGGCCAAGTCCCCCGGACACAAGACGAGGCCGGAGGCCTGGTATATACACGTAGCTAAGAAGAGCTCATCCAGACTGGGAACGGTGTGCCAGCAGCCGCGGTAACATCACCACAACGTATTCGGTCACAAATTGATCGGAGGGAGAAATCGTCCGCAGGATCTCAAACTTTAACTAAGGACTAGTACTACATAGGCTCGAGAAGAGCTACCGTTTGCAGGGTCGCCGGGTACCGCTTAACCATAAAAGATCCACTCAGGTAGCCGTCCAGTTTCCTCTGAAATGATGGGGCGAGAAACACGGCTGGGCGTTATACGAGTGCTTTAGAATATGAGGAGAGACAGGGGTATATTCAAGG'
spike_18S = 'CTTCGTTATCGTCACGGAGAGAACTGCCTTTAGCGATCTGTCTAGAGACCCGCCAAATATAAGCCTTGGGGTTATTAACAAAGACGCCAAACACTAGTGAATATGACAGTGAAGGGGTGTGAGATAGCTTTACTGGGTGAGAAAACACTCGTTAAAAAGAATTAGACCGGATAATCCCCGAGGGGCCGTAGGCATGGACTTGTCGTTGCCACCGAGCATAGCGGTTTCGAAATAGCCGAGATGGGCACTGGCGAATTAACCCACTGGTTTATATGGATCCGATGGGTTCACTTAATAAGCTCGTACCAGGGATGAATAAAGCGTTACGAGAATTATAAACATGGAGTTCCTATTGATTTGAGGTTAATACCGAACGGGAACATTTGTCGATCATGCTTCACATAGAGT'

# For 16S - identify spike ASVs
spike_16S_ix = agrep(spike_16S, asv_16S)
spiketab_16S = seqtab_16S[spike_16S_ix, ]
dim(spiketab_16S)
#taxa_16S[spike_16S_ix,]

# For 18S - identify spike ASVs
spike_18S_ix = agrep(spike_18S, asv_18S)
spiketab_18S = seqtab_18S[spike_18S_ix, ]


## Exclude spike reads and Metazoa from taxa 
#For 16S
ix_taxa_16S =  setdiff(1:nrow(taxa_16S), spike_16S_ix)
taxa_16S = taxa_16S[ix_taxa_16S,]
seqtab_16S = seqtab_16S[ix_taxa_16S,]
#For 18S
ix_taxa_18S =  setdiff(1:nrow(taxa_18S), spike_18S_ix)
ix_taxa_18S = intersect(ix_taxa_18S, which(taxa_18S[,3] != 'Metazoa'))
ix_taxa_18S = union(ix_taxa_18S, which(is.na(taxa_18S[,3])))
taxa_18S = taxa_18S[ix_taxa_18S,]
seqtab_18S = seqtab_18S[ix_taxa_18S,]

# For 16S
norm_seqtab_16S = seqtab_16S
for (i in 1:ncol(seqtab_16S)) {
  norm_seqtab_16S[,i] = seqtab_16S[,i]/sum(seqtab_16S[,i])
}

# For 18S
norm_seqtab_18S = seqtab_18S
for (i in 1:ncol(seqtab_18S)) {
  norm_seqtab_18S[,i] = seqtab_18S[,i]/sum(seqtab_18S[,i])
}

# Rarefying reads per sample to the same counts
# For 16S
m=min(colSums(seqtab_16S))  # should we use alfat/beta diversity to select m?
r_seqtab_16S = rrarefy(x = t(seqtab_16S), sample = m)
colSums(r_seqtab_16S)
r_seqtab_16S = t(r_seqtab_16S)
# which(colSums(seqtab) == min(colSums(seqtab)))

r_norm_seqtab_16S = r_seqtab_16S
for (i in 1:ncol(r_seqtab_16S)) {
  r_norm_seqtab_16S[,i] = r_seqtab_16S[,i]/sum(r_seqtab_16S[,i])
}

# For 18S
m=min(colSums(seqtab_18S))
r_seqtab_18S = rrarefy(x = t(seqtab_18S), sample = m)
colSums(r_seqtab_18S)
r_seqtab_18S = t(r_seqtab_18S)
# which(colSums(seqtab) == min(colSums(seqtab)))

r_norm_seqtab_18S = r_seqtab_18S
for (i in 1:ncol(r_seqtab_18S)) {
  r_norm_seqtab_18S[,i] = r_seqtab_18S[,i]/sum(r_seqtab_18S[,i])
}


dir.create(argv$o)

export_table<-function(dfcount,dftax,file_name ) {
  #dfcount=norm_seqtab_16S;dftax=taxa_16S;file_name="norm_seqtab_16S.tsv" 
dataset_VAE <- data.frame(dfcount)
taxa_VAE<-data.frame(dftax)

taxa_VAE_tmp = taxa_VAE %>% filter(!is.na(Order))

dataset_VAE_temp <- dataset_VAE[row.names(dataset_VAE) %in% row.names(taxa_VAE_tmp),]

new_row_names=rep("", nrow(dataset_VAE_temp))
tokeep=rep(0, nrow(dataset_VAE_temp))
for ( i in 1:nrow(taxa_VAE_tmp)) {
  tempotext=row.names(taxa_VAE_tmp)[i]
  if (tempotext %in% row.names(dataset_VAE_temp)) {
    for (c in 1:ncol(taxa_VAE_tmp)) {
      if (!is.na(taxa_VAE_tmp[i,c])) tempotext=paste(tempotext,taxa_VAE_tmp[i,c],sep="_")
    }
    new_row_names[i]  <-tempotext
    tokeep[i]=i
  } 


}


new_row_names_clean=new_row_names[tokeep[tokeep>0]]

names(dataset_VAE_temp)<-gsub("X","",names(dataset_VAE_temp))
row.names(dataset_VAE_temp) <-new_row_names_clean


dataset_VAE_final <- data.frame(Samples=names(dataset_VAE_temp),    # Append new column to front of data
                                    data.frame(t(dataset_VAE_temp)))

fwrite(dataset_VAE_final, file=paste(argv$o, file_name, sep="/"), quote=FALSE, sep='\t', row.names = F)


fwrite(dataset_VAE_temp, file=paste(argv$o, paste0("T_",file_name), sep="/"), quote=FALSE, sep='\t', row.names = F)


return(dataset_VAE_final)
}


DF_norm_16S=export_table(norm_seqtab_16S,taxa_16S,"norm_seqtab_16S.tsv" )


DF_r_norm_16S=export_table(r_norm_seqtab_16S,taxa_16S,"r_norm_seqtab_16S.tsv" )



#Average of Temp,Salinity, Lat, Lon per station
#abiotics=metadata 

abiotics=metadata[,c(5,16,17,28)]

rest_metadata=metadata[,c(39,40,41,44,45,46,49:87)]
  

rest_metadata=rest_metadata %>% mutate_if(is.character, function(x) as.numeric(x))

ab_factors=c("Temperature",  "oxygen",  "ammonium",  "phosphate", "phosphoru",  "nitrogen",  "silicate",  "humic_substance", "chlorophyll")

for (patern in ab_factors) {
  item=paste0(patern,"_average")
  abiotics[[item]]=apply(rest_metadata[,grep(patern,names(rest_metadata),ignore.case = T)],1, function(x) mean(x,
    na.rm=T))
}

abiotics$salinity_average=as.numeric(metadata$salinity_average)

#abiotics=abiotics[,c(5,16,17,28,40,45,50,53,56,59,62,65,68,71,74,77,80,83,86)]


names(abiotics)[1] <- "Samples"

abiotics$Samples<-factor(abiotics$Samples, levels = unique(abiotics$Samples))
#str(abiotics)
abiotics=abiotics %>% mutate_if(is.character, function(x) as.numeric(x))

abiotics_means=aggregate(abiotics, list(abiotics$Samples), function(x) mean(x,na.rm = T))
abiotics_means$Samples=abiotics_means$Group.1
abiotics_means=abiotics_means[,!names(abiotics_means) %in% c("Group.1")]
fwrite(abiotics_means, file=paste(argv$o, "abiotic_avg_factors.tsv", sep="/"), quote=FALSE, sep='\t', row.names = F)

####



missing_data_stations=abiotics_means$Samples[! abiotics_means$Samples %in% DF_r_norm_16S$Samples]
### Anders code ends here

DF=merge(DF_r_norm_16S,abiotics_means, by="Samples")
DF=DF[,!names(DF) %in% c("Samples")]

F_analysis_with_CV<-function(MDF, target, myNtree, dirout) {
  prefi=paste0("RF_",target)
  #MDF=DF; target="Salinity"
  dir.create(dirout)
  set.seed(123)

  patern=paste("ASV",target,  sep="|")
  X0=MDF[,grep(patern,names(MDF),ignore.case = T)]
  X0=na.omit(X0)
  
  full_name=names(MDF)[grep(target,names(MDF),ignore.case = T)]
  names(X0)[grep(full_name,names(X0) )]<-"Response"
  
  df_split <- initial_split(X0, strata = Response)
  X <- training(df_split)
  df_test <- testing(df_split)
  


 # try(bestmtry <- tuneRF(X[, !names(X) %in% "Response"], X$Response,ntreeTry = myNtree, stepFactor = 1.2, improve = 0.01, trace=T, plot= T), silent = TRUE)
#  if (exists("bestmtry")) {
#    Mtry=bestmtry[match(min(bestmtry[,2]), bestmtry[,2]),1]
    
#    rf2 <- randomForest(X[,!names(X) %in% "Response"],X$Response, mtry=Mtry,ntree = myNtree, keep.forest=TRUE, importance=TRUE )
#  } else {rf2 <- randomForest(X[,!names(X) %in% "Response"],X$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE ) }
 
  rf2 <- randomForest(X[,!names(X) %in% "Response"],X$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE )   
  result<-rfcv(X[!names(X) %in% "Response"], X$Response)
  
  NV<-as.integer(names(which.min(result$error.cv)))
  total_features=length(names(X))-1
  if (NV >= nrow(X) & NV < total_features ) {  
    imp<-as.data.frame(importance(rf2))
    imp<-imp[order(imp[[1]], decreasing = TRUE), ]
    imp<-imp[1:NV,]
    
    selected_seqs=c("Response",row.names(imp))
    
    Xs<-X[, (names(X) %in% selected_seqs)]
    
 #   try(bestmtry2 <- tuneRF(Xs[, !names(Xs) %in% "Response"], Xs$Response,ntreeTry = myNtree,stepFactor = 1.2, improve = 0.01, trace=T, plot= T), silent = TRUE)
#    if (exists("bestmtry2")) {
#      Mtry2=bestmtry2[match(min(bestmtry2[,2]), bestmtry2[,2]),1]
      
#      rf <- randomForest(Xs[,!names(Xs) %in% "Response"],Xs$Response, mtry=Mtry2,ntree = myNtree, keep.forest=TRUE, importance=TRUE )
    rf <- randomForest(Xs[,!names(Xs) %in% "Response"],Xs$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE )
    
      
    } else {
      cat(paste("INFO: RF CV", prefi, "No sample has been removed", sep=" - "), "\n")
      rf<-rf2
      #dataMatrix_sub<-MDF
    }
    
#  }
  
  #rf <- randomForest(X[,!names(X) %in% "Response"],X$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE )
  
  seq_import<-as.data.frame(varImpPlot(rf))
  
  seq_import$Type="ASV"
  #seq_import$Type[grep("ASV.*",row.names(seq_import))]<-"ASV"
  seq_import$Type[grep("Archaea",row.names(seq_import))]<-"Archaea"
  seq_import$Type[grep("Bacteria",row.names(seq_import))]<-"Bacteria"
  
  newnames=c()
  vec=row.names(seq_import)
  for (i in 1:nrow(seq_import)) { 
    if (length(grep("ASV_",vec)) != 0) {
      m=strsplit(vec[i], split = "\\.")[[1]]
      n=paste(m[1],m[length(m)-1], sep=" ")
    } else {
      n=vec[i]
    }
    newnames<-c(newnames,n)
  }
  
  row.names(seq_import)<-newnames
  
  seq_import=seq_import[order(seq_import$`%IncMSE`, decreasing = TRUE),]
  #SEL=row.names(seq_import)
  
  #re=list(SELEC=SEL, data=Xs, RF=rf)
  #Exporting resutls
  sink(paste(dirout,paste(prefi,"txt",sep="."), sep = "/"))
  print(rf)
  print("Importance")
  print(seq_import)
  sink()
  ###
  
  #Creating figures
  seq_import=seq_import[(seq_import$`%IncMSE` > 4 | seq_import$`%IncMSE` < -4),]
  gC <- ggplot(data = seq_import, aes(x = reorder(row.names(seq_import),`%IncMSE`), y = `%IncMSE`)) +
    geom_bar(stat = "identity",aes(fill=Type)) + ggtitle("Random Forest importance - %IncMSE > |4|") +
    theme_bw()+       
    guides(fill=guide_legend(title=""))+ theme(legend.title = element_text(size=12)) +
    scale_color_manual(values = c("#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"))+
    theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(),axis.text.x = element_text(size=6)) +
    coord_flip()
  
  ggsave(paste(dirout,paste(prefi,"pdf",sep="."),sep = "/"), gC, width = 82, height = 42, units = "cm")
  
  predicted_tr<-predict(rf, X, type="response")
  ggDF=data.frame(Actual=X$Response,Predicted=predicted_tr)
  wTr<-ggplot(data=ggDF, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm )+stat_cor(label.y = max(X$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("train",prefi,"pdf",sep="."),sep = "/"), wTr, width = 42, height = 42, units = "cm")
  
  
  predicted_te<-predict(rf, df_test[,grep("ASV",names(df_test))], type="response")
  ggDFtest=data.frame(Actual=df_test$Response,Predicted=predicted_te)
  wTe<-ggplot(data=ggDFtest, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm )+stat_cor(label.y = max(df_test$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("test",prefi,"pdf",sep="."),sep = "/"), wTe, width = 42, height = 42, units = "cm")
  
  
  #return(re)
}

Dout=paste(argv$o,"ML_average_CV",sep="/")

ab_factors=c("Salinity", "Temperature",  "oxygen",  "ammonium",  "phosphate", "phosphoru",  "nitrogen",  "silicate",  "humic_substance", "chlorophyll")

for (fc in ab_factors) {
  cat("RF with ",fc, " \n")
  F_analysis_with_CV(DF, fc, 1000, Dout) 
}


#######
#xgb_analysis<-function(MDF, target, max_depth,nrounds_max,dirout) {  
xgb_analysis<-function(MDF, target, dirout) {
  #MDF=DF; target="salinity"; dirout=Dout2
  dir.create(dirout, recursive = T)
  patern=paste("ASV",target,  sep="|")
  
  
  X0=MDF[,grep(patern,names(MDF),ignore.case = T)]
  X0=na.omit(X0)
  
  full_name=names(MDF)[grep(target,names(MDF),ignore.case = T)]
  names(X0)[grep(full_name,names(X0) )]<-"Response"
  
  set.seed(123)
  
  df_split <- initial_split(X0, strata = Response)
  X <- training(df_split)
  df_test <- testing(df_split)
  
  
  prefi=paste0("xgb_",target)
  dtrain <- xgb.DMatrix(data = as.matrix(X[,!names(X) %in% "Response"]), label=X$Response)
  dtest <- xgb.DMatrix(data = as.matrix(df_test[,!names(df_test) %in% "Response"]), label=df_test$Response)
  
  df_split2 <- initial_split(X, strata = Response)
  X2 <- training(df_split2)
  df_valid <- testing(df_split2)
  dtrain2 <- xgb.DMatrix(data = as.matrix(X2[,!names(X2) %in% "Response"]), label=X2$Response)
  dvalid <- xgb.DMatrix(data = as.matrix(df_valid[,!names(df_valid) %in% "Response"]), label=df_valid$Response)
  
  watchlist <- list(train=dtrain2, val=dvalid)
  #watchlist <- list(train=dtrain)
  ###
  params <- list(objective = "reg:squarederror", gamma=5) #booster = "gbtree", 
  xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, verbose=F, early.stop.round = 20, maximize = F)
  #bst <- xgb.train(params = params, data = dtrain, nrounds = xgbcv$best_iteration, watchlist =watchlist , print.every.n = 10, early.stop.round = 10, maximize = F , eval_metric = "error")
  #list(val=dtest,train=dtrain)
  ###
  
  bst <- xgb.train(data=dtrain, nrounds=xgbcv$best_iteration, watchlist=watchlist, eval.metric = "error", #eval.metric = "logloss", 
                   objective = "reg:squarederror", verbose=0)
  
  
  predtr <- predict(bst, dtrain)
  ggDF=data.frame(Actual=X$Response,Predicted=predtr)
  wTr<-ggplot(data=ggDF, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm )+stat_cor(label.y = max(X$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("train",prefi,"pdf",sep="."),sep = "/"), wTr, width = 42, height = 42, units = "cm")
  
  
  
  pred <- predict(bst, dtest)
  
  ggDFtest=data.frame(Actual=df_test$Response,Predicted= pred)
  wTe<-ggplot(data=ggDFtest, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm )+stat_cor(label.y = max(df_test$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("test",prefi,"pdf",sep="."),sep = "/"), wTe, width = 42, height = 42, units = "cm")
  
  
  importance_matrix <- xgb.importance(model = bst)
  
  importance_matrix2 <- importance_matrix[order(importance_matrix$Gain),]
  
  sink(paste(dirout,paste(prefi,"txt",sep="."), sep = "/"))
  print(importance_matrix2)
  sink()
  

  xgb_ggplot <- xgb.ggplot.importance(importance_matrix = importance_matrix[1:60]) +theme_bw() #,rel_to_first = TRUE, xlab = "Relative importance")
  
  ggsave(paste(dirout,paste(prefi,"pdf",sep="."),sep = "/"), xgb_ggplot, width = 82, height = 42, units = "cm")
  gr<-xgb.plot.tree(model = bst,render=FALSE) #, fname=paste(dirout,paste(prefi,"tree2","pdf",sep="."),sep = "/"))
  
  export_graph(gr, paste(dirout,paste(prefi,"tree","pdf",sep="."),sep = "/"), width=1500, height=1900)
}

Dout2=paste(argv$o,"ML_avg","xgb",sep="/")

ab_factors=c("Salinity", "Temperature",  "oxygen",  "ammonium",  "phosphate", "phosphoru",  "nitrogen",  "silicate",  "humic_substance", "chlorophyll")
#fc="oxygen"
for (fc in ab_factors) {
  cat("xgb with ",fc, " \n")
#xgb_analysis(DF, fc, 7, 7,Dout2)
  xgb_analysis(DF, fc, Dout2)
  }

# save model to binary local file
#xgb.save(bst, "xgboost.model")

# load binary model to R
#bst2 <- xgb.load("xgboost.model")
#pred2 <- predict(bst2, test$data)



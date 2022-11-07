# Functions needed for sc and Bcr data matching



#==============================================================================#
# sub-function-1: cloneMax_consensusCountMax
# input: 
#   - df_bcr: a bcr data
#   - df_barcode: cell_id or cluster_cb used for collapse clone
#   - default use vj_junction_aa for grouping of clonal type
#==============================================================================#

cloneMax_consensusCountMax<- function(df_bcr){
  # heavy version for nonUniq_lvD2_bcr
  uniq_clusterCb <- unique(df_bcr$cluster_cb)
  #Con.df <- data.frame(matrix(NA, length(uniq_clusterCb),length(colnames(df_bcr))+1))
  #colnames(Con.df) <- c("uniq_cluster",colnames(df_bcr))
  #Con.df$uniq_cluster <- uniq_clusterCb
  
  finalBcr <- list()
  for(y in seq_len(length(uniq_clusterCb))){
    cluster.i <- uniq_clusterCb[y]
    location.i <- which(cluster.i==df_bcr$cluster_cb)
    data.y<- df_bcr[location.i,]
    cloneCount <- data.y %>% group_by(vj_junction_aa )%>% 
      summarise(cb_no=n()) %>% arrange(desc(cb_no))
    
    cloneCountMax <- cloneCount[which.max(cloneCount$cb_no),]
    data.y.Max <- data.y %>% filter(vj_junction_aa %in% cloneCountMax$vj_junction_aa)
    data.y.Max_consensusCountMax <- data.y.Max[which.max(data.y.Max$consensus_count),]
    
    finalBcr[[y]]<- data.y.Max_consensusCountMax
    
  }
  final<- do.call(rbind, finalBcr)
  return(final)
}


###################################
# sub-function-2
###################################
cloneMax_consensusCountMax_cel_id<- function(df_bcr){
  # heavy version for nonUniq_lvD2_bcr
  uniq_cel_id <- unique(df_bcr$cel_id)
  #Con.df <- data.frame(matrix(NA, length(uniq_cel_id),length(colnames(df_bcr))+1))
  #colnames(Con.df) <- c("uniq_cluster",colnames(df_bcr))
  #Con.df$uniq_cluster <- uniq_cel_id
  
  finalBcr <- list()
  for(y in seq_len(length(uniq_cel_id))){
    cluster.i <- uniq_cel_id[y]
    location.i <- which(cluster.i==df_bcr$cel_id)
    data.y<- df_bcr[location.i,]
    cloneCount <- data.y %>% group_by(vj_junction_aa )%>% 
      summarise(cb_no=n()) %>% arrange(desc(cb_no))
    
    cloneCountMax <- cloneCount[which.max(cloneCount$cb_no),]
    data.y.Max <- data.y %>% filter(vj_junction_aa %in% cloneCountMax$vj_junction_aa)
    data.y.Max_consensusCountMax <- data.y.Max[which.max(data.y.Max$consensus_count),]
    
    finalBcr[[y]]<- data.y.Max_consensusCountMax
    
  }
  final<- do.call(rbind, finalBcr)
  return(final)
}


###################################
# sub-function-3
###################################
# https://stackoverflow.com/questions/5721883/agrep-only-return-best-matches

# define the "closestMatch" function to calculate the distance and select the 
# - lv distance <=2(similarity>=0.83)
# - 
closestMatch= function(string,stringVector,sim_threshold){
  similarity = levenshteinSim(string, stringVector)
  lv_match <- data.frame(stringVector,similarity)
  
  #  using "filter(similarity==max(similarity) )" to select multiple max value
  closest<-lv_match %>% 
    filter(similarity>sim_threshold) %>%
    mutate(cluster_cb=string)%>% mutate(cell_id=stringVector)%>% 
    filter(similarity==max(similarity) ) %>% 
    select(cell_id,similarity,cluster_cb)
}
 



##==============================================================================##
#    MAIN FUNCTION-1: LV correction for heavy chain 
##==============================================================================##
  

lv_repairBcrHeavy <- function(file_List,df_igblastBcrPath, df_TRUST4BcrPath, df_scRdsPath,outputPathPlot, 
                               outputPathSummary){
 
    # allow one bp mismatch and calculate the threshold
    threshold<- round((12-1)/12, digits = 1)
    final <- list()
    for(i in seq_along(file_List)){
      # convert the file name to fit the rds files(from xiaoli) which have additional "2B" in the file name 
      fileName<- file_List[i]
      
      if(fileName=="HD1500"){
        fileNameRds<- "HD2B1500"
      }else if(fileName=="HD3600"){
        fileNameRds<-"HD2B3600"
      }else if(fileName=="HD960"){
        fileNameRds<-"HD2B960"
      }else {
        fileNameRds<- fileName
      }
      
      
      ###############
      # READ BCR DATA (igblast)
      ###############
      
      bcr_db<- read_tsv(file=paste0(df_igblastBcrPath,fileName,"_noPearIgblastBarcodeAirr_db-pass.tsv"))
      bcr_db$cell_id= sapply(strsplit(bcr_db$sequence_id,"_"), "[",1)  
      dim(bcr_db)
      bcr_db<-within(bcr_db, 
                     vj_junction_aa<-paste0(str_split(v_call,",",simplify = T)[,1],"_",str_split(j_call,",",simplify = T)[,1],"_",junction_aa) )
      
      colnames(bcr_db)
      bcr_db$vj_junction_aa[1:4]
      # View(bcr_db)
      
      
      ###############
      # READ BCR DATA (TRUST4 *barcode_airr.tsv)
      ###############
      bcr_db_TRUST<- read_tsv(file=paste0(df_TRUST4BcrPath,"/",fileName,"HL_barcode_airr.tsv"))
      bcr_db_TRUST_select <- bcr_db_TRUST %>% select(sequence_id,consensus_count)
      
      dim(bcr_db_TRUST_select)# slightly more bcr than igblast 
      
      
      ###############
      # add consensus count (TRUST4 *barcode_airr.tsv) to bcr(igblast)
      ###############
      
      bcr <- left_join(x=bcr_db,y=bcr_db_TRUST_select,by=c("sequence_id"))
      dim(bcr)
      dim(bcr_db_TRUST_select)
      dim(bcr_db)
      #View(bcr)
      
      heavy_bcr <- bcr[grep("IGH",bcr$v_call ,invert= F),]
      
      ###############
      # Read GEX data
      ###############
      if(fileName=="HD3600"|fileName=="HD1500"|fileName=="HD960"){
        gex_db <- readRDS(paste0(df_scRdsPath,fileNameRds,".rds"))
        sc_cell_id<- sapply(strsplit(Cells(gex_db),"_"), "[",2)  
        
        gex_cell_id <- data.frame(sc_cell_id,sc_cell_id)
        colnames(gex_cell_id) <- c("sc_cell_id","repeat")
        dim(gex_cell_id)
        head(gex_cell_id)
        
      }else{
        gex_db <- readRDS(paste0(df_scRdsPath,fileNameRds,".rds"))
        sc_cell_id<- sapply(strsplit(Cells(gex_db),"_"), "[",3)  
        
        gex_cell_id <- data.frame(sc_cell_id,sc_cell_id)
        colnames(gex_cell_id) <- c("sc_cell_id","repeat")
        dim(gex_cell_id)
        head(gex_cell_id)
      }
      
      
      ################################################################################
      #     1. -  separate sc cell_id in and not in BCR heavy 
      ################################################################################
      ###############
      # sc  not in bcr
      ###############
      not_gex_cell_id <- gex_cell_id %>%  filter(! sc_cell_id %in% heavy_bcr$cell_id) 
      dim(not_gex_cell_id)
      dim(gex_cell_id)
      
      ################################################################################
      #     2.-  separate bcr cell_id in and not in sc heavy 
      ################################################################################
      ###############
      # bcr in sc
      ###############
      matched_bcr_db <- heavy_bcr%>%  filter( cell_id %in% gex_cell_id$sc_cell_id) 
      
      ###############
      # bcr not in sc
      ###############
      not_bcr_db <- heavy_bcr %>%  filter(! cell_id %in% gex_cell_id$sc_cell_id) 
      dim(not_bcr_db)
      dim(matched_bcr_db)
      dim(heavy_bcr)
      
      
      ################################################################################
      #     3.-  match  sc cluster_cb to not_in bcr cell_id  allowing lv distance=2
      ################################################################################
       
      ##########################
      #   lv calculation
      ##########################
      #-------------------------------------------------------------------------------
      # https://stackoverflow.com/questions/5721883/agrep-only-return-best-matches
      
      lv<-list()
      for(y in seq_len(length(not_gex_cell_id$sc_cell_id))){
        not_sc_cell_id <- not_gex_cell_id$sc_cell_id[y]
        
        # define the "closestMatch" function to calculate the distance and select the 
        # - lv distance <=2(similarity>=0.83)
        # - 
        closestMatch= function(string,stringVector,sim_threshold){
          similarity = levenshteinSim(string, stringVector)
          lv_match <- data.frame(stringVector,similarity)
          
          #  using "filter(similarity==max(similarity) )" to selet multipe max value
          closest<-lv_match %>% 
            filter(similarity>sim_threshold) %>%
            mutate(cluster_cb=string)%>% mutate(cell_id=stringVector)%>% 
            filter(similarity==max(similarity) ) %>% 
            select(cell_id,similarity,cluster_cb)
        }
        
        # similarity=0.83333 is equal to lv=2 distance(2 mutation/insertion/deletion)
        lv_closest <- closestMatch(not_sc_cell_id,not_bcr_db$cell_id,sim_threshold=threshold)
        lv[[y]]<- lv_closest
      }
      
      
      lv2Index <- do.call(rbind,lv)
      
      dim(lv2Index)
      #View(finalBcr)
      head(lv2Index)
      length(unique(lv2Index$cluster_cb))  # 503 
      ################################################################################
      #     4.-  inner_join the lv2 bcr_sc index with not_in bcr 
      ################################################################################
      not_bcr_db_lv2 <- inner_join(x=not_bcr_db, y=lv2Index,by=c("cell_id"))
      
      dim(not_bcr_db_lv2)
      dim(not_bcr_db)
      length(unique(not_bcr_db_lv2$cell_id))   # 1203
      length(unique(not_bcr_db_lv2$cluster_cb))  # 503
      
      ##################
      # func: view the conflict clone type with same lv cell_id
      # heavy version (borrow from pipeline2 code)
      ##################
       pdf(file = paste0(outputPathPlot,fileName,"_ld1_Heavy_conflictCloneDistribution.pdf"))
      
      p_uniquClone_lvD2 <- not_bcr_db_lv2 %>% group_by(cluster_cb, vj_junction_aa )%>%
        summarise(cb_no=n()) %>% arrange(desc(cb_no))
      
      p<- p_uniquClone_lvD2 %>% group_by(cluster_cb) %>%  summarise(diffClone_no=n()) %>%
        arrange(desc(diffClone_no)) %>%
        ggplot(aes(x=diffClone_no))+
        geom_histogram(binwidth=1,aes(y=..density..))+
        geom_density(alpha=.2, fill="#FF6666") +
        labs(title=fileName ,x="heavy chain different cloneType  per cluster cb", y = "Density")
      print(p)
      dev.off()
      
      
      ##################
      # func: extract unique clone
      
      ##################
      lvD2_bcr <-not_bcr_db_lv2 %>% group_by(cluster_cb, vj_junction_aa )%>%
        summarise(cb_no=n()) %>% arrange(desc(cb_no))
      # extract unique cluster_cb and vj_junction_aa
      uniq_cluster<- lvD2_bcr %>% group_by(cluster_cb) %>% summarise(conflict_no=n()) %>%
        arrange(desc(conflict_no)) %>% filter(conflict_no==1)
      # retrive the full bcr infor
      uniq_lvD2_bcr <- not_bcr_db_lv2 %>% filter(cluster_cb %in% uniq_cluster$cluster_cb)
      
      # collapse the unique clone which have same cluster_cb and same vj_junctiona_aa but different seq_id
      uniq_lvD2_bcr1<-cloneMax_consensusCountMax(df_bcr=uniq_lvD2_bcr)
      dim(uniq_lvD2_bcr1)
      length(unique(uniq_lvD2_bcr1$cell_id))
      ##################
      # Extract   non-unique clone
      ##################
      # group cluster_cb with vj_junction_aa
      lvD2_bcr <- not_bcr_db_lv2 %>% group_by(cluster_cb, vj_junction_aa )%>%
        summarise(cb_no=n()) %>% arrange(desc(cb_no))
      # extract non unique cluster_cb and vj_junction_aa
      nonUniq_cluster<- lvD2_bcr %>% group_by(cluster_cb) %>% summarise(conflict_no=n()) %>%
        arrange(desc(conflict_no)) %>% filter(conflict_no>1)
      # retrive the full bcr infor
      nonUniq_lvD2_bcr <- not_bcr_db_lv2 %>% filter(cluster_cb %in% nonUniq_cluster$cluster_cb)
 
      #########
      # collapse by choosing the max cb no. , if multiple max, choose the 1st max
      #########
      finalBcr1<-cloneMax_consensusCountMax(nonUniq_lvD2_bcr)
      dim(finalBcr1)
      length(unique(finalBcr1$cell_id))
      
      #########
      # rbind the matched bcr with unique lv2 bcr  and  non unique maxBcrCountMaxConsensusCount_collapsed bcr 
      ######### 
      matched_bcr_db$cluster_cb<-matched_bcr_db$cell_id
      matched_bcr_db$similarity<-"1"
      
      dim(matched_bcr_db)
      length(unique(matched_bcr_db$cell_id))
      
      colnames(uniq_lvD2_bcr)
      colnames(finalBcr1)
      
      # combine the lv2 bcr with unique or non unique bcr for the same cluster_cb
      lv2_bcr <- rbind(uniq_lvD2_bcr,finalBcr1)
      dim(lv2_bcr)
      length(unique(lv2_bcr$cell_id))
      length(unique(lv2_bcr$cluster_cb))
      #-------------------------------------------------------------------------------
      # Note  different cluster_cb might matches to the same cell_id 
      # to collapse the redundent cluster_cb, first group by cell_id >> select the 
      # one with highest similarity >> choose the one with max consensus count 
      
      data2 <-list()
      for(y in seq_len(length(lv2_bcr$cell_id))){
        
        uniq_cell_id_list <- unique(lv2_bcr$cell_id)
        cell_id.i <- uniq_cell_id_list[y]
        location.i <- which(cell_id.i==lv2_bcr$cell_id)
        data.y<- lv2_bcr[location.i,]
        
        uniqMaxSim <-data.y %>% filter(similarity==max(similarity) ) %>% filter(consensus_count==max(consensus_count)) 
        data2[[y]]<- uniqMaxSim
      }
      finalBcr2 <- do.call(rbind, data2)
      
      dim(finalBcr2)
      length(unique(finalBcr2$cell_id))
      #---------------------------------------------------------------------------
      test <- finalBcr2 %>%group_by(cell_id) %>% summarise(n=n()) %>% 
        arrange(desc(n)) %>% head(5)
      test1<-finalBcr2%>% filter(cell_id %in%test$cell_id) %>% 
        select(sequence_id,similarity,consensus_count,cell_id,cluster_cb)
      test1
      
      #-------------------------------------------------------------------------------
      finalCombinedBcr <- rbind(matched_bcr_db,finalBcr2)
      
      dim(finalCombinedBcr)
      length(finalCombinedBcr$cell_id)
      length(finalCombinedBcr$cluster_cb)
      # extract the sequence_id and cluster_cb and make a new variable as sequence_id_unique
      finalCombinedBcr_pip11 <- finalCombinedBcr %>% 
        select(sequence_id,cell_id,cluster_cb)
      head(finalCombinedBcr_pip11)
      
      final[[i]]<- finalCombinedBcr
      # Write to pip11 results
      write_tsv(x=finalCombinedBcr, file =paste0(outputPathSummary,"/",fileNameRds,"_lv1CorrectHeavyBcr.tsv") )

    }
  
  
  final_lv2<-setNames(final, nm = fileList)
  return(final_lv2)
  
  }
 

  

##==============================================================================##
#    MAIN FUNCTION-2: LV correction for light chain 
##==============================================================================##

 
lv_repairBcrLight <- function(file_List,df_igblastBcrPath, df_TRUST4BcrPath, df_scRdsPath,outputPathPlot,
                              outputPathSummary ){
  
    threshold<- round((12-1)/12, digits = 1)
    final <- list()
    for(i in seq_along(file_List)){
      # convert the file name to fit the rds files(from xiaoli) which have additional "2B" in the file name 
      fileName<- file_List[i]
      
      if(fileName=="HD1500"){
        fileNameRds<- "HD2B1500"
      }else if(fileName=="HD3600"){
        fileNameRds<-"HD2B3600"
      }else if(fileName=="HD960"){
        fileNameRds<-"HD2B960"
      }else {
        fileNameRds<- fileName
      }
      
      
      ###############
      # READ BCR DATA (igblast)
      ###############
      
      bcr_db<- read_tsv(file=paste0(df_igblastBcrPath,fileName,"_noPearIgblastBarcodeAirr_db-pass.tsv"))
      bcr_db$cell_id= sapply(strsplit(bcr_db$sequence_id,"_"), "[",1)  
      dim(bcr_db)
      bcr_db<-within(bcr_db, 
                     vj_junction_aa<-paste0(str_split(v_call,",",simplify = T)[,1],"_",str_split(j_call,",",simplify = T)[,1],"_",junction_aa) )
      
      colnames(bcr_db)
      bcr_db$vj_junction_aa[1:4]
      # View(bcr_db)
      
      
      ###############
      # READ BCR DATA (TRUST4 *barcode_airr.tsv)
      ###############
      bcr_db_TRUST<- read_tsv(file=paste0(df_TRUST4BcrPath,"/",fileName,"HL_barcode_airr.tsv"))
      bcr_db_TRUST_select <- bcr_db_TRUST %>% select(sequence_id,consensus_count)
      
      dim(bcr_db_TRUST_select)# slightly more bcr than igblast 
      
      
      ###############
      # add consensus count (TRUST4 *barcode_airr.tsv) to bcr(igblast)
      # 
      ###############
      
      bcr <- left_join(x=bcr_db,y=bcr_db_TRUST_select,by=c("sequence_id"))
      dim(bcr)
      dim(bcr_db_TRUST_select)
      dim(bcr_db)
      #View(bcr)
      
      #heavy_bcr <- bcr[grep("IGH",bcr$v_call ,invert= F),]
      light_bcr <- bcr[grep("IGH",bcr$v_call ,invert= T),]
      
      ###############
      # Read GEX data
      ###############
      if(fileName=="HD3600"|fileName=="HD1500"|fileName=="HD960"){
        gex_db <- readRDS(paste0(df_scRdsPath,fileNameRds,".rds"))
        sc_cell_id<- sapply(strsplit(Cells(gex_db),"_"), "[",2)  
        
        gex_cell_id <- data.frame(sc_cell_id,sc_cell_id)
        colnames(gex_cell_id) <- c("sc_cell_id","repeat")
        dim(gex_cell_id)
        head(gex_cell_id)
        
      }else{
        gex_db <- readRDS(paste0(df_scRdsPath,fileNameRds,".rds"))
        sc_cell_id<- sapply(strsplit(Cells(gex_db),"_"), "[",3)  
        
        gex_cell_id <- data.frame(sc_cell_id,sc_cell_id)
        colnames(gex_cell_id) <- c("sc_cell_id","repeat")
        dim(gex_cell_id)
        head(gex_cell_id)
      }
      
      
      ################################################################################
      #     1. -  separate sc cell_id in and not in BCR heavy 
      ################################################################################
      ###############
      # sc  not in bcr
      ###############
      not_gex_cell_id <- gex_cell_id %>%  filter(! sc_cell_id %in% light_bcr$cell_id) 
      dim(not_gex_cell_id)
      dim(gex_cell_id)
      
      ################################################################################
      #     2.-  separate bcr cell_id in and not in sc light 
      ################################################################################
      ###############
      # bcr in sc
      ###############
      matched_bcr_db <- light_bcr%>%  filter( cell_id %in% gex_cell_id$sc_cell_id) 
      
      ###############
      # bcr not in sc
      ###############
      not_bcr_db <- light_bcr %>%  filter(! cell_id %in% gex_cell_id$sc_cell_id) 
      dim(not_bcr_db)
      dim(matched_bcr_db)
      dim(light_bcr)
      
      
      ################################################################################
      #     3.-  match  sc cluster_cb to not_in bcr cell_id  allowing lv distance=2
      ################################################################################
      
      ##########################
      # lv calculation
      ##########################
      #-------------------------------------------------------------------------------
      # https://stackoverflow.com/questions/5721883/agrep-only-return-best-matches
      
      lv<-list()
      for(y in seq_len(length(not_gex_cell_id$sc_cell_id))){
        not_sc_cell_id <- not_gex_cell_id$sc_cell_id[y]
        
        # define the "closestMatch" function to calculate the distance and select the 
        # - lv distance <=2(similarity>=0.83)
        # - 
        closestMatch= function(string,stringVector,sim_threshold){
          similarity = levenshteinSim(string, stringVector)
          lv_match <- data.frame(stringVector,similarity)
          
          #  using "filter(similarity==max(similarity) )" to selet multipe max value
          closest<-lv_match %>% 
            filter(similarity>sim_threshold) %>%
            mutate(cluster_cb=string)%>% mutate(cell_id=stringVector)%>% 
            filter(similarity==max(similarity) ) %>% 
            select(cell_id,similarity,cluster_cb)
        }
        
        # similarity=0.83333 is equal to lv=2 distance(2 mutation/insertion/deletion)
        lv_closest <- closestMatch(not_sc_cell_id,not_bcr_db$cell_id,sim_threshold=threshold)
        lv[[y]]<- lv_closest
      }
      
      
      lv2Index <- do.call(rbind,lv)
      
      dim(lv2Index)
      #View(finalBcr)
      head(lv2Index)
      length(unique(lv2Index$cluster_cb))  # 503 
      ################################################################################
      #     4.-  inner_join the lv2 bcr_sc index with not_in bcr 
      ################################################################################
      not_bcr_db_lv2 <- inner_join(x=not_bcr_db, y=lv2Index,by=c("cell_id"))
      
      dim(not_bcr_db_lv2)
      dim(not_bcr_db)
      length(unique(not_bcr_db_lv2$cell_id))   # 1203
      length(unique(not_bcr_db_lv2$cluster_cb))  # 503
      
      ##################
      # func: view the conflict clone type with same lv cell_id
      # heavy version (borrow from pipeline2 code)
      ##################
      pdf(file = paste0(outputPathPlot,fileName,"_ld1_","Light_conflictCloneDistribution.pdf"))
      
      p_uniquClone_lvD2 <- not_bcr_db_lv2 %>% group_by(cluster_cb, vj_junction_aa )%>%
        summarise(cb_no=n()) %>% arrange(desc(cb_no))
      
      p<- p_uniquClone_lvD2 %>% group_by(cluster_cb) %>%  summarise(diffClone_no=n()) %>%
        arrange(desc(diffClone_no)) %>%
        ggplot(aes(x=diffClone_no))+
        geom_histogram(binwidth=1,aes(y=..density..))+
        geom_density(alpha=.2, fill="#FF6666") +
        labs(title=fileName ,x="Light chain different cloneType  per cluster cb", y = "Density")
      print(p)
      dev.off()
      
      
      ##################
      # func: extract unique clone
      
      ##################
      lvD2_bcr <-not_bcr_db_lv2 %>% group_by(cluster_cb, vj_junction_aa )%>%
        summarise(cb_no=n()) %>% arrange(desc(cb_no))
      # extract unique cluster_cb and vj_junction_aa
      uniq_cluster<- lvD2_bcr %>% group_by(cluster_cb) %>% summarise(conflict_no=n()) %>%
        arrange(desc(conflict_no)) %>% filter(conflict_no==1)
      # retrive the full bcr infor
      uniq_lvD2_bcr <- not_bcr_db_lv2 %>% filter(cluster_cb %in% uniq_cluster$cluster_cb)
      
      # collapse the unique clone which have same cluster_cb and same vj_junctiona_aa but different seq_id
      uniq_lvD2_bcr1<-cloneMax_consensusCountMax(df_bcr=uniq_lvD2_bcr)
      dim(uniq_lvD2_bcr1)
      length(unique(uniq_lvD2_bcr1$cell_id))
      ##################
      # Extract and sankey plot non-unique clone
      ##################
      # group cluster_cb with vj_junction_aa
      lvD2_bcr <- not_bcr_db_lv2 %>% group_by(cluster_cb, vj_junction_aa )%>%
        summarise(cb_no=n()) %>% arrange(desc(cb_no))
      # extract non unique cluster_cb and vj_junction_aa
      nonUniq_cluster<- lvD2_bcr %>% group_by(cluster_cb) %>% summarise(conflict_no=n()) %>%
        arrange(desc(conflict_no)) %>% filter(conflict_no>1)
      # retrive the full bcr infor
      nonUniq_lvD2_bcr <- not_bcr_db_lv2 %>% filter(cluster_cb %in% nonUniq_cluster$cluster_cb)
      
      
     
      
      #########
      # collapse by choosing the max cb no. , if multiple max, choose the 1st max
      #########
      finalBcr1<-cloneMax_consensusCountMax(nonUniq_lvD2_bcr)
      dim(finalBcr1)
      length(unique(finalBcr1$cell_id))
      #########
      # rbind the matched bcr with unique lv2 bcr  and  non unique maxBcrCountMaxConsensusCount_collapsed bcr 
      ######### 
      matched_bcr_db$cluster_cb<-matched_bcr_db$cell_id
      matched_bcr_db$similarity<-"1"
      
      dim(matched_bcr_db)
      length(unique(matched_bcr_db$cell_id))
      
      colnames(uniq_lvD2_bcr)
      colnames(finalBcr1)
      
      # combine the lv2 bcr with unique or non unique bcr for the same cluster_cb
      lv2_bcr <- rbind(uniq_lvD2_bcr,finalBcr1)
      dim(lv2_bcr)
      length(unique(lv2_bcr$cell_id))
      length(unique(lv2_bcr$cluster_cb))
      #-------------------------------------------------------------------------------
      # Note  different cluster_cb might matches to the same cell_id 
      # to collapse the redundent cluster_cb, we need to group by cell_id and select the 
      # most similarity then max consensus count 
      
      data2 <-list()
      for(y in seq_len(length(lv2_bcr$cell_id))){
        
        uniq_cell_id_list <- unique(lv2_bcr$cell_id)
        cell_id.i <- uniq_cell_id_list[y]
        location.i <- which(cell_id.i==lv2_bcr$cell_id)
        data.y<- lv2_bcr[location.i,]
        
        uniqMaxSim <-data.y %>% filter(similarity==max(similarity) ) %>% filter(consensus_count==max(consensus_count)) 
        data2[[y]]<- uniqMaxSim
      }
      finalBcr2 <- do.call(rbind, data2)
      
      dim(finalBcr2)
      length(unique(finalBcr2$cell_id))
      #---------------------------------------------------------------------------
      test <- finalBcr2 %>%group_by(cell_id) %>% summarise(n=n()) %>% 
        arrange(desc(n)) %>% head(5)
      test1<-finalBcr2%>% filter(cell_id %in%test$cell_id) %>% 
        select(sequence_id,similarity,consensus_count,cell_id,cluster_cb)
      test1
      
      #-------------------------------------------------------------------------------
      finalCombinedBcr <- rbind(matched_bcr_db,finalBcr2)
      
      dim(finalCombinedBcr)
      length(finalCombinedBcr$cell_id)
      length(finalCombinedBcr$cluster_cb)
      # extract the sequence_id and cluster_cb and make a new variable as sequence_id_unique
      finalCombinedBcr_pip11 <- finalCombinedBcr %>% 
        select(sequence_id,cell_id,cluster_cb)
      head(finalCombinedBcr_pip11)
      
      final[[i]]<- finalCombinedBcr
      # Write to pip11 results
      write_tsv(x=finalCombinedBcr, file =paste0(outputPathSummary,"/",fileNameRds,"_lv1CorrectLightBcr.tsv") )
       
    }
    
 
  final_lv2<-setNames(final, nm = fileList)
  return(final_lv2)
  
}

 
























##############################################################################
# match bcr data with sc data (allowing one mismatch)
##############################################################################

summary_scMatchToBcrHL<- function(df_fileList,df_bcrInputPath,rdsInputPath,outputPath,outputSummary){
  finalMatchB<- list()
  finalMatchB_lv1<- list()
  for(i in seq_along(df_fileList)){
    name.i<- df_fileList[i]
    
    
    # Read BCR H data
    # convert junction_aa to CDR3aa by removing one flanking from both end 
    
    bcr_db_H <- read_tsv(paste0(df_bcrInputPath, name.i,"_lv1CorrectHeavyBcr.tsv"))
    
    bcr_db_H$v_call <- sapply(str_split(bcr_db_H$v_call, pattern = "\\*"),"[",1)
    bcr_db_H$j_call <- sapply(str_split(bcr_db_H$j_call, pattern = "\\*"),"[",1)
    bcr_db_H<-within(bcr_db_H, 
                     cdr3aa<-substr(junction_aa,2,nchar(junction_aa)-1) )
    
    # Read BCR L data
    
    bcr_db_L <- read_tsv(paste0(df_bcrInputPath, name.i,"_lv1CorrectLightBcr.tsv"))
    bcr_db_L$v_call <- sapply(str_split(bcr_db_L$v_call, pattern = "\\*"),"[",1)
    bcr_db_L$j_call <- sapply(str_split(bcr_db_L$j_call, pattern = "\\*"),"[",1)
    bcr_db_L<-within(bcr_db_L, 
                     cdr3aa<-substr(junction_aa,2,nchar(junction_aa)-1) )
    
    # Read GEX data
    gex_db <- readRDS(paste0(rdsInputPath, name.i,".rds"))
    
    
    ## format cell IDs
    if(name.i=="p36T1"){
      sample_id<-"p36_T1"
      sampleName<-"p36D4"
      dpi<- 4
      
    }else if(name.i=="p36T2"){
      sample_id<-"p36_T2"
      sampleName<-"p36D45"
      dpi<- 45
    }else if(name.i=="p36T3"){
      sample_id<-"p36_T3"
      sampleName<-"p36D92"
      dpi<-92
    }else if(name.i=="p36T4"){
      sample_id<-"p36_T4"
      sampleName<-"p36D169"
      dpi<-169
    }else if(name.i=="p45T2"){
      sample_id<-"p45_T2"
      sampleName<-"p45D24"
      dpi<-24
    }else if(name.i=="HD2B1500"){
      sample_id<-"HD2B1500"
      sampleName<-"HD1"
      dpi<-0
    }else if(name.i=="HD2B3600"){
      sample_id<-"HD2B3600"
      sampleName<-"HD2"
      dpi<-0
    }else if(name.i=="HD2B960"){
      sample_id<-"HD2B960"
      sampleName<-"HD1_stD2"
      dpi<-0
    }else {
      print("now only support p36 time points")
    }
    
    bcr_db_H$sample<- sample_id
    bcr_db_H$day<- dpi
    
    bcr_db_L$sample<- sample_id
    bcr_db_L$day<- dpi
    # make dpi as day_f (factor)
    bcr_db_H$day_f<- factor(dpi, levels=c(0,4,24,45,92,169))
    bcr_db_L$day_f<- factor(dpi, levels=c(0,4,24,45,92,169))
    
    # NOTE: cell_id becomes the cluster_cb after lv1 correction
    bcr_db_H$cell_id_unique = paste0(bcr_db_H$sample, "_", bcr_db_H$cell_id)
    bcr_db_H$cluster_cb_unique = paste0(bcr_db_H$sample, "_", bcr_db_H$cluster_cb)
    bcr_db_H$cell_id_unique[1:4]
    
    bcr_db_L$cell_id_unique = paste0(bcr_db_L$sample, "_", bcr_db_L$cell_id)
    bcr_db_L$cluster_cb_unique = paste0(bcr_db_L$sample, "_", bcr_db_L$cluster_cb)
    bcr_db_L$cell_id_unique[1:4]
    
    #########################################################################
    # overlap between sc B cells and BCR  
    #########################################################################
    
    #####################
    # matching after lv1 correction
    #####################
    # match index to find the position of the GEX cells in the BCR data
    match.index_lv1_H = match(Cells(gex_db), bcr_db_H$cluster_cb_unique)
    match.index_lv1_L= match(Cells(gex_db), bcr_db_L$cluster_cb_unique)
    
    #####################
    # Add BCR h and L v,d.j and CDR3 aa to sc data as meta data
    #####################
    v_H <- bcr_db_H$v_call
    d_H <- bcr_db_H$d_call
    j_H <- bcr_db_H$j_call
    cdr3aa_H<- bcr_db_H$cdr3aa
    
    v_L <- bcr_db_L$v_call 
    j_L <- bcr_db_L$j_call
    cdr3aa_L<- bcr_db_L$cdr3aa
    
    
    #############################################################################
    # lv1 based matching 
    #############################################################################
    
    gex_db$v_H_lv1= unlist(lapply(match.index_lv1_H, function(x){ifelse(is.na(x),NA, v_H[x])}))
    gex_db$d_H_lv1= unlist(lapply(match.index_lv1_H, function(x){ifelse(is.na(x),NA, d_H[x])}))
    gex_db$j_H_lv1= unlist(lapply(match.index_lv1_H, function(x){ifelse(is.na(x),NA, j_H[x])}))
    gex_db$cdr3aa_H_lv1= unlist(lapply(match.index_lv1_H, function(x){ifelse(is.na(x),NA, cdr3aa_H[x])}))
    
    
    gex_db$v_L_lv1= unlist(lapply(match.index_lv1_L, function(x){ifelse(is.na(x),NA, v_L[x])})) 
    gex_db$j_L_lv1= unlist(lapply(match.index_lv1_L, function(x){ifelse(is.na(x),NA, j_L[x])}))
    gex_db$cdr3aa_L_lv1= unlist(lapply(match.index_lv1_L, function(x){ifelse(is.na(x),NA, cdr3aa_L[x])}))
    
    # creat a new data frame with label.main and Single.R.label
    hl_pair_db_lv1<- data.frame(gex_db$v_H_lv1,gex_db$v_L_lv1,gex_db$SingleR.labels)
    hl_pair_db1_lv1 <- hl_pair_db_lv1 %>% mutate(v_hl=paste0(gex_db.v_H_lv1,gex_db.v_L_lv1))
    hl_pair_db1_lv1$v_hl<- gsub("NA", "",hl_pair_db1_lv1$v_hl) 
    
    # SAVE the sc with BCR integration with lv1 correcte match
    saveRDS(object =gex_db, file = paste0(outputSummary,sampleName,"sc_BCR_integration_lv1.rds"))
    
    
    hl_pair_db1_lv1[is.na(hl_pair_db1_lv1)] <- ""
    head(hl_pair_db1_lv1)
    
    hl_pair_lv1 <- hl_pair_db1_lv1 %>%  
      mutate(pair=ifelse(nchar(v_hl)==0,"ND",ifelse(nchar(gex_db.v_L_lv1)==0,"H",ifelse(nchar(gex_db.v_H_lv1)==0,"L","HL"))))
    
    head(hl_pair_lv1)
    
    hl_pair_lv1$pair <- factor(hl_pair_lv1$pair ,levels = c("ND","L","H","HL"))
    hl_pair_lv1 <- rownames_to_column(hl_pair_lv1)
    hl_pair_lv1 <- dplyr::rename(hl_pair_lv1, unique_id=rowname)
    hl_pair_lv1_sum<- hl_pair_lv1 %>% group_by(gex_db.SingleR.labels,pair) %>%
      dplyr::summarise(count=n()) %>% mutate(percent=round(count/sum(count),digits =3)*100)
    head(hl_pair_lv1_sum)
    
    
    #############################################################################
    # extract pure B cells
    #############################################################################
    
    ###########
    # subset b cells
    ###########
    B.fine<- grep("B cells|Plasmablasts",unique(hl_pair_lv1_sum$gex_db.SingleR.labels),value = TRUE)
    hl_pair_lv1_sum_pureB <- hl_pair_lv1_sum %>% filter(gex_db.SingleR.labels %in% B.fine )
    hl_pair_lv1_sum_pureB 
    
    # total b cells
    hl_pair_lv1$sample<- sampleName
    
    hl_pair_lv1_pureB <- hl_pair_lv1 %>% filter(gex_db.SingleR.labels %in% B.fine ) %>% group_by(pair,sample) %>% 
      dplyr::summarise(count=n()) %>% ungroup() %>%mutate(percent=round(count/sum(count),digits = 3)*100)
    head(hl_pair_lv1_pureB)
    
    finalMatchB_lv1[[i]]<- hl_pair_lv1_pureB
    
  }
  
  finalMatchB_lv1 <- do.call(rbind,finalMatchB_lv1)
  finalMatchB_lv1$sample <- factor(finalMatchB_lv1$sample,
                                   levels = c("HD2","HD1","HD1_stD2","p36D4","p36D45","p36D169"))
  
  # save the data for matching plot
  write_tsv(x=finalMatchB_lv1, file =paste0(outputSummary,"fig5b_sourceData.tsv") )
  
  ###########
  # barplot
  ###########
  
  colours<- c( "ND"="grey", "L"="cadetblue1",  "H"="mediumorchid1", "HL"="green4")
  #show_col(colours)
  
  #--------------------------------------
  # https://stackoverflow.com/questions/70791096/x-axis-labels-cut-off-in-ggplot-when-rotating
  #--------------------------------------
  margin_spacer <- function(x) {
    # where x is the column in your dataset
    left_length <- nchar(levels(factor(x)))[1]
    if (left_length > 8) {
      return( left_length +5)
    }
    else
      return(0)
  }
  
  
  p_all_pairing_fill_lv1_no2 <- ggplot(finalMatchB_lv1, aes(fill=pair, y=percent, x=sample)) + 
    geom_bar(position="stack", stat="identity") + 
    ggtitle(paste0(sampleName," B cell sc with matched BCR(lv1)")) +
    geom_text(aes( y=percent, x=sample,
                   label = paste0(round(percent,digits = 1),"%")),size = 4, position = position_stack(vjust = 0.5))+ theme(
                     plot.title = element_text(color="black", size=8, face="bold"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust=1,size=10, face="bold"),
                     axis.text.y = element_text( size=10, face="bold"),
                     plot.margin = margin(l = 0 + margin_spacer(finalMatchB_lv1$sample),r = 2,t=2)
                   )  + scale_fill_manual(values = colours)
  
  print(p_all_pairing_fill_lv1_no2)
  
  ggsave(filename=paste0(outputPath,"LV1_all_pairing_combine_pt2.png"), width=12, height=10,
         plot=p_all_pairing_fill_lv1_no2,units="cm",dpi=300)
  
  
  
}





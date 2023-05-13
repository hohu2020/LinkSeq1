

#==============================================================================#
#  function: combineSc_Bcr_singleR
#==============================================================================#

combineSc_Bcr_singleR<- function(file_list, output_summary,inputBcr,inputSc){
  all_sc<- list()
  for(i in seq_along(fileList)){

    ###############################
    # singleR assignment for sc data
    ###############################
    ############
    # read Sc data
    ############

    name.i<- fileList[i]
    name.i
    gex_db <- readRDS(paste0(rdsInputPath,name.i,".rds"))

    # load ref data for singleR prediction
    ref_Monaco <- MonacoImmuneData()
    minPercent.mt<- 15
    gex_db<- subset(gex_db, subset = nFeature_RNA > 150 & nFeature_RNA < 2000& percent.mt < minPercent.mt )
    # https://www.biostars.org/p/9479638/
    # predict gex

    pred.gex_db<-SingleR(GetAssayData(gex_db, assay = "RNA", slot = "data"), ref = ref_Monaco, labels = ref_Monaco$label.main)

    pred.gex_db_fine<-SingleR(GetAssayData(gex_db, assay = "RNA", slot = "data"), ref = ref_Monaco, labels = ref_Monaco$label.fine)
    # assign the label.main to p36D1 meta data
    gex_db$label.main<- pred.gex_db$labels
    gex_db$SingleR.labels <- pred.gex_db_fine$labels


    if( name.i=="p36T1"){
      sample_id<-"p36_T1"
      sampleName<-"p36D4"
      dpi<- 4
    }else if( name.i=="p36T2"){
      sample_id<-"p36_T2"
      dpi<- 45
      sampleName<-"p36D45"
    }else if( name.i=="p36T4"){
      sample_id<-"p36_T4"
      dpi<- 169
      sampleName<-"p36D169"
    }else if( name.i=="p45T2"){
      sample_id<-"p45_T2"
      dpi<- 24
      sampleName<-"p45D24"
    }else if( name.i=="p45T5"){
      sample_id<-"p45_T5"
      dpi<- 164
      sampleName<-"p45D164"
    }else if( name.i=="HD2B1500"){
      sample_id<-"HD2B1500"
      dpi<- 0
      sampleName<-"HD1"
    }else if( name.i=="HD2B3600"){
      sample_id<-"HD2B3600"
      dpi<- 0
      sampleName<-"HD2"
    }else if( name.i=="HD2B960"){
      sample_id<-"HD2B960"
      dpi<- 0
      sampleName<-"HD1_stD2"
    }else if( name.i=="p16T1B"){
      sample_id<-"p16T1B"
      dpi<- 9
      sampleName<-"p16D9"
    }else if( name.i=="p18T1"){
      sample_id<-"p18T1A08"
      dpi<- 13
      sampleName<-"p18D13"
    }else {
      print("now only support COVID and HD bcr data")
    }


    gex_db$sample<- sampleName
    gex_db$id<- sample_id
    gex_db$orig.ident<- NULL
    Idents(gex_db)
    str(gex_db)



    if(name.i=="p45T5"|name.i=="p45T2"|name.i=="p16T1B"|name.i=="p18T1"){
      # no bcr to match sc
      table(gex_db$label.main)
      table(gex_db$SingleR.labels)

      # save the singleR annotated file
      write_rds(x=gex_db,file =paste0(output_summary,name.i,"_singleR_AnnoSC_HX.rds") )


    }else{
      ###############################
      # add BCR lv1 heavy V,CDR3,J to sc meta data
      ###############################

      ############
      # read BCR
      ############
      # Read BCR H data

      bcr_db_H <- read_tsv(paste0(bcrInputPath, name.i,"_lv1CorrectHeavyBcr.tsv"))
      bcr_db_H$v_call <- sapply(str_split(bcr_db_H$v_call, pattern = "\\*"),"[",1)
      bcr_db_H$j_call <- sapply(str_split(bcr_db_H$j_call, pattern = "\\*"),"[",1)
      bcr_db_H<-within(bcr_db_H,
                       cdr3aa<-substr(junction_aa,2,nchar(junction_aa)-1) )

      # Read BCR L data
      bcr_db_L <- read_tsv(paste0(bcrInputPath, name.i,"_lv1CorrectLightBcr.tsv"))
      bcr_db_L$v_call <- sapply(str_split(bcr_db_L$v_call, pattern = "\\*"),"[",1)
      bcr_db_L$j_call <- sapply(str_split(bcr_db_L$j_call, pattern = "\\*"),"[",1)
      bcr_db_L<-within(bcr_db_L,
                       cdr3aa<-substr(junction_aa,2,nchar(junction_aa)-1) )
      # Add meta infor to BCR
      bcr_db_H$sample<- sampleName
      bcr_db_H$id<- sample_id
      bcr_db_H$day<- dpi

      bcr_db_L$sample<- sampleName
      bcr_db_L$id<- sample_id
      bcr_db_L$day<- dpi

      ############
      # Match barcode
      ############
      # NOTE: cell_id becomes the cluster_cb after lv1 correction
      bcr_db_H$cell_id_unique = paste0(bcr_db_H$id, "_", bcr_db_H$cell_id)
      bcr_db_H$cluster_cb_unique = paste0(bcr_db_H$id, "_", bcr_db_H$cluster_cb)
      bcr_db_H$cell_id_unique[1:4]

      bcr_db_L$cell_id_unique = paste0(bcr_db_L$id, "_", bcr_db_L$cell_id)
      bcr_db_L$cluster_cb_unique = paste0(bcr_db_L$id, "_", bcr_db_L$cluster_cb)
      bcr_db_L$cell_id_unique[1:4]
      ###
      # LV1 based matching
      ###
      # match index to find the position of the GEX cells in the BCR data
      match.index_lv1_H = match(Cells(gex_db), bcr_db_H$cluster_cb_unique)
      match.index_lv1_L= match(Cells(gex_db), bcr_db_L$cluster_cb_unique)

      ############
      # Add BCR h and L v,d.j and CDR3 aa to sc data as meta data
      ############
      v_H <- bcr_db_H$v_call
      d_H <- bcr_db_H$d_call
      j_H <- bcr_db_H$j_call
      cdr3aa_H<- bcr_db_H$cdr3aa
      cluster_cb_unique_H <-bcr_db_H$cluster_cb_unique

      v_L <- bcr_db_L$v_call
      j_L <- bcr_db_L$j_call
      cdr3aa_L<- bcr_db_L$cdr3aa

      ############
      # lv1 based matching
      ############

      gex_db$v_H_lv1= unlist(lapply(match.index_lv1_H, function(x){ifelse(is.na(x),NA, v_H[x])}))
      gex_db$d_H_lv1= unlist(lapply(match.index_lv1_H, function(x){ifelse(is.na(x),NA, d_H[x])}))
      gex_db$j_H_lv1= unlist(lapply(match.index_lv1_H, function(x){ifelse(is.na(x),NA, j_H[x])}))
      gex_db$cdr3aa_H_lv1= unlist(lapply(match.index_lv1_H, function(x){ifelse(is.na(x),NA, cdr3aa_H[x])}))
      gex_db$cluster_cb_unique= unlist(lapply(match.index_lv1_H, function(x){ifelse(is.na(x),NA, cluster_cb_unique_H[x])}))


      gex_db$v_L_lv1= unlist(lapply(match.index_lv1_L, function(x){ifelse(is.na(x),NA, v_L[x])}))
      gex_db$j_L_lv1= unlist(lapply(match.index_lv1_L, function(x){ifelse(is.na(x),NA, j_L[x])}))
      gex_db$cdr3aa_L_lv1= unlist(lapply(match.index_lv1_L, function(x){ifelse(is.na(x),NA, cdr3aa_L[x])}))

      ############
      # save the singleR annotated file
      ############
      write_rds(x=gex_db,file =paste0(output_summary,name.i,"_singleR_AnnoSC_HX.rds") )

    }

    all_sc[[i]]<- gex_db

  }
  all_sc <- setNames(all_sc,nm=newName)
  sapply(all_sc, dim)
  return(all_sc)

}








#==============================================================================#
# general function for dotplot
# input :
#         -sc data()
#         -Dot plot for features
#==============================================================================#

# FUNCTION
customMarkers_bySampleByCluster<- function(df,markers.to.plot,outputPlotPath,threshold.cell){

  ######################
  # data formatting: extract pure b cells
  ######################
  My_Theme_Dot = theme(
    plot.title = element_blank(),
    axis.title.x =element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle =90,   hjust=1, size=20, face="bold"),
    axis.text.y =  element_text( size=20, face="bold"))

  #Idents(immune.combinedB) <- immune.combinedB$label.main
  immune.combined_test <- df

  Idents(immune.combined_test)<- immune.combined_test$SingleR.labels

  immune.combined_fine_B <- subset(immune.combined_test,
                                   idents=grep(pattern = "B cells|Plasmablasts", x = immune.combined_test$SingleR.labels,value = T))
  Idents(immune.combined_fine_B)<- immune.combined_fine_B$sample

  #all_main_colours<- c( "HD1"="#b2df8a",  "HD2"="#762a83", "HD1_stD2"="#08519c", "p36D4"="#e31a1c","p36D45"="#9ecae1","p36D169"="#ff7f00","p45D24"="#c51b7d","p45D164"="gray95")

  all_main_colours<- c( "#b2df8a", "#762a83", "#08519c", "#e31a1c","#9ecae1","#ff7f00","#c51b7d","gray95","red","green","blue")


  immune.combined_fine_B$sample <- factor(immune.combined_fine_B$sample ,
                                          levels =c("HD2","HD1","HD1_stD2","p16D9","p18D13","p36D4","p36D45","p36D169","p45D24","p45D164")  )

  #levels(immune.combined_fine_B$sample)<- c("HD1","HD2","HD1_stD2","p36D4","p36D45","p36D169","p45D24","p45D164")
  #Idents(immune.combined_fine_B) <- factor(Idents(immune.combined_fine_B),  levels = c("Exhausted B cells","Switched memory B cells","Non-switched memory B cells","Naive B cells"))

  ###########################################
  # PLOT BY SAMPLE
  ###########################################
  ######################
  # dot plot by samples
  ######################
  dimplot<- DotPlot(immune.combined_fine_B, assay="RNA",features = markers.to.plot, dot.scale = 16, group.by = "sample", cluster.idents=T) +
    RotatedAxis()
  dimplot
  #------------------------------------------
  # https://www.biostars.org/p/446743/
  #------------------------------------------
  dimplot$data$id <- factor(x = dimplot$data$id, levels = c("HD2", "HD1", "HD1_stD2","p16D9","p18D13","p36D4","p36D45","p36D169" , "p45D24", "p45D164"))
  dimplot_flip<- dimplot +  coord_flip() & My_Theme_Dot
  dimplot_flip


}












#==============================================================================#
# heatmap_VJ
# df_fileList: file list ; df_bcrInputPath: lv1_corrected BCR; rdsInputPath: QCed sc data
# type: "heavy"/"light"/"HL"
#==============================================================================#
heatmap_VJ <- function(df_fileList,df_bcrInputPath,rdsInputPath,outputPath,type,vjSummary){
  final<- list()
  for(i in df_fileList){

    # Read GEX data
    gex_db <- readRDS(paste0(rdsInputPath,i,".rds"))

    # load ref data for singleR prediction
    ref_Monaco <- MonacoImmuneData()
    minPercent.mt<- 15
    gex_db<- subset(gex_db, subset = nFeature_RNA > 150 & nFeature_RNA < 2000& percent.mt < minPercent.mt )
    # https://www.biostars.org/p/9479638/
    # predict gex

    pred.gex_db<-SingleR(GetAssayData(gex_db, assay = "RNA", slot = "data"), ref = ref_Monaco, labels = ref_Monaco$label.main)

    pred.gex_db_fine<-SingleR(GetAssayData(gex_db, assay = "RNA", slot = "data"), ref = ref_Monaco, labels = ref_Monaco$label.fine)
    # assign the label.main to p36D1 meta data
    gex_db$label.main<- pred.gex_db$labels
    gex_db$SingleR.labels <- pred.gex_db_fine$labels

    if(type=="heavy"){
      # Read BCR data
      bcr_db <- read_tsv(paste0(df_bcrInputPath,i,"_lv1CorrectHeavyBcr.tsv"))
      ## Standardize cell IDs

      if(i=="p36T1"){
        sample_id<-"p36_T1"
      }else if(i=="p36T2"){
        sample_id<-"p36_T2"
      }else if(i=="p36T3"){
        sample_id<-"p36_T3"
      }else if(i=="p36T4"){
        sample_id<-"p36_T4"
      }else if(i=="p45T2"){
        sample_id<-"p45_T2"
      }else if(i=="p36T4ep300"){
        sample_id<-"p36_T4"
      }else if(i=="p36T4ep500nt"){
        sample_id<-"p36_T4"
      }else if(i=="HD2B960"){
        sample_id<-"HD2B960"
      }else if(i=="HD2B3600"){
        sample_id<-"HD2B3600"
      }else if(i=="HD2B1500"){
        sample_id<-"HD2B1500"
      }else {
        print("now only support p36 time points")
      }

      bcr_db$sample<- sample_id
      # Make cell IDs in BCR match those in Seurat Object
      # bcr_db$cell_id = sapply(strsplit(bcr_db$sequence_id,"_"), "[",1)

      # NOTE: cell_id becomes the cluster_cb after lv2 correction
      bcr_db$cell_id_unique = paste0(bcr_db$sample, "_", bcr_db$cell_id)
      bcr_db$cluster_cb_unique = paste0(bcr_db$sample, "_", bcr_db$cluster_cb)

      bcr_db$cell_id_unique[1:4]

      # subset heavy bcr for match to have unique cell_id(because paired light chain will have the same cell_id)
      # Note: lv2 correction already select the heavy chain
      # bcr_db <- bcr_db[grep("IGH",bcr_db$v_call ,invert= F),]

      #####################
      # matching without lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index = match(Cells(gex_db), bcr_db$cell_id_unique)

      #####################
      # matching after lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index_lv2 = match(Cells(gex_db), bcr_db$cluster_cb_unique)


      # In this data, not all cells are B cells.
      # What proportion of cells don’t have BCRs?
      scWithOutMatchedBcr<- NULL
      scWithOutMatchedBcr$scMissBcr <- mean(!is.na(match.index))
      scWithOutMatchedBcr$scMissBcr_lv2 <- mean(!is.na(match.index_lv2))
      scWithOutMatchedBcr$id <-i

      final[[i]]<- scWithOutMatchedBcr


      #------------------------------------------------------------------------
      ### Find the BCR cells in the GEX data

      # Match indices to find the position of the BCR cells in the GEX data
      # Different from finding the position of the GEX cells in the BCR data!

      match.index = match(bcr_db$cluster_cb_unique, Cells(gex_db))
      # What proportion of BCRs don’t have GEX information?
      mean(is.na(match.index))
      mean(!is.na(match.index))
      ### Transfer GEX annotations into the BCR data
      # Add annotations to BCR data
      cell.annotation = as.character(gex_db@meta.data$SingleR.labels)
      bcr_db$gex_annotation= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.annotation[x])}))
      bcr_db$gex_annotation[1:5]

      # Add UMAP coordinates to BCR data
      #umap1 = gex_db@reductions$umap@cell.embeddings[,1]
      #umap2 = gex_db@reductions$umap@cell.embeddings[,2]
      #bcr_db$gex_umap1= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap1[x])}))
      #bcr_db$gex_umap2= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap2[x])}))
      bcr_db[1:5,] %>%
        select(cell_id,gex_annotation)



      ### Remove cells without GEX data
      # Remove cells that didn’t match
      bcr_db = filter(bcr_db, !is.na(gex_annotation))


      ################################################################
      # UMAP cluster labeling need to be customized for individual sample
      # !!!!!!!!!!!!!!!!!!!!!!!!!
      ################################################################
      ### Ensure information transferred from Seurat object
      col_anno = c("Exhausted B cells"="dodgerblue2", "Naive B cells"="firebrick2", "Non-switched memory B cells"="seagreen","Plasmablasts"="darkgoldenrod2", "Switched memory B cells"="black")
      # Plot UMAP from bcr_db
      #pdf(paste0(outputPath,i, "umap_from_bcr_db.pdf"))

      #bcr_umap <- ggplot(bcr_db) +geom_point(aes(x = gex_umap1, y = gex_umap2, color = gex_annotation)) +scale_colour_manual(values=col_anno) +theme_bw()

      #print(bcr_umap)
      #dev.off()
      #----------------------------------------------------------------------
      # bcr_dbH <-bcr_db[grep("IGH",bcr_db$v_call,invert = F),] # It's already Heavy chain after lv2 correction


      bcr_hm1 <- bcr_db  %>% select(sequence_id,v_call,j_call,gex_annotation)
      # only use the 1st v gene hit
      bcr_hm1$v_call <- sapply(str_split(bcr_hm1$v_call, pattern = "-"),"[",1)
      bcr_hm1$j_call <- sapply(str_split(bcr_hm1$j_call, pattern = "\\*"),"[",1)
      head(bcr_hm1)
      bcr_hm1$v_call <- factor(bcr_hm1$v_call, levels = c("IGHV1", "IGHV2", "IGHV3", "IGHV4", "IGHV5", "IGHV6", "IGHV7" ))
      bcr_hm1$j_call <- factor(bcr_hm1$j_call, levels = c("IGHJ1", "IGHJ2", "IGHJ3", "IGHJ4", "IGHJ5", "IGHJ6"))



      #----------
      # https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
      # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
      freq <- table(bcr_hm1$v_call,bcr_hm1$j_call)
      freq

      # write the IGHV AND IGHJ frequency
      cellNumber= dim(bcr_hm1)
      totalNumber=length(gex_db$label.main)
      IGHV_Count=rowSums(freq)
      IGJ_Count=colSums(freq)
      i.sum <- list(totalNumber,cellNumber,IGHV_Count, IGJ_Count)

      lapply(i.sum, function(x) write.table( data.frame(x), file =paste0(vjSummary,i,"_IGH_VJ_countSummary.csv") , append= T, sep=',' ))

      # scale the row or column before the heatmap ploting
      cal_z_score <- function(x){(x - mean(x)) / sd(x)}

      data_subset_norm <- t(apply(freq, 1, cal_z_score))

      #################
      # prepare row and column annotation bar plot
      #################
      # row and column data for barplot
      v_h<-as.data.frame(table(bcr_hm1$v_call))
      j_h<-as.data.frame(table(bcr_hm1$j_call))
      # topN for labeling

      # topN < group_by(vj) %>% summarise(n=n()) %>%  arrange(desc(n)) %>%  top_n()
      topN <- rownames_to_column(v_h) %>% arrange(desc(Freq)) %>%top_n(10,v_h$Freq)
      topN_nb<- as.numeric(topN$rowname)


      # barplot annotation
      v_r_ha = rowAnnotation(
        V_count = anno_barplot(v_h$Freq),foo = anno_mark(at = topN_nb, labels = topN$Var1), annotation_name_rot = 0)
      v_r_ha1 = rowAnnotation(
        V_count = anno_barplot(v_h$Freq))

      v_c_ha = HeatmapAnnotation(J_count = anno_barplot(j_h$Freq),which = "col")



      #-------------------------------------------------------------------------
      # heatmap


      #pheatmap(data_subset_norm)
      # use complex heat map to remove the cluster of row and column
      # Turn off the clustering fucntion


      #--------------------------------------------------------------------
      # adjust legend size
      # https://github.com/jokergoo/ComplexHeatmap/issues/3
      col_fun = colorRamp2(c(-5, 0, 4), c("green", "white", "red"))
      # lgd = list(title = "Freq", col_fun = col_fun,legend_gp = gpar(fontsize=20), labels_gp = gpar(font=200))
      lgd = list(title = "Freq", col_fun = col_fun,labels_gp = gpar(fontsize=14))




      #----------------------------------------------------------------------



      col_fun1a = colorRamp2(c(0, max(v_h$Freq)/100, max(v_h$Freq)/4), c(  "white", "orange","red"))
      ph1a<- Heatmap(freq, column_title = paste0(i," heavy VJ gene usage"),
                     name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun1a,
                     row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                     top_annotation = v_c_ha, right_annotation =v_r_ha1,
                     row_names_side = "left",heatmap_legend_param =lgd,
                     width = ncol(freq)*unit(10, "mm"),
                     height = nrow(freq)*unit(10, "mm")) # turn off row clustering

      draw(ph1a)

      #-------------------------------------------------------------------------


      # save as png file
      png(filename=paste0(outputPath,i,"_",type,"_main_VJ.png"),width=16,height=16,units="cm",res=300)
      draw(ph1a)
      dev.off()



    }else if(type=="light"){
      # Read BCR data
      bcr_db <- read_tsv(paste0(df_bcrInputPath,i,"_lv1CorrectLightBcr.tsv"))
      ## Standardize cell IDs

      if(i=="p36T1"){
        sample_id<-"p36_T1"
      }else if(i=="p36T2"){
        sample_id<-"p36_T2"
      }else if(i=="p36T3"){
        sample_id<-"p36_T3"
      }else if(i=="p36T4"){
        sample_id<-"p36_T4"
      }else if(i=="p45T2"){
        sample_id<-"p45_T2"
      }else if(i=="p36T4ep300"){
        sample_id<-"p36_T4"
      }else if(i=="p36T4ep500nt"){
        sample_id<-"p36_T4"
      }else if(i=="HD2B960"){
        sample_id<-"HD2B960"
      }else if(i=="HD2B3600"){
        sample_id<-"HD2B3600"
      }else if(i=="HD2B1500"){
        sample_id<-"HD2B1500"
      }else {
        print("now only support p36 time points")
      }

      bcr_db$sample<- sample_id
      # Make cell IDs in BCR match those in Seurat Object
      # bcr_db$cell_id = sapply(strsplit(bcr_db$sequence_id,"_"), "[",1)

      # NOTE: cell_id becomes the cluster_cb after lv2 correction
      bcr_db$cell_id_unique = paste0(bcr_db$sample, "_", bcr_db$cell_id)
      bcr_db$cluster_cb_unique = paste0(bcr_db$sample, "_", bcr_db$cluster_cb)

      bcr_db$cell_id_unique[1:4]

      # subset heavy bcr for match to have unique cell_id(because paired light chain will have the same cell_id)
      # Note: lv2 correction already select the heavy chain
      # bcr_db <- bcr_db[grep("IGH",bcr_db$v_call ,invert= F),]

      #####################
      # matching without lv2 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index = match(Cells(gex_db), bcr_db$cell_id_unique)

      #####################
      # matching after lv2 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index_lv2 = match(Cells(gex_db), bcr_db$cluster_cb_unique)


      # In this data, not all cells are B cells.
      # What proportion of cells don’t have BCRs?
      scWithOutMatchedBcr<- NULL
      scWithOutMatchedBcr$scMissBcr <- mean(!is.na(match.index))
      scWithOutMatchedBcr$scMissBcr_lv2 <- mean(!is.na(match.index_lv2))
      scWithOutMatchedBcr$id <-i

      final[[i]]<- scWithOutMatchedBcr


      #------------------------------------------------------------------------
      ### Find the BCR cells in the GEX data

      # Match indices to find the position of the BCR cells in the GEX data
      # Different from finding the position of the GEX cells in the BCR data!

      match.index = match(bcr_db$cluster_cb_unique, Cells(gex_db))
      # What proportion of BCRs don’t have GEX information?
      mean(is.na(match.index))
      mean(!is.na(match.index))
      ### Transfer GEX annotations into the BCR data
      # Add annotations to BCR data
      cell.annotation = as.character(gex_db@meta.data$SingleR.labels)
      bcr_db$gex_annotation= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.annotation[x])}))
      bcr_db$gex_annotation[1:5]

      # Add UMAP coordinates to BCR data
      #umap1 = gex_db@reductions$umap@cell.embeddings[,1]
      #umap2 = gex_db@reductions$umap@cell.embeddings[,2]
      #bcr_db$gex_umap1= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap1[x])}))
      #bcr_db$gex_umap2= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap2[x])}))
      bcr_db[1:5,] %>%
        select(cell_id,gex_annotation)



      ### Remove cells without GEX data
      # Remove cells that didn’t match
      bcr_db = filter(bcr_db, !is.na(gex_annotation))


      ################################################################
      # UMAP cluster labeling need to be customized for individual sample
      # !!!!!!!!!!!!!!!!!!!!!!!!!
      ################################################################
      ### Ensure information transferred from Seurat object
      col_anno = c("Exhausted B cells"="dodgerblue2", "Naive B cells"="firebrick2", "Non-switched memory B cells"="seagreen","Plasmablasts"="darkgoldenrod2", "Switched memory B cells"="black")
      # Plot UMAP from bcr_db
      #pdf(paste0(outputPath,i, "umap_from_bcr_db.pdf"))

      #bcr_umap <- ggplot(bcr_db) +geom_point(aes(x = gex_umap1, y = gex_umap2, color = gex_annotation)) +scale_colour_manual(values=col_anno) +theme_bw()

      #print(bcr_umap)
      #dev.off()
      #----------------------------------------------------------------------

      # bcr_dbH <-bcr_db[grep("IGH",bcr_db$v_call,invert = F),] # It's already Heavy chain after lv2 correction


      bcr_hm1 <- bcr_db  %>% select(sequence_id,v_call,j_call,gex_annotation)
      # only use the 1st v gene hit
      bcr_hm1$v_call <- sapply(str_split(bcr_hm1$v_call, pattern = "-"),"[",1)

      # replace D in all light chain k v
      bcr_hm1$v_call <-  str_remove(bcr_hm1$v_call,pattern = "D")
      bcr_hm1$j_call <- sapply(str_split(bcr_hm1$j_call, pattern = "\\*"),"[",1)
      head(bcr_hm1)
      # add level to the V and J gene
      bcr_hm1$v_call <- factor(bcr_hm1$v_call, levels = c("IGKV1","IGKV2","IGKV3",  "IGKV4",  "IGKV5",  "IGKV6", "IGKV7", "IGLV1",  "IGLV2",  "IGLV3",  "IGLV4",  "IGLV5",  "IGLV6",  "IGLV7",  "IGLV8",  "IGLV9","IGLV10" ,"IGLV11" ))
      bcr_hm1$j_call <- factor( bcr_hm1$j_call , levels = c("IGKJ1", "IGKJ2", "IGKJ3", "IGKJ4", "IGKJ5", "IGLJ1", "IGLJ2", "IGLJ3","IGLJ4","IGLJ5","IGLJ6", "IGLJ7" ) )


      freq <- table(bcr_hm1$v_call,bcr_hm1$j_call)

      # write the IGHV AND IGHJ frequency
      totalNumber=length(gex_db$label.main)
      cellNumber= dim(bcr_hm1)
      IGV_Count=rowSums(freq)
      IGJ_Count=colSums(freq)
      i.sum <- list(totalNumber,cellNumber,IGV_Count, IGJ_Count)

      lapply(i.sum, function(x) write.table( data.frame(x), file =paste0(vjSummary,i,"_light_VJ_countSummary.csv") , append= T, sep=',' ))
      #----------
      # https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
      # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/


      # scale the row or column before the heatmap ploting
      cal_z_score <- function(x){(x - mean(x)) / sd(x)}
      data_subset_norm <- t(apply(freq, 1, cal_z_score))


      #################
      # prepare row and column annotation bar plot
      #################
      # row and column data for barplot
      v_l<-as.data.frame(table(bcr_hm1$v_call))
      j_l<-as.data.frame(table(bcr_hm1$j_call))
      head(v_l)
      # topN for labeling

      # topN < group_by(vj) %>% summarise(n=n()) %>%  arrange(desc(n)) %>%  top_n()
      topN <- rownames_to_column(v_l) %>% arrange(desc(Freq)) %>%top_n(10,v_l$Freq)
      topN_nb<- as.numeric(topN$rowname)
      head(topN_nb)

      # barplot annotation
      v_r_ha = rowAnnotation(
        V_count = anno_barplot(v_l$Freq),foo = anno_mark(at = topN_nb, labels = topN$Var1), annotation_name_rot = 270)

      v_r_ha1 = rowAnnotation(
        V_count = anno_barplot(v_l$Freq) )

      v_c_ha = HeatmapAnnotation(J_count = anno_barplot(j_l$Freq),which = "col")

      # heatmap

      #pheatmap(data_subset_norm)
      # use complex heat map to remove the cluster of row and column
      # Turn off the clustering fucntion
      col_fun = colorRamp2(c(-5, 0, 4), c("green", "white", "red"))
      lgd = list(title = "Freq", col_fun = col_fun ,labels_gp = gpar(fontsize=14))



      # for rotated light heatmap
      # top annotation
      #v_c_ha1 = HeatmapAnnotation(J_count = anno_barplot(j_l$Freq),which = "col")
      v_c_ha_a = HeatmapAnnotation(V_count = anno_barplot(v_l$Freq),which = "col")
      # side annotation
      #v_r_ha1 = rowAnnotation(V_count = anno_barplot(v_l$Freq) )
      v_r_ha_a = rowAnnotation(J_count = anno_barplot(j_l$Freq) )
      # plot
      t_freq <- t(freq)
      col_fun1a = colorRamp2(c(0, sum(v_l$Freq)/500, max(v_l$Freq)/5), c(  "white", "orange","red"))
      ph1b<- Heatmap(t_freq, column_title = paste0(i," light VJ gene usage"),
                     name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun1a,
                     row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                     top_annotation = v_c_ha_a, right_annotation =v_r_ha_a,
                     row_names_side = "left",heatmap_legend_param =lgd,
                     width = ncol(t_freq)*unit(5, "mm"),
                     height = nrow(t_freq)*unit(5, "mm")) # turn off row clustering

      draw(ph1b)


      #-------------------------------------------------------------------------



      png(filename=paste0(outputPath,i,"_",type,"_main_VJ.png"),width=20,height=16,units="cm",res=300)
      draw(ph1b)
      dev.off()

    }else if(type=="HL"){

      #########################################################################
      # type="HL"
      #########################################################################

      # Read BCR data
      bcr_dbH <- read_tsv(paste0(df_bcrInputPath,i,"_lv1CorrectHeavyBcr.tsv"))
      bcr_dbL <- read_tsv(paste0(df_bcrInputPath,i,"_lv1CorrectLightBcr.tsv"))
      ## Standardize cell IDs

      if(i=="p36T1"){
        sample_id<-"p36_T1"
      }else if(i=="p36T2"){
        sample_id<-"p36_T2"
      }else if(i=="p36T3"){
        sample_id<-"p36_T3"
      }else if(i=="p36T4"){
        sample_id<-"p36_T4"
      }else if(i=="p45T2"){
        sample_id<-"p45_T2"
      }else if(i=="p36T4ep300"){
        sample_id<-"p36_T4"
      }else if(i=="p36T4ep500nt"){
        sample_id<-"p36_T4"
      }else if(i=="HD2B960"){
        sample_id<-"HD2B960"
      }else if(i=="HD2B3600"){
        sample_id<-"HD2B3600"
      }else if(i=="HD2B1500"){
        sample_id<-"HD2B1500"
      }else {
        print("now only support p36 time points")
      }

      bcr_dbH$sample<- sample_id
      bcr_dbL$sample<- sample_id
      dim(bcr_dbH)
      dim(bcr_dbL)

      # Make cell IDs in BCR match those in Seurat Object
      # bcr_db$cell_id = sapply(strsplit(bcr_db$sequence_id,"_"), "[",1)

      # NOTE: cell_id becomes the cluster_cb after lv2 correction
      bcr_dbH$cell_id_unique = paste0(bcr_dbH$sample, "_", bcr_dbH$cell_id)
      bcr_dbH$cluster_cb_unique = paste0(bcr_dbH$sample, "_", bcr_dbH$cluster_cb)

      bcr_dbL$cell_id_unique = paste0(bcr_dbL$sample, "_", bcr_dbL$cell_id)
      bcr_dbL$cluster_cb_unique = paste0(bcr_dbL$sample, "_", bcr_dbL$cluster_cb)
      bcr_dbL$cell_id_unique[1:4]

      # subset heavy bcr for match to have unique cell_id(because paired light chain will have the same cell_id)
      # Note: lv2 correction already select the heavy chain
      # bcr_db <- bcr_db[grep("IGH",bcr_db$v_call ,invert= F),]

      #####################
      # matching without lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.indexH = match(Cells(gex_db), bcr_dbH$cell_id_unique)

      #####################
      # matching after lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index_lv2H = match(Cells(gex_db), bcr_dbH$cluster_cb_unique)


      # In this data, not all cells are B cells.
      # What proportion of cells don’t have BCRs?
      scWithOutMatchedBcrH<- NULL
      scWithOutMatchedBcrH$scMatchBcr <- mean(!is.na(match.indexH))
      scWithOutMatchedBcrH$scMatchBcr_lv2 <- mean(!is.na(match.index_lv2H))
      scWithOutMatchedBcrH$id <-i
      scWithOutMatchedBcrH

      final[[i]]<- scWithOutMatchedBcrH


      #------------------------------------------------------------------------
      ### Find the BCR cells in the GEX data

      # Match indices to find the position of the BCR cells in the GEX data
      # Different from finding the position of the GEX cells in the BCR data!

      match.indexH = match(bcr_dbH$cluster_cb_unique, Cells(gex_db))
      # What proportion of BCRs don’t have GEX information?
      mean(is.na(match.indexH))
      mean(!is.na(match.indexH))
      ### Transfer GEX annotations into the BCR data
      # Add annotations to BCR data
      cell.annotation = as.character(gex_db@meta.data$SingleR.labels)
      bcr_dbH$gex_annotation= unlist(lapply(match.indexH,function(x){ifelse(is.na(x),NA, cell.annotation[x])}))
      bcr_dbH$gex_annotation[1:5]

      # Add UMAP coordinates to BCR data
      #umap1 = gex_db@reductions$umap@cell.embeddings[,1]
      #umap2 = gex_db@reductions$umap@cell.embeddings[,2]
      #bcr_db$gex_umap1= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap1[x])}))
      #bcr_db$gex_umap2= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap2[x])}))
      bcr_dbH[1:5,] %>%
        select(cell_id_unique,gex_annotation)



      ### Remove cells without GEX data
      # Remove cells that didn’t match
      bcr_dbH = filter(bcr_dbH, !is.na(gex_annotation))


      ################################################################
      # UMAP cluster labeling need to be customized for individual sample
      # !!!!!!!!!!!!!!!!!!!!!!!!!
      ################################################################
      ### Ensure information transferred from Seurat object
      col_anno = c("Exhausted B cells"="dodgerblue2", "Naive B cells"="firebrick2", "Non-switched memory B cells"="seagreen","Plasmablasts"="darkgoldenrod2", "Switched memory B cells"="black")
      # Plot UMAP from bcr_db
      #pdf(paste0(outputPath,i, "umap_from_bcr_db.pdf"))

      #bcr_umap <- ggplot(bcr_db) +geom_point(aes(x = gex_umap1, y = gex_umap2, color = gex_annotation)) +scale_colour_manual(values=col_anno) +theme_bw()

      #print(bcr_umap)
      #dev.off()
      #----------------------------------------------------------------------
      # Format the V J gene usage for Heavy
      # only use the 1st v gene hit
      bcr_dbH$v_call <- sapply(str_split(bcr_dbH$v_call, pattern = "-"),"[",1)
      bcr_dbH$j_call <- sapply(str_split(bcr_dbH$j_call, pattern = "\\*"),"[",1)
      head(bcr_dbH[1:4])



      # Format the V J gene usage for light
      # only use the 1st v gene hit
      bcr_dbL$v_call <- sapply(str_split(bcr_dbL$v_call, pattern = "-"),"[",1)

      # replace D in all light chain k v
      bcr_dbL$v_call <-  str_remove(bcr_dbL$v_call,pattern = "D")

      bcr_dbL$j_call <- sapply(str_split(bcr_dbL$j_call, pattern = "\\*"),"[",1)
      head(bcr_dbL)

      # combine H L BCR
      bcr_dbH_sel <- bcr_dbH %>%
        select(sequence_id,v_call,j_call,gex_annotation,cluster_cb_unique) %>%
        mutate(vj=paste0(v_call,"_",j_call))

      bcr_dbL_sel <- bcr_dbL %>% select(sequence_id,v_call,j_call,cluster_cb_unique)%>%
        mutate(vj=paste0(v_call,"_",j_call))

      bcr_HL <- inner_join(x=bcr_dbH_sel,y=bcr_dbL_sel,
                           by=c("cluster_cb_unique"),suffix = c("_heavy","_light"))

      # add level to the HV and LV
      bcr_HL$v_call_light <- factor(bcr_HL$v_call_light, levels = c("IGKV1","IGKV2","IGKV3",  "IGKV4",  "IGKV5",  "IGKV6", "IGKV7", "IGLV1",  "IGLV2",  "IGLV3",  "IGLV4",  "IGLV5",  "IGLV6",  "IGLV7",  "IGLV8",  "IGLV9","IGLV10" ,"IGLV11" ))

      bcr_HL$v_call_heavy <- factor(bcr_HL$v_call_heavy, levels=c("IGHV1", "IGHV2", "IGHV3", "IGHV4", "IGHV5", "IGHV6", "IGHV7"))

      #########################################################################
      # part-1. HL VJ  heatmap
      #########################################################################

      bcr_HL_sel<- bcr_HL %>% select(vj_heavy, vj_light,gex_annotation) %>% mutate(HL_vj=paste0(vj_heavy,"_",vj_light))

      dim(bcr_HL_sel)


      #####################################
      # 1. HL VJ  heatmap
      # highlight top N based on pure heavy vj count
      #####################################
      # https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
      # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
      freq <- table(bcr_HL_sel$vj_heavy,bcr_HL_sel$vj_light)
      head(freq)

      # scale the row or column before the heatmap ploting
      cal_z_score <- function(x){
        (x - mean(x)) / sd(x)
      }

      data_subset_norm <- t(apply(freq, 1, cal_z_score))
      data_subset_norm

      #################
      # prepare row and column annotation bar plot
      #################
      # row and column data for barplot
      vj_H<-as.data.frame(table(bcr_HL_sel$vj_heavy))
      vj_L<-as.data.frame(table(bcr_HL_sel$vj_light))
      vj_HL<- as.data.frame(table(bcr_HL_sel$HL_vj))
      # topN for labeling

      # topN based on heavy chain vj
      topN <- rownames_to_column(vj_H) %>% arrange(desc(Freq)) %>%top_n(5,vj_H$Freq)
      topN_nb<- as.numeric(topN$rowname)


      # barplot annotation
      v_r_ha = rowAnnotation(
        H_VJc = anno_barplot(vj_H$Freq),foo = anno_mark(at = topN_nb, labels = topN$Var1), annotation_name_rot = 0)

      v_r_ha1 = rowAnnotation(H_VJc = anno_barplot(vj_H$Freq) )

      v_c_ha = HeatmapAnnotation(L_VJc = anno_barplot(vj_L$Freq),which = "col")
      #-------------------------------------------------------------------------
      # heatmap
      pdf(file=paste0(outputPath,i,"_lv1_HL_mainVJ_usageHeatmap.pdf"))

      #pheatmap(data_subset_norm)
      # use complex heat map to remove the cluster of row and column
      # Turn off the clustering fucntion
      col_fun = colorRamp2(c(-3, 0, 5), c("green", "white", "red"))
      lgd = list(title = "Freq", col_fun = col_fun ,labels_gp = gpar(fontsize=14))

      ph1<- Heatmap(data_subset_norm, column_title = paste0(i," heavy VJ vs light VJ gene usage(scaled)"),
                    name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                    top_annotation = v_c_ha, right_annotation =v_r_ha,
                    row_names_side = "left", heatmap_legend_param = lgd,
                    show_column_names = T, show_row_names = T,
                    width = ncol(data_subset_norm)*unit(2.5, "mm"),
                    height = nrow(data_subset_norm)*unit(2.5, "mm")) # turn off row clustering

      draw(ph1)

      col_fun1a = colorRamp2(c(0, max(vj_H$Freq)/100, max(vj_H$Freq)/4), c(  "white", "orange","red"))
      ph1a<- Heatmap(freq, column_title = paste0(i,"  heavy VJ vs light VJ gene usage"),
                     name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun1a,
                     row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                     top_annotation = v_c_ha, right_annotation =v_r_ha1,
                     row_names_side = "left",heatmap_legend_param =lgd,
                     width = ncol(freq)*unit(2.5, "mm"),
                     height = nrow(freq)*unit(2.5, "mm")) # turn off row clustering

      draw(ph1a)



      #####################################
      # 2. HL VJ  heatmap
      # top N based on paired HL vj count but only label the heavy VJ because, same heavy vj has different light vj
      #####################################

      # topN based on heavy and light paired  vj
      row_vjH <- rownames_to_column(vj_H)

      # combine heavy vj row order to the bcr data
      order_bcr_HL_sel<- inner_join(x=row_vjH, y=bcr_HL_sel, by=c("Var1"="vj_heavy"))
      head(order_bcr_HL_sel)

      # extract top n paired heavy and light vj
      topN_hl <- rownames_to_column(vj_HL) %>% arrange(desc(Freq)) %>%top_n(5,vj_HL$Freq)

      # extract the row number and corresponding heavy vj genes
      topN_hl_order<- order_bcr_HL_sel %>% filter(HL_vj%in%topN_hl$Var1)  %>%
        select(rowname,Var1)  %>% unique()
      topN_hl_order

      # extract the row number
      topN_nb_hl_order<- as.numeric(topN_hl_order$rowname)
      topN_nb_hl_order


      # project to heavy vj data for corresponding row number
      HL_vj_r_ha = rowAnnotation(
        H_vj_count = anno_barplot(vj_H$Freq),foo = anno_mark(at = topN_nb_hl_order, labels = topN_hl_order$Var1))


      lgd = list(title = "Freq", col_fun = col_fun ,labels_gp = gpar(fontsize=14))

      # plot
      ph2<- Heatmap(data_subset_norm, column_title = paste0(i," heavy VJ vs light VJ gene usage(scaled) \n topN paired VJ genes"),
                    name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                    top_annotation = v_c_ha, right_annotation =HL_vj_r_ha,
                    row_names_side = "left", heatmap_legend_param = lgd,
                    show_column_names = FALSE, show_row_names = FALSE,
                    width = ncol(data_subset_norm)*unit(2.5, "mm"),
                    height = nrow(data_subset_norm)*unit(2.5, "mm")) # turn off row clustering

      draw(ph2)

      col_fun1a = colorRamp2(c(0, max(vj_H$Freq)/500, max(vj_H$Freq)/4), c(  "white", "orange","red"))
      ph2a<- Heatmap(freq, column_title = paste0(i," heavy VJ vs light VJ gene usage(1/500-orange)"),
                     name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun1a,
                     row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                     top_annotation = v_c_ha, right_annotation =v_r_ha1,
                     row_names_side = "left",heatmap_legend_param =lgd,
                     width = ncol(freq)*unit(2.5, "mm"),
                     height = nrow(freq)*unit(2.5, "mm")) # turn off row clustering

      draw(ph2a)

      dev.off()


      #########################################################################
      # Part-2: heavy and light v gene heatmap
      #########################################################################
      bcr_HL_V_sel<- bcr_HL %>% select(v_call_heavy, v_call_light,gex_annotation)  %>%
        mutate(HL_v=paste0(v_call_heavy,"_",v_call_light))

      head(bcr_HL_V_sel)




      # https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
      # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
      freq_v <- table(bcr_HL_V_sel$v_call_heavy, bcr_HL_V_sel$v_call_light)
      freq_v

      # scale the row or column before the heatmap ploting
      cal_z_score <- function(x){
        (x - mean(x)) / sd(x)
      }

      data_subset_norm_v <- t(apply(freq_v, 1, cal_z_score))
      data_subset_norm_v


      #####################################
      # 1 HL V  heatmap
      # top N based on paired HL vj count but only label the heavy VJ because, same heavy vj has different light vj
      #####################################

      # use complex heat map to remove the cluster of row and column
      # Turn off the clustering fucntion
      col_fun = colorRamp2(c(-3, 0, 5), c("green", "white", "red"))
      #################
      # prepare row and column annotation bar plot
      #################
      # row and column data for barplot
      v_h<-as.data.frame(table(bcr_HL_V_sel$v_call_heavy))
      v_l<-as.data.frame(table(bcr_HL_V_sel$v_call_light))
      # topN for labeling

      # topN < group_by(vj) %>% summarise(n=n()) %>%  arrange(desc(n)) %>%  top_n()
      topN <- rownames_to_column(v_h) %>% arrange(desc(Freq)) %>%top_n(10,v_h$Freq)
      topN_nb<- as.numeric(topN$rowname)
      topN_nb

      # barplot annotation
      v_r_ha = rowAnnotation(
        H_Vc = anno_barplot(v_h$Freq),foo = anno_mark(at = topN_nb, labels = topN$Var1), annotation_name_rot = 00)

      v_r_ha1 = rowAnnotation(H_Vc = anno_barplot(v_h$Freq) )

      v_c_ha = HeatmapAnnotation(L_Vc = anno_barplot(v_l$Freq),which = "col")
      lgd = list(title = "Freq", col_fun = col_fun ,labels_gp = gpar(fontsize=14))


      ph1<- Heatmap(data_subset_norm_v, column_title = paste0(i," HL v gene usage(scaled)"),
                    name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                    top_annotation = v_c_ha, right_annotation =v_r_ha,
                    row_names_side = "left", heatmap_legend_param = lgd,
                    width = ncol(data_subset_norm_v)*unit(3, "mm"),
                    height = nrow(data_subset_norm_v)*unit(3, "mm")) # turn off row clustering
      draw(ph1)


      col_fun1a = colorRamp2(c(0, max(v_h$Freq)/500, max(v_h$Freq)/4), c(  "white", "orange","red"))
      ph1a<- Heatmap(freq_v, column_title = paste0(i," HL v gene usage"),
                     name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun1a,
                     row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                     top_annotation = v_c_ha, right_annotation =v_r_ha1,
                     row_names_side = "left",heatmap_legend_param =lgd,
                     width = ncol(freq_v)*unit(3, "mm"),
                     height = nrow(freq_v)*unit(3, "mm")) # turn off row clustering

      draw(ph1a)
      #----------------------------------------------------------------------

      pdf(file=paste0(outputPath,i,"_lv1_HL_mainV_usageHeatmap.pdf"))

      draw(ph1)
      draw(ph1a)


      #####################################
      # 2 HL V  heatmap
      # top N based on paired HL v count but only label the heavy V because, same heavy v has different light v pairs
      #####################################

      # Lable top HL paired vj

      # topN based on heavy and light paired  V genes
      row_vH <- rownames_to_column(v_h)
      head(row_vH)

      # combine heavy vj row order to the bcr data
      order_bcr_HL_v<- inner_join(x=row_vH, y=bcr_HL_V_sel, by=c("Var1"="v_call_heavy"))
      head(order_bcr_HL_v)

      # extract top n paired heavy and light v
      v_HL<- as.data.frame(table(bcr_HL_V_sel$HL_v))
      topN_hl_v <- rownames_to_column(v_HL) %>% arrange(desc(Freq)) %>%top_n(20,v_HL$Freq)

      # extract the row number and corresponding heavy vj genes
      topN_hl_v_order<- order_bcr_HL_v %>% filter(HL_v%in%topN_hl_v$Var1)  %>%
        select(rowname,Var1)  %>% unique()
      topN_hl_v_order

      # extract the row number
      topN_nb_hl_v_order<- as.numeric(topN_hl_v_order$rowname)
      topN_nb_hl_v_order


      # project to heavy vj data for corresponding row number

      HL_v_r_ha = rowAnnotation(
        H_Vc = anno_barplot(v_h$Freq),foo = anno_mark(at = topN_nb_hl_v_order, labels = topN_hl_v_order$Var1), annotation_name_rot = 0)
      lgd = list(title = "Freq", col_fun = col_fun ,labels_gp = gpar(fontsize=14))


      # plot
      ph2<- Heatmap(data_subset_norm_v, column_title = paste0(i," HL V gene usage \n topN paired HL V genes"),
                    name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                    top_annotation = v_c_ha, right_annotation =HL_v_r_ha,
                    row_names_side = "left",heatmap_legend_param = lgd,
                    width = ncol(data_subset_norm_v)*unit(3, "mm"),
                    height = nrow(data_subset_norm_v)*unit(3, "mm")) # turn off row clustering

      draw(ph2)

      col_fun2a = colorRamp2(c(0, max(v_h$Freq)/500, max(v_h$Freq)/4), c(  "white", "orange","red"))
      ph2a<- Heatmap(freq_v, column_title = paste0(i," HL V gene usage"),
                     name = i, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun2a,
                     row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16),
                     top_annotation = v_c_ha, right_annotation =v_r_ha1,
                     row_names_side = "left",heatmap_legend_param =lgd,
                     width = ncol(freq_v)*unit(3, "mm"),
                     height = nrow(freq_v)*unit(3, "mm")) # turn off row clustering

      draw(ph2a)

      #-------------------------------------------------------------------------
      dev.off()


    }

  }
  summary<- do.call(rbind.data.frame,final)
  return(summary)
}










#==============================================================================#
#  function: biased_VHL
# Heatmap plot for the H and L v gene usage
# df_fileList: file list ; df_bcrInputPath: lv1_corrected BCR; rdsInputPath: QCed sc data
# type: "heavy"/"light"/"HL"
#==============================================================================#


biased_VHL <- function(df_fileList,df_bcrInputPath,rdsInputPath,outputPlot,outputSummary,type){
  final<- list()
  for(i in df_fileList){

    # Read GEX data
    gex_db <- readRDS(paste0(rdsInputPath,i,".rds"))

    if(type=="heavy"){
      # Read BCR data
      bcr_db <- read_tsv(paste0(df_bcrInputPath,i,"_lv1CorrectHeavyBcr.tsv"))
      ## Standardize cell IDs

      if(i=="p36T1"){
        sample_id<-"p36_T1"
        sample_NC<-"p36D4"
      }else if(i=="p36T2"){
        sample_id<-"p36_T2"
        sample_NC<-"p36D45"
      }else if(i=="p36T4"){
        sample_id<-"p36_T4"
        sample_NC<-"p36D169"
      }else if(i=="HD2B960"){
        sample_id<-"HD2B960"
        sample_NC<-"HD1_stD2"
      }else if(i=="HD2B3600"){
        sample_id<-"HD2B3600"
        sample_NC<-"HD2"
      }else if(i=="HD2B1500"){
        sample_id<-"HD2B1500"
        sample_NC<-"HD1"
      }else {
        print("now only support p36&HD time points")
      }

      bcr_db$sample<- sample_id
      bcr_db$sampleNC<- sample_NC

      # Make cell IDs in BCR match those in Seurat Object
      # bcr_db$cell_id = sapply(strsplit(bcr_db$sequence_id,"_"), "[",1)

      # NOTE: cell_id becomes the cluster_cb after lv1 correction
      bcr_db$cell_id_unique = paste0(bcr_db$sample, "_", bcr_db$cell_id)
      bcr_db$cluster_cb_unique = paste0(bcr_db$sample, "_", bcr_db$cluster_cb)

      bcr_db$cell_id_unique[1:4]

      # subset heavy bcr for match to have unique cell_id(because paired light chain will have the same cell_id)
      # Note: lv1 correction already select the heavy chain
      # bcr_db <- bcr_db[grep("IGH",bcr_db$v_call ,invert= F),]

      #####################
      # matching without lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index = match(Cells(gex_db), bcr_db$cell_id_unique)

      #####################
      # matching after lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index_lv1 = match(Cells(gex_db), bcr_db$cluster_cb_unique)


      # In this data, not all cells are B cells.
      # What proportion of cells don’t have BCRs?
      scWithOutMatchedBcr<- NULL
      scWithOutMatchedBcr$scMissBcr <- mean(!is.na(match.index))
      scWithOutMatchedBcr$scMissBcr_lv1 <- mean(!is.na(match.index_lv1))
      scWithOutMatchedBcr$id <-i

      final[[i]]<- scWithOutMatchedBcr


      #------------------------------------------------------------------------
      ### Find the BCR cells in the GEX data

      # Match indices to find the position of the BCR cells in the GEX data
      # Different from finding the position of the GEX cells in the BCR data!

      match.index = match(bcr_db$cluster_cb_unique, Cells(gex_db))
      # What proportion of BCRs don’t have GEX information?
      mean(is.na(match.index))
      mean(!is.na(match.index))
      ### Transfer GEX annotations into the BCR data
      # Add annotations to BCR data
      cell.annotation = as.character(gex_db@meta.data$SingleR.labels)
      bcr_db$gex_annotation= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.annotation[x])}))
      bcr_db$gex_annotation[1:5]

      # Add UMAP coordinates to BCR data
      #umap1 = gex_db@reductions$umap@cell.embeddings[,1]
      #umap2 = gex_db@reductions$umap@cell.embeddings[,2]
      #bcr_db$gex_umap1= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap1[x])}))
      #bcr_db$gex_umap2= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap2[x])}))
      bcr_db[1:5,] %>%
        select(cell_id,gex_annotation)



      ### Remove cells without GEX data
      # Remove cells that didn’t match
      bcr_db = filter(bcr_db, !is.na(gex_annotation))
      col_anno = c("Exhausted B cells"="dodgerblue2", "Naive B cells"="firebrick2", "Non-switched memory B cells"="seagreen","Plasmablasts"="darkgoldenrod2", "Switched memory B cells"="black")
      bcr_hm <- bcr_db  %>% select(sequence_id,v_call,j_call,gex_annotation)
      # only use the 1st v gene hit
      bcr_hm$v_call <- sapply(str_split(bcr_hm$v_call, pattern = ","),"[",1)
      bcr_hm$j_call <- sapply(str_split(bcr_hm$j_call, pattern = ","),"[",1)
      head(bcr_hm)

      bcr_hm1<- bcr_hm
      # further remove the * after to reduce the V diversity
      bcr_hm1$v_call <- sapply(str_split(bcr_hm$v_call, pattern = "\\*"),"[",1)
      bcr_hm1$j_call <- sapply(str_split(bcr_hm$j_call, pattern = "\\*"),"[",1)
      bcr_hm1 <- bcr_hm1 %>% mutate(vj=paste0(v_call,"_", j_call))
      head(bcr_hm1)

      #----------
      # https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
      # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
      freq <- table(bcr_hm1$v_call,bcr_hm1$j_call)
      freq

      # scale the row or column before the heatmap ploting
      cal_z_score <- function(x){
        (x - mean(x)) / sd(x)
      }

      data_subset_norm <- t(apply(freq, 1, cal_z_score))
      data_subset_norm

      #################
      # prepare row and column annotation bar plot
      #################
      # row and column data for barplot
      v_h<-as.data.frame(table(bcr_hm1$v_call))
      j_h<-as.data.frame(table(bcr_hm1$j_call))
      # topN for labeling

      # topN < group_by(vj) %>% summarise(n=n()) %>%  arrange(desc(n)) %>%  top_n()
      topN <- rownames_to_column(v_h) %>% arrange(desc(Freq)) %>%top_n(10,v_h$Freq)
      topN_nb<- as.numeric(topN$rowname)


      # barplot annotation
      v_r_ha = rowAnnotation(
        V_count = anno_barplot(v_h$Freq),foo = anno_mark(at = topN_nb, labels = topN$Var1), annotation_name_rot = 270)

      v_c_ha = HeatmapAnnotation(J_count = anno_barplot(j_h$Freq),which = "col")



      #-------------------------------------------------------------------------
      # heatmap
      pdf(file=paste0(outputPlot,sample_NC,"_lv1_heavy_vj_usageHeatmap.pdf"))

      #pheatmap(data_subset_norm)
      # use complex heat map to remove the cluster of row and column
      # Turn off the clustering fucntion
      col_fun = colorRamp2(c(-3, 0, 5), c("green", "white", "red"))
      lgd = list(title = "Freq", col_fun = col_fun)

      ph1<- Heatmap(data_subset_norm, column_title = paste0(sample_NC," heavy VJ gene usage"),
                    name = sample_NC, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 8),
                    top_annotation = v_c_ha, right_annotation =v_r_ha,
                    row_names_side = "left",heatmap_legend_param =lgd) # turn off row clustering

      draw(ph1)

      #-------------------------------------------------------------------------

      dev.off()


    }else if(type=="light"){
      # Read BCR data
      bcr_db <- read_tsv(paste0(df_bcrInputPath,i,"_lv1CorrectLightBcr.tsv"))
      ## Standardize cell IDs

      if(i=="p36T1"){
        sample_id<-"p36_T1"
        sample_NC<-"p36D4"
      }else if(i=="p36T2"){
        sample_id<-"p36_T2"
        sample_NC<-"p36D45"
      }else if(i=="p36T4"){
        sample_id<-"p36_T4"
        sample_NC<-"p36D169"
      }else if(i=="HD2B960"){
        sample_id<-"HD2B960"
        sample_NC<-"HD1_stD2"
      }else if(i=="HD2B3600"){
        sample_id<-"HD2B3600"
        sample_NC<-"HD2"
      }else if(i=="HD2B1500"){
        sample_id<-"HD2B1500"
        sample_NC<-"HD1"
      }else {
        print("now only support p36&HD time points")
      }

      bcr_db$sample<- sample_id
      bcr_db$sampleNC<- sample_NC
      # Make cell IDs in BCR match those in Seurat Object
      # bcr_db$cell_id = sapply(strsplit(bcr_db$sequence_id,"_"), "[",1)

      # NOTE: cell_id becomes the cluster_cb after lv1 correction
      bcr_db$cell_id_unique = paste0(bcr_db$sample, "_", bcr_db$cell_id)
      bcr_db$cluster_cb_unique = paste0(bcr_db$sample, "_", bcr_db$cluster_cb)

      bcr_db$cell_id_unique[1:4]

      # subset heavy bcr for match to have unique cell_id(because paired light chain will have the same cell_id)
      # Note: lv1 correction already select the heavy chain
      # bcr_db <- bcr_db[grep("IGH",bcr_db$v_call ,invert= F),]

      #####################
      # matching without lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index = match(Cells(gex_db), bcr_db$cell_id_unique)

      #####################
      # matching after lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index_lv1 = match(Cells(gex_db), bcr_db$cluster_cb_unique)


      # In this data, not all cells are B cells.
      # What proportion of cells don’t have BCRs?
      scWithOutMatchedBcr<- NULL
      scWithOutMatchedBcr$scMissBcr <- mean(!is.na(match.index))
      scWithOutMatchedBcr$scMissBcr_lv1 <- mean(!is.na(match.index_lv1))
      scWithOutMatchedBcr$id <-i

      final[[i]]<- scWithOutMatchedBcr


      #------------------------------------------------------------------------
      ### Find the BCR cells in the GEX data

      # Match indices to find the position of the BCR cells in the GEX data
      # Different from finding the position of the GEX cells in the BCR data!

      match.index = match(bcr_db$cluster_cb_unique, Cells(gex_db))
      # What proportion of BCRs don’t have GEX information?
      mean(is.na(match.index))
      mean(!is.na(match.index))
      ### Transfer GEX annotations into the BCR data
      # Add annotations to BCR data
      cell.annotation = as.character(gex_db@meta.data$SingleR.labels)
      bcr_db$gex_annotation= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.annotation[x])}))
      bcr_db$gex_annotation[1:5]

      # Add UMAP coordinates to BCR data
      #umap1 = gex_db@reductions$umap@cell.embeddings[,1]
      #umap2 = gex_db@reductions$umap@cell.embeddings[,2]
      #bcr_db$gex_umap1= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap1[x])}))
      #bcr_db$gex_umap2= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap2[x])}))
      bcr_db[1:5,] %>%
        select(cell_id,gex_annotation)



      ### Remove cells without GEX data
      # Remove cells that didn’t match
      bcr_db = filter(bcr_db, !is.na(gex_annotation))


      ################################################################
      # UMAP cluster labeling need to be customized for individual sample
      # !!!!!!!!!!!!!!!!!!!!!!!!!
      ################################################################
      ### Ensure information transferred from Seurat object
      col_anno = c("Exhausted B cells"="dodgerblue2", "Naive B cells"="firebrick2", "Non-switched memory B cells"="seagreen","Plasmablasts"="darkgoldenrod2", "Switched memory B cells"="black")
      # Plot UMAP from bcr_db
      #pdf(paste0(outputPlot,i, "umap_from_bcr_db.pdf"))

      #bcr_umap <- ggplot(bcr_db) +geom_point(aes(x = gex_umap1, y = gex_umap2, color = gex_annotation)) +scale_colour_manual(values=col_anno) +theme_bw()

      #print(bcr_umap)
      #dev.off()
      #----------------------------------------------------------------------

      # bcr_dbH <-bcr_db[grep("IGH",bcr_db$v_call,invert = F),] # It's already Heavy chain after lv1 correction


      bcr_hm <- bcr_db  %>% select(sequence_id,v_call,j_call,gex_annotation)
      # only use the 1st v gene hit
      bcr_hm$v_call <- sapply(str_split(bcr_hm$v_call, pattern = ","),"[",1)
      bcr_hm$j_call <- sapply(str_split(bcr_hm$j_call, pattern = ","),"[",1)
      head(bcr_hm)

      bcr_hm1<- bcr_hm
      # further remove the * after to reduce the V diversity
      bcr_hm1$v_call <- sapply(str_split(bcr_hm$v_call, pattern = "\\*"),"[",1)
      bcr_hm1$j_call <- sapply(str_split(bcr_hm$j_call, pattern = "\\*"),"[",1)
      bcr_hm1 <- bcr_hm1 %>% mutate(vj=paste0(v_call,"_", j_call))
      head(bcr_hm1)


      #----------
      # https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
      # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
      freq <- table(bcr_hm1$v_call,bcr_hm1$j_call)
      freq

      # scale the row or column before the heatmap ploting
      cal_z_score <- function(x){
        (x - mean(x)) / sd(x)
      }

      data_subset_norm <- t(apply(freq, 1, cal_z_score))
      data_subset_norm

      #################
      # prepare row and column annotation bar plot
      #################
      # row and column data for barplot
      v_l<-as.data.frame(table(bcr_hm1$v_call))
      j_l<-as.data.frame(table(bcr_hm1$j_call))
      head(v_l)
      # topN for labeling

      # topN < group_by(vj) %>% summarise(n=n()) %>%  arrange(desc(n)) %>%  top_n()
      topN <- rownames_to_column(v_l) %>% arrange(desc(Freq)) %>%top_n(10,v_l$Freq)
      topN_nb<- as.numeric(topN$rowname)
      head(topN_nb)

      # barplot annotation
      v_r_ha = rowAnnotation(
        V_count = anno_barplot(v_l$Freq),foo = anno_mark(at = topN_nb, labels = topN$Var1), annotation_name_rot = 270)

      v_c_ha = HeatmapAnnotation(J_count = anno_barplot(j_l$Freq),which = "col")

      # heatmap
      pdf(file=paste0(outputPlot,sample_NC,"_lv1_light_vj_usageHeatmap.pdf"))
      #-------------------------------------------------------------------------

      #pheatmap(data_subset_norm)
      # use complex heat map to remove the cluster of row and column
      # Turn off the clustering fucntion
      col_fun = colorRamp2(c(-3, 0, 5), c("green", "white", "red"))
      lgd = list(title = "Freq", col_fun = col_fun)

      ph1<- Heatmap(data_subset_norm, column_title = paste0(sample_NC," light VJ gene usage"),
                    name = sample_NC, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 8),
                    top_annotation = v_c_ha, right_annotation =v_r_ha,
                    row_names_side = "left",
                    heatmap_legend_param = lgd) # turn off row clustering

      draw(ph1)

      #-------------------------------------------------------------------------

      dev.off()

    }else if(type=="HL"){

      #########################################################################
      # type="HL"
      #########################################################################

      # Read BCR data
      bcr_dbH <- read_tsv(paste0(df_bcrInputPath,i,"_lv1CorrectHeavyBcr.tsv"))
      bcr_dbL <- read_tsv(paste0(df_bcrInputPath,i,"_lv1CorrectLightBcr.tsv"))
      ## Standardize cell IDs

      if(i=="p36T1"){
        sample_id<-"p36_T1"
        sample_NC<-"p36D4"
      }else if(i=="p36T2"){
        sample_id<-"p36_T2"
        sample_NC<-"p36D45"
      }else if(i=="p36T4"){
        sample_id<-"p36_T4"
        sample_NC<-"p36D169"
      }else if(i=="HD2B960"){
        sample_id<-"HD2B960"
        sample_NC<-"HD1_stD2"
      }else if(i=="HD2B3600"){
        sample_id<-"HD2B3600"
        sample_NC<-"HD2"
      }else if(i=="HD2B1500"){
        sample_id<-"HD2B1500"
        sample_NC<-"HD1"
      }else {
        print("now only support p36&HD time points")
      }

      bcr_dbH$sampleNC<- sample_NC
      bcr_dbL$sampleNC<- sample_NC

      bcr_dbH$sample<- sample_id
      bcr_dbL$sample<- sample_id
      dim(bcr_dbH)
      dim(bcr_dbL)

      # Make cell IDs in BCR match those in Seurat Object
      # bcr_db$cell_id = sapply(strsplit(bcr_db$sequence_id,"_"), "[",1)

      # NOTE: cell_id becomes the cluster_cb after lv1 correction
      bcr_dbH$cell_id_unique = paste0(bcr_dbH$sample, "_", bcr_dbH$cell_id)
      bcr_dbH$cluster_cb_unique = paste0(bcr_dbH$sample, "_", bcr_dbH$cluster_cb)

      bcr_dbL$cell_id_unique = paste0(bcr_dbL$sample, "_", bcr_dbL$cell_id)
      bcr_dbL$cluster_cb_unique = paste0(bcr_dbL$sample, "_", bcr_dbL$cluster_cb)
      bcr_dbL$cell_id_unique[1:4]

      # subset heavy bcr for match to have unique cell_id(because paired light chain will have the same cell_id)
      # Note: lv1 correction already select the heavy chain
      # bcr_db <- bcr_db[grep("IGH",bcr_db$v_call ,invert= F),]

      #####################
      # matching without lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.indexH = match(Cells(gex_db), bcr_dbH$cell_id_unique)

      #####################
      # matching after lv1 correction
      #####################
      # match index to find the position of the GEX cells in the BCR data
      match.index_lv1H = match(Cells(gex_db), bcr_dbH$cluster_cb_unique)


      # In this data, not all cells are B cells.
      # What proportion of cells don’t have BCRs?
      scWithOutMatchedBcrH<- NULL
      scWithOutMatchedBcrH$scMatchBcr <- mean(!is.na(match.indexH))
      scWithOutMatchedBcrH$scMatchBcr_lv1 <- mean(!is.na(match.index_lv1H))
      scWithOutMatchedBcrH$id <-i
      scWithOutMatchedBcrH

      final[[i]]<- scWithOutMatchedBcrH


      #------------------------------------------------------------------------
      ### Find the BCR cells in the GEX data

      # Match indices to find the position of the BCR cells in the GEX data
      # Different from finding the position of the GEX cells in the BCR data!

      match.indexH = match(bcr_dbH$cluster_cb_unique, Cells(gex_db))
      # What proportion of BCRs don’t have GEX information?
      mean(is.na(match.indexH))
      mean(!is.na(match.indexH))
      ### Transfer GEX annotations into the BCR data
      # Add annotations to BCR data
      cell.annotation = as.character(gex_db@meta.data$SingleR.labels)
      bcr_dbH$gex_annotation= unlist(lapply(match.indexH,function(x){ifelse(is.na(x),NA, cell.annotation[x])}))
      bcr_dbH$gex_annotation[1:5]

      # Add UMAP coordinates to BCR data
      #umap1 = gex_db@reductions$umap@cell.embeddings[,1]
      #umap2 = gex_db@reductions$umap@cell.embeddings[,2]
      #bcr_db$gex_umap1= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap1[x])}))
      #bcr_db$gex_umap2= unlist(lapply(match.index, function(x){ifelse(is.na(x),NA, umap2[x])}))
      bcr_dbH[1:5,] %>%
        select(cell_id_unique,gex_annotation)



      ### Remove cells without GEX data
      # Remove cells that didn’t match
      bcr_dbH = filter(bcr_dbH, !is.na(gex_annotation))


      ################################################################
      # UMAP cluster labeling need to be customized for individual sample
      # !!!!!!!!!!!!!!!!!!!!!!!!!
      ################################################################
      ### Ensure information transferred from Seurat object
      col_anno = c("Exhausted B cells"="dodgerblue2", "Naive B cells"="firebrick2", "Non-switched memory B cells"="seagreen","Plasmablasts"="darkgoldenrod2", "Switched memory B cells"="black")
      # Plot UMAP from bcr_db
      #pdf(paste0(outputPlot,i, "umap_from_bcr_db.pdf"))

      #bcr_umap <- ggplot(bcr_db) +geom_point(aes(x = gex_umap1, y = gex_umap2, color = gex_annotation)) +scale_colour_manual(values=col_anno) +theme_bw()

      #print(bcr_umap)
      #dev.off()
      #----------------------------------------------------------------------
      # Format the V J gene usage for Heavy
      # only use the 1st v gene hit
      bcr_dbH$v_call <- sapply(str_split(bcr_dbH$v_call, pattern = ","),"[",1)
      bcr_dbH$j_call <- sapply(str_split(bcr_dbH$j_call, pattern = ","),"[",1)
      head(bcr_dbH[1:4])

      bcr_dbH1<- bcr_dbH
      # further remove the * after to reduce the V diversity
      bcr_dbH1$v_call <- sapply(str_split(bcr_dbH$v_call, pattern = "\\*"),"[",1)
      bcr_dbH1$j_call <- sapply(str_split(bcr_dbH$j_call, pattern = "\\*"),"[",1)

      head(bcr_dbH1)
      # Format the V J gene usage for light
      # only use the 1st v gene hit
      bcr_dbL$v_call <- sapply(str_split(bcr_dbL$v_call, pattern = ","),"[",1)
      bcr_dbL$j_call <- sapply(str_split(bcr_dbL$j_call, pattern = ","),"[",1)
      head(bcr_dbL)

      bcr_dbL1<- bcr_dbL
      # further remove the * after to reduce the V diversity
      bcr_dbL1$v_call <- sapply(str_split(bcr_dbL$v_call, pattern = "\\*"),"[",1)
      bcr_dbL1$j_call <- sapply(str_split(bcr_dbL$j_call, pattern = "\\*"),"[",1)

      head(bcr_dbL1)

      # combine H L BCR
      bcr_dbH_sel <- bcr_dbH1 %>%
        select(sequence_id,v_call,j_call,gex_annotation,cluster_cb_unique) %>%
        mutate(vj=paste0(v_call,"_",j_call))

      bcr_dbL_sel <- bcr_dbL1 %>% select(sequence_id,v_call,j_call,cluster_cb_unique)%>%
        mutate(vj=paste0(v_call,"_",j_call))

      bcr_HL <- inner_join(x=bcr_dbH_sel,y=bcr_dbL_sel,
                           by=c("cluster_cb_unique"),suffix = c("_heavy","_light"))
      # save the integrated BCR data and ready for plotting based on this data
      write_tsv(x=bcr_HL, file = paste0(outputSummary,sample_NC,"_iBcr.tsv"))

      #########################################################################
      # part-1. HL VJ  heatmap
      #########################################################################

      bcr_HL_sel<- bcr_HL %>% select(vj_heavy, vj_light,gex_annotation) %>% mutate(HL_vj=paste0(vj_heavy,"_",vj_light))

      dim(bcr_HL_sel)


      #####################################
      # 1. HL VJ  heatmap
      # top N based on pure heavy vj count
      #####################################
      # https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
      # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
      freq <- table(bcr_HL_sel$vj_heavy,bcr_HL_sel$vj_light)
      head(freq)

      # scale the row or column before the heatmap ploting
      cal_z_score <- function(x){
        (x - mean(x)) / sd(x)
      }

      data_subset_norm <- t(apply(freq, 1, cal_z_score))
      data_subset_norm

      #################
      # prepare row and column annotation bar plot
      #################
      # row and column data for barplot
      vj_H<-as.data.frame(table(bcr_HL_sel$vj_heavy))
      vj_L<-as.data.frame(table(bcr_HL_sel$vj_light))
      vj_HL<- as.data.frame(table(bcr_HL_sel$HL_vj))
      # topN for labeling

      # topN based on heavy chain vj
      topN <- rownames_to_column(vj_H) %>% arrange(desc(Freq)) %>%top_n(10,vj_H$Freq)
      topN_nb<- as.numeric(topN$rowname)


      # barplot annotation
      v_r_ha = rowAnnotation(
        H_vj_count = anno_barplot(vj_H$Freq),foo = anno_mark(at = topN_nb, labels = topN$Var1), annotation_name_rot = 270)
      v_c_ha = HeatmapAnnotation(L_vj_count = anno_barplot(vj_L$Freq),which = "col")
      #-------------------------------------------------------------------------
      # heatmap
      pdf(file=paste0(outputPlot,sample_NC,"_lv1_HL_vj_usageHeatmap.pdf"))

      #pheatmap(data_subset_norm)
      # use complex heat map to remove the cluster of row and column
      # Turn off the clustering fucntion
      col_fun = colorRamp2(c(-3, 0, 5), c("green", "white", "red"))
      lgd = list(title = "Freq", col_fun = col_fun)

      ph1<- Heatmap(data_subset_norm, column_title = paste0(sample_NC," heavy VJ vs light VJ gene usage"),
                    name = sample_NC, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 8),
                    top_annotation = v_c_ha, right_annotation =v_r_ha,
                    row_names_side = "left", heatmap_legend_param = lgd,
                    show_column_names = FALSE, show_row_names = FALSE) # turn off row clustering

      draw(ph1)


      #####################################
      # 2. HL VJ  heatmap
      # top N based on paired HL vj count but only label the heavy VJ because, same heavy vj has different light vj
      #####################################

      # topN based on heavy and light paired  vj
      row_vjH <- rownames_to_column(vj_H)

      # combine heavy vj row order to the bcr data
      order_bcr_HL_sel<- inner_join(x=row_vjH, y=bcr_HL_sel, by=c("Var1"="vj_heavy"))
      head(order_bcr_HL_sel)

      # extract top n paired heavy and light vj
      topN_hl <- rownames_to_column(vj_HL) %>% arrange(desc(Freq)) %>%top_n(20,vj_HL$Freq)

      # extract the row number and corresponding heavy vj genes
      topN_hl_order<- order_bcr_HL_sel %>% filter(HL_vj%in%topN_hl$Var1)  %>%
        select(rowname,Var1)  %>% unique()
      topN_hl_order

      # extract the row number
      topN_nb_hl_order<- as.numeric(topN_hl_order$rowname)
      topN_nb_hl_order


      # project to heavy vj data for corresponding row number
      HL_vj_r_ha = rowAnnotation(
        H_vj_count = anno_barplot(vj_H$Freq),foo = anno_mark(at = topN_nb_hl_order, labels = topN_hl_order$Var1))
      lgd = list(title = "Freq", col_fun = col_fun)

      # plot
      ph2<- Heatmap(data_subset_norm, column_title = paste0(sample_NC," heavy VJ vs light VJ gene usage \n topN paired VJ genes"),
                    name = sample_NC, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 8),
                    top_annotation = v_c_ha, right_annotation =HL_vj_r_ha,
                    row_names_side = "left", heatmap_legend_param = lgd,
                    show_column_names = FALSE, show_row_names = FALSE) # turn off row clustering

      draw(ph2)
      dev.off()


      #########################################################################
      # Part-2: heavy and light v gene heatmap
      #########################################################################
      bcr_HL_V_sel<- bcr_HL %>% select(v_call_heavy, v_call_light,gex_annotation)  %>%
        mutate(HL_v=paste0(v_call_heavy,"_",v_call_light))

      head(bcr_HL_V_sel)


      # https://www.quora.com/How-do-I-get-a-frequency-count-based-on-two-columns-variables-in-an-R-dataframe
      # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
      freq_v <- table(bcr_HL_V_sel$v_call_heavy, bcr_HL_V_sel$v_call_light)
      freq_v

      # scale the row or column before the heatmap ploting
      cal_z_score <- function(x){
        (x - mean(x)) / sd(x)
      }

      data_subset_norm_v <- t(apply(freq_v, 1, cal_z_score))
      data_subset_norm_v


      #####################################
      # 1 HL V  heatmap
      # top N based on paired HL vj count but only label the heavy VJ because, same heavy vj has different light vj
      #####################################

      # use complex heat map to remove the cluster of row and column
      # Turn off the clustering fucntion
      col_fun = colorRamp2(c(-3, 0, 5), c("green", "white", "red"))
      #################
      # prepare row and column annotation bar plot
      #################
      # row and column data for barplot
      v_h<-as.data.frame(table(bcr_HL_V_sel$v_call_heavy))
      v_l<-as.data.frame(table(bcr_HL_V_sel$v_call_light))
      # topN for labeling

      # topN < group_by(vj) %>% summarise(n=n()) %>%  arrange(desc(n)) %>%  top_n()
      topN <- rownames_to_column(v_h) %>% arrange(desc(Freq)) %>%top_n(10,v_h$Freq)
      topN_nb<- as.numeric(topN$rowname)
      topN_nb

      # barplot annotation
      v_r_ha = rowAnnotation(
        H_V_count = anno_barplot(v_h$Freq),foo = anno_mark(at = topN_nb, labels = topN$Var1), annotation_name_rot = 270)

      v_c_ha = HeatmapAnnotation(L_V_count = anno_barplot(v_l$Freq),which = "col")
      lgd = list(title = "Freq", col_fun = col_fun)


      ph1<- Heatmap(data_subset_norm_v, column_title = paste0(sample_NC," HL v gene usage"),
                    name = sample_NC, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 8),
                    top_annotation = v_c_ha, right_annotation =v_r_ha,
                    row_names_side = "left", heatmap_legend_param = lgd) # turn off row clustering
      draw(ph1)

      # save the plot for manuscript
      png(filename=paste0(outputPlot,sample_NC,"_lv1_HL_v_usageHeatmap.png"),width=20,height=20,units="cm",res=300)
      draw(ph1)
      dev.off()
      #----------------------------------------------------------------------

      pdf(file=paste0(outputPlot,sample_NC,"_lv1_HL_v_usageHeatmap.pdf"))

      draw(ph1)


      #####################################
      # 2 HL V  heatmap
      # top N based on paired HL v count but only label the heavy V because, same heavy v has different light v pairs
      #####################################

      # Lable top HL paired vj

      # topN based on heavy and light paired  V genes
      row_vH <- rownames_to_column(v_h)
      head(row_vH)

      # combine heavy vj row order to the bcr data
      order_bcr_HL_v<- inner_join(x=row_vH, y=bcr_HL_V_sel, by=c("Var1"="v_call_heavy"))
      head(order_bcr_HL_v)

      # extract top n paired heavy and light v
      v_HL<- as.data.frame(table(bcr_HL_V_sel$HL_v))
      topN_hl_v <- rownames_to_column(v_HL) %>% arrange(desc(Freq)) %>%top_n(20,v_HL$Freq)

      # extract the row number and corresponding heavy vj genes
      topN_hl_v_order<- order_bcr_HL_v %>% filter(HL_v%in%topN_hl_v$Var1)  %>%
        select(rowname,Var1)  %>% unique()
      topN_hl_v_order

      # extract the row number
      topN_nb_hl_v_order<- as.numeric(topN_hl_v_order$rowname)
      topN_nb_hl_v_order


      # project to heavy vj data for corresponding row number

      HL_v_r_ha = rowAnnotation(
        v = anno_barplot(v_h$Freq),foo = anno_mark(at = topN_nb_hl_v_order, labels = topN_hl_v_order$Var1), annotation_name_rot = 270)
      lgd = list(title = "Freq", col_fun = col_fun)


      # plot
      ph2<- Heatmap(data_subset_norm_v, column_title = paste0(sample_NC," HL V gene usage \n topN paired HL V genes"),
                    name = sample_NC, cluster_rows = FALSE,cluster_columns = FALSE,col=col_fun,
                    row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 8),
                    top_annotation = v_c_ha, right_annotation =HL_v_r_ha,
                    row_names_side = "left",heatmap_legend_param = lgd) # turn off row clustering

      draw(ph2)

      #-------------------------------------------------------------------------
      dev.off()


    }

  }
  summary<- do.call(rbind.data.frame,final)
  return(summary)
}















#==============================================================================#
# gex by sample (fig-4 g)
#  
#==============================================================================#



combine_sc_multiTime_isoPositive<- function(df_fileList,rdsInput_Path,outputPath_summary,outputPath_plot){
  final<- list()
  final_wide<-list()
  final_wide_ct<- list()
  final_gex <- list()
  for(i in  df_fileList ){
    
    ## Standardize cell IDs
    if(i=="p36T1"){
      sample_id<-"p36D4"
      dpi<- 1
      
    }else if(i=="p36T2"){
      sample_id<-"p36D45"
      dpi<- 42
    }else if(i=="p36T3"){
      sample_id<-"p36D92"
      dpi<-89
    }else if(i=="p36T4"){
      sample_id<-"p36D169"
      dpi<-166
    }else if(i=="p45T2"){
      sample_id<-"p45D24"
      dpi<-16
    }else if(i=="p45T5"){
      sample_id<-"p45D164"
      dpi<-16
    }else if(i=="HD2B1500"){
      sample_id<-"HD1"
      dpi<-0
    }else if(i=="HD2B3600"){
      sample_id<-"HD2"
      dpi<-0
    }else if(i=="HD2B960"){
      sample_id<-"HD1_stD2"
      dpi<-0
    }else if( i=="p16T1B"){
      sample_id<-"p16D9"
      dpi<- 2
    }else if(i=="p18T1"){
      sample_id<-"p18D13"
      dpi<- 1
    }else {
      print("now only support p36 time points")
    }
    
    
    # Read GEX data
    gex_db <- readRDS(paste0(outputSummary,i,"_singleR_AnnoSC_HX.rds"))
    
    #----------------------------------------------------------------------------
    # testing code
    # Add heavy and light chain isotype as features to the dataset
    
    
    gex_db[["percent.IgA1.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGHA1")
    gex_db[["percent.IgA2.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGHA2")
    
    
    gex_db[["percent.IgM.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGHM$")
    
    
    gex_db[["percent.IgG1.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGHG1")
    gex_db[["percent.IgG2.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGHG2")
    gex_db[["percent.IgG3.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGHG3")
    gex_db[["percent.IgG4.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGHG4")
    
    
    gex_db[["percent.IgD.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGHD$")
    
    gex_db[["percent.IgE.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGHE$")
    
    
    gex_db[["percent.K.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGKC")
    
    gex_db[["percent.L.iso"]] <- PercentageFeatureSet(gex_db, pattern = "^IGLC")
    
    count_matrics <- as.data.frame(gex_db@assays$RNA@counts)
    # Extract antibody heavy  chain isotype infor
    #new_matrics<- count_matrics %>% rownames_to_column()
    # new_matrics[1:5,1:5]
    tdb<- t <- t(count_matrics)
    tdb[1:4,1:4]
    tdb_s<- tdb[,c( "IGHA1","IGHA2",   "IGHE" ,   "IGHG4",   "IGHG2"  ,   "IGHG1" ,  "IGHG3",   "IGHD",    "IGHM")]
    tdb_s_df <- tdb_s %>% as.data.frame()  %>% rownames_to_column()
    head(tdb_s_df )
    
    
    ##############################
    # add gex_annotation to the iso plot
    ##############################
    match.index = match(tdb_s_df$rowname, Cells(gex_db))
    
    # What proportion of BCRs don’t have GEX information?
    mean(is.na(match.index))
    mean(!is.na(match.index))
    
    ### Transfer GEX annotations into the BCR data
    
    
    # Add annotations to BCR data
    cell.annotation = as.character(gex_db@meta.data$SingleR.labels)
    tdb_s_df$gex_annotation= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.annotation[x])}))
    tdb_s_df$gex_annotation[1:5]
    
    ##############################
    # add isotype percentage to the iso plot
    ##############################
    
    # Add annotations to BCR data
    # igM
    cell.igM.pct = as.character(gex_db@meta.data$percent.IgM.iso)
    tdb_s_df$igM.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igM.pct[x])}))
    tdb_s_df$igM.pct[1:5]
    tdb_s_df[1:5,]
    # IgD
    cell.igD.pct = as.character(gex_db@meta.data$percent.IgD.iso)
    tdb_s_df$igD.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igD.pct[x])}))
    tdb_s_df[1:5,]
    
    # IgG
    cell.igG1.pct = as.character(gex_db@meta.data$percent.IgG1.iso)
    tdb_s_df$igG1.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igG1.pct[x])}))
    
    cell.igG2.pct = as.character(gex_db@meta.data$percent.IgG2.iso)
    tdb_s_df$igG2.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igG2.pct[x])}))
    
    cell.igG3.pct = as.character(gex_db@meta.data$percent.IgG3.iso)
    tdb_s_df$igG3.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igG3.pct[x])}))
    
    cell.igG4.pct = as.character(gex_db@meta.data$percent.IgG4.iso)
    tdb_s_df$igG4.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igG4.pct[x])}))
    tdb_s_df[1:5,]
    
    
    # IgA
    cell.igA1.pct = as.character(gex_db@meta.data$percent.IgA1.iso)
    tdb_s_df$igA1.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igA1.pct[x])}))
    
    cell.igA2.pct = as.character(gex_db@meta.data$percent.IgA2.iso)
    tdb_s_df$igA2.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igA2.pct[x])}))
    tdb_s_df[1:5,]
    
    
    # IgE
    cell.igE.pct = as.character(gex_db@meta.data$percent.IgE.iso)
    tdb_s_df$igE.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igE.pct[x])}))
    tdb_s_df[1:5,]
    
    # Igk
    #cell.igK.pct = as.character(gex_db@meta.data$percent.K.iso)
    #tdb_s_df$igK.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igK.pct[x])}))
    #tdb_s_df[1:5,]
    
    
    # Ig_lambda
    #cell.igLambda.pct = as.character(gex_db@meta.data$percent.L.iso)
    #tdb_s_df$igL.pct= unlist(lapply(match.index,function(x){ifelse(is.na(x),NA, cell.igLambda.pct[x])}))
    #tdb_s_df[1:5,]
    
    # only select b cells
    tdb_s_df_b <- tdb_s_df[grep(pattern = "B cells|Plasmablasts",tdb_s_df$gex_annotation ),]
    dim(tdb_s_df)  # for p36t1:  4407  cells
    dim(tdb_s_df_b)  # for p36t1:  3155  are B  cells
    tdb_s_df_b[1:4,1:4]
    head(tdb_s_df_b)
    table( tdb_s_df_b$gex_annotation)
    # select gex and rowname for gex plot by sample
    tdb_s_df_b_gexBySample<- tdb_s_df_b %>% select(rowname,gex_annotation)
    tdb_s_df_b_gexBySample$sample<- sample_id
    head(tdb_s_df_b_gexBySample)
    colnames(tdb_s_df_b_gexBySample)[1]<- "cell_id_uniq"
    
    #----------------------------------------------------------------------------
    # for count 
    wide_ct <- tdb_s_df_b[,c(1:11)]
    head(wide_ct)
    names(wide_ct)[names(wide_ct)=="rowname"]<- "cell_id_unique"
    
    tdb_long_count <- tdb_s_df_b[,c(1:11)] %>% pivot_longer(cols =starts_with("IGH") , names_to = "iso",values_to = "count") %>%
      filter(count>0) %>% mutate(cell_id_unique=rowname) %>% mutate(cell_id_unique_iso=paste0(cell_id_unique,"_",iso))%>%
      select(cell_id_unique_iso,cell_id_unique,gex_annotation,iso,count)
    head(tdb_long_count)
    dim(tdb_long_count)  # for p36t1:  1985 out of 3155 B cells have isotype information
    
    table(tdb_long_count$iso)
    
    # for percentage
    tdb_long_percent_1<- tdb_s_df_b[,c(1,11:20)]
    colnames(tdb_long_percent_1)<- c("cell_id_unique_pt","gex_annotation","IGHM","IGHD","IGHG1" ,"IGHG2","IGHG3",
                                     "IGHG4","IGHA1", "IGHA2" ,"IGHE" )
    wide_pt<- tdb_long_percent_1
    head(wide_pt)
    #wide_ptCt
    
    wide_pt$sample<- sample_id
    
    tdb_long_percent<- tdb_long_percent_1 %>% pivot_longer(cols =IGHM:IGHE , names_to = "iso_pt",values_to = "percent")%>%
      filter(percent>0) %>% mutate(gex_annotation_pt=gex_annotation)%>%
      mutate(cell_id_unique_iso=paste0(cell_id_unique_pt,"_",iso_pt))%>%
      select(cell_id_unique_iso,cell_id_unique_pt,gex_annotation_pt,iso_pt,percent)
    
    t2 <-head(tdb_long_percent,10)
    dim(tdb_long_percent)
    head(tdb_long_percent)
    table(tdb_long_percent$iso_pt)
    
    tdb_long <- inner_join(x=tdb_long_count,y=tdb_long_percent,by=c("cell_id_unique_iso")) %>%
      select(cell_id_unique,gex_annotation,iso,count,percent)
    dim(tdb_long)
    head(tdb_long,5)
    
    #----------------------------------------------------------------------------
    # compare combined two data pt and count
    #table(ifelse(tdb_long$gex_annotation==tdb_long$gex_annotation_pt, "yes","no"))
    #table(ifelse(tdb_long$cell_id_unique==tdb_long$cell_id_unique_pt, "yes","no"))
    
    # add sample name
    tdb_long$sample<- sample_id
    head(tdb_long)
    table(tdb_long$gex_annotation)
    
    final[[i]]<- tdb_long
    final_wide[[i]] <- wide_pt
    final_wide_ct[[i]] <-wide_ct
    final_gex[[i]]<-tdb_s_df_b_gexBySample
    
    
    
  }
  # combine individual samples
  combined_data <- do.call(rbind,final)
  combined_data_wide<- do.call(rbind,final_wide)
  combined_data_wideCt<- do.call(rbind,final_wide_ct)
  combined_data_gex <- do.call(rbind,final_gex)
  # add group infor to combined_data_gex
  combined_data_gex$group <- sapply(substr(combined_data_gex$sample,1,2), "[",1)
  # add factor order to samples
  combined_data$sample <- factor(combined_data$sample,
                                 levels = c( "HD2","HD1","HD1_stD2", "p16D9","p18D13","p36D4","p36D45","p36D169","p45D24","p45D164"))
  combined_data_wide$sample <- factor(combined_data_wide$sample,
                                      levels = c("HD2","HD1","HD1_stD2", "p16D9","p18D13","p36D4","p36D45","p36D169","p45D24","p45D164" ))
  combined_data_wideCt$sample <- factor(combined_data_wide$sample,
                                        levels = c("HD2","HD1","HD1_stD2", "p16D9","p18D13","p36D4","p36D45","p36D169","p45D24","p45D164"))
  combined_data_gex$sample <- factor(combined_data_gex$sample,
                                     levels = c("HD2","HD1","HD1_stD2", "p16D9","p18D13","p36D4","p36D45","p36D169","p45D24","p45D164"))
  
  write_tsv(combined_data_wide,file = paste0(outputPath_summary,"combined_iso_widePt.tsv"))
  write_tsv(combined_data_wideCt,file = paste0(outputPath_summary,"combined_data_wideCt.tsv"))
  
  ##########################
  # barplot of gex by sample 
  ##########################
  colours<- c( "Exhausted B cells"="red1",  "Naive B cells"="royalblue2", "Non-switched memory B cells"="seagreen2", "Plasmablasts"="grey20","Switched memory B cells"="yellow2")
  
  
  
  My_Theme = theme(
    panel.border = element_rect(colour="black", fill = NA),
    strip.text = element_text(size = 4.5),
    panel.background = element_blank(),
    plot.title = element_text(color="black", size=8, face="bold"),
    axis.title.x = element_text(color="black", size=12),
    axis.title.y = element_text(color="black", size=12),
    axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size=10, face="bold"),
    axis.text.y = element_text( size=10, face="bold"))
  
  
  # stack plot 
  p_gex_1 <- combined_data_gex %>% group_by(gex_annotation,sample,  cell_id_uniq)%>%
    dplyr::summarise(count=n())  %>%ungroup() %>%
    ggplot( aes(fill=gex_annotation, y=count, x=sample)) +
    geom_bar(position="stack", stat="identity") +
    labs(title="B cell types by sample time ",x ="Samples", y = "Counts") +  
    scale_fill_manual(values = colours) +My_Theme 
  p_gex_1
  
  
  
  
  #---------------------------------------------------------------
  # percnetage plot 
  p_gex_pt1 <- combined_data_gex %>% group_by(gex_annotation,sample,  cell_id_uniq,group)%>%
    dplyr::summarise(count=n())  %>%ungroup() %>%
    ggplot( aes(fill=gex_annotation, y=count, x=sample)) +
    geom_bar(position="fill", stat = "identity") +labs(title="B cell types by sample time ",
                                                       x ="Samples", y = "Percentage") +  scale_fill_manual(values = colours)  +My_Theme
  p_gex_pt1
  
  ggsave(filename=paste0(outputPath_plot,"gex_bySample.png"),
         plot=p_gex_pt1, width=length(unique(combined_data_gex$sample))*2, height=12, units="cm",dpi = 300)
  
  
  
  return(combined_data)
}



###################################### Base processing ##################################################

##### LogQuant Transformation #####
LogQuantTrans <- function(df, annotation, basic, logBase = 10) {
  # Does a log (the base can be selected) quantile normalization of the data frame given
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations
  # logbase = base for the log normalization, 10 is chosen by default
  countsCols <- c((length(basic)+1):(length(c(basic,rownames(annotation)))))
  dfLog <- df
  dfLog[,countsCols] <- log(df[,countsCols]+1, logBase)
  dfLogQuan <- dfLog
  dfLogQuan[,countsCols] <- normalize.quantiles(as.matrix(dfLog[,countsCols]))
  dfLogQuan <- as.data.frame(dfLogQuan)
  colnames(dfLogQuan) <- colnames(dfLog)
  rownames(dfLogQuan) <- rownames(dfLog)
  return(dfLogQuan)
}

### Transformation to RPKMs ###
RPKMs <- function(df, annotation, basic) {
  # Normalises the given dataframe to RPKMs
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations
  basic$Length <- basic$End - basic$Start
  basic$Strand <- NULL # Needed to keep the number of columns the same
  countsCols <- c((length(basic)+1):(length(c(basic,rownames(annotation)))))
  frag_lengths_kb <- basic$Length / 1000
  Total_counts_per_M <- colSums(df[,countsCols]) / 1000000
  df[,countsCols] <- mapply('/', df[,countsCols], Total_counts_per_M)
  df[,countsCols] <- df[,countsCols] / frag_lengths_kb
  return(df)
}

##### Filter peaks #####
peaks_found <- function(df, annotation, name = "All", cut_off = 5){
  # Filters out peaks in less than a third of the population
  # If a peak has an read count > cut_off in 1/3 or more of a population it is considered to be found
  # df = counts table
  # annotation = annotation table
  # name = AgeCelltype that is looked at or "All" (default) if all samples should be considered
  # cut_off = Peaks above this cutoff in 1/3 or more of a population are considered found; default = 5
  countsCols <- c((length(basic)+1):(length(c(basic,rownames(annotation)))))
  df[,countsCols] <- lapply(df[,countsCols], as.numeric)
  data <- df[,countsCols]
  data[] <- 0
  data[df[,countsCols] > cut_off] <- 1
  if (name == "All"){
    n <- nrow(annotation) / 3
    peaks_found <- which(rowSums(data) > n)
  } else {
    index <- which(annotation$AgeCelltype == name)
    n <- length(index) / 3
    peaks_found <- which(rowSums(data[,index]) > n)
  }
  return(peaks_found)
}

###################################### Differential analysis ##################################################

##### Differential analysis with DESeq2 #####
DESeq2 <- function(df, annotation, basic, GroupCol = "Cell_type"){
  # Main function for differential analysis
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations
  # GroupCol = column of the annotation df specifying the groups to compare; default = "Cell_type"
  countsCols <- c((length(basic)+1):(length(c(basic,rownames(annotation)))))
  data <- as.matrix(df[,countsCols])
  data <- floor(data + 0.5) # DEseq takes integer values only but we have 1/2 reads due to HOMERs processing read-pairs
  dds <- DESeqDataSetFromMatrix(countData = data, colData = annotation, design = formula(paste("~", GroupCol)))
  dds <- DESeq(dds)
  return(dds)
}

Pairwise_differential <- function(My_DESeq, df, group1, group2, Groupcol = "Cell_type", p_max){
  # Pairwise comparisons for differential analysis
  # My_DESeq = DESeq object returned by running the DESeq2() function
  # df = dataframe to add columns with differential analysis results to
  # group1 = first level of the GroupCol to compare
  # group2 = second level of the GroupCol to compare
  # GroupCol = column of the annotation df specifying the compared groups; default = "Cell_type"
  # p_max = pvalue cutoff used
  temp <- results(My_DESeq, contrast = c(Groupcol,group1,group2), alpha = p_max)
  print(summary(temp))
  df$temp_log2FoldChange <- temp$log2FoldChange
  df$temp_pvalue <- temp$pvalue
  df$temp_padj <- temp$padj
  colnames(df)[(ncol(df)-2):ncol(df)] <- c(paste0("log2FoldChange_",group1,"_vs_",group2),
                                           paste0("pvalue_",group1,"_vs_",group2),
                                           paste0("padj_",group1,"_vs_",group2))
  return(df)
}

find_up_down <- function(df, group1, group2, p_max, log_FC_min = 0, dir = "up"){
  # Find significantly up/down regulated regions
  # df = dataframe containing a column with the differential analysis results from the Pairwise_differential() function
  # group1 = first level of the GroupCol to compare
  # group2 = second level of the GroupCol to compare
  # p_max = pvalue cutoff used
  # log_FC_min = Fold change cutoff used; default=0
  # dir = Determines if only upregulated ("up", default) regions, downregulated ("down") regions or both ("both") should be returned
  padj <- df[,colnames(df) == paste0("padj_",group1,"_vs_",group2)]
  df_2 <- df[((!is.na(padj)) & (padj < p_max)),]
  padj_2 <- df_2[,colnames(df) == paste0("padj_",group1,"_vs_",group2)]
  df_2 <- df_2[order(padj_2),]
  
  logFC <- df_2[,colnames(df_2) == paste0("log2FoldChange_",group1,"_vs_",group2)]
  
  if(dir == "up"){
    df_return <- df_2[logFC > log_FC_min,]
  } else if(dir == "down"){
    df_return <- df_2[logFC < (-log_FC_min),]
  } else if (dir == "both"){
    df_return <- df_2[(logFC < (-log_FC_min)) | logFC > log_FC_min,]
  } else {
    df_return <- NULL
  }
  return(df_return)
}

###################################### Plotting ##################################################

##### Boxplot of distribution #####
boxplt <- function(df, annotation, basic, plotName) {
  # Creates a boxplot with all the counts columns
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations
  # plotName = string containing title of plot
  countsCols <- c((length(basic)+1):(length(c(basic,rownames(annotation)))))
  par(oma = c(6,0,0,0))
  boxplot(df[,countsCols], las = 2, cex.axis = 0.3, main = paste(plotName,"- Boxplot",sep = " "), names = annotation$Short_name)
}

##### PCA plot #####
PCAplt <- function(df, annotation, basic, groupVector, groupVector_shape = NULL, plotName,
                   colorVector = NULL, showFilename = TRUE, pc1 = 1, pc2 = 2) {
  # Creates a PCA plot with all the counts columns
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations
  # groupVector = vector that defines how the points are coloured (usually a column of the annotation df)
  # groupVector_shape = vector that defines how the points are shaped, can be left empty
  # plotName = string containing title of plot
  # colorVector = vector defining the colours used in the PCA plot, can be left empty
  # showFilename = should sample names be displayed next to the points; default = TRUE
  # pc1 = PC used for the x-axes; default = 1
  # pc2 = PC used for the y-axes; default = 2
  # return_pca = should the pca data be returned; default = FALSE
  countsCols <- c((length(basic)+1):(length(c(basic,rownames(annotation)))))
  data <- df[,countsCols]
  names(data) <- annotation$Short_name
  row.names(data) <- df$Location
  groupVector <- as.factor(groupVector)
  groupVector <- droplevels(groupVector)

  par(oma = c(0,0,0,0))
  par(mar = c(4,4,4,10) + 0.1)
  par(xpd = TRUE)
  pca <- prcomp(t(data), scale. = TRUE)
  
  p <- ggplot(data = as.data.frame(pca$x),
              aes_string(x = colnames(pca$x)[pc1], y = colnames(pca$x)[pc2], colour = groupVector)) +
        theme_classic() +
        labs(title = paste(plotName,"\n n =",nrow(data),"- PCA",sep = " "),
             x = paste("PC",pc1," - ",(round(summary(pca)[[6]][2,pc1],2)*100),"%",sep = ""),
             y = paste("PC",pc2," - ",(round(summary(pca)[[6]][2,pc2],2)*100),"%",sep = ""))
  
  if(!is.null(groupVector_shape)) {
    groupVector_shape <- as.factor(groupVector_shape)
    groupVector_shape <- droplevels(groupVector_shape)
    shape_order <- c(1,15,17,6,0,16,2,8)
    p <- p + geom_point(aes(shape = groupVector_shape), size = 2, stroke = 1) +
             scale_shape_manual(values = shape_order[1:length(groupVector_shape)])
  } else {
    p <- p + geom_point() 
  }
  if(showFilename == TRUE) {p <- p + geom_text_repel(aes(label = colnames(data)))}
  if(!is.null(colorVector)) {p <- p + scale_color_manual(name = "", drop = FALSE, values = colorVector)}
  print(p)
}

##### Venn diagram #####
vennplt <- function(df1, df2, df3 = NA, df4 = NA, categories, n_cat = 2){
  # Creates venn diagrams with 2 or 3 categories
  # df1, df2, df3, and df4 = vectors with the regions of interest to be overlapped; df3 and df4 only need to be specified if n_cat = 3 or 4
  # categories = vector with names of the groups to be overlapped
  # n_cat = number of groups to be overlapped; 2,3 or 4 are possible; default = 2
  if (n_cat == 2) {
    overlaps <- c(length(df1), length(df2), length(intersect(df1,df2)))
    
    venn.plot <- draw.pairwise.venn(
      area1 = overlaps[1], area2 = overlaps[2], cross.area = overlaps[3],
      category = categories,
      fill = c("dodgerblue", "darkorange1"),
      cat.col = c("dodgerblue", "darkorange1"),
      margin = 0.08,
      ind = TRUE
    )  
  } else if (n_cat == 3) {
    overlaps <- c(length(df1), length(df2), length(df3),
                  length(intersect(df1,df2)), length(intersect(df1,df3)), length(intersect(df2,df3)),
                  length(intersect(intersect(df1,df2),df3)))
    
    venn.plot <- draw.triple.venn(
      area1 = overlaps[1], area2 = overlaps[2], area3 = overlaps[3],
      n12 = overlaps[4], n13 = overlaps[5], n23 = overlaps[6], n123 = overlaps[7],
      category = categories,
      fill = c("dodgerblue", "darkorange1", "seagreen3"),
      cat.col = c("dodgerblue", "darkorange1", "seagreen3"),
      margin = 0.08,
      ind = TRUE
    )   
  } else if (n_cat == 4) {
    overlaps <- c(length(df1), length(df2), length(df3), length(df4),
                  length(intersect(df1,df2)), length(intersect(df1,df3)), length(intersect(df1,df4)), length(intersect(df2,df3)), length(intersect(df2,df4)), length(intersect(df3,df4)), 
                  length(intersect(intersect(df1,df2),df3)), length(intersect(intersect(df1,df2),df4)), length(intersect(intersect(df1,df3),df4)), length(intersect(intersect(df2,df3),df4)),
                  length(intersect(intersect(intersect(df1,df2),df3),df4)))
    venn.plot <- draw.quad.venn(
      area1 = overlaps[1], area2 = overlaps[2], area3 = overlaps[3], area4 = overlaps[4], 
      n12 = overlaps[5], n13 = overlaps[6], n14 = overlaps[7], n23 = overlaps[8], n24 = overlaps[9], n34 = overlaps[10], n123 = overlaps[11], n124 = overlaps[12], n134 = overlaps[13], n234 = overlaps[14], n1234 = overlaps[15],
      category = categories,
      fill = c("dodgerblue", "darkorange1", "seagreen3", "grey"),
      cat.col = c("dodgerblue", "darkorange1", "seagreen3", "grey"),
      margin = 0.08,
      ind = TRUE
    )
  } else {
    stop("VENN diagrams can only be made for 2, 3 or 4 populations")
  }
  
  return(venn.plot)
}

##### Correlation heatmap #####
Correlation_Heatmap <- function(df, annotation, basic, plotName, clustered = TRUE, cor = "spearman", colours_Celltypes, colours_Age, rownames = FALSE){
  # Plots a heatmap showing the correlation between samples
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations
  # plotName = string containing title of plot
  # clustered = Defining if the heatmap should be clustered; default = TRUE
  # cor = Method used to calculate the correlation; "spearman" or "pearson"; default = "spearman"
  # colours_Celltypes and colours_Age = vectors defining the colours used for annotation
  # rownames = should row names be shown; default = FALSE
  countsCols <- c((length(basic)+1):(length(c(basic,rownames(annotation)))))
  myorder <- order(annotation$Cell_type, annotation$Age) # needed if the heatmap is not clustered
  
  plot_df <- cor(df[countsCols[myorder]], df[countsCols[myorder]], method = cor)
  rownames(plot_df) <- annotation$Short_name[myorder]
  colnames(plot_df) <- annotation$Short_name[myorder]
  
  group_annotation <- data.frame(Population = as.character(annotation$Cell_type[myorder]),
                                 Age = as.character(annotation$Age[myorder]))
  row.names(group_annotation) <- colnames(plot_df)
  names(colours_Celltypes) <- unique(annotation$Cell_type[myorder])
  names(colours_Age) <- unique(annotation$Age[myorder])
  my_colours <- list(Population = colours_Celltypes,
                     Age = colours_Age)
  
  sampleDists <- as.dist(1 - plot_df)
  print(pheatmap(plot_df, cluster_cols = clustered, cluster_rows = clustered, show_colnames = FALSE,
                 show_rownames = rownames, main = plotName, annotation_col = group_annotation,
                 annotation_row=group_annotation, annotation_colors=my_colours, clustering_distance_rows = sampleDists,
                 clustering_distance_cols = sampleDists))
}

##### Row normalisation of heatmap #####
Row_norm <- function(df, annotation, basic){
  # Row normalises the input dataframe
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations  
  countsCols <- c((length(basic)+1):(length(c(basic,row.names(annotation)))))
  df_rownorm <- df
  df_rownorm[,countsCols] <- t(apply(df[,countsCols], 1, function(x)(x-min(x))/(max(x)-min(x))))
  return(df_rownorm)
}

##### Heatmap of accessibility #####
Heat_map <- function(df, annotation, basic, plotName, xcluster = FALSE,rownames = TRUE,ycluster = TRUE,
                     colours_Celltypes, colours_Age, nclust = NA, my_row_order = NA, order = "Celltype", fixed_height = FALSE) {
  # This function creates a heatmap with the expression of each peak in each sample
  # Samples and peaks are clustered based on hierarchical clustering
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations
  # plotName = string containing title of plot
  # xcluster = should samples be clustered; default=FALSE
  # rownames = should rownames be displayed; default=FALSE
  # ycluster = should the regions be clustered; default =TRUE
  # colours_Celltypes and colours_Age = vectors defining the colours used for annotation
  # nclust = number of kmeans clusters to make, by default no clusters are made
  # my_row_order = order to display the rows of the heatmap in (can be left unspecified)
  # order = Column of the annotation df that the samples should be ordered by; default = "Cell_type"
  # fixed_height = should the squares of the heatmap have a fixed height; default = FALSE
  countsCols <- c((length(basic)+1):(length(c(basic,rownames(annotation)))))

  if(order == "Celltype"){
    myorder <- order(annotation$Cell_type, annotation$Age) # needed if the heatmap is not clustered
  } else if (order == "Age"){
    myorder <- order(annotation$Age, annotation$Cell_type) # needed if the heatmap is not clustered
  }

  if (is.na(my_row_order)[1]){
    my_row_order <- c(1:nrow(df))
    cluster_rows <- T
  } else {
    cluster_rows <- F
  }
  
  group_annotation <- data.frame(Population = as.character(annotation$Cell_type[myorder]),
                                 Age = as.character(annotation$Age[myorder]))
  row.names(group_annotation) <- colnames(df[,countsCols[myorder]])
  names(colours_Celltypes) <- unique(annotation$Cell_type[myorder])
  names(colours_Age) <- unique(annotation$Age[myorder])
  my_colours <- list(Population = colours_Celltypes,
                     Age = colours_Age)
  
  if (fixed_height == FALSE){
    pheatmap(df[my_row_order,countsCols[myorder]], cluster_cols = xcluster, cluster_rows = ycluster,
           main = plotName, show_colnames = FALSE, annotation_col = group_annotation,
           show_rownames = rownames, labels_row = as.character(df$Gene.Name[my_row_order]), fontsize_row = 2,
           annotation_colors = my_colours, clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
           kmeans_k = nclust, drop_levels = TRUE)
  } else {
    pheatmap(df[my_row_order,countsCols[myorder]], cluster_cols = xcluster, cluster_rows = ycluster,
             main = plotName, show_colnames = FALSE, annotation_col = group_annotation,
             show_rownames = rownames, labels_row = as.character(df$Gene.Name[my_row_order]), fontsize_row = 2,
             annotation_colors = my_colours, clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
             kmeans_k = nclust, drop_levels = TRUE, cellheight = 0.3, cellwidth = 6)
  }
}

##### Heatmap of accessibility with peak annotation #####
Heat_map_with_row_annotation <- function(df, annotation, basic, row_annotation, plotName, xcluster = FALSE, rownames = TRUE,
                                         ycluster = TRUE, colours_Celltypes, colours_Age, my_row_order = NA, order = "Celltype") {
  # This function creates a heatmap with the expression of each peak in each sample
  # Samples and peaks are clustered based on hierarchical clustering
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations
  # row_annotation = dataframe containingthe information to annotate the peaks
  # plotName = string containing title of plot
  # xcluster = should samples be clustered; default=FALSE
  # rownames = should rownames be displayed; default=FALSE
  # ycluster = should the regions be clustered; default=TRUE
  # colours_Celltypes and colours_Age = Vectors defining the colours used for annotation
  # my_row_order = order to display the rows of the heatmap in (can be left unspecified)
  # order = column of the annotation df that the samples should be ordered by; default = "Cell_type"
  countsCols <- c((length(basic)+1):(length(c(basic,rownames(annotation)))))

  if(order == "Celltype"){
    myorder <- order(annotation$Cell_type, annotation$Age) # needed if the heatmap is not clustered
  } else if (order == "Age"){
    myorder <- order(annotation$Age, annotation$Cell_type) # needed if the heatmap is not clustered
  }
  
  if (sum(is.na(my_row_order)) > 0){
    my_row_order <- c(1:nrow(df))
    cluster_rows <- T
  } else {
    cluster_rows <- F
  }
  
  group_annotation <- data.frame(Population = as.character(annotation$Cell_type[myorder]),
                                 Age = as.character(annotation$Age[myorder]))
  rownames(group_annotation) <- colnames(df[,countsCols[myorder]])
  names(colours_Celltypes) <- unique(annotation$Cell_type[myorder])
  names(colours_Age) <- unique(annotation$Age[myorder])
  my_colours <- list(Population = colours_Celltypes,
                     Age = colours_Age,
                     Diff_1 = c("TRUE" = "black", "FALSE" = "white"),
                     Diff_2 = c("TRUE" = "black", "FALSE" = "white"))

  row_annotation <- as.data.frame(row_annotation)
  rownames(row_annotation) <- rownames(df)
  colnames(row_annotation) <- c("Diff_1","Diff_2")
  row_annotation$Diff_1 <- as.factor(row_annotation$Diff_1)
  row_annotation$Diff_2 <- as.factor(row_annotation$Diff_2)

  pheatmap(df[my_row_order,countsCols[myorder]], cluster_cols = xcluster, cluster_rows = ycluster,
           main = plotName, show_colnames = FALSE, annotation_col = group_annotation,
           show_rownames = rownames, labels_row = as.character(df$Gene.Name[my_row_order]), fontsize_row = 2,
           annotation_colors = my_colours, clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
           annotation_row = row_annotation, drop_levels = TRUE)
}

##### Divide heatmap into clusters #####
divide_into_clust <- function(my_heatmap, k, df){
  # Divides a heatmap into n clusters
  # my_heatmap = heatmap returned by the Heat_map() function
  # k = number of clusters to divide the heatmap into
  # df = dataframe that the heatmap is based on
  my_heatmap_ord <- my_heatmap$tree_row$order
  
  pos_list_unordered <- list()
  locs <- rep(NA, k)
  for (i in 1:k){
    pos_list_unordered[[paste0("clust",i)]] <- which(cutree(my_heatmap$tree_row, k=k) == i)
    locs[i] <- df[pos_list_unordered[[paste0("clust",i)]],] %>% pull(Location) %>% .[1]
  }
  cluster_order <- df[my_heatmap_ord,] %>% pull(Location) %>% {match(locs,.)} %>% order() # Returns the clusters in order top down
  
  pos_list <- list()
  for (i in 1:k){
    pos_list[[paste0("clust",i)]] <- pos_list_unordered[[paste0("clust",cluster_order[i])]]
  }
  
  return(pos_list)
}

##### Box and Whiskers plot #####
box_and_whiskers <- function(df, annotation, basic, ymin = 0.7, ymax = 1.3, plot_name = "", groups = "AgeCelltype"){
  # df = counts table
  # annotation = annotation table
  # basic = array for annotations
  # ymin = lower cutoff of the plot; default = 0.7
  # ymax = upper cutoff of the plot; default = 1.3
  # plotName = string containing title of plot
  # groups = column of the annotation df that the samples should be grouped by ; default = "AgeCelltype" 
  temp <- rowMeans(df[,c((length(basic)+1):(length(c(basic,rownames(annotation)))))])
  df_norm <- df
  df_norm[,c((length(basic)+1):(length(c(basic,rownames(annotation)))))] <- df_norm[,c((length(basic)+1):(length(c(basic,rownames(annotation)))))]/temp

  my_df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(my_df) <- c("Med_temp", "Group", "Peak")
  
  if (groups == "AgeCelltype"){
    for (i in annotation$AgeCelltype) {
      Med_temp <- rowMedians(as.matrix(df_norm[,length(basic)+which(annotation$AgeCelltype==i)]))
      my_df_temp <- data.frame(Med_temp)
      my_df_temp$Group <- i
      my_df_temp$Peak <- row.names(my_df_temp)
      my_df <- rbind(my_df,my_df_temp)
    }
    
    p <- ggplot(my_df, aes(x = Group, y = Med_temp)) +
      geom_boxplot(outlier.shape = NA, color = "grey30") +
      labs(title = plot_name, x = "", y = "Median normalized ATAC signal") +
      theme_classic() +
      ylim(ymin, ymax)
  } else if (groups == "Age"){
    for (i in annotation$Age) {
      Med_temp <- rowMedians(as.matrix(df_norm[,length(basic)+which(annotation$Age==i)]))
      my_df_temp <- data.frame(Med_temp)
      my_df_temp$Group <- i
      my_df_temp$Peak <- row.names(my_df_temp)
      my_df <- rbind(my_df,my_df_temp)
    }
    my_df$Group <- factor(my_df$Group, levels=c("Young","Adult","Old"))
    
    p <- ggplot(my_df, aes(x = Group, y = Med_temp)) +
      geom_boxplot(outlier.shape = NA, color = "grey30", fill = colours_Age) +
      labs(title = plot_name, x = "", y = "Median normalized ATAC signal") +
      theme_classic() +
      ylim(ymin, ymax)
  }

  print(p)
}

###################################### Gene ontology analysis ##################################################

##### GREAT analysis #####
GREAT_analysis <- function(INTERESTING_GREAT_GO, input_file_names, skip = FALSE, topn = 5, sig_cut_off = 0.05){
  # INTERESTING_GREAT_GO = Ontologies to look at
  # input_file_names = vector of bed files to be included
  # skip = should the initial step be skipped and only plotting done; default=FALSE
  # topn = maximum number of enriched terms to return for each region list
  # sig_cut_off = maximum adjusted p-value for a term to be considered significantly enriched
  input_path <- "../output_files/ATACseq_analysis/"
  output_path <- "../output_files/GREAT/"
  
  diff_genes <- gsub("|.bed", "", input_file_names)
  
  if (skip == FALSE){
    for (i in 1:length(input_file_names)){
      input_file_path <- paste0(input_path, input_file_names[i])
      peaks <- read.table(input_file_path, row.names = NULL, header = F, sep = "\t", comment.char = "", quote = "\"")
      peaks <- peaks[,1:4]
      colnames(peaks) <- c("Chr", "Start", "End", "Location")
      
      job <- submitGreatJob(peaks, species = "mm10", rule = "basalPlusExt", version = "4.0.4")
      
      for (GREAT_GO in INTERESTING_GREAT_GO){
        tb = getEnrichmentTables(job,ontology = GREAT_GO)
        df_specific <- as.data.frame(tb[GREAT_GO])
        colnames(df_specific) <- gsub(paste0(gsub(" ", ".", GREAT_GO), "."), "", colnames(df_specific))
        write.table(df_specific, paste0(output_path,GREAT_GO,"_",diff_genes[i],".txt"), quote = F, row.names = F, col.names = T, sep = "\t")
      }
    }
  }
  
  Fold_enrichment_cutoff <- 2
  p <- list()
  
  diff_genes_short <- gsub("Diff_merged_Y_O_main_", "", diff_genes) # To make column names shorter
  diff_genes_short <- gsub("Oneg_Ypos_subtracted_MACS_p0.01_", "", diff_genes_short) # To make column names shorter

  for (GREAT_GO in INTERESTING_GREAT_GO){
    df_go <- NULL
    
    for (i in 1:length(diff_genes)){
      df_specific <- read.table(paste0(output_path, GREAT_GO,"_",diff_genes[i],".txt"),
                                row.names = NULL, header = T, sep = "\t", comment.char = "", quote = "\"")
      df_specific <- df_specific[1:min(which(df_specific$Binom_Adjp_BH==1)-1),c("name","Binom_Adjp_BH","Binom_Fold_Enrichment","Binom_Raw_PValue","Hyper_Adjp_BH")]
      colnames(df_specific) <- c("name",
                                 paste0("Binom_Adjp_BH_",diff_genes_short[i]),
                                 paste0("Binom_Fold_Enrichment_",diff_genes_short[i]),
                                 paste0("Binom_Raw_PValue_",diff_genes_short[i]),
                                 paste0("Hyper_Adjp_BH_",diff_genes_short[i]))
      if (!is.null(df_go)){
        df_go <- merge(df_go, df_specific, all = TRUE)
        df_go[is.na(df_go)] <- 1 
      }else{
        df_go <- df_specific
      }
    }
    
    colnames(df_go) <- gsub("Binom_Adjp_BH_", "", colnames(df_go))
    topn_go_per_cat <- NULL
    for (i in seq(2,ncol(df_go),4)){
      df_go <- df_go[order(df_go[,i+2]),]
      topn_go_per_cat_temp <- df_go[df_go[,i] < sig_cut_off,]
      topn_go_per_cat_temp <- topn_go_per_cat_temp[topn_go_per_cat_temp[,(i+3)] < sig_cut_off,]
      topn_go_per_cat_temp <- topn_go_per_cat_temp[topn_go_per_cat_temp[,(i+1)] > Fold_enrichment_cutoff,]
      topn_go_per_cat_temp <- topn_go_per_cat_temp[1:min(topn,nrow(topn_go_per_cat_temp)),]
      topn_go_per_cat <- rbind(topn_go_per_cat, topn_go_per_cat_temp)
    }
    topn_go_per_cat <- unique(topn_go_per_cat)
    topn_go_per_cat <- topn_go_per_cat$name
    
    df_go[,seq(3,ncol(df_go),4)] <- NULL
    df_go[,seq(3,ncol(df_go),3)] <- NULL
    df_go[,seq(3,ncol(df_go),2)] <- NULL
    
    df_plot <- -log10(df_go[df_go$name %in% topn_go_per_cat,-1])
    names <- df_go[df_go$name %in% topn_go_per_cat,1]
    df_plot[df_plot==Inf] <- 300 # Otherwise it would be Inf
    row.names(df_plot) <- names
    
    mat <- as.matrix(df_plot)
    
    quantile_breaks <- function(xs, n = 10) {
      breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
      breaks[!duplicated(breaks)]
    }
    
    mat_breaks <- quantile_breaks(mat, n = 30)
    
    p[[GREAT_GO]] <- pheatmap(df_plot, show_rownames = T,
                              main = paste0("Top ",topn," pvalues per category \n with cut-off ",sig_cut_off," ",GREAT_GO),
                              cluster_cols = FALSE, cluster_rows = TRUE, treeheight_row = 0, fontsize_col = 9, fontsize_row = 12,
                              color = colorRampPalette(brewer.pal(9, "Blues"))(length(mat_breaks) - 1), breaks = mat_breaks)
    }
  return(p)
}

###################################### Motif enrichment ##################################################

##### Motif enrichment #####
split_name <- function(df){
  # Reformats the output from HOMER for Motif enrichment
  # Used by the motif_plots() function
  # df = dataframe returned by HOMER after motif enrichment
  df$motif_name <- sapply(strsplit(as.character(df$Motif.Name),"\\/"),`[`, 1)
  df$motif_name_short <- sapply(strsplit(as.character(df$motif_name),"\\("),`[`, 1)
  df$motif_family <- sapply(strsplit(as.character(df$motif_name),"\\("),`[`, 2)
  df$motif_family <- sapply(strsplit(as.character(df$motif_family),"\\)"),`[`, 1)
  df$experiment <- sapply(strsplit(as.character(df$Motif.Name),"\\/"),`[`, 2)
  df$database <- sapply(strsplit(as.character(df$Motif.Name),"\\/"),`[`, 3)
  df$temp <- sapply(strsplit(as.character(df$X..of.Target.Sequences.with.Motif),"%"),`[`, 1)
  df$temp2 <- sapply(strsplit(as.character(df$X..of.Background.Sequences.with.Motif),"%"),`[`, 1)
  df$Fold_enrichment <- as.numeric(df$temp) / as.numeric(df$temp2)
  df$logq <- df$q.value..Benjamini.
  df$logq[df$logq == 0] <- 0.0001
  df$logq <- -(log10(df$logq))
  df$temp <- NULL
  df$temp2 <- NULL
  return(df)
}

motif_plots <- function(df_list, logp_min = 30, logp_max = 100000000, TFs = c(), logq_min = 2, logq = FALSE){
  # Makes a plot to show the enriched motifs in different peak sets
  # df_list = List of dataframes from the knownResults.txt files output by HOMER
  # logp_min = Minimal log p-value for a TF to be included; default = 30
  # logp_max = Limit to the colour scale at this maximal log p-value (by default not limited)
  # TFs = if specified only these TFs will be included in the plot
  # logq_min = Minimal log q-value for a TF to be included; default = 2; only used if logq=TRUE is specified
  # logq = Use log q-value instead of log p-value to select TFs; default = FALSE
  
  n <- length(df_list)
  
  for (i in 1:n){
    df_list[[i]] <- split_name(df_list[[i]])
    df_list[[i]] <- df_list[[i]][order(df_list[[i]]$Motif.Name),]
  }
  
  df_logp <- df_list[[1]][,c("motif_name","motif_name_short","motif_family","Motif.Name","Consensus")]
  df_logq <- df_list[[1]][,c("motif_name","motif_name_short","motif_family","Motif.Name","Consensus")]
  df_FC <- df_list[[1]][,c("motif_name","motif_name_short","motif_family","Motif.Name","Consensus")]
  df_Target_Percent <- df_list[[1]][,c("motif_name","motif_name_short","motif_family","Motif.Name","Consensus")]
  
  for (i in 1:n){
    df_logp[,(i+5)] <- -df_list[[i]]$Log.P.value
    df_logq[,(i+5)] <- df_list[[i]]$logq
    df_FC[,(i+5)] <- df_list[[i]]$Fold_enrichment
    df_Target_Percent[,(i+5)] <- df_list[[i]]$X..of.Target.Sequences.with.Motif %>%
      as.character() %>% strsplit("%") %>% sapply(`[`, 1) %>% as.numeric()
  }
  
  if(length(TFs) > 0) {
    index <- which(df_logq$motif_name_short %in% TFs) # Plot selected TFs
  } else {
    index <- c()
    for (i in 1:n){
      if (logq == FALSE){
        index_temp <- which(df_logp[,(5+i)] >= logp_min) # Plot all TFs with logp > logp_min
        index <- unique(c(index,index_temp))
      } else {
        index_temp <- which(df_logq[,(5+i)] >= logq_min) # Plot all TFs with logq > logq_min
        index <- unique(c(index,index_temp))        
      }
    }
  }
  
  df_logp_sig <- df_logp[index,]
  df_FC_sig <- df_FC[index,]
  df_logq_sig <- df_logq[index,]
  df_Target_Percent_sig <- df_Target_Percent[index,]
  
  df_logp_sig_melted <- melt(df_logp_sig)
  df_FC_sig_melted <- melt(df_FC_sig)
  df_logq_sig_melted <- melt(df_logq_sig)
  df_Target_Percent_sig_melted <- melt(df_Target_Percent_sig)
  df_Target_Percent_sig_melted$logp <- df_logp_sig_melted$value
  df_Target_Percent_sig_melted$logq <- df_logq_sig_melted$value
  
  df_Target_Percent_sig_melted$logp[df_Target_Percent_sig_melted$logp > logp_max] <- logp_max # For logp values > logp_max use logp_max instead
  
  quantile_breaks <- function(xs, n = 30) {
    breaks <- quantile(xs, probs=seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  
  mat_breaks <- quantile_breaks(df_Target_Percent_sig_melted$logp)
  
  if(length(TFs) > 0) {
    df_Target_Percent_sig_melted$motif_name_short <- factor(df_Target_Percent_sig_melted$motif_name_short, levels = TFs)
  } else {
    df_Target_Percent_sig_melted$motif_name_short <- factor(df_Target_Percent_sig_melted$motif_name, levels = unique(df_Target_Percent_sig_melted$motif_name[order(df_Target_Percent_sig_melted$motif_family)]))
  }
  
  p <- ggplot(df_Target_Percent_sig_melted, aes(x = motif_name_short, y = variable)) +
    geom_point( aes(size = value,color = logp)) +
    theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    scale_colour_gradient(low = "white", high = "#073763") +
    scale_size_area(max_size = 4) +
    labs(color = "-ln(p-value)", size = "Peaks with Motif (%)") +
    scale_y_discrete(labels = as.character(c(1:n))) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  if (logq == FALSE){
    p <- p + ggtitle(paste0("TF families with -ln(p-value)>", logp_min))
  } else {
    p <- p + ggtitle(paste0("TF families with -ln(q-value)>", logq_min))
  }
  
  print(p)
}
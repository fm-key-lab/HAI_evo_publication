## Load libraries
library(ape) ## to modify tree file
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggtree)
library(ggstar) ## to plot stars as shape
library(ggnewscale)  # for multiple fill scales in ggplot
library(RColorBrewer)
# library(treeio)

## NOTE:
# One need to let at first run refbased plottings through to get all the first observational dates of that species and then the denovo dates!
use_ml_tree = TRUE
analysis_types <- c("refbased", "denovo")
outpath <- '/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/phylogenetic_trees'

for (analysis_type in analysis_types){
  ## correct sources which have been collapsed before
  corrected_source <- c(
                          'P07_MIBI101' = 'Pleural empyema',
                          'P07_MIBI103' = 'Pleural empyema',
                          'P07_MIBI87' = 'Pleural empyema',
                          'P07_MIBI88' = 'Pleural empyema',
                          'P07_MIBI89' = 'Pleural empyema',
                          'P07_MIBI91' = 'Pleural empyema',
                          'P07_MIBI93' = 'Pleural empyema',
                          'P07_MIBI95' = 'Pleural empyema',
                          'P07_MIBI97' = 'Pleural empyema',
                          'P07_MIBI99' = 'Pleural empyema'
  )
  if (analysis_type == "refbased"){
    analysis_path <- "/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/refbased_2024_08/analysis"
    species_observation_range <- c()
    all_STs <- c("23", "45", "58", "95", "97", "117", "131", "4260", "16602") ## manually defined upfront after running it once before so fixed colors per ST are used
    palette_colors <- brewer.pal(n = 12, name = "Set3")
    na_color <- palette_colors[9]
    palette_colors <- c(palette_colors[2:8], palette_colors[10:12])[1:length(all_STs)] ## remove first color which is similar to location color!
    color_map <- setNames(palette_colors, all_STs)
    color_map <- c(color_map, 'NA' = na_color) ## add NA value
    
  }else if(analysis_type == "denovo"){
    analysis_path <- "/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2022/denovo_2024_08/analysis/backup"
    species_observation_range <- list(
      "P07 A. pittii" = c("2021-07-27", "2021-08-17"),
      "P07 B. thetaiotaomicron" = c("2021-08-10", "2021-08-26"),
      "P07 E. coli" = c("2021-06-29", "2021-07-20"),
      "P07 E. hormaechei" = c("2021-07-27", "2021-09-21"),
      "P07 P. aeruginosa" = c("2021-07-16", "2021-08-17"),
      "P10 K. michiganensis" = c("2021-07-26", "2021-08-02"),
      "P21 E. coli" = c("2021-09-02", "2021-10-19"),
      "P21 P. mirabilis" = c("2021-09-02", "2021-11-02"),
      "P21 S. aureus" = c("2021-09-02", "2021-09-08"))
    additional_indels <- list(
      "P07 E. hormaechei" = list(
        'ampD[93]'= c('P07_1047', 'P07_1540', 'P07_1049', 'P07_1545', 'P07_1553', 'P07_1046'),
        'ampD[130]'= c('P07_1201'),
        'ampD[160]'= c('P07_1149')
    ))
    remove_singletons <- TRUE ## For plotting the MGE and indel heatmaps, set to true if also all singletons should be plotted
  }
  species_paths <- list.files(path = analysis_path, pattern = "^P00.*", full.names = TRUE)
  
  ## NOTE: FROM S aureus the blood culture positive strain will be removed due to contamination!
  remove_contamination_isolate <- list("P0021_Saureus_ASM1342v1" = c("P0021_TP12_210915_D_Blood_ST5__P21_MIBI03", "P0021_TP12_210915_D_Blood_ST5__P21_MIBI03_S2"))
  
  for (species_path in sort(species_paths)){
    patient_species <- gsub('.*/', '', species_path)
    patient <- strsplit(patient_species, '_')[[1]][1]
    patient <- paste0(substring(patient, 1,1), substring(patient, 4))
    if (analysis_type == "refbased"){
      species_label <- strsplit(patient_species, '_')[[1]][2]
      species <- paste0(substring(species_label, 1,1), '. ', substring(species_label, 2))
      species <- gsub('YT3$', '', species) ## applies to Ehorm only 
    }else if(analysis_type == "denovo"){
      species_label <- strsplit(patient_species, '_')[[1]][3]
      species <- paste0(substring(species_label, 1,1), '. ', substring(species_label, 2))
    }
    print(paste('Processing', species))

    most_recent_analysis_path <- sort(list.files(path = species_path, pattern = "", full.names = TRUE))
    if (length(most_recent_analysis_path) == 0) {
      next
    }else{most_recent_analysis_path <- most_recent_analysis_path[[length(most_recent_analysis_path)]]}
    ## Load tree
    if (use_ml_tree == TRUE){
      #tree_files <- list.files(path = most_recent_analysis_path, pattern = ".*.bestTree", full.names = TRUE)
      tree_files <- list.files(path = most_recent_analysis_path, pattern = ".*.bestTreeCollapsed", full.names = TRUE)
      tree_type = 'MLtree'
    }else{
      tree_files <- list.files(path = most_recent_analysis_path, pattern = ".*_TP0.nwk.tree", full.names = TRUE)
      tree_type = 'PStree'
    }
    if (length(tree_files) == 0) {
      next
    }else{tree <- read.tree(tree_files[length(tree_files)])}

    ## Load data to rescale tree based on number of core SNVs
    num_pos_coregenome_file <- paste0(most_recent_analysis_path, "/num_pos_coregenome.csv")
    # Rescaling of ML tree incorrect due to various substitutions/branch!!
    if (file.exists(num_pos_coregenome_file)){
      num_pos_coregenome <- read.csv(num_pos_coregenome_file)
      tree$edge.length = tree$edge.length * num_pos_coregenome$number_core_pos
    }
    
    ## Load metadata 
    ## e.g. MGE, Indels, in silico AMR prediction or AMR measurements
    indel_table_file <- paste0(most_recent_analysis_path, "/indel_table.csv")
    
    mge_most_recent_path <- strsplit(most_recent_analysis_path, "/")[[1]]
    dirname <- mge_most_recent_path[length(mge_most_recent_path)]
    mge_most_recent_path[length(mge_most_recent_path)] <- paste0('2025_08_09_', substr(dirname, 12, nchar(dirname)))
    mge_most_recent_path <- paste(mge_most_recent_path, collapse = "/")
    #mge_table_file <- paste0(mge_most_recent_path, "/mge_analysis/mobile_genetic_elements_summary.csv")
    mge_table_file <- paste0(mge_most_recent_path, "/mge_analysis_ml_tree_based/mobile_genetic_elements_cov_breadth.csv")
    
    ## ensure tree is rooted
    outgroup_tip <- grep("^OG", tree$tip.label, value = TRUE)
    rerooted_tree <- root(tree, outgroup = outgroup_tip, resolve.root = TRUE)
    rerooted_tree <- rotateConstr(rerooted_tree, sort(rerooted_tree$tip.label, decreasing = TRUE)) ## sort internal branches by sampleid
    
    ## remove reference from tree
    tips_to_remove <- grep("_ref$", rerooted_tree$tip.label, value = TRUE)
    tips_to_remove <- append(tips_to_remove, outgroup_tip)
    if (patient_species %in% names(remove_contamination_isolate)){
      tips_to_remove <- append(tips_to_remove, remove_contamination_isolate[[patient_species]])
    }
    rerooted_tree_wo_ref <- drop.tip(rerooted_tree, tips_to_remove)
    
    ## Extract metadata from tree
    tip_metadata <- strsplit(rerooted_tree_wo_ref$tip.label, "_")
    tip_sampleid <- strsplit(rerooted_tree_wo_ref$tip.label, "__")
    
    if (analysis_type == "refbased"){
      sampleid <- sapply(tip_sampleid, function(x) x[2])
      date <- sapply(tip_metadata, function(x) x[3])
      date <- as.Date(date, format = "%y%m%d")
      min_species_date <- min(date, na.rm = TRUE)
      max_species_date <- max(date, na.rm = TRUE)
      species_observation_range[[paste(patient, species, collapse = '_')]] <- c(min_species_date, max_species_date)
      days_since_obs = as.numeric(date - min_species_date)
      STs <- sapply(tip_metadata, function(str_list) {
        str_contains_ST <- grepl('ST', str_list)
        if (any(str_contains_ST)) {
          st = str_list[str_contains_ST]
          st = ifelse(grepl("^STNF", st), "NA", st)
          return(sub('\\?$', 'x', sub('^ST', '', st))) ## remove low confidence and ST strings
          } else {return(NA)}
        })
      if (grepl("Ehormaechei", patient_species)){
        STs <- rep(NA, length(tip_metadata)) ## has not STs@
        }
      ST_variant <- grepl('x$', STs, perl = TRUE) ## & !grepl('NFx$', STs, perl = TRUE) ## grep true STs (not NF) and check if they are an allele of the ST
      STs <- sub('x$', '', STs) ## remove the x at the end which is saved in ST_variant
      indel_table_t <- data.frame()
      mge_table_t <- data.frame()
    }else if(analysis_type == "denovo"){
      sampleid <- sapply(tip_sampleid, function(x) x[2])
      date <- sapply(tip_metadata, function(x) x[3])
      date <- as.Date(date, format = "%Y%m%d")
      refbased_pat_spec <- paste(patient, strsplit(species, '-')[[1]][1], collapse = '_')
      min_species_date <- as.Date(species_observation_range[[refbased_pat_spec]][1], format = "%Y-%m-%d")
      max_species_date <- as.Date(species_observation_range[[refbased_pat_spec]][2], format = "%Y-%m-%d")
      days_since_obs = as.numeric(date - min_species_date)
      STs <- rep(NA, length(tip_metadata))
      ST_variant <- rep(NA, length(tip_metadata))
      
      ########################################
      ## read in indel table if available      
      ########################################
      if (file.exists(indel_table_file)){
        indel_table <- read.csv(indel_table_file)
        indel_table_metadata_cols <- c('alt_start', 'alt_stop', 'alt_stop_pos', 'chr', 'chr_boundary', 'distance1', 'distance2', 'frameshift', 'gene', 'gene1', 'gene2', 'gene_num', 'gene_num_global', 'indel_alt', 'indel_anc', 'indel_pos_end', 'indel_pos_start', 'indel_ref', 'indel_size_gp', 'loc1', 'loc2', 'locustag', 'locustag1', 'locustag2', 'max_indel_GL_diff_to_ref', 'mut', 'mut_translation', 'note', 'num_isolates', 'orthologtag', 'orthologtag1', 'orthologtag2', 'pos', 'product', 'product1', 'product2', 'protein_id', 'sequence', 'strand', 'translation', 'type')
        
        # Unnest all columns to get indel_table with one allele per row
        for (col in c('alt_stop', 'indel_alt', 'indel_pos_start', 'indel_pos_end')) {
          indel_table[[col]] <- indel_table[[col]] %>%
            str_replace_all("\\[|\\]|'|,|\\(|\\)", "") %>%
            str_replace_all("np.int64", "") %>%
            str_replace_all("\\.", "")
          if (col %in% c('indel_pos_start', 'indel_pos_end')) {
            indel_table[[col]] <- ifelse(is.na(indel_table[[col]]), "-1", as.character(indel_table[[col]]))
          } else {
            indel_table[[col]] <- ifelse(is.na(indel_table[[col]]), "", indel_table[[col]])
          }
          indel_table[[col]] <- str_split(indel_table[[col]], " ")
        }
        indel_table <- indel_table %>% unnest(cols = all_of(c('alt_stop', 'indel_alt', 'indel_pos_start', 'indel_pos_end')))
        
        ## get aa start position
        indel_table$indel_pos_aa_start <- indel_table %>%
          select(all_of(c('indel_pos_start', 'indel_pos_end'))) %>%
          mutate(across(everything(), as.integer)) %>%
          mutate(across(everything(), ~ ifelse(is.na(.), -1, .))) %>%
          pmap_int(~ floor(min(c(...), na.rm = TRUE) / 3))
        indel_table$indel_pos_aa_start[indel_table$indel_pos_aa_start == -1] <- ""
        ## check if all columns of interest are present
        for (col in c('gene', 'locustag', 'gene1', 'locustag1', 'gene2', 'locustag2')){
          if (!(col %in% colnames(indel_table))){
            indel_table[[col]] <- ''
          }
        }
        indel_table <- indel_table %>%
          mutate(gene_id = ifelse(gene != '.', gene, gsub('.*_', 'LT', locustag)),
                 gene_id1 = ifelse(gene1 != '.', gene1, gsub('.*_', 'LT', locustag1)),
                 gene_id2 = ifelse(gene2 != '.', gene2, gsub('.*_', 'LT', locustag2)),
                 mut_name = ifelse(gene_id != '', paste0(gene_id, '[',  indel_pos_aa_start, ']'), paste0(gene_id1, '|', gene_id2))) %>%
          select(-c('gene_id', 'gene_id1', 'gene_id2'))
        indel_table_metadata_cols <- c(indel_table_metadata_cols, 'indel_pos_aa_start', 'mut_name')
        smpl_cols <- keep(colnames(indel_table), ~ !.x %in% indel_table_metadata_cols && !grepl('_GCA_|_GCF_|_ASM_', .x))
        
        indel_table_t <- indel_table %>%
          mutate(across(all_of(smpl_cols), ~ .x == indel_alt)) %>%
          select(all_of(c('mut_name', smpl_cols))) %>%
          group_by(mut_name) %>%
          summarize(across(all_of(smpl_cols), ~ any(.x))) %>%
          t() %>%
          as.data.frame()
        colnames(indel_table_t) <- indel_table_t[1, ]
        indel_table_t <- indel_table_t[-1, ,drop = FALSE]
        ## add additionally identified indels
        if (refbased_pat_spec %in% names(additional_indels)){
          indels_to_add = additional_indels[[refbased_pat_spec]]
          for (mut in names(indels_to_add)){
            indel_table_t[[mut]] <- rownames(indel_table_t) %in% indels_to_add[[mut]]
          }
        }
        indel_table_t <- indel_table_t[, sort(colnames(indel_table_t)), drop = FALSE]
        indel_table_t[] <- lapply(indel_table_t, as.logical)
        if (isTRUE(remove_singletons)){
          no_singleton_indels <- colnames(indel_table_t)[colSums(indel_table_t) > 1]
          print(paste('Removing', length(colnames(indel_table_t)) - length(no_singleton_indels),'/',length(colnames(indel_table_t)), 'singleton indels'))
          indel_table_t <- select(indel_table_t, all_of(no_singleton_indels))
        }
      }else{
        print(paste('No Indel table available for', species))
        indel_table_t <- data.frame()
      }
      
      ########################################
      ## read in mge_table if available
      ########################################
      if (file.exists(mge_table_file)){
        mge_table <- read.csv(mge_table_file)
        ## read in summary file to get chromosome and positional coordinates
        mge_table_summary <- read.csv(gsub('mobile_genetic_elements_cov_breadth.csv', 'mobile_genetic_elements_summary.csv', mge_table_file))
        mge_table_summary$COORDINATES <- sapply(mge_table_summary$COORDINATES, function(x){
          chr_s <- strsplit(x, ':')[[1]][1]
          p_s <- strsplit(strsplit(x, ':')[[1]][2], '-')[[1]][1]
          p_s <- round(as.numeric(p_s)/1000, 0)
          p_e <- strsplit(x, '-')[[1]]
          p_e <- round(as.numeric( p_e[length(p_e)] )/1000, 0)
          chr_e <- strsplit(x, ';')[[1]]
          if (length(chr_e)==1){
            paste0(chr_s, ':', p_s, '-', p_e)
          }else{
            chr_e <- strsplit(chr_e[length(chr_e)], ':')[[1]][1]
            paste0(chr_s, ':', p_s, '-', chr_e, ':', p_e)
          }
        })
        mge_table <- left_join(mge_table, mge_table_summary[c('TYPE', 'START', 'END', 'LENGTH', 'COORDINATES')], by = c('TYPE', 'START', 'END', 'LENGTH'))
        ## select only deletions
        mge_table <- mge_table[mge_table$TYPE == 0, ] 
        ## subset to isolates of interest and coordinate column only 
        mge_smpl_id_cols <- keep(colnames(mge_table), ~ !.x %in% c("TYPE", "START", "END", "LENGTH", "COORDINATES") && !grepl('_GCA_|_GCF_|_ASM_', .x))
        mge_table_t <- mge_table[c('COORDINATES', mge_smpl_id_cols)] %>%
          t() %>%
          as.data.frame()
        colnames(mge_table_t) <- mge_table_t[1, ]
        mge_table_t <- mge_table_t[-1, ,drop = FALSE]
        mge_table_t <- mge_table_t[, sort(colnames(mge_table_t)), drop = FALSE]
        mge_table_t[] <- lapply(mge_table_t, as.numeric)
        if (isTRUE(remove_singletons)){
          ## get only those mges which are at least in 2 isolates covered ≥ 95% or ≤ 75%
          no_singleton_mges <- colnames(mge_table_t)[(colSums(mge_table_t >= 0.95) > 1) & (colSums(mge_table_t <= 0.75) > 1)] ## This checks for isolates which are at least covered ≥2 isolates and deleted in ≥2 at the same time! 
          print(paste('Removing', length(colnames(mge_table_t)) - length(no_singleton_mges),'/',length(colnames(mge_table_t)), 'MGEs due to too low breadth difference amongst samples'))
          mge_table_t <- select(mge_table_t, all_of(no_singleton_mges))
        }
      }else{
        print(paste('No MGE table available for', species))
        mge_table_t <- data.frame()
      }
      
    }else{
      sampleid <- NA
      indel_table_t <- data.frame()
      mge_table_t <- data.frame()
    }
    
    tree_annotation_heatmap_df <- tibble(
      tip_label = rerooted_tree_wo_ref$tip.label,
      patient = sapply(tip_metadata, function(x) x[1]),
      timepoint = sapply(tip_metadata, function(x) x[2]),
      date = date,
      sampletype = sapply(tip_metadata, function(x) x[4]),
      location = sapply(tip_metadata, function(x) x[5]),
      sampleid = sampleid,
      days_since_obs = days_since_obs,
      ST = STs,
      ST_variant = ST_variant
    ) %>%
      mutate(sampleid = if_else(str_detect(tip_label, "^OG_"), 
                                gsub("_P.*_TP0$|^OG_|seq", "", tip_label), 
                                sampleid)
             ) %>%
      as.data.frame()
    tree_annotation_heatmap_df$location <- ifelse(tree_annotation_heatmap_df$sampleid %in% names(corrected_source), corrected_source[tree_annotation_heatmap_df$sampleid], tree_annotation_heatmap_df$location)
    
    ## set everything for outgroup to NA
    rows_to_na <- grepl("^OG_", tree_annotation_heatmap_df$tip_label)
    cols_to_na <- setdiff(names(tree_annotation_heatmap_df), c("sampleid", "tip_label"))
    tree_annotation_heatmap_df[rows_to_na, cols_to_na] <- NA
    
    rownames(tree_annotation_heatmap_df) <- tree_annotation_heatmap_df$sampleid
    
    ## replace tree tip labels with sampelid only 
    rerooted_tree_wo_ref$tip.label <- tree_annotation_heatmap_df$sampleid[match(rerooted_tree_wo_ref$tip.label, tree_annotation_heatmap_df$tip_label)] 

    ## get nodes to highlight 
    tips_to_highlight_branch <- c()
    if (length(tips_to_highlight_branch) > 0){
      highlight_clade <- MRCA(rerooted_tree_wo_ref, tips_to_highlight_branch)
    }
    
    ## get nodes to annotate
    if (grepl("Ehormaechei", patient_species)){
      gt_to_isolate_file <- "/Users/ad_lemur/Documents/Data/Drylab/mf_2020_hap/muller_plots/P07_Ehormaechei-c1/P07_Ehormaechei_genotype_to_isolates.csv"
      gt_to_isolates <- read.csv(gt_to_isolate_file, header=FALSE, row.names = NULL, sep = ',')
      
      muts_to_annotate <- c('rpiR[T138P]', 'iscR[G60D]', 'gyrA[S83Y]', 'gyrA[S83F]', 'LT12530[G275W]', 'fimZ[F126L]', 'citG[L68Q]', 'gyrA[D87G]', 
                            "dacA[T243I]", "gyrA[G81C]", "narJ[E78*]", "narL[E49*]", "romA|acrR")
      
      tips_to_mark_mut <- list()
      for (mut in muts_to_annotate){
        isolates_with_mut <- gt_to_isolates[gt_to_isolates$V1 == mut, 'V2']
        tips_to_mark_mut[[mut]] <- strsplit(isolates_with_mut, ";")[[1]]
      }
      
    }else{tips_to_mark_mut <- c()}
    
    ## generate for heatmaps a color scale
    location_colors <- c("Nasal" = "#72e0ba",
                         "Oral" = "#41c49e",
                         "Skin" = "#1b9e77",
                         "Rectal" = "#146252",
                         "Gastric" = "#0f3f33",
                         "Lung" = "#eb4824", 
                         "Blood" = "#eb4824",
                         "Pleural empyema" = "#EB785E")
    mandatory_sites_for_legend <- c("Nasal", "Oral", "Skin", "Rectal")
    
    ####################
    ## plot initial tree 
    ####################
    max_tree_depth <- max(node.depth.edgelength(rerooted_tree_wo_ref))
    number_tips <- length(rerooted_tree_wo_ref$tip.label)
    ## heatmap specific variables
    colnames_offset <- number_tips/250
    width = 0.1
    width_offset = max_tree_depth*(width*1.1)
    plot_offset = max_tree_depth*0.025
    
    print(species)
    print(max_tree_depth)
    
    ## plot tree
    treescale_width = NULL
    treescale_x <- max_tree_depth + plot_offset
    if ((max_tree_depth >= 0.95) & (max_tree_depth <= 20)){
      treescale_width <- 1
      treescale_x <- treescale_x - treescale_width
    }else if ((max_tree_depth > 20) & (max_tree_depth <= 60)){
      treescale_width <- 5
      treescale_x <- treescale_x - treescale_width
    }else if((max_tree_depth > 60) & (max_tree_depth <= 250)){
      treescale_width <- 20
      treescale_x <- treescale_x - treescale_width
    }else if(max_tree_depth > 250){
      treescale_width <- 100
      treescale_x <- treescale_x - treescale_width
    }
    ## set ylim adjustment to allow plotting of legend to the pdf 
    ylim_adj <- 5 ## set default ylim adjustment
    if (length(tree$tip.label) < 60){
      ylim_adj <- (60+2*ylim_adj-length(tree$tip.label))/2
    }
    
    tree_plt <- ggtree(rerooted_tree_wo_ref, size = 0.8) +
      geom_tiplab(align = TRUE, linetype = "dotted", size = 0, linesize = 0.2) +
      geom_rootedge(max_tree_depth/50, size = 0.8) + 
      geom_treescale(x = treescale_x, y = -4, width = treescale_width) + ## add tree scale
      xlim(-max_tree_depth/25, max_tree_depth*1.5) + ## expand on right of plot
      ylim(-ylim_adj, length(tree$tip.label)+ylim_adj) + 
      coord_cartesian(clip = "off") + ## do not clip stuff outside of plotting area
      ggtitle(paste(patient, species))
      
    ####################
    ## annotate tree 
    ####################
    tree_plt_data <- tree_plt$data
    
    ## get nodes to mark via a star
    if (length(tips_to_mark_mut) > 0){
      for (mut in names(tips_to_mark_mut)){
        ## subset list of isolates just to those which are present in tree (required to avoid errors); Note: If wrong tree loaded, then some mutations might not be annotated!!
        isolates_on_node = tips_to_mark_mut[[mut]][tips_to_mark_mut[[mut]] %in% rerooted_tree_wo_ref$tip.label]
        mark_mut_node <- MRCA(rerooted_tree_wo_ref, isolates_on_node)
        # Find midpoint of branch given the node
        derived <- tree_plt_data[tree_plt_data$node == mark_mut_node, ]
        ancestral <- tree_plt_data[tree_plt_data$node == derived$parent, ]
        mid_x <- mean(c(ancestral$x, derived$x))
        ## add linebreak to label
        mut_w_linebreak <- gsub('\\[', '\n\\[', mut)
        ## mark mutations
        tree_plt <- tree_plt +
          annotate('text', x = mid_x, y = min(derived$y)+0.5, label = mut_w_linebreak, lineheight = 0.8, vjust = 0)
      }
    }
    
    ## generate heatmap dataframes 
    level_order <- union(mandatory_sites_for_legend, unique(na.omit(tree_annotation_heatmap_df$location)))
    loc_df <- data.frame(Location = tree_annotation_heatmap_df$location) %>%
      `rownames<-`(tree_annotation_heatmap_df$sampleid)
    
    date_df <- data.frame(Days = tree_annotation_heatmap_df$days_since_obs) %>%
      `rownames<-`(tree_annotation_heatmap_df$sampleid)
    
    if (!all(is.na(tree_annotation_heatmap_df$ST))){
      st_df <- data.frame(ST = tree_annotation_heatmap_df$ST) %>%
        `rownames<-`(tree_annotation_heatmap_df$sampleid) %>%
        mutate(ST = ifelse(is.na(ST), 'NA', ST),
               ST = factor(ST, levels = names(color_map)))
      st_variant_df <- tree_annotation_heatmap_df[(tree_annotation_heatmap_df$ST_variant == TRUE) & !is.na(tree_annotation_heatmap_df$ST_variant) , c('sampleid', 'ST_variant')]
      st_variant_df <- left_join(st_variant_df, select(tree_plt$data, c('label', 'x', 'y')), by = c('sampleid' = 'label'))
    }else(st_df <- data.frame())
    
    ## add location heatmap 
    tree_plt <- gheatmap(tree_plt, loc_df, 
                         offset = plot_offset, width = width,
                         color = 'black',
                         colnames_angle=90, colnames_position = "top", 
                         colnames_offset_y = colnames_offset, hjust = 0) +
      scale_fill_manual(values = location_colors, breaks = level_order, na.value="white", name="Location", drop = FALSE) 
    current_width_offset <- plot_offset + width_offset
      
    ## add new scale and days since observation heatmap 
    tree_plt <- tree_plt + new_scale_fill()
    tree_plt <- gheatmap(tree_plt, date_df, 
                         offset=current_width_offset, width=width,
                         color = 'black',
                         colnames_angle=90, colnames_position = "top", 
                         colnames_offset_y = colnames_offset, hjust = 0) +
      scale_fill_gradient(low = "white", high = "black", na.value="white",
                          limits = c(0, max_species_date-min_species_date), 
                          name="Days since\nobservation")
    current_width_offset <- current_width_offset + width_offset
    
    if (any(!is.na(st_df$ST))) {
      ## add new scale and ST heatmap 
      tree_plt <- tree_plt + new_scale_fill() +
        xlim(0, max_tree_depth + current_width_offset + width_offset*1.5)
      tree_plt <- gheatmap(tree_plt, st_df, 
                           offset=current_width_offset, width=width,
                           color = 'black',
                           colnames_angle=90, colnames_position = "top", 
                           colnames_offset_y = colnames_offset, hjust = 0) +
        scale_fill_manual(values = color_map, name="Sequence type") +
        geom_point(data = st_variant_df, aes(y = y), x = max_tree_depth+plot_offset+width_offset*2+width*0.94, shape = 21, fill = "black", size = 2)
      current_width_offset <- current_width_offset + width_offset
    }
    
    if ((ncol(indel_table_t) > 0)){
      num_indels <- ncol(indel_table_t)
      tree_plt <- tree_plt + new_scale_fill() +
        xlim(0, max_tree_depth + current_width_offset + width_offset/3*num_indels + width_offset*1.5)
      tree_plt <- gheatmap(tree_plt, indel_table_t, 
                           offset=current_width_offset + ((max_tree_depth*width)/3), width=width/3*num_indels, ## NOTE: The alignment is done arbitrarily as the offset changes somehow when using multicolumn heatmaps!
                           color = 'black',
                           colnames_angle=90, colnames_position = "top", 
                           colnames_offset_y = colnames_offset, hjust = 0) +
        scale_fill_manual(values = c('TRUE' = '#615952', 'FALSE' = 'white'), labels = c("TRUE" = "present", "FALSE" = "absent"), name="Indel") ## dark blue hexcolor: #34344C
      current_width_offset <- current_width_offset + width_offset/3*num_indels
    }
    
    if ((ncol(mge_table_t) > 0)){
      num_mges <- ncol(mge_table_t)
      tree_plt <- tree_plt + new_scale_fill() +
        xlim(0, max_tree_depth + current_width_offset + width_offset/3*num_mges + width_offset*1.5)
      tree_plt <- gheatmap(tree_plt, mge_table_t, 
                           offset=current_width_offset + (max_tree_depth*width/3) - (width_offset - max_tree_depth*width)/2, width=width/3*num_mges, ## NOTE: The alignment is done arbitrarily as the offset changes somehow when using multicolumn heatmaps!
                           color = 'black',
                           colnames_angle=90, colnames_position = "top", 
                           colnames_offset_y = colnames_offset, hjust = 0) +
        scale_fill_gradient(low = "white", high = "#5F758B", na.value="white",
                            limits = c(0, 1), 
                            name="MGE coverage\nbreadth")
      current_width_offset <- current_width_offset + width_offset/3*num_mges
    }
    
    ## save plot
    plot_width <- (8 / (max_tree_depth + plot_offset)) * (max_tree_depth + current_width_offset) ## account for new width of plot after adding metadata; default width = 8
    plot_height <- (number_tips+2*ylim_adj)/10 ## account for length of tree (number tip labels + y limit offset added during tree generation (ylim()))
    
    pdf(paste0(outpath, "/", patient, '_', species_label, '_', analysis_type, '_' , tree_type, '.pdf'), 
        width = plot_width, height = plot_height, paper = "special", onefile = F)
    print(tree_plt)
    dev.off()
  }
}

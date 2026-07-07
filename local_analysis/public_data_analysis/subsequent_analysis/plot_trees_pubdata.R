## Load libraries
library(ape) ## to modify tree file
library(tidyverse)
library(Biostrings)
library(stringr)
library(ggplot2)
library(ggtree)
library(ggstar) ## to plot stars as shape
library(ggnewscale)  # for multiple fill scales in ggplot
library(RColorBrewer)
# library(treeio)

## Functions
get_descendant_tips <- function(tree, node) {
  # Find all descendants
  get_desc <- function(n) {
    children <- tree$edge[tree$edge[, 1] == n, 2]
    if (length(children) == 0) return(integer(0))
    c(children, unlist(lapply(children, get_desc)))
  }
  
  descendants <- get_desc(node)
  
  # Return tip labels only
  tip_ids <- descendants[descendants <= length(tree$tip.label)]
  tree$tip.label[tip_ids]
}


## Okabe ito extended by 9th color and with grey
country_colors = list('USA'=	'#d55e00',
                     'Canada'=	'#da711e',
                     'Mexico'=	'#df843c',
                     'Guadeloupe'=	'#e49659',
                     'Brazil'=	'#e9a977',
                     'Argentina'=	'#eebc95',
                     'Colombia'=	'#f2cfb2',
                     'France'=	'#009e73',
                     'Spain'=	'#1ea983',
                     'Lithuania'=	'#3cb594',
                     'Croatia'=	'#59c0a4',
                     'Slovakia'=	'#77cbb4',
                     'Greece'=	'#95d7c5',
                     'Poland'=	'#b2e2d5',
                     'China'=	'#0072b2',
                     'Japan'=	'#2486bd',
                     'Pakistan'=	'#4799c8',
                     'Russia'=	'#6badd2',
                     'Philippines'=	'#8fc1dd',
                     'Singapore'=	'#b2d5e8',
                     'Australia'=	'#faf7c6',
                     'South Africa'=	'#cc79a7',
                     'Tanzania'=	'#f0d7e5' 
) 

## colors based on Clustal X amino acid colors and adjusted from there
nt_colors = list('A' = '#9ecf7d',
                 'T' = '#d8666c',
                 'C' = '#7393cb',
                 'G' = '#f0bf4f',
                 'N' = 'grey')
## colors extended via glasbey from the tol_muted palette
study_colors = c('#332288','#88ccee','#44aa99','#117733','#999933','#ddcc77','#cc6677','#882255','#aa4499','#dddddd','#183d2d','#5d79d7','#108271','#7d2d00','#ba9696','#d26dba','#35457d','#145561','#9e654d','#494da6','#8a9eff','#5df7ff','#006541','#8a4586','#dfa2c6','#bec6a6','#044100','#8ab669','#20bedf','#dba224','#aa5d00','#41285d','#9a7586','#be8261','#492841','#106d6d','#9e6db6','#ffd74d','#92716d','#e3dfae','#9a4118','#7d7992','#4d5586','#3d399e','#0071ae','#554520','#ce96ca','#8e8639','#65ceb2','#187959')
nonsyn_syn_color = list('Nonsynonymous' = '#454545',
                        'Synonymous' = '#cfcfcf')


## NOTE:
# One need to let at first run refbased plottings through to get all the first observational dates of that species and then the denovo dates!
outpath <- '/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/figures/phylogenetic_trees'

analysis_path <- "/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/public_data/pathogendetection_2025_12_respiratory/case_Pub_data_resp_P07_Ehormaechei-c1_240825/2026_05_28_P07_Ehormaechei-c1_240825_gubbinsdefault_ambigALL/"
metadata_path <- "/Users/ad_lemur/Nextcloud/keylab/projects/mf_2020_hap/labbook/data_analysis/2023/public_data/pathogendetection_2025_12_respiratory/metadata/"
species_observation_range <- c()

patient_species <- gsub('.*/', '', analysis_path)
patient <- 'Public data'
species <- 'Ehormaechei-c1_240825'

## Load tree
tree_files <- list.files(path = paste0(analysis_path, '/raxml_bootstrapped_tree1000/'), pattern = ".*_ML_GTRG_snp_msa_all_1000bs_with_bootstraps.raxml.supportTBE$", full.names = TRUE)
tree <- read.tree(tree_files[length(tree_files)])

## ensure tree is rooted
outgroup_tip <- grep("^Equasihorm", tree$tip.label, value = TRUE)
rerooted_tree <- root(tree, outgroup = outgroup_tip, resolve.root = TRUE)
rerooted_tree <- rotateConstr(rerooted_tree, sort(rerooted_tree$tip.label, decreasing = TRUE)) ## sort internal branches by sampleid

## remove reference from tree
tips_to_remove <- grep("_ref$", rerooted_tree$tip.label, value = TRUE)
tips_to_remove <- append(tips_to_remove, outgroup_tip)
# rerooted_tree_wo_ref <- rerooted_tree
rerooted_tree_wo_ref <- drop.tip(rerooted_tree, tips_to_remove)

## Load metadata
tree_annotation_heatmap_df <- read.csv(paste0(metadata_path, 'metadata_samples.csv'))
rename_dict <- c(
  "Sample" = "tip_label" 
)
names(tree_annotation_heatmap_df) <- ifelse(
  names(tree_annotation_heatmap_df) %in% names(rename_dict),
  rename_dict[names(tree_annotation_heatmap_df)],
  names(tree_annotation_heatmap_df)
)

## Load mappings of bioprojects to study names
bioproject_to_study <- read.csv(paste0(metadata_path, 'HAI_pubdata_bioprojects.csv'), skip = 1)
bioproject_to_study$study_name <- paste(bioproject_to_study$First_Author, bioproject_to_study$Year)
bioproject_to_study$study_name <- ifelse(bioproject_to_study$study_name == ' NA', bioproject_to_study$BioProject, bioproject_to_study$study_name)

tree_annotation_heatmap_df <- left_join(tree_annotation_heatmap_df, bioproject_to_study %>% select('BioProject', 'study_name'), by = 'BioProject', relationship = 'many-to-many')


## set everything for outgroup to NA
rows_to_na <- grepl("^Equasihorm", tree_annotation_heatmap_df$tip_label)

## get nodes to highlight 
tips_to_highlight_branch <- c()
if (length(tips_to_highlight_branch) > 0){
  highlight_clade <- MRCA(rerooted_tree_wo_ref, tips_to_highlight_branch)
}

## get nodes to annotate
snp_table_file <- paste0(analysis_path, "/snp_table.csv")
snp_table <- read.csv(snp_table_file, header=TRUE, row.names = NULL, sep = ',')

annotation_mutation_file <- list.files(path = analysis_path, pattern = ".*annotation_mutations_dedup_repolarized.csv", full.names = TRUE) 
annotation_mutation <- read.csv(annotation_mutation_file[length(annotation_mutation_file)], header=TRUE, row.names = NULL, sep = ',')
annotation_mutation <- annotation_mutation[, c('chr', 'pos', 'nt_pos', 'AA_gt_ref', 'AA_gt_anc', 'AA_gt_alt')]
snp_table <- left_join(snp_table, annotation_mutation, by = c('chr', 'pos', 'nt_pos'))

## get ancestral seq from treetime reconstruction for proper rooting
anc_reconstruction_fasta <- paste0(analysis_path, "/node_sequence_inference/ancestral_sequences.fasta")
anc_seq <- readDNAStringSet(anc_reconstruction_fasta)$NODE_0000001
anc_seq <- unlist(strsplit(as.character(anc_seq), split = ""))
snp_table$nt_anc <- anc_seq

snp_table_fimZ <- snp_table[(snp_table$gene == 'fimZ'), ]
rownames(snp_table_fimZ) <- NULL
metadata_colnames <- c(colnames(snp_table_fimZ)[1:15], c('AA_gt_ref', 'AA_gt_anc', 'AA_gt_alt'))
isolate_colnames <- colnames(snp_table_fimZ)[!(colnames(snp_table_fimZ) %in% metadata_colnames)]

mut_df <- snp_table_fimZ
mut_df[isolate_colnames] <- lapply(mut_df[isolate_colnames], function(col){ifelse((col == mut_df$nt_anc), '', col)})
mut_df[isolate_colnames] <- lapply(mut_df[isolate_colnames], function(col){ifelse((col == '?'), 'N', col)})
rownames(mut_df) <- ifelse(mut_df$AA_gt_anc != mut_df$AA_gt_alt,
                          paste0(mut_df$AA_gt_anc, mut_df$aa_pos, mut_df$AA_gt_alt),
                          paste0(mut_df$AA_gt_anc, mut_df$aa_pos))
mut_df <- mut_df %>% 
  select(all_of(isolate_colnames)) %>%
  t() 

colnames(mut_df)[colnames(mut_df) == "P149PP"] <- "P149" ## replace multiallelic site of a syn mutation

####################
## plot initial tree 
####################
max_tree_depth <- max(node.depth.edgelength(rerooted_tree_wo_ref))
number_tips <- length(rerooted_tree_wo_ref$tip.label)
## heatmap specific variables
colnames_offset <- number_tips/250
width = 0.1
width_offset = max_tree_depth*(width*1.1)
plot_offset = 0
open_angle = 20
tree_rotation = 90+open_angle ## degree how much the tree is rotated
tree_scale_loc = round((360-tree_rotation)*number_tips/(360-open_angle), 0)

## plot tree
tree_plt <- ggtree(rerooted_tree_wo_ref, linewidth = 0.8, layout="fan", open.angle = open_angle) +
  geom_tiplab(align = TRUE, linetype = "dotted", size = 0, linesize = 0.2) +
  geom_rootedge(max_tree_depth/50, size = 0.8) + 
  geom_treescale(x = max_tree_depth*1.65, y = tree_scale_loc) + ## add tree scale
  geom_point2(aes(subset = !isTip, fill = as.numeric(label)), shape = 21, color = 'black', size = 2) +
  scale_fill_gradient(low = "white", high = "black", name = "Transfer bootstrap\nexpectation value") + ## add bootstrap values
  ggtitle(paste(patient, species))

  
####################
## annotate tree 
####################
tree_plt_data <- tree_plt$data



## generate heatmap dataframes 
country_df <- data.frame('Country' = tree_annotation_heatmap_df$country) %>%
  `rownames<-`(tree_annotation_heatmap_df$tip_label)

study_df <- data.frame('Study' = tree_annotation_heatmap_df$study_name) %>%
  `rownames<-`(tree_annotation_heatmap_df$tip_label)


tree_plt <- rotate_tree(tree_plt, angle = tree_rotation)
  
## add country heatmap 
tree_plt <- tree_plt + new_scale_fill()
tree_plt <- gheatmap(tree_plt, country_df, 
                     offset=plot_offset, width=width,
                     color = 'black', 
                     colnames_angle=0, colnames_position = "top", 
                     colnames_offset_y = colnames_offset, hjust = 1) +
  scale_fill_manual(values = scales::alpha(country_colors, 0.9), breaks = names(country_colors), name="Country", drop = FALSE) 

## add study heatmap 
tree_plt <- tree_plt + new_scale_fill()
tree_plt <- gheatmap(tree_plt, study_df, 
                     offset = plot_offset+width_offset, width = width,
                     color = 'black',
                     colnames_angle=0, colnames_position = "top", 
                     colnames_offset_y = colnames_offset, hjust = 1) +
  scale_fill_manual(values = scales::alpha(study_colors, 0.75), name="BioProject ID", drop = FALSE)

## add mutation heatmap 
tree_plt <- tree_plt + new_scale_fill() +
  xlim(0, max_tree_depth + plot_offset+width_offset + width_offset/3*length(colnames(mut_df)) + width_offset*2.3)
tree_plt <- gheatmap(tree_plt, mut_df, 
                     offset = plot_offset+width_offset*2.3, width = width/3*length(colnames(mut_df)),
                     color = 'black',
                     colnames_angle=0, colnames_position = "top", 
                     colnames_offset_y = colnames_offset, hjust = 1) +
  scale_fill_manual(values = nt_colors, breaks = names(nt_colors), name="fimZ SNVs", drop = FALSE, na.value = 'white')

print(tree_plt)

## save plot
pdf(paste0(outpath, "/", gsub(' ', '_', patient), '-', species, '_pathogendetection_resp.pdf'), 
    width = 18, height = 15, paper = "special", onefile = F)
print(tree_plt)
dev.off()
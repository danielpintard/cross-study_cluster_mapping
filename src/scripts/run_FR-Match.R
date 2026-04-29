# config - accept args, set wd and load libraries
print("Loading Libraries...")
library(FRmatch)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("At least one argument must be supplied", call. = False)
}

working_dir <- args[1] # should be src/scripts/
fr_match_dirs_path <- args [2] # base path to tmpdir containing FR-Match file dirs for each dataset
results_base_dir <- args[3]

#Ensure proper working dir
print(sprintf("Working out of: '%s'", working_dir))
setwd(working_dir)

dirs <- list.dirs(path=fr_match_dirs_path, full.names = TRUE) # list of dirs containing data files to create sce_objects
folder_names <- basename(dirs)
data_ids <- sub("_FRMatch_files$", "", folder_names)
dir_id_map <- setNames(dirs, data_ids)

# FUTURE:
# NOTE: consider creating conditional operations for whether or not fscore and cluster_order info will be included in cluster matching run
create_sce_object <- function(dir_path) {
    mtx <- fread(paste(dir_path, "matrix.csv"))
    cluster_info <- fread(paste(dir_path, "clusters.csv"))
    markers_info <- fread(paste(dir_path, "markers.csv"))
    fscores <- fread(paste(dir_path, "fscores.csv"))
    dendro_order <- fread(paste(dir_path, "dendrogram.csv"))
    unique_markers <- unique(markers_info$markers)
    
    # format + explode markers and associated cluster information
    NSF_formatted_mkrs <- markers_info %>% select(clusterName, markers) %>% mutate(markers = gsub("\\[|\\]|\\'| ", "", markers)) %>%
        separate_rows(markers, sep = ",") %>% as.data.frame()
    
    sce.obj <- make_data_object(
        dat = mtx,
        tab = cluster_info,
        markers = unique_markers,
        cluster_marker_info = NSF_formatted_mkrs,
        f_score = fscores,
        cluster_order = dendro_order
    )

    sce.obj.norm <- normalization(sce.obj)

    print(sprintf("Data read from %s, SCE object created...", dir_path))
    return(sce.obj.norm)

run_FRMatch <- function(sce.query.norm, sce.ref.norm, query_id, ref_id, res_base_dir) {

    ############################################ RUNNING FR TEST ############################################
    print(sprintf("Running FRMatch test in reference: %s -> query: %s mapping direction", ref_id, query_id))
    rst.query_to_ref <- FRmatch(sce.query = sce.query.norm, sce.ref = sce.ref.norm, subsamp.size = 15) 
    print(sprintf("Running FRMatch test in query: %s -> reference: %s mapping direction", query_id, ref_id))
    rst.ref_to_query <- FRmatch(sce.query = sce.ref.norm, sce.ref = sce.query.norm, subsamp.size = 15)

    ############################################ PLOTTING AND STORING RESULTS ############################################
    #### plot where E1/reference is x axis and E2/query data is y axis ####
    #plot one-way and two-way match results
    match_matrix <- plot_bi_FRmatch(
        rst.ref_to_query, rst.query_to_ref, 
        name.E1 = ref_id, name.E2 = query_id, 
        filename = sprintf("%s/%s_to_%s_mapping/%s_to_%s_bidir_plot_E1E2.png",res_base_dir, ref_id, query_id), return.value = TRUE
        ) 

    # plot two-way results only
    plot_bi_FRmatch(
        rst.ref_to_query, rst.query_to_ref, 
        name.E1 = ref_id, name.E2 = query_id, two.way.only = TRUE, 
        filename = sprintf("%s/%s_to_%s_mapping/%s_to_%s_bidir_plot(two-way)_E1E2.png", res_base_dir, ref_id, query_id)
        ) 
    
    match_df <- get_values(match_matrix)
    map_filename <- sprintf("%s/%s_to_%s_mapping/s%_to_%s_FRMatch_mapping_results.csv", res_base_dir, ref_id, query_id)
    write.csv(match_df, file = map_filename, row.names = FALSE)
    print(sprintf("Saved mapping results to %s", map_filename))

    #### plot where E1/query is x axis and E2/reference data is y axis ####
    #plot one-way and two-way match results
    plot_bi_FRmatch(
        rst.query_to_ref, rst.ref_to_query,
        name.E1 = query_id, name.E2 = ref_id, filename = sprintf("%s/%s_to_%s_mapping/%s_to_%s_matching_E2E1.png", res_base_dir, 
        ref_id, query_id, # since we want these plots to still get saved in the same results file
        query_id, ref_id)
        ) 
    # plot two-way results only
    plot_bi_FRmatch(
        rst.query_to_ref, rst.ref_to_query,
        name.E1 = query_id, name.E2 = ref_id, two.way.only = TRUE, 
        filename = sprintf("%s/%s_to_%s_mapping/%s_to_%s_matching_(two-way)_E2E1.png", res_base_dir, 
        ref_id, query_id, # since we want these plots to still get saved in the same results file
        query_id, ref_id)
        ) 
}

# using create_sce_object, loop through dirs and run FR-Match on non-duplicative, pairwise combinations of datasets being mapped onto each other
# NOTE: need to estimate amount of memory necessary for these to be stored
# alternatively, these could be written to and pulled from lscratch too
adatas <- list() # initialize empty array/vector ? \_("/)_/ idk how to in R
for (id in names(dir_id_map)) {
    dirpath <- dir_id_map[[id]] # access dirpath from named dir_id_map vector

    gc()

    adata <- create_sce_object(dirpath)
    adatas[[id]] <- adata # store ID: sce_obj in adatas vector
}

dataset_ids <- names(adatas) # get keys/names/ids from named vector
num_ds <- length(dataset_ids)

for (i in 1:(num_ds - 1)) {
    for (j in (i + 1):num_datasets) {
        id1 <- dataset_ids[i]
        id2 <- dataset_ids[j]

        sce1 <- adatas[[id1]]
        sce2 <- adatas[[id2]]
        print(sprintf("=== Running FR-Match on: %s and %s ===", id1, id2))
        runFRMatch(sce.query.norm = sce1, sce.ref.norm = sce2, query_id = id1, ref_id = id2)
    }
}

# #################################################################################################
# # scratch ground start #
# dir <- "/data/pintardde/class_marker_pipeline/results/HLCA_Core_results/tables/"
# markers_info <- fread(paste0(dir, "combined_markers_eval_results.csv"))

# sub("_[^_]+$", "", "HLCA_Core_results")
# # scratch ground end #
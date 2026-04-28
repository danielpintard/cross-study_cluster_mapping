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
fr_match_dirs_path <- args [2]

print(sprintf("Working out of: '%s'", working_dir))
setwd(working_dir)

dirs <- list.dirs(path=fr_match_dirs_path, full.names = TRUE) # list of dirs containing data files to create sce_objects

# scratch ground start #
dir <- "/data/pintardde/class_marker_pipeline/results/HLCA_Core_results/tables/"
markers_info <- fread(paste0(dir, "combined_markers_eval_results.csv"))

sub("_[^_]+$", "", "HLCA_Core_results")
# scratch ground end #

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
    return list("standard_sce"=sce.obj, "normalized_sce"=sce.obj.norm) # [ ] TODO: check if this is even possible to return more than one object in R functions

run_FRMatch <- function(sce.query, sce.ref, query_id, ref_id) {
    print(sprintf("Running FRMatch test in reference: %s -> query: %s mapping direction", ref_id, query_id))
    query_to_ref <- FRMatch(sce.query = sce.query, sce.ref = sce.ref, subsamp.size = 15) 
    print(sprintf("Running FRMatch test in query: %s -> reference: %s mapping direction", query_id, ref_id))
    ref_to_query <- FRMatch(sce.query = sce.ref, sce.ref = sce.query, subsamp.size = 15)
}

# using create_sce_object, loop through dirs and run FR-Match on non-duplicative, pairwise combinations of datasets being mapped onto each other
# adatas <- () # initialize empty list ? \_("/)_/ idk


# TODO: W.I.P. 
"""
I think the logic should switch to looping through dirs to create 
sce_objs to store and then we retreive those stored objs from mem to run FRMatch
"""
# for (dir1 in dirs) {
#     data_id1 = sub("_[^_]+$", "", dir1) # regex to ID last underscore and chars behind it and remove it
#     for (dir2 in dirs) {
#         if data_id1 == sub("_[^_]+$", "", dir2) {
#             # if the datasets are identical, skip
#             next
#         }
#         else {
#             # if they are not identical, continue
#             data_id2 = sub("_[^_]+$", "", dir2)

#         }
#     }
# } 

# define where data is stored and create of list dirs

# read in files and create sce objects when needed (calculate conversion from folder to file)


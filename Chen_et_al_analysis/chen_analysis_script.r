#
# use perviously published data from Chen et al. paper
# https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2020.00591/full

#===================================================
# analysis from the Chen et al using data downloaded from Adaptive
#===================================================
# install the replicateFest
#library(pak)
#pak::pak("OncologyQS/replicateFest")

library(replicateFest)
library(tools)


# reproducing analysis from the paper
# exclude 5 samples
excludeSamp = c("TSTLAEQVAW_2", "TSTLAEQVAW_3", "TSTLSEQVAW_3",
                "TSNLQEQIGW_2", "TSNLQEQIAW_2")

# specify path to tsv files downloaded from Adaptive
# https://clients.adaptivebiotech.com/pub/chan-2020-fi
inputDir = "."
# list paths to files with data
files = list.files(inputDir, full.names = T,
                   pattern = "tsv", recursive = TRUE)
# remove samples that were excluded from analysis according to the paper
files = setdiff(files, sapply(excludeSamp,grep,files, value = T ))

filenames = file_path_sans_ext(basename(files))
sampAnnot = splitFileName(filenames)

# specify cross-reactive conditions
xrCond = setdiff(sampAnnot$condition, c("AY9", "CEF","uncultured"))

# run all clones in a patient and time point and return the results
res = runExperiment(files,
                    peptides = sampAnnot$condition,
                    "NoPeptide",
                    fdrThr = 0.01,
                    orThr = 5,
                    nReads = 50,
                    percentThr = 0,
                    xrCond = xrCond,
                    excludeCond = "uncultured",
                    outputFile = "Chen_data-output_xr_FDR01_OR5.xlsx",
                    saveToFile = T)

# significant clones reported in the paper for Gag and Nef
clonesGag = c("CASSLDPGANTEAFF", "CASSPGVGNTEAFF", "CASSPRQAGLVTQYF")
clonesNef = c("CASSLDLRTFTYEQYF", "CASSLERVGYNEQFF",
              "CASSLLAGGSLDEQFF", "CASSPRWGDAGELFF",
              "CAWETGVRDGYTF", "CAISLMGTEAFF")

# intersections of published clones and cross-reactive clones
intersect(res$cross_reactive$clone,clonesGag)
intersect(res$cross_reactive$clone,clonesNef)
# intersection of publised clones and expanded clones
intersect(res$ref_comparison_only$clone,clonesGag)
intersect(res$ref_comparison_only$clone,clonesNef)
# cross-reactive clones that were not reported in the paper
setdiff(res$cross_reactive$clone, c(clonesGag,clonesNef))
# missed cross-reactive clones
setdiff(c(clonesGag,clonesNef), res$cross_reactive$clone)

# save results
saveRDS(res, file = "Chen_data-output_xr_FDR01_OR5.rds")


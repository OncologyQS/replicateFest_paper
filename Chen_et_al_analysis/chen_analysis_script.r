#
# use perviously published data from Kellie's paper
# https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2020.00591/full

#===================================================
# repeat analysis from the Chen et al using data downloaded from Adaptive
#===================================================

library(replicateFest)
library(tools)


# try to reproduce analysis
# exclude 5 samples
excludeSamp = c("TSTLAEQVAW_2", "TSTLAEQVAW_3", "TSTLSEQVAW_3",
                "TSNLQEQIGW_2", "TSNLQEQIAW_2")

inputDir = "C:/Users/ldanilo1/OneDrive - Johns Hopkins/JHU/Manafest/2025-10-09_test_data_for_package/from_Adaptive"
#inputDir = "C:/Users/Luda/OneDrive - Johns Hopkins/JHU/Manafest/2025-10-09_test_data_for_package/from_Adaptive"
# list paths to files with data
files = list.files(inputDir, full.names = T,
                   pattern = "tsv", recursive = TRUE)
# remove samples that were excluded from analysis according to paper
files = setdiff(files, sapply(excludeSamp,grep,files, value = T ))

filenames = file_path_sans_ext(basename(files))
sampAnnot = splitFileName(filenames)

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

#=========================
# create a matrix for Leslie to simulate
# run analysis with FDR < 0.05 to remove more expanded clones
res_05 = runExperiment(files,
                    peptides = sampAnnot$condition,
                    "NoPeptide",
                    fdrThr = 0.05,
                    orThr = 5,
                    nReads = 50,
                    percentThr = 0,
                    xrCond = xrCond,
                    excludeCond = "uncultured",
                    saveToFile = F)

# take clones to test that are not expanded for 10 peptides and no peptide
# get mergedData
dat = readMergeSave(files)
save(dat$mergedData,dat$ntData,file="Chen_inputData.rda")
mergedData = dat$mergedData
# get clones to test
clonesToTest = replicateFest::getClonesToTest(mergedData) # 371
length(clonesToTest)
# exclude expanded
clonesToSimulate = setdiff(clonesToTest,
                           res_05$ref_comparison_only$clone) #350

# restrict to samples that go to the simulations
samp = unique(sampAnnot$condition)[c(3:4,6:14)]
files = sampAnnot[(which(sampAnnot$condition %in% samp)),"file"]
clonesToTest = replicateFest::getClonesToTest(mergedData[files]) # 76
length(clonesToTest)
clonesToSimulate = setdiff(clonesToTest,
                           res_05$ref_comparison_only$clone) #63
length(clonesToSimulate)

# create matrix with counts for 10 peptide conditions and no peptide
# use only peptides with all 3 replicates
mat = getFreqOrCount(clonesToSimulate,
                     mergedData = mergedData,
                     samp = ,
                     colSuf = "",
                     returnFreq = F)
saveRDS(mat, file = "countMatrixForSim.rds")

#===================================

analysisRes = fitModelSet(clonesToTest,
                          mergedData,
                          sampAnnot$condition,
                          excludeCond = "uncultured",
                          refSamp = "NoPeptide",
                          c.corr=1)
rownames(analysisRes) = analysisRes$clone


tab = getExpanded(analysisRes,mergedData,
                      refSamp = "NoPeptide",
                      fdrThr = 0.01,
                      orThr = 5)
# second best results
screen_scndBest = analysisRes[,grepl("second",colnames(analysisRes))]

resTable = createResTableReplicates(analysisRes,
                                    mergedData,
                                    refSamp = "NoPeptide",
                                    fdrThr = 0.01,
                                    orThr = 5)

posCloneRep = getPositiveClonesReplicates(analysisRes,
                                          mergedData,
                                          "NoPeptide")
# There are no positive clones, because no clones met
# the significance threshold for the second best condition by FDR

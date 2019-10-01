# First get presence / absence data for the CNVs
cyp6.pa.calls <- read.table('../CNV_stats/tables_for_phase2_paper/cyp6aap.csv', header = T, row.names = 1)
colnames(cyp6.pa.calls) <- paste('Cyp6', colnames(cyp6.pa.calls), sep = '_')
cyp9k1.pa.calls <- read.table('../CNV_stats/tables_for_phase2_paper/cyp9k1.csv', header = T, row.names = 1)
colnames(cyp9k1.pa.calls) <- paste('Cyp9k1', colnames(cyp9k1.pa.calls), sep = '_')
# For the other two regions, we need to build the tables. 
load('/home/eric/Liverpool/CNV_v2/counting_output_v4_3R/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/Gste2_analysis_shrunk_data.Rdata')
# The +0 converts the logical values to numeric
gste.pa.calls <- (read.based.gst.duplications >= 1) + 0
load('/home/eric/Liverpool/CNV_v2/counting_output_v4_3R/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/CYP6M2-Z1_analysis_shrunk_data.Rdata')
cyp6mz.pa.calls <- (read.based.cyp.duplications >= 1) + 0
# Remove all objects that aren't those tables
rm(list = ls()[!grepl('pa.calls', ls())])

# Join the tables
pa.calls <- do.call(cbind, list(cyp6.pa.calls, cyp6mz.pa.calls, cyp9k1.pa.calls, gste.pa.calls))


# Next lets get a table of coverage calls. Where there is an NA, we replace the coverage calls by 1, because
# we know the dup is at least present
cyp6.coverage.calls <- read.table('/home/eric/Manuscripts/GSTE_new/Supplementary/cyp6_description_report_190201/Cyp6_coverage_calls.csv', header = T, sep = ',', row.names = 1)
colnames(cyp6.coverage.calls) <- paste('Cyp6', colnames(cyp6.coverage.calls), sep = '_')
cyp6.coverage.calls[is.na(cyp6.coverage.calls)] <- 1
#
cyp6mz.coverage.calls <- read.table('/home/eric/Manuscripts/GSTE_new/Supplementary/cyp6m2-z1_description_report_190201/Cyp6m2-z1_coverage_calls.csv', header = T, sep = ',', row.names = 1)
colnames(cyp6mz.coverage.calls) <- paste('Cyp6mz', colnames(cyp6mz.coverage.calls), sep = '_')
cyp6mz.coverage.calls[is.na(cyp6mz.coverage.calls)] <- 1
#
cyp9k1.coverage.calls <- read.table('/home/eric/Manuscripts/GSTE_new/Supplementary/cyp9k1_description_report_190201/Cyp9k1_coverage_calls.csv', header = T, sep = ',', row.names = 1)
colnames(cyp9k1.coverage.calls) <- paste('Cyp9k1', colnames(cyp9k1.coverage.calls), sep = '_')
cyp9k1.coverage.calls[is.na(cyp9k1.coverage.calls)] <- 1
#
gste.coverage.calls <- read.table('/home/eric/Manuscripts/GSTE_new/Supplementary/GSTE_description_report_190201/Gste_coverage_calls.csv', header = T, sep = ',', row.names = 1)
colnames(gste.coverage.calls) <- paste('Gst', colnames(gste.coverage.calls), sep = '_')
gste.coverage.calls[is.na(gste.coverage.calls)] <- 1

coverage.calls <- do.call(cbind, list(cyp6.coverage.calls, cyp6mz.coverage.calls, cyp9k1.coverage.calls, gste.coverage.calls))
rownames(coverage.calls) <- sub('\\\\', '', rownames(coverage.calls))


# Next just get the tables of hmm-based CNVs. 
load('../CNV_stats/CNV_stats.Rdata')
# Need to turn the list of samples carrying each CNV into a table
filtered.sample.names <- rownames(meta.reduced)
hmm.calls.list <- lapply(duplications.by.cluster.allchrom, function(x) filtered.sample.names %in% rownames(x))
hmm.calls <- do.call(cbind, hmm.calls.list) + 0
rownames(hmm.calls) <- filtered.sample.names

# Get the CNVs that passed population frequency filtering
hmm.nonsingle.calls <- hmm.calls[,rownames(subset(all.CNV.ranges.allchrom, singleton == F))]
hmm.filterpass.calls <- hmm.calls[,rownames(subset(all.CNV.ranges.allchrom, goodfreq))]

# Get the colors for the different populations:
colourscheme <- unlist(list('AOcol'= rgb(0.69439448188332953, 0.070034602810354785, 0.092318341048324815),
                            'BFcol'= rgb(0.98357554884517895, 0.4127950837799147, 0.28835064675293715),
                            'BFgam'= rgb(0.57960786223411564, 0.77019609212875362, 0.87372549772262575),
                            'CIcol'= rgb(0.98823529481887817, 0.62614381313323975, 0.50849674145380652),
                            'CMgam'= rgb(0.090196083486080159, 0.39294118285179136, 0.67058825492858887),
                            'FRgam'= rgb(0.47320263584454852, 0.43267974257469177, 0.69934642314910889),
                            'GAgam'= rgb(0.21568628648916882, 0.62875819206237793, 0.3333333432674408),
                            'GHcol'= rgb(0.89019608497619629, 0.18562091638644537, 0.15294117977221808),
                            'GHgam'= rgb(0.2909804046154022, 0.59450982809066777, 0.78901962041854856),
                            'GM'= rgb(0.939607846736908, 0.47137255668640132, 0.094901964068412781),
                            'GNcol'= rgb(0.99358708227381987, 0.83234141714432663, 0.76249136363758763),
                            'GNgam'= rgb(0.81411765813827519, 0.8839215755462646, 0.94980392456054685),
                            'GQgam'= rgb(0.7764706015586853, 0.77908498048782349, 0.88235294818878174),
                            'GW'= rgb(0.99607843160629272, 0.73490197658538814, 0.28000001013278963),
                            'KE'= rgb(0.58608230025160546, 0.58608230025160546, 0.58608230025160546),
                            'UGgam'= rgb(0.6810457706451416, 0.871895432472229, 0.65620917081832886)))

allcolours <- colourscheme[as.character(meta.reduced$population)]

pa.pca <- prcomp(pa.calls)
coverage.pca <- prcomp(coverage.calls)
hmm.pca <- prcomp(hmm.calls)
hmm.nonsingle.pca <- prcomp(hmm.nonsingle.calls)
hmm.filterpass.pca <- prcomp(hmm.filterpass.calls)

# Now let's try to do a pca on these. 
png('PCA_comparison.png', width = 720)
par(mfrow = c(2,3), mar = c(3,3,2.5,1), mgp = c(1.5,0.3,0), tcl = -0.3, cex = 0.9)
plot(pa.pca$x[,1], pa.pca$x[,2], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 1', ylab = 'PC 2', main = 'CNVs from discordant reads')
plot(coverage.pca$x[,1], coverage.pca$x[,2], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 1', ylab = 'PC 2', main = 'CNVs from discordant\nreads with copy-number')
plot(hmm.pca$x[,1], hmm.pca$x[,2], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 1', ylab = 'PC 2', main = 'HMM-based CNVs')
plot(hmm.nonsingle.pca$x[,1], hmm.nonsingle.pca$x[,2], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 1', ylab = 'PC 2', main = 'HMM-based CNVs non-singletons')
plot(hmm.filterpass.pca$x[,1], hmm.filterpass.pca$x[,2], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 1', ylab = 'PC 2', main = 'HMM-based CNVs > 5% freq')
plot(c(0,1), c(0,1), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
legend(-0.05, 1, names(colourscheme)[1:8], col = colourscheme[1:8], pch = 19, lty = 0, bty = 'n', cex = 1.5)
legend(0.5, 1, names(colourscheme[9:16]), col = colourscheme[9:16], pch = 19, lty = 0, bty = 'n', cex = 1.5)
dev.off()

# Since the HMM-based CNVs are the only ones that look reasonable, let's plot those with more PCs
png('PCA_detailed.png')
par(mfrow = c(2,2), mar = c(3,3,2.5,1), mgp = c(1.5,0.3,0), tcl = -0.3)
plot(hmm.nonsingle.pca$x[,1], hmm.nonsingle.pca$x[,2], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 1', ylab = 'PC 2', main = 'HMM-based CNVs non-singletons')
plot(hmm.nonsingle.pca$x[,3], hmm.nonsingle.pca$x[,4], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 3', ylab = 'PC 4', main = 'HMM-based CNVs non-singletons')
plot(hmm.nonsingle.pca$x[,5], hmm.nonsingle.pca$x[,6], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 5', ylab = 'PC 6', main = 'HMM-based CNVs non-singletons')
plot(c(0,1), c(0,1), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
legend(-0.05, 1, names(colourscheme)[1:8], col = colourscheme[1:8], pch = 19, lty = 0, bty = 'n', cex = 1.5)
legend(0.5, 1, names(colourscheme[9:16]), col = colourscheme[9:16], pch = 19, lty = 0, bty = 'n', cex = 1.5)
dev.off()

png('PCA_filterpass_detailed.png')
par(mfrow = c(2,2), mar = c(3,3,2.5,1), mgp = c(1.5,0.3,0), tcl = -0.3)
plot(hmm.filterpass.pca$x[,1], hmm.filterpass.pca$x[,2], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 1', ylab = 'PC 2', main = 'HMM-based CNVs > 5% freq')
plot(hmm.filterpass.pca$x[,3], hmm.filterpass.pca$x[,4], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 3', ylab = 'PC 4', main = 'HMM-based CNVs > 5% freq')
plot(hmm.filterpass.pca$x[,5], hmm.filterpass.pca$x[,6], col = allcolours, pch = 19, cex = 0.8, xlab = 'PC 5', ylab = 'PC 6', main = 'HMM-based CNVs > 5% freq')
plot(c(0,1), c(0,1), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
legend(-0.05, 1, names(colourscheme)[1:8], col = colourscheme[1:8], pch = 19, lty = 0, bty = 'n', cex = 1.5)
legend(0.5, 1, names(colourscheme[9:16]), col = colourscheme[9:16], pch = 19, lty = 0, bty = 'n', cex = 1.5)
dev.off()

# We have decided to go with the non-singleton HMM-based CNVs, with more PCs, and with a barplot of 
# proportion of variance explained
pca.variance <- hmm.nonsingle.pca$sdev^2
prop.pca.variance <- pca.variance / sum(pca.variance)

png('PCA_nonsingle_full.png', width = 720)
par(mfrow = c(2,3), mar = c(3,3,1,1), mgp = c(1.5,0.3,0), tcl = -0.3, cex = 0.9)
for (i in 1:4){
	plot(hmm.filterpass.pca$x[,i*2-1], hmm.filterpass.pca$x[,i*2], col = allcolours, pch = 19, cex = 0.8, xlab = paste('PC', i*2-1), ylab = paste('PC', i*2))
}
plot(c(0,1), c(0,1), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
legend(-0.05, 1, names(colourscheme)[1:8], col = colourscheme[1:8], pch = 19, lty = 0, bty = 'n', cex = 1.5)
legend(0.5, 1, names(colourscheme[9:16]), col = colourscheme[9:16], pch = 19, lty = 0, bty = 'n', cex = 1.5)
barplot(prop.pca.variance[1:10]*100, names.arg = 1:10, border = NA, ylab = 'Variance explained (%)', xlab = 'Principal component', cex.names = 0.94)
dev.off()

png('PCA_nonsingle_full_nolegend.png', width = 720)
par(mfrow = c(2,3), mar = c(3,3,1,1), mgp = c(1.5,0.3,0), tcl = -0.3, cex = 0.9)
for (i in 1:5){
	plot(hmm.filterpass.pca$x[,i*2-1], hmm.filterpass.pca$x[,i*2], col = allcolours, pch = 19, cex = 0.8, xlab = paste('PC', i*2-1), ylab = paste('PC', i*2))
}
barplot(prop.pca.variance[1:10]*100, names.arg = 1:10, border = NA, ylab = 'Variance explained (%)', xlab = 'Principal component', cex.names = 0.94)
dev.off()

# Now write the different HMM call matrices to file
write.table(hmm.calls, 'HMM_calls.csv', sep = '\t', col.names = NA)
write.table(hmm.pca$x, 'HMM_PCA.csv', sep = '\t', col.names = NA)
write.table(hmm.nonsingle.calls, 'HMM_nonsingle_calls.csv', sep = '\t', col.names = NA)
write.table(hmm.nonsingle.pca$x, 'HMM_nonsingle_PCA.csv', sep = '\t', col.names = NA)
write.table(hmm.filterpass.calls, 'HMM_filterpass_calls.csv', sep = '\t', col.names = NA)
write.table(hmm.filterpass.pca$x, 'HMM_filterpass_PCA.csv', sep = '\t', col.names = NA)

save.image('PCA.Rdata')

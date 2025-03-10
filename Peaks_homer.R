# Analyse scATACseq and get a list of DE peaks from condition 1 vs condition 2 in a specific cluster of interest

# Homer bg peak set:
  # This depends on the number of cells and DE peaks you have
  # Options can be all the open and closed peaks from you comparision
  # All peaks in the assay
  # A 'resonable' number of peaks randomly selected from your assay, re-iterate this by resampling and choose motifs that consistently come up

# Determing a resobale P value fro homer

# Extract 'open' and 'closed' peaks
neutr_open_mutant <- rownames(markers.neutr[markers.neutr$avg_logFC < -0.2, ])

neutr_closed_mutant <- rownames(markers.neutr[markers.neutr$avg_logFC > 0.2, ])

# Combine together

all_peaks_DE <- rbind(neutr_open_mutant, neutr_closed_mutant)


# Randomly shuffle the open and closed peaks over several iterations

all_peaks_DE$ran1 <- as.data.frame(sample(all_peaks_DE$group))[,1]
all_peaks_DE$ran2 <- as.data.frame(sample(all_peaks_DE$group))[,1]
all_peaks_DE$ran3 <- as.data.frame(sample(all_peaks_DE$group))[,1]
all_peaks_DE$ran4 <- as.data.frame(sample(all_peaks_DE$group))[,1]
all_peaks_DE$ran5 <- as.data.frame(sample(all_peaks_DE$group))[,1]
all_peaks_DE$ran6 <- as.data.frame(sample(all_peaks_DE$group))[,1]


all_peaks_DE$chr <- matrix(unlist(strsplit(row.names(all_peaks_DE),"-")),ncol=3,byrow=T)[,1]
all_peaks_DE$start <- matrix(unlist(strsplit(row.names(all_peaks_DE),"-")),ncol=3,byrow=T)[,2]
all_peaks_DE$end <- matrix(unlist(strsplit(row.names(all_peaks_DE),"-")),ncol=3,byrow=T)[,3]

# Now save file ready for input to homer

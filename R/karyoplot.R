# karyoplot_v1.R

plot_karyo_vcf <- function(tumor_vcf_list, chr, chr_lengths){

  library(karyoploteR)

  # Set up the plot with three data panels
  kp <- plotKaryotype(plot.type = 4, genome = "hg19", chromosomes = chr, main = paste("All Variants Signals in Simulation"))

  # Define track heights and offsets
  tr.i <- 1/4
  tr.o <- 1/3

  # Define labels based on plot_type
  labels = c("DP", "ALT_AD", "BAF")

  # Loop through three SNP types: Paternal Hetero, Maternal Hetero, and Homo
  for (tn in 0:2) {

    # Set plotting value and max value
    vcf_df <- tumor_vcf_list[[chr]]
    y_val <- vcf_df[,labels[tn+1]]
    max_val <- max(y_val, na.rm = TRUE)
    # Set data panel background and axis for each SNP type
    kpDataBackground(kp, data.panel = 1, r0 = tr.o*tn, r1 = tr.o*tn + tr.i)
    kpAxis(kp, ymin = 0, ymax = max_val, r0=tr.o*tn, r1=tr.o*tn + tr.i - 0.03, col="gray50", cex=0.5)
    kpText(kp, chr=chr, x = round(chr_lengths[chr]/2), y=0.5, col="black",
           r0=tr.o*tn + tr.i+0.01, r1=tr.o*tn + tr.i+0.06, labels=labels[tn+1], cex=1)

    # Plot points for the combined SNP data based on plot_type
    kpPoints(kp, chr = vcf_df$CHROM, x = as.integer(vcf_df$POS), y = y_val, col = "black",
             r0 = tr.o*tn, r1 = tr.o*tn + tr.i - 0.03, cex = 0.3, ymin = 0,
             ymax = max_val, data.panel = 1)
  }
}

plot_karyo_snp <- function(segment_snp_list, tumor_vcf_list, chr, chr_lengths, plot_type = "BAF") {
  library(karyoploteR)

  # Validate plot_type
  if (!plot_type %in% c("BAF", "ALT_AD")) {
    stop("Invalid plot_type. Choose either 'BAF' or 'ALT_AD'.")
  }

  # Set up the plot with three data panels
  kp <- plotKaryotype(plot.type = 4, genome = "hg19", chromosomes = chr, main = paste("Ground Truth SNP", plot_type, "in Simulation"))

  # Define track heights and offsets
  tr.i <- 1/4
  tr.o <- 1/3

  # Define labels based on plot_type
  base_name <- c(" of Paternal Hetero SNP", " of Maternal Hetero SNP", " of Homo SNP")
  labels <- c(paste0(plot_type, base_name))

  # Loop through three SNP types: Paternal Hetero, Maternal Hetero, and Homo
  snp_types <- c("paternal_hetero", "maternal_hetero", "homo")

  for (tn in 0:2) {

    # Combine SNP data with VCF data based on the SNP type
    if (snp_types[tn+1] == "paternal_hetero") {
      combined_df <- inner_join(segment_snp_list[[chr]]$paternal$hetero, tumor_vcf_list[[chr]], by = c("POS", "REF", "ALT"))
    } else if (snp_types[tn+1] == "maternal_hetero") {
      combined_df <- inner_join(segment_snp_list[[chr]]$maternal$hetero, tumor_vcf_list[[chr]], by = c("POS", "REF", "ALT"))
    } else { # homo
      combined_df <- inner_join(segment_snp_list[[chr]]$maternal$homo, tumor_vcf_list[[chr]], by = c("POS", "REF", "ALT"))
    }

    # Set plotting value and max value
    y_val <- combined_df[, plot_type]
    max_val <- max(y_val, na.rm = TRUE)

    # Set data panel background and axis for each SNP type
    kpDataBackground(kp, data.panel = 1, r0 = tr.o*tn, r1 = tr.o*tn + tr.i)
    kpAxis(kp, ymin = 0, ymax = max_val, r0=tr.o*tn, r1=tr.o*tn + tr.i - 0.03, col="gray50", cex=0.5)
    kpText(kp, chr=chr, x = round(chr_lengths[chr]/2), y=0.5, col="black",
           r0=tr.o*tn + tr.i+0.01, r1=tr.o*tn + tr.i+0.06, labels=labels[tn+1], cex=1)

    # Plot points for the combined SNP data based on plot_type
    kpPoints(kp, chr = combined_df$CHROM.y, x = combined_df$POS, y = y_val, col = "black",
             r0 = tr.o*tn, r1 = tr.o*tn + tr.i - 0.03, cex = 0.3, ymin = 0,
             ymax = max_val, data.panel = 1)
  }
}


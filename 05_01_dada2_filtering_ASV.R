
library(dada2)
library(data.table)
library(dplyr)
library(ShortRead)  # for streaming / sampling fastq

# -------------------------
# User-configurable params
# -------------------------
metadata_file <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/stats_ad_qc_tr_final_meta_final.csv"
fastq_dir_16s <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_16s_trimmed_dada2"
fastq_dir_its <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_its_trimmed_dada2"
filtered_dir_16s <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_16s_dada2_filtered"
filtered_dir_its <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/fastq_pr_its_dada2_filtered"
results_dir <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/dada2_results"
summary_file <- "/media/jwchoi/ssd2/projects/microbiome/16S_1/metadata/dada2_ASV_results.csv"

dir.create(filtered_dir_16s, recursive=TRUE, showWarnings=FALSE)
dir.create(filtered_dir_its, recursive=TRUE, showWarnings=FALSE)
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

# Pipeline behavior params
BATCH_SIZE <- 50                # number of samples processed per batch (tune for memory)
MAX_READS_PER_SAMPLE <- 500000  # if a filtered sample has more reads than this, subsample to this cap
POOL_MODE <- "pseudo"           # "FALSE", "pseudo", or "TRUE" for dada pooling
NBases_LEARN <- 1e8             # nbases for learnErrors (increase for better model if memory allows)
MAXEE <- c(1,1)                 # filterAndTrim maxEE (for PE; SE uses first element)
TRUNCQ <- 5
MIN_ABUNDANCE_REL <- 1e-6       # relative abundance cutoff (fraction of total reads)
MIN_ABUNDANCE_ABS <- 5          # absolute reads cutoff
PREVALENCE_MIN <- 2             # minimum number of samples ASV must appear in
VERBOSE <- TRUE

# -------------------------
# Helper functions
# -------------------------
logmsg <- function(...) {
  if (VERBOSE) cat(sprintf(...), "\n")
}

# Reservoir subsampling of a FASTQ file using ShortRead streaming
# - infile: input fastq (can be gz)
# - outfile: path to write sampled fastq (gz)
# - cap: desired number of reads to retain
# This performs random reservoir sampling without loading entire file into memory.
subsample_fastq_reservoir <- function(infile, outfile, cap) {
  # Return original infile path if cap is NA or negative (no subsampling)
  if (is.null(cap) || cap <= 0) return(infile)
  rs <- FastqStreamer(infile, n=1e5)  # chunk size
  reservoir <- list()   # store ShortReadQ objects as list of cycles
  total_seen <- 0L
  set.seed(12345)  # reproducible sampling
  
  while(length(fq <- yield(rs))) {
    # fq is ShortReadQ object
    chunk_n <- length(fq)
    for (i in seq_len(chunk_n)) {
      total_seen <- total_seen + 1L
      if (length(reservoir) < cap) {
        reservoir[[length(reservoir) + 1L]] <- fq[i]
      } else {
        # replace with probability cap/total_seen
        j <- sample.int(total_seen, 1)
        if (j <= cap) reservoir[[j]] <- fq[i]
      }
    }
  }
  close(rs)
  
  # combine reservoir into ShortReadQ and write out gzipped fastq
  if (length(reservoir) == 0) {
    stop("No reads found during subsampling.")
  }
  combined <- do.call(c, reservoir)
  writeFastq(combined, outfile, compress=TRUE)
  return(outfile)
}

# Utility: safe list.files for matching patterns
safe_list_files <- function(dir, pattern) {
  if (!dir.exists(dir)) return(character(0))
  list.files(dir, pattern=pattern, full.names=TRUE, recursive=TRUE)
}

# Get fastq files for sample (improved matching)
get_fastq_files <- function(sample_id, layout, base_dir) {
  search_dir <- file.path(base_dir, tolower(layout))
  if (!dir.exists(search_dir)) {
    warning("Search dir missing: ", search_dir)
    return(NULL)
  }
  
  if (layout == "PE") {
    fpat <- paste0(sample_id, ".*_1|", sample_id, ".*_R1|", sample_id, "_1\\.")
    rpat <- paste0(sample_id, ".*_2|", sample_id, ".*_R2|", sample_id, "_2\\.")
    ffiles <- safe_list_files(search_dir, fpat)
    rfiles <- safe_list_files(search_dir, rpat)
    if (length(ffiles) >= 1 && length(rfiles) >= 1) {
      # choose first pair that seems to match (user can refine if needed)
      return(list(forward=ffiles[1], reverse=rfiles[1]))
    }
  } else {
    spat <- paste0(sample_id, ".*_s|", sample_id, ".*_SE|", sample_id, ".*\\.fastq")
    sfiles <- safe_list_files(search_dir, spat)
    if (length(sfiles) >= 1) return(list(single=sfiles[1]))
  }
  return(NULL)
}

# -------------------------
# Read metadata
# -------------------------
logmsg("Reading metadata...")
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
metadata$SRR <- trimws(metadata$SRR)
metadata$BioProject <- trimws(metadata$BioProject)
metadata$Seq <- trimws(metadata$Seq)
metadata$Layout <- trimws(metadata$Layout)

# optional trunc lengths in metadata
if (!("Trunc1" %in% names(metadata))) metadata$Trunc1 <- NA
if (!("Trunc2" %in% names(metadata))) metadata$Trunc2 <- NA

summary_df <- data.frame(
  BioProject=character(), Seq=character(), Layout=character(),
  Sample_no=integer(), Total_read_no=integer(), ASV_no=integer(),
  stringsAsFactors = FALSE
)

# -------------------------
# processing per bioproject / seq type
# -------------------------
process_bioproject <- function(bioproject, seq_type) {
  logmsg("Processing: %s %s", bioproject, seq_type)
  bp_metadata <- metadata[metadata$BioProject == bioproject & metadata$Seq == seq_type, ]
  if (nrow(bp_metadata) == 0) {
    logmsg("No samples for %s %s", bioproject, seq_type)
    return()
  }
  
  layout <- unique(bp_metadata$Layout)[1]
  if (seq_type == "16S") {
    base_dir <- fastq_dir_16s; filtered_base_dir <- filtered_dir_16s
  } else {
    base_dir <- fastq_dir_its; filtered_base_dir <- filtered_dir_its
  }
  filtered_project_dir <- file.path(filtered_base_dir, tolower(layout), bioproject)
  dir.create(filtered_project_dir, recursive=TRUE, showWarnings=FALSE)
  
  # collect sample files and truncation parameters
  sample_list <- list()
  for (i in 1:nrow(bp_metadata)) {
    row <- bp_metadata[i, ]
    sid <- row$SRR
    files <- get_fastq_files(sid, layout, base_dir)
    if (is.null(files)) {
      warning("Files not found for sample: ", sid)
      next
    }
    sample_list[[sid]] <- list(files=files,
                               trunc1 = ifelse(is.na(row$Trunc1), NA, as.numeric(row$Trunc1)),
                               trunc2 = ifelse(is.na(row$Trunc2), NA, as.numeric(row$Trunc2)))
  }
  if (length(sample_list) == 0) {
    logmsg("No valid files for %s %s", bioproject, seq_type); return()
  }
  logmsg("Found %d samples", length(sample_list))
  
  # Filtering step (filterAndTrim) - create filtered files and optionally subsample very deep ones
  filtered_files <- list()
  for (sid in names(sample_list)) {
    info <- sample_list[[sid]]
    if (layout == "PE") {
      fnF <- info$files$forward; fnR <- info$files$reverse
      outF <- file.path(filtered_project_dir, paste0(sid, "_1_filt.fastq.gz"))
      outR <- file.path(filtered_project_dir, paste0(sid, "_2_filt.fastq.gz"))
      truncF <- ifelse(is.na(info$trunc1), NULL, info$trunc1)
      truncR <- ifelse(is.na(info$trunc2), NULL, info$trunc2)
      truncVec <- if (is.null(truncF) || is.null(truncR)) NULL else c(truncF, truncR)
      # call filterAndTrim with truncLen if present
      if (is.null(truncVec)) {
        out <- tryCatch({
          filterAndTrim(fnF, outF, fnR, outR, maxN=0, maxEE=MAXEE, truncQ=TRUNCQ,
                        rm.phix=TRUE, compress=TRUE, multithread=FALSE)
        }, error=function(e) {
          warning("filterAndTrim failed for ", sid, ": ", e$message)
          return(NULL)
        })
      } else {
        out <- tryCatch({
          filterAndTrim(fnF, outF, fnR, outR, truncLen = truncVec, maxN=0, maxEE=MAXEE, truncQ=TRUNCQ,
                        rm.phix=TRUE, compress=TRUE, multithread=FALSE)
        }, error=function(e) { warning("filterAndTrim failed for ", sid, ": ", e$message); return(NULL) })
      }
      if (!is.null(out) && out[1,"reads.out"] > 0) {
        # optional: check reads.out and subsample if too many
        reads_out <- out[1,"reads.out"]
        if (!is.na(reads_out) && reads_out > MAX_READS_PER_SAMPLE) {
          logmsg("Sample %s has %d reads after filtering -> subsampling to %d", sid, reads_out, MAX_READS_PER_SAMPLE)
          # subsample each read file
          subsampledF <- file.path(filtered_project_dir, paste0(sid, "_1_filt_sub.fastq.gz"))
          subsampledR <- file.path(filtered_project_dir, paste0(sid, "_2_filt_sub.fastq.gz"))
          subsampledF <- subsample_fastq_reservoir(outF, subsampledF, MAX_READS_PER_SAMPLE)
          subsampledR <- subsample_fastq_reservoir(outR, subsampledR, MAX_READS_PER_SAMPLE)
          filtered_files[[sid]] <- list(forward=subsampledF, reverse=subsampledR, reads_out=MAX_READS_PER_SAMPLE)
        } else {
          filtered_files[[sid]] <- list(forward=outF, reverse=outR, reads_out=reads_out)
        }
      } else {
        warning("No reads passed filtering for ", sid)
      }
    } else {
      fn <- info$files$single
      outS <- file.path(filtered_project_dir, paste0(sid, "_s_filt.fastq.gz"))
      truncF <- ifelse(is.na(info$trunc1), NULL, info$trunc1)
      if (is.null(truncF)) {
        out <- tryCatch({
          filterAndTrim(fn, outS, maxN=0, maxEE=MAXEE[1], truncQ=TRUNCQ, rm.phix=TRUE, compress=TRUE, multithread=FALSE)
        }, error=function(e) { warning("filterAndTrim failed for ", sid, ": ", e$message); return(NULL) })
      } else {
        out <- tryCatch({
          filterAndTrim(fn, outS, truncLen=truncF, maxN=0, maxEE=MAXEE[1], truncQ=TRUNCQ, rm.phix=TRUE, compress=TRUE, multithread=FALSE)
        }, error=function(e) { warning("filterAndTrim failed for ", sid, ": ", e$message); return(NULL) })
      }
      if (!is.null(out) && out[1,"reads.out"] > 0) {
        reads_out <- out[1,"reads.out"]
        if (!is.na(reads_out) && reads_out > MAX_READS_PER_SAMPLE) {
          logmsg("Sample %s has %d reads after filtering -> subsampling to %d", sid, reads_out, MAX_READS_PER_SAMPLE)
          subsampled <- file.path(filtered_project_dir, paste0(sid, "_s_filt_sub.fastq.gz"))
          subsampled <- subsample_fastq_reservoir(outS, subsampled, MAX_READS_PER_SAMPLE)
          filtered_files[[sid]] <- list(single=subsampled, reads_out=MAX_READS_PER_SAMPLE)
        } else {
          filtered_files[[sid]] <- list(single=outS, reads_out=reads_out)
        }
      } else {
        warning("No reads passed filtering for ", sid)
      }
    }
  } # end for each sample
  
  if (length(filtered_files) == 0) {
    logmsg("No samples passed filtering for %s %s", bioproject, seq_type); return()
  }
  
  logmsg("Samples passed filtering: %d", length(filtered_files))
  
  # --- Learn error rates ---
  logmsg("Learning error rates using nbases = %g ...", NBases_LEARN)
  # build a vector of files for learning; use up to first BATCH_SIZE or all if small
  if (layout == "PE") {
    learnFs <- sapply(filtered_files, function(x) x$forward)
    learnRs <- sapply(filtered_files, function(x) x$reverse)
    # use learnErrors with nbases parameter
    errF <- learnErrors(learnFs, nbases=NBases_LEARN, multithread=FALSE)
    errR <- learnErrors(learnRs, nbases=NBases_LEARN, multithread=FALSE)
  } else {
    learnFs <- sapply(filtered_files, function(x) x$single)
    errF <- learnErrors(learnFs, nbases=NBases_LEARN, multithread=FALSE)
    errR <- NULL
  }
  
  # process in batches to limit memory; create cumulative seqtab
  sample_ids <- names(filtered_files)
  cum_seqtab <- NULL
  total_reads_all <- 0
  
  for (batch_start in seq(1, length(sample_ids), by=BATCH_SIZE)) {
    batch_ids <- sample_ids[batch_start : min(batch_start + BATCH_SIZE - 1, length(sample_ids))]
    logmsg("Processing batch %d - %d (%d samples)", batch_start, batch_start+length(batch_ids)-1, length(batch_ids))
    
    # per-batch containers
    derepFs <- vector("list", length(batch_ids))
    dadaFs <- vector("list", length(batch_ids))
    names(derepFs) <- names(dadaFs) <- batch_ids
    if (layout == "PE") {
      derepRs <- vector("list", length(batch_ids))
      dadaRs <- vector("list", length(batch_ids))
      names(derepRs) <- names(dadaRs) <- batch_ids
    }
    
    # derep + dada per sample
    for (sid in batch_ids) {
      tryCatch({
        if (layout == "PE") {
          ffile <- filtered_files[[sid]]$forward
          rfile <- filtered_files[[sid]]$reverse
          drpF <- derepFastq(ffile); drpF$uniques <- drpF$uniques  # ensure object
          drpR <- derepFastq(rfile); drpR$uniques <- drpR$uniques
          derepFs[[sid]] <- drpF
          derepRs[[sid]] <- drpR
          dadaFs[[sid]] <- dada(drpF, err=errF, multithread=FALSE, pool=POOL_MODE)
          dadaRs[[sid]] <- dada(drpR, err=errR, multithread=FALSE, pool=POOL_MODE)
        } else {
          sfile <- filtered_files[[sid]]$single
          drp <- derepFastq(sfile)
          derepFs[[sid]] <- drp
          dadaFs[[sid]] <- dada(drp, err=errF, multithread=FALSE, pool=POOL_MODE)
        }
      }, error=function(e) {
        warning("DADA failed for sample ", sid, ": ", e$message)
      })
    } # end per-sample loop
    
    # mergePairs and makeSequenceTable for batch
    if (layout == "PE") {
      # prepare derep lists for mergePairs (must be in same order as dada lists)
      drpF_list <- derepFs[batch_ids]
      drpR_list <- derepRs[batch_ids]
      dfs <- dadaFs[batch_ids]; drs <- dadaRs[batch_ids]
      mergers <- mergePairs(dfs, drpF_list, drs, drpR_list, verbose=TRUE)
      batch_seqtab <- makeSequenceTable(mergers)
    } else {
      batch_seqtab <- makeSequenceTable(dadaFs[batch_ids])
    }
    
    # accumulate total reads
    total_reads_all <- total_reads_all + sum(batch_seqtab)
    logmsg("Batch seqtab dims: %d samples x %d seqs ; batch reads: %d", nrow(batch_seqtab), ncol(batch_seqtab), sum(batch_seqtab))
    
    # merge to cumulative seqtab (if exists)
    if (is.null(cum_seqtab)) {
      cum_seqtab <- batch_seqtab
    } else {
      cum_seqtab <- mergeSequenceTables(cum_seqtab, batch_seqtab)
    }
    
    # clean up batch objects
    rm(derepFs, derepRs, dadaFs, dadaRs, batch_seqtab, mergers)
    gc()
  } # end batches
  
  logmsg("Total reads across all batches (pre-chimera): %d", total_reads_all)
  # Remove chimeras on cumulative seqtab
  logmsg("Removing chimeras using method='consensus' ...")
  seqtab.nochim <- removeBimeraDenovo(cum_seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
  
  # Aggressive but reasonable abundance filtering
  asv_counts <- colSums(seqtab.nochim)
  total_reads <- sum(seqtab.nochim)
  min_reads_relative <- max(MIN_ABUNDANCE_ABS, round(total_reads * MIN_ABUNDANCE_REL))
  min_abundance <- max(min_reads_relative, MIN_ABUNDANCE_ABS)
  logmsg("Applying abundance filter: min_abundance = %d (total reads = %d)", min_abundance, total_reads)
  seqtab.filtered <- seqtab.nochim[, asv_counts >= min_abundance, drop=FALSE]
  
  # Prevalence filtering
  asv_prevalence <- colSums(seqtab.filtered > 0)
  seqtab.final <- seqtab.filtered[, asv_prevalence >= PREVALENCE_MIN, drop=FALSE]
  
  # Save results for this project
  result_prefix <- paste(tolower(seq_type), layout, bioproject, sep = "_")
  asv_df <- as.data.frame(seqtab.final)
  asv_df$Sample <- rownames(asv_df)
  asv_df <- asv_df[, c("Sample", setdiff(names(asv_df), "Sample"))]
  csv_file <- file.path(results_dir, paste0(result_prefix, ".csv"))
  write.csv(asv_df, csv_file, row.names = FALSE)
  rds_file <- file.path(results_dir, paste0(result_prefix, ".rds"))
  saveRDS(seqtab.final, rds_file)
  
  # update summary
  new_row <- data.frame(BioProject=bioproject, Seq=seq_type, Layout=layout,
                        Sample_no=nrow(seqtab.final), Total_read_no=total_reads,
                        ASV_no=ncol(seqtab.final), stringsAsFactors = FALSE)
  summary_df <<- rbind(summary_df, new_row)
  write.csv(summary_df, summary_file, row.names = FALSE)
  
  logmsg("Completed %s %s: %d samples, %d ASVs", bioproject, seq_type, nrow(seqtab.final), ncol(seqtab.final))
  # final cleanup
  rm(cum_seqtab, seqtab.nochim, seqtab.filtered, seqtab.final)
  gc()
}

# -------------------------
# Main loop
# -------------------------
logmsg("Starting DADA2 pipeline (memory-aware) ...")
unique_combos <- unique(metadata[, c("BioProject", "Seq")])
for (i in seq_len(nrow(unique_combos))) {
  bioproject <- unique_combos$BioProject[i]
  seq_type <- unique_combos$Seq[i]
  tryCatch({
    process_bioproject(bioproject, seq_type)
  }, error=function(e) {
    warning("Processing failed for ", bioproject, " ", seq_type, ": ", e$message)
  })
}

logmsg("DADA2 run finished. Summary file saved: %s", summary_file)
print(summary_df)

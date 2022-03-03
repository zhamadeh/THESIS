###################################################################################
##################     Workflow Part 1: Collect initial metrics   #################
###################################################################################

suppressMessages(suppressPackageStartupMessages(library(BiocManager)))
suppressMessages(suppressPackageStartupMessages(BiocManager::install("bamsignals")))
suppressMessages(suppressPackageStartupMessages(BiocManager::install("GenomicAlignments")))
suppressMessages(suppressPackageStartupMessages(library(IRanges)))
suppressMessages(suppressPackageStartupMessages(library(bamsignals)))
suppressMessages(suppressPackageStartupMessages(library(GenomeInfoDb)))
suppressMessages(suppressPackageStartupMessages(library(GenomicRanges)))
suppressMessages(suppressPackageStartupMessages(library(GenomicAlignments)))
suppressMessages(suppressPackageStartupMessages(library(Rsamtools)))



qc.spikiness <- function(counts) {
	if (is.null(counts)) {
		return(NA)
	}
	counts <- as.vector(counts)
	sum.counts <- sum(counts)
	spikiness <- sum(abs(diff(counts))) / sum.counts
	return(spikiness)
}

bamToGRanges <- function(bamfile, bamindex=bamfile,chromosomes=NULL,pairedEndReads=FALSE,remove.duplicate.reads=FALSE,min.mapq=10,max.fragment.width=1000,blacklist=NULL,what='mapq') {

	## Input checks
	if (!is.null(blacklist)) {
		if ( !(is.character(blacklist) | class(blacklist)=='GRanges') ) {
			stop("'blacklist' has to be either a bed(.gz) file or a GRanges object")
		}
	}

	## Check if bamindex exists
	bamindex.raw <- sub('\\.bai$', '', bamindex)
	bamindex <- paste0(bamindex.raw,'.bai')
	if (!file.exists(bamindex)) {
		ptm <- startTimedMessage("Making bam-index file ...")
		bamindex.own <- Rsamtools::indexBam(bamfile)
		warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
		bamindex <- bamindex.own
		stopTimedMessage(ptm)
	}
	chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if non of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('The specified chromosomes ', chrstring, ' do not exist in the data. Pay attention to the naming convention in your data, e.g. "chr1" or "1".')
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
	}

	## Import the file into GRanges
	ptm <- startTimedMessage("Reading file ",basename(bamfile)," ...")
	gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
	if (!remove.duplicate.reads) {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter = min.mapq))
		} else {
			data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter = min.mapq))
		}
	} else {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter = min.mapq, flag=scanBamFlag(isDuplicate=FALSE)))
		} else {
			data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(which=range(gr), what=what, mapqFilter = min.mapq, flag=scanBamFlag(isDuplicate=FALSE)))
		}
	}
	stopTimedMessage(ptm)

	if (length(data.raw) == 0) {
		if (pairedEndReads) {
			stop(paste0("No reads imported. Does your file really contain paired end reads? Try with 'pairedEndReads=FALSE'"))
		}
		stop(paste0('No reads imported! Check your BAM-file ', bamfile))
	}

	## Convert to GRanges and filter
	if (pairedEndReads) {
		ptm <- startTimedMessage("Converting to GRanges ...")
		data <- GenomicAlignments::granges(data.raw, use.mcols = TRUE, on.discordant.seqnames='drop') # treat as one fragment
		stopTimedMessage(ptm)

		ptm <- startTimedMessage("Filtering reads ...")
		# if (!is.na(min.mapq)) {
		# 	mapq.first <- mcols(GenomicAlignments::first(data.raw))$mapq
		# 	mapq.last <- mcols(GenomicAlignments::last(data.raw))$mapq
		# 	mapq.mask <- mapq.first >= min.mapq & mapq.last >= min.mapq
		# 	if (any(is.na(mapq.mask))) {
		# 		warning(paste0(bamfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NA' to keep all reads."))
		# 	}
		# 	data <- data[which(mapq.mask)]
		# }
		# Filter out too long fragments
		data <- data[width(data)<=max.fragment.width]
		stopTimedMessage(ptm)
	} else {
		ptm <- startTimedMessage("Converting to GRanges ...")
		data <- GenomicAlignments::granges(data.raw, use.mcols = TRUE) # treat as one fragment
		stopTimedMessage(ptm)

		ptm <- startTimedMessage("Filtering reads ...")
		# if (!is.na(min.mapq)) {
		# 	if (any(is.na(mcols(data)$mapq))) {
		# 		warning(paste0(bamfile,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NA' to keep all reads."))
		# 		mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
		# 	}
		# 	data <- data[mcols(data)$mapq >= min.mapq]
		# }
		# Filter out too long fragments
		data <- data[width(data)<=max.fragment.width]
		stopTimedMessage(ptm)
	}

	## Exclude reads falling into blacklisted regions
	if (!is.null(blacklist)) {
		ptm <- startTimedMessage("Filtering blacklisted regions ...")
		if (is.character(blacklist)) {
			if (grepl('^chr', seqlevels(data)[1])) {
				chromosome.format <- 'UCSC'
			} else {
				chromosome.format <- 'NCBI'
			}
			black <- importBed(blacklist, skip=0, chromosome.format=chromosome.format)
		} else if (class(blacklist)=='GRanges') {
			black <- blacklist
		} else {
			stop("'blacklist' has to be either a bed(.gz) file or a GRanges object")
		}
		overlaps <- findOverlaps(data, black)
		idx <- setdiff(1:length(data), S4Vectors::queryHits(overlaps))
		data <- data[idx]
		stopTimedMessage(ptm)
	}
	return(data)
}

fixedWidthBins <- function(bamfile=NULL, assembly=NULL, chrom.lengths=NULL, chromosome.format, binsizes=1e6, stepsizes=NULL, chromosomes=NULL) {

	### Check user input ###
	if (length(binsizes) == 0) {
		return(list())
	}
	if (is.null(bamfile) & is.null(assembly) & is.null(chrom.lengths)) {
		stop("Please specify either a 'bamfile', 'assembly' or 'chrom.lengths'")
	}
	if (is.null(bamfile) & is.null(chrom.lengths)) {
		trigger.error <- chromosome.format
	}
	if (!is.null(stepsizes)) {
		if (length(stepsizes) != length(binsizes)) {
			stop("Need one element in 'stepsizes' for each element in 'binsizes'.")
		}
		if (any(binsizes < stepsizes)) {
			stop("'stepsizes' must be smaller/equal than 'binsizes'")
		}
	}

	### Get chromosome lengths ###
	if (!is.null(bamfile)) {
		ptm <- startTimedMessage(paste0("Reading header from ", bamfile, " ..."))
		chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
		stopTimedMessage(ptm)
	} else if (!is.null(assembly)) {
		if (is.character(assembly)) {
			ptm <- startTimedMessage("Fetching chromosome lengths from UCSC ...")
			df <- GenomeInfoDb::getChromInfoFromUCSC(assembly)
			stopTimedMessage(ptm)
		} else if (is.data.frame(assembly)) {
			df <- assembly
		} else {
			stop("Unknown assembly")
		}
		chrom.lengths <- df$size
		if (chromosome.format=='UCSC') {
		} else if (chromosome.format=='NCBI') {
			df$chrom = sub('^chr', '', df$chrom)
		}
		names(chrom.lengths) <- df$chrom
		chrom.lengths <- chrom.lengths[!is.na(names(chrom.lengths))]
		chrom.lengths <- chrom.lengths[!is.na(chrom.lengths)]
	} else if (!is.null(chrom.lengths)) {
		chrom.lengths <- chrom.lengths[!is.na(names(chrom.lengths))]
		chrom.lengths <- chrom.lengths[!is.na(chrom.lengths)]
	}
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	## Stop if none of the specified chromosomes exist
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('Could not find length information for any of the specified chromosomes: ', chrstring, '. Pay attention to the naming convention in your data, e.g. "chr1" or "1".')
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning('Could not find length information for the following chromosomes: ', diffs)
	}

	### Making fixed-width bins ###
	bins.list <- list()
	for (ibinsize in 1:length(binsizes)) {
		binsize <- binsizes[ibinsize]
		ptm <- startTimedMessage("Making fixed-width bins for bin size ", binsize, " ...")
		chrom.lengths.floor <- floor(chrom.lengths / binsize) * binsize
		bins <- unlist(GenomicRanges::tileGenome(chrom.lengths.floor[chroms2use], tilewidth=binsize), use.names=FALSE)
		bins <- bins[end(bins) > 0] # no chromosomes that are smaller than binsize
		if (any(width(bins)!=binsize)) {
			stop("tileGenome failed")
		}
		# seqlengths(bins) <- as.integer(chrom.lengths[names(seqlengths(bins))])
		seqlengths(bins) <- chrom.lengths[chroms2use]
		if (!is.null(stepsizes)) {
			shift.bp <- 0
			stepsize <- stepsizes[ibinsize]
			bins.list.step <- GRangesList()
			while (shift.bp < binsize) {
				bins.list.step[[as.character(shift.bp)]] <- suppressWarnings( trim(shift(bins, shift.bp)) )
				shift.bp <- stepsize + shift.bp
			}
			bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE), '_stepsize_', format(stepsize, scientific=TRUE, trim=TRUE))]] <- bins.list.step
		} else {
			bins.list[[paste0('binsize_', format(binsize, scientific=TRUE, trim=TRUE))]] <- bins
		}

		skipped.chroms <- setdiff(seqlevels(bins), as.character(unique(seqnames(bins))))
		if (length(skipped.chroms)>0) {
			warning("The following chromosomes were skipped because they are smaller than binsize ", binsize, ": ", paste0(skipped.chroms, collapse=', '))
		}
		stopTimedMessage(ptm)

	}

	return(bins.list)

}

startTimedMessage <- function(...) {

	x <- paste0(..., collapse='')
	message(x, appendLF=FALSE)
	ptm <- proc.time()
	return(ptm)

}

stopTimedMessage <- function(ptm) {

	time <- proc.time() - ptm
	message(" ", round(time[3],2), "s")

}

transCoord <- function(gr) {
  cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
  cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
  names(cum.seqlengths.0) <- GenomeInfoDb::seqlevels(gr)
  gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  return(gr)
}

collectLibraryStats <- function(folder){
	met=data.frame(file=c(),coverage=c(),spikiness=c(),evenness.mean=c(),evenness.med=c())
	#file=list.files(folder,full.names = T,pattern="\\.bam$")[1]
	for (file in list.files(folder,full.names = T,pattern="\\.bam$")){
		bamindex=file
		bamindex.raw <- sub('\\.bai$', '', bamindex)
		bamindex <- paste0(bamindex.raw,'.bai')
		if (!file.exists(bamindex)) {
			ptm <- startTimedMessage("Making bam-index file ...")
			bamindex.own <- Rsamtools::indexBam(file)
			warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
			bamindex <- bamindex.own
			stopTimedMessage(ptm)
		}
		chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(file))
		chrom.lengths.df <- data.frame(chromosome=names(chrom.lengths), length=chrom.lengths)
		bins <- fixedWidthBins(chrom.lengths=chrom.lengths,  binsizes=1e+06, stepsizes=1e+06)
		counts.ss <- bamsignals::bamCount(file, bins[[1]][[1]], mapqual=10, paired.end="filter", tlenFilter=c(0,1000), verbose=FALSE, ss=TRUE, filteredFlag=1024)
		pcounts <- counts.ss[1,]
		mcounts <- counts.ss[2,]
		counts <- pcounts + mcounts
		readsperbin <- round(sum(as.numeric(counts)) / length(counts), 2)
		countmatrix <- matrix(c(counts,mcounts,pcounts), ncol=3)
		colnames(countmatrix) <- c('counts','mcounts','pcounts')
		mcols(bins[[1]][[1]]) <- as(countmatrix, 'DataFrame')

		data=bamToGRanges(file)

		genome.length <- sum(as.numeric(seqlengths(data)))
		data.strand <- data
		strand(data.strand) <- '*'
		coverage <- sum(as.numeric(width(data.strand))) / genome.length
		#genome.covered <- sum(as.numeric(width(reduce(data.strand)))) / genome.length

		binned.reads = bins[[1]][[1]]

		med.reads.per.mb=median(binned.reads$counts)
		mean.reads.per.mb=mean(binned.reads$counts)

		binned.reads$var.med = (abs(binned.reads$counts - med.reads.per.mb)) /
										sd(binned.reads$counts)

		binned.reads$var.mean = (abs(binned.reads$counts - mean.reads.per.mb)) /
										sd(binned.reads$counts)

		evenness.med = sum(binned.reads$var.med) / (3* length(binned.reads))
		evenness.mean = sum(binned.reads$var.mean) / (3* length(binned.reads))

		spikiness=qc.spikiness(bins[[1]][[1]]$counts)

		print(paste0("Coverage:",coverage))
		#print(paste0("Genome.covered:",genome.covered))
		print(paste0("Spikiness",spikiness))
		print(paste0("Evenness.med",evenness.med))
		print(paste0("Evenness.mean",evenness.mean))

		row=data.frame(file=basename(file),coverage=coverage,spikiness=spikiness,evenness.mean=evenness.mean,evenness.med=evenness.med)
		met=rbind(met,row)
	}
	
	write.table(met,"thesis.bam_feb23-2022-metrics.txt",sep="\t",quote=F,row.names = F,col.names = T)
}

collectLibraryStats("THESIS/FEB23_2022/GOOD/")

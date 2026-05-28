
# function to run the Battenberg pipeline
# heavily modified by Oriol Pich and Kerstin Thol to run on txWGs data
# and to run for multiple solutions, including multiple solutions from PURPLE and
# manually QC'd TRAXERx421 WES solutions

#' Run the Battenberg pipeline
#'
#' @param analysis The mode of Battenberg copy number analysis to be undertaken: 'paired' for tumour-normal pair, 'cell_line' for Cell line tumour-only and 'germline' for germline CNV of normal sample (Default: 'paired')
#' @param samplename Sample identifier (tumour or germline), this is used as a prefix for the output files. If allele counts are supplied separately, they are expected to have this identifier as prefix.
#' @param normalname Matched normal identifier, this is used as a prefix for the output files. If allele counts are supplied separately, they are expected to have this identifier as prefix.
#' @param sample_data_file A BAM or CEL file for the sample
#' @param normal_data_file A BAM or CEL file for the normal-pair (paired analysis)
#' @param imputeinfofile Full path to a Battenberg impute info file with pointers to Impute2 reference data
#' @param g1000prefix Full prefix path to 1000 Genomes SNP loci data, as part of the Battenberg reference data
#' @param problemloci Full path to a problem loci file that contains SNP loci that should be filtered out
#' @param gccorrectprefix Full prefix path to GC content files, as part of the Battenberg reference data, not required for SNP6 data (Default: NULL)
#' @param repliccorrectprefix Full prefix path to replication timing files, as part of the Battenberg reference data, not required for SNP6 data (Default: NULL)
#' @param g1000allelesprefix Full prefix path to 1000 Genomes SNP alleles data, as part of the Battenberg reference data, not required for SNP6 data (Default: NA)
#' @param ismale A boolean set to TRUE if the donor is male, set to FALSE if female, not required for SNP6 data (Default: NA)
#' @param data_type String that contains either wgs or snp6 depending on the supplied input data (Default: wgs)
#' @param impute_exe Pointer to the Impute2 executable (Default: impute2, i.e. expected in $PATH)
#' @param allelecounter_exe Pointer to the alleleCounter executable (Default: alleleCounter, i.e. expected in $PATH)
#' @param nthreads The number of concurrent processes to use while running the Battenberg pipeline (Default: 8)
#' @param platform_gamma Platform scaling factor, suggestions are set to 1 for wgs and to 0.55 for snp6 (Default: 1)
#' @param phasing_gamma Gamma parameter used when correcting phasing mistakes (Default: 1)
#' @param segmentation_gamma The gamma parameter controls the size of the penalty of starting a new segment during segmentation. It is therefore the key parameter for controlling the number of segments (Default: 10)
#' @param segmentation_gamma_multisample The gamma parameter controls the size of the penalty of starting a new segment during mutlisample segmentation. It is the key parameter for controlling the number of segments (Default: 10)
#' @param segmentation_kmin Kmin represents the minimum number of probes/SNPs that a segment should consist of (Default: 3)
#' @param phasing_kmin Kmin used when correcting for phasing mistakes (Default: 3)
#' @param clonality_dist_metric  Distance metric to use when choosing purity/ploidy combinations (Default: 0)
#' @param ascat_dist_metric Distance metric to use when choosing purity/ploidy combinations (Default: 1)
#' @param min_ploidy Minimum ploidy to be considered (Default: 1.6)
#' @param max_ploidy Maximum ploidy to be considered (Default: 4.8)
#' @param min_rho Minimum purity to be considered (Default: 0.1)
#' @param max_rho Maximum purity to be considered (Default: 1.0)
#' @param min_goodness Minimum goodness of fit required for a purity/ploidy combination to be accepted as a solution (Default: 0.63)
#' @param uninformative_BAF_threshold The threshold beyond which BAF becomes uninformative (Default: 0.51)
#' @param min_normal_depth Minimum depth required in the matched normal for a SNP to be considered as part of the wgs analysis (Default: 10)
#' @param min_base_qual Minimum base quality required for a read to be counted when allele counting (Default: 20)
#' @param min_map_qual Minimum mapping quality required for a read to be counted when allele counting (Default: 35)
#' @param max_allowed_state The maximum CN state allowed (Default 250)
#' @param cn_upper_limit Maximum number of copy number that can be called (Default 1000)
#' @param calc_seg_baf_option Sets way to calculate BAF per segment: 1=mean, 2=median, 3=ifelse median==0 | 1, mean, median (Default (paired): 3, cell_line & germline: 1)
#' @param skip_allele_counting Provide TRUE when allele counting can be skipped (i.e. its already done) (Default: FALSE)
#' @param skip_preprocessing Provide TRUE when preprocessing is already complete (Default: FALSE)
#' @param skip_phasing  Provide TRUE when phasing is already complete (Default: FALSE)
#' @param skip_segmentation Provide TRUE when segmentation and multisample rephasing are already complete, i.e. BAFsegmented.txt files exist (Default: FALSE)
#' @param usebeagle Should use beagle5 instead of impute2 Default: FALSE
#' @param beaglejar Full path to Beagle java jar file Default: NA
#' @param beagleref.template Full path template to Beagle reference files where the chromosome is replaced by 'CHROMNAME' Default: NA
#' @param beagleplink.template Full path template to Beagle plink files where the chromosome is replaced by 'CHROMNAME' Default: NA
#' @param beaglemaxmem Integer Beagle max heap size in Gb  Default: 10
#' @param beaglenthreads Integer number of threads used by beagle5 Default:1
#' @param beaglewindow Integer size of the genomic window for beagle5 (cM) Default:40
#' @param beagleoverlap Integer size of the overlap between windows beagle5 Default:4
#' @param javajre Path to the Java JRE executable, only required for haplotype reconstruction with Beagle (default java, i.e. in $PATH)
#' @param snp6_reference_info_file Reference files for the SNP6 pipeline only (Default: NA)
#' @param apt.probeset.genotype.exe Helper tool for extracting data from CEL files, SNP6 pipeline only (Default: apt-probeset-genotype)
#' @param apt.probeset.summarize.exe  Helper tool for extracting data from CEL files, SNP6 pipeline only (Default: apt-probeset-summarize)
#' @param norm.geno.clust.exe  Helper tool for extracting data from CEL files, SNP6 pipeline only (Default: normalize_affy_geno_cluster.pl)
#' @param birdseed_report_file Sex inference output file, SNP6 pipeline only (Default: birdseed.report.txt)
#' @param heterozygousFilter Legacy option to set a heterozygous SNP filter, SNP6 pipeline only (Default: "none")
#' @param prior_breakpoints_file A two column file with prior breakpoints to be used during segmentation (Default: NULL)
#' @param genomebuild Genome build upon which the 1000G SNP coordinates were obtained (Default: hg19; options: "hg19" or "hg38")  
#' @param externalhaplotypefile Vcf containing externally obtained haplotype blocks (Default: NA)
#' @param write_battenberg_phasing Write the Battenberg phasing results as vcf to disk, e.g. for multisample cases (Default: TRUE)
#' @param multisample_maxlag Maximal number of upstream SNPs used in the multisample haplotyping to inform the haplotype at another SNP (Default: 100)
#' @param multisample_relative_weight_balanced Relative weight to give to haplotype info from a sample without allelic imbalance in the region (Default: 0.25)
#' @param enhanced_grid_search Should use multi-start, parallelized and multi-approach grid search (Default: FALSE)
#' @param cn_confidence_level Main confidence level to use for the bootstrapped confidence intervals on nMajor and nMinor. Set to FALSE, NULL, NA or 0 to skip those intervals entirely (Default 0.95)
#' @param test_mode Enable test mode for faster execution on production data. Overrides noperms, impute_region_size, and enhanced_grid_search with fast defaults unless explicitly set (Default: FALSE)
#' @param noperms Number of bootstrap permutations for subclonal copy number confidence intervals (Default: 1000; test_mode default: 10)
#' @param impute_region_size Size of genomic windows in bp used by IMPUTE2 for phasing. Larger values mean fewer chunks and faster runtime at the cost of more RAM (Default: 5000000; test_mode default: 20000000)
#' @param test_chromosomes Character vector of chromosome names to restrict the analysis to, e.g. c("1", "10", "21"). NULL means use all chromosomes from the impute info file (Default: NULL)
#' @param organize_output Organize output files into subdirectories (results/, intermediate/, plots/, logs/) after processing. On restart with any skip flag, files are flattened back to the working directory first (Default: TRUE)
#' @author sd11, jdemeul, Naser Ansari-Pour, Julio Cesar Cortes Rios
#' @export
battenberg = function(analysis="paired",
                      samplename,
                      normalname,
                      sample_data_file,
                      normal_data_file,
                      imputeinfofile,
                      g1000prefix,
                      problemloci,
                      gccorrectprefix=NULL,
                      repliccorrectprefix=NULL,
                      g1000allelesprefix=NA,
                      ismale=NA,
                      data_type="wgs",
                      impute_exe="impute2",
                      allelecounter_exe="alleleCounter",
                      nthreads=8,
                      platform_gamma=1,
                      phasing_gamma=1,
                      segmentation_gamma=10,
                      segmentation_kmin=3,
                      phasing_kmin=1,
                      clonality_dist_metric=0,
                      ascat_dist_metric=1,
                      min_ploidy=1.6,
                      max_ploidy=4.8,
                      min_rho=0.1,
                      max_rho=1.0,
                      min_goodness=0.63,
                      cn_confidence_level=0.95,
                      uninformative_BAF_threshold=0.51,
                      min_normal_depth=10,
                      min_base_qual=20,
                      min_map_qual=35,
                      max_allowed_state=250,
                      cn_upper_limit=1000,
                      calc_seg_baf_option=3,
                      skip_allele_counting=F,
                      skip_preprocessing=F,
                      skip_phasing=F,
                      skip_segmentation=F,
                      externalhaplotypefile = NA,
                      usebeagle=FALSE,
                      beaglejar=NA,
                      beagleref.template=NA,
                      beagleplink.template=NA,
                      beaglemaxmem=10,
                      beaglenthreads=1,
                      beaglewindow=40,
                      beagleoverlap=4,
                      javajre="java",
                      write_battenberg_phasing = T,
                      multisample_relative_weight_balanced = 0.25,
                      multisample_maxlag = 90,
                      segmentation_gamma_multisample = 5,
                      snp6_reference_info_file=NA,
                      apt.probeset.genotype.exe="apt-probeset-genotype",
                      apt.probeset.summarize.exe="apt-probeset-summarize",
                      norm.geno.clust.exe="normalize_affy_geno_cluster.pl",
                      birdseed_report_file="birdseed.report.txt",
                      heterozygousFilter="none",
                      prior_breakpoints_file=NULL,
                      genomebuild="hg19",
                      chrom_coord_file=NULL,
		                  enhanced_grid_search = F,
                      purple_path=NULL,
                      WES_solutions=NULL,
                      test_mode=FALSE,
                      noperms=1000,
                      impute_region_size=5000000,
                      test_chromosomes=NULL,
                      organize_output=TRUE,
                      seed=as.integer(1),
                      debug_parallel=FALSE) {
  requireNamespace("parallel")
  requireNamespace("R.utils")
  libs <- .libPaths()

  # --- Test mode overrides ---
  if (test_mode) {
    message("WARNING: Running in TEST MODE — results are NOT suitable for production use. ",
            "Overriding noperms, impute_region_size, and enhanced_grid_search for speed.")
    # Only override parameters the user didn't explicitly supply
    if (missing(noperms)) noperms <- 10
    if (missing(impute_region_size)) impute_region_size <- 20000000
    if (missing(enhanced_grid_search)) enhanced_grid_search <- TRUE
    if (!is.null(test_chromosomes)) {
      print(paste0("Chromosome subsetting active: restricting to chromosomes ", paste(test_chromosomes, collapse=", ")))
    }
  } else {
    test_chromosomes <- NULL
  }

  # --- Flatten organized subdirectories for restart compatibility ---
  any_skip <- any(skip_preprocessing, skip_phasing, skip_segmentation, skip_allele_counting)
  if (organize_output && any_skip) {
    flatten_organized_output()
    unzip_all_files()
  }

  pipeline_start_time <- Sys.time()
  format_duration <- function(seconds) {
    if (is.null(seconds) || is.na(seconds)) {
      return("NA")
    }
    secs <- as.integer(round(seconds))
    hh <- secs %/% 3600
    mm <- (secs %% 3600) %/% 60
    ss <- secs %% 60
    sprintf("%02d:%02d:%02d", hh, mm, ss)
  }
  emit_stage_event <- function(stage, status, stage_start = NULL, extra_info = "") {
    now <- Sys.time()
    total_elapsed <- as.numeric(difftime(now, pipeline_start_time, units = "secs"))
    stage_elapsed <- if (is.null(stage_start)) NA_real_ else as.numeric(difftime(now, stage_start, units = "secs"))
    token <- paste0("BATTENBERG_STAGE_", toupper(status), ":", toupper(stage))
    msg <- paste0(
      "[", format(now, "%Y-%m-%d %H:%M:%S %Z"), "] ",
      token,
      " stage=", stage,
      " status=", status,
      " stage_elapsed=", format_duration(stage_elapsed),
      " total_elapsed=", format_duration(total_elapsed),
      if (nchar(extra_info) > 0) paste0(" ", extra_info) else ""
    )
    cat(msg, "\n")
  }

  emit_stage_event(
    stage = "pipeline",
    status = "start",
    stage_start = pipeline_start_time,
    extra_info = paste0("samples=", length(samplename), " analysis=", analysis, " data_type=", data_type)
  )

  # Write a run-parameters log at the very start so failed runs can be inspected
  run_log_file <- paste0(samplename[1], "_battenberg_run_params.log")
  run_log_lines <- c(
    paste0("# Battenberg run parameters"),
    paste0("# Written: ", format(pipeline_start_time, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("# Working directory: ", getwd()),
    paste0("# R version: ", R.version$version.string),
    paste0("# Platform: ", R.version$platform),
    "",
    paste0("analysis = ", analysis),
    paste0("samplename = ", paste(samplename, collapse=", ")),
    paste0("normalname = ", paste(normalname, collapse=", ")),
    paste0("sample_data_file = ", paste(sample_data_file, collapse=", ")),
    paste0("normal_data_file = ", paste(normal_data_file, collapse=", ")),
    paste0("imputeinfofile = ", imputeinfofile),
    paste0("g1000prefix = ", g1000prefix),
    paste0("g1000allelesprefix = ", paste(g1000allelesprefix, collapse=", ")),
    paste0("problemloci = ", problemloci),
    paste0("gccorrectprefix = ", if (is.null(gccorrectprefix)) "NULL" else gccorrectprefix),
    paste0("repliccorrectprefix = ", if (is.null(repliccorrectprefix)) "NULL" else repliccorrectprefix),
    paste0("ismale = ", ismale),
    paste0("data_type = ", data_type),
    paste0("genomebuild = ", genomebuild),
    paste0("impute_exe = ", impute_exe),
    paste0("allelecounter_exe = ", allelecounter_exe),
    paste0("nthreads = ", nthreads),
    paste0("seed = ", seed),
    paste0("platform_gamma = ", platform_gamma),
    paste0("phasing_gamma = ", phasing_gamma),
    paste0("segmentation_gamma = ", segmentation_gamma),
    paste0("segmentation_kmin = ", segmentation_kmin),
    paste0("phasing_kmin = ", phasing_kmin),
    paste0("clonality_dist_metric = ", clonality_dist_metric),
    paste0("ascat_dist_metric = ", ascat_dist_metric),
    paste0("min_ploidy = ", min_ploidy),
    paste0("max_ploidy = ", max_ploidy),
    paste0("min_rho = ", min_rho),
    paste0("max_rho = ", max_rho),
    paste0("min_goodness = ", min_goodness),
    paste0("cn_confidence_level = ", cn_confidence_level),
    paste0("uninformative_BAF_threshold = ", uninformative_BAF_threshold),
    paste0("min_normal_depth = ", min_normal_depth),
    paste0("min_base_qual = ", min_base_qual),
    paste0("min_map_qual = ", min_map_qual),
    paste0("max_allowed_state = ", max_allowed_state),
    paste0("cn_upper_limit = ", cn_upper_limit),
    paste0("calc_seg_baf_option = ", calc_seg_baf_option),
    paste0("skip_allele_counting = ", paste(skip_allele_counting, collapse=", ")),
    paste0("skip_preprocessing = ", paste(skip_preprocessing, collapse=", ")),
    paste0("skip_phasing = ", paste(skip_phasing, collapse=", ")),
    paste0("skip_segmentation = ", paste(skip_segmentation, collapse=", ")),
    paste0("externalhaplotypefile = ", paste(externalhaplotypefile, collapse=", ")),
    paste0("usebeagle = ", usebeagle),
    paste0("beaglejar = ", if (is.na(beaglejar)) "NA" else beaglejar),
    paste0("beagleref.template = ", if (is.na(beagleref.template)) "NA" else beagleref.template),
    paste0("beagleplink.template = ", if (is.na(beagleplink.template)) "NA" else beagleplink.template),
    paste0("beaglemaxmem = ", beaglemaxmem),
    paste0("beaglenthreads = ", beaglenthreads),
    paste0("beaglewindow = ", beaglewindow),
    paste0("beagleoverlap = ", beagleoverlap),
    paste0("javajre = ", javajre),
    paste0("write_battenberg_phasing = ", write_battenberg_phasing),
    paste0("multisample_relative_weight_balanced = ", multisample_relative_weight_balanced),
    paste0("multisample_maxlag = ", multisample_maxlag),
    paste0("segmentation_gamma_multisample = ", segmentation_gamma_multisample),
    paste0("snp6_reference_info_file = ", if (is.na(snp6_reference_info_file)) "NA" else snp6_reference_info_file),
    paste0("apt.probeset.genotype.exe = ", apt.probeset.genotype.exe),
    paste0("apt.probeset.summarize.exe = ", apt.probeset.summarize.exe),
    paste0("norm.geno.clust.exe = ", norm.geno.clust.exe),
    paste0("birdseed_report_file = ", birdseed_report_file),
    paste0("heterozygousFilter = ", heterozygousFilter),
    paste0("prior_breakpoints_file = ", if (is.null(prior_breakpoints_file)) "NULL" else prior_breakpoints_file),
    paste0("chrom_coord_file = ", if (is.null(chrom_coord_file)) "NULL" else chrom_coord_file),
    paste0("enhanced_grid_search = ", enhanced_grid_search),
    paste0("purple_path = ", if (is.null(purple_path)) "NULL" else purple_path),
    paste0("WES_solutions = ", if (is.null(WES_solutions)) "NULL" else WES_solutions),
    paste0("test_mode = ", test_mode),
    paste0("noperms = ", noperms),
    paste0("impute_region_size = ", impute_region_size),
    paste0("test_chromosomes = ", if (is.null(test_chromosomes)) "NULL" else paste(test_chromosomes, collapse=", ")),
    paste0("organize_output = ", organize_output),
    paste0("debug_parallel = ", debug_parallel)
  )
  writeLines(run_log_lines, con=run_log_file)

  if (analysis == "cell_line"){
    calc_seg_baf_option=1
    phasing_gamma=1
    phasing_kmin=2
    segmentation_gamma=20
    segmentation_kmin=3
    # no matched normal required, but we  are generating normal counts which have this name coded
    normalname = paste0(samplename, "_normal")
    # other cell_line specific parameter values
    min_ploidy=min_ploidy
    max_ploidy=max_ploidy
    min_rho=0.99
    max_rho=1.01
  }
  if (analysis == "germline"){
    calc_seg_baf_option=1
    phasing_gamma=3
    phasing_kmin=1
    segmentation_gamma=3
    segmentation_kmin=3
    # no matched normal required, but we  are generating normal counts which have this name coded
    normalname = paste0(samplename, "_normal")
    min_ploidy=1.5
    max_ploidy=2.5
    min_rho=0.99
    max_rho=1.01
  }
  
  if (data_type=="wgs" & is.na(ismale)) {
    stop("Please provide a boolean denominator whether this sample represents a male donor")
  }
  
  if (data_type=="wgs" & is.na(g1000allelesprefix)) {
    stop("Please provide a path to 1000 Genomes allele reference files")
  }
  
  if (data_type=="wgs" & is.null(gccorrectprefix)) {
    stop("Please provide a path to GC content reference files")
  }
  
  if (data_type=="wgs" && !file.exists(problemloci)) {
    stop("Please provide a path to a problematic loci file")
  }
  
  if (!file.exists(imputeinfofile)) {
    stop("Please provide a path to an impute info file")
  }
  
  # check whether the impute_info.txt file contains correct paths
  check.imputeinfofile(imputeinfofile = imputeinfofile, is.male = ismale, usebeagle = usebeagle)
  
  # check whether multisample case
  nsamples <- length(samplename)
  if (nsamples > 1) {
    if (length(skip_allele_counting) < nsamples) {
      skip_allele_counting = rep(skip_allele_counting[1], nsamples)
    }
    if (length(skip_preprocessing) < nsamples) {
      skip_preprocessing = rep(skip_preprocessing[1], nsamples)
    }
    if (length(skip_phasing) < nsamples) {
      skip_phasing = rep(skip_phasing[1], nsamples)
    }
  }
  
  if (data_type=="wgs" | data_type=="WGS") {
    if (nsamples > 1) {
      print(paste0("Running Battenberg in multisample mode on ", nsamples, " samples: ", paste0(samplename, collapse = ", ")))
    }
    chrom_names = get.chrom.names(imputeinfofile, ismale, analysis=analysis)
  } else if (data_type=="snp6" | data_type=="SNP6") {
    if (nsamples > 1) {
      stop(paste0("Battenberg multisample mode has not been tested with SNP6 data"))
    }
    chrom_names = get.chrom.names(imputeinfofile, TRUE)
    logr_file = paste(samplename, "_mutantLogR.tab", sep="")
    allelecounts_file = NULL
  }
  print(chrom_names) 

  # Apply chromosome subsetting if requested
  if (!is.null(test_chromosomes)) {
    missing_chroms <- setdiff(test_chromosomes, chrom_names)
    if (length(missing_chroms) > 0) {
      warning("Requested test chromosomes not found in impute info: ", paste(missing_chroms, collapse=", "))
    }
    chrom_names <- chrom_names[chrom_names %in% test_chromosomes]
    if (length(chrom_names) == 0) {
      stop("No valid chromosomes remain after applying test_chromosomes filter")
    }
    print(paste0("After test_chromosomes subsetting: ", paste(chrom_names, collapse=", ")))
  }

  print(nsamples)

  preprocessing_stage_start <- Sys.time()
  phasing_stage_start <- Sys.time()
  single_seg_stage_start <- Sys.time()

  preprocessing_target <- sum(!skip_preprocessing)
  preprocessing_done <- 0
  phasing_target <- sum(!skip_phasing)
  phasing_done <- 0
  single_seg_target <- nsamples
  single_seg_done <- 0

  for (sampleidx in 1:nsamples) {
    if (!skip_preprocessing[sampleidx]) {
      if (data_type=="wgs" | data_type=="WGS") {
        # Setup for parallel computing
        clp_setup = setup_cluster_for_debugging(nthreads=nthreads, libs=libs, enable_debug=debug_parallel)
        clp = clp_setup$cluster
        
        if (analysis == "paired"){
          
          if (is.null(normalname)|is.na(normalname)){
            stop("No normal sample is specified for 'paired analysis' - a normal paired BAM is required")
            }
          prepare_wgs(chrom_names=chrom_names,
                      tumourbam=sample_data_file[sampleidx],
                      normalbam=normal_data_file,
                      tumourname=samplename[sampleidx],
                      normalname=normalname,
                      g1000allelesprefix=g1000allelesprefix,
                      g1000prefix=g1000prefix,
                      gccorrectprefix=gccorrectprefix,
                      repliccorrectprefix=repliccorrectprefix,
                      min_base_qual=min_base_qual,
                      min_map_qual=min_map_qual,
                      allelecounter_exe=allelecounter_exe,
                      min_normal_depth=min_normal_depth,
                      nthreads=nthreads,
                      skip_allele_counting=skip_allele_counting[sampleidx],
                      skip_allele_counting_normal = (sampleidx > 1),
                      seed=seed)
          
        } else if (analysis == "cell_line") {
          prepare_wgs_cell_line(chrom_names=chrom_names,
                                chrom_coord=chrom_coord_file,
                                tumourbam=sample_data_file,
                                tumourname=samplename,
                                g1000lociprefix=g1000prefix,
                                g1000allelesprefix=g1000allelesprefix, 
                                gamma_ivd=1e5,
                                kmin_ivd=50,
                                centromere_noise_seg_size=1e6,
                                centromere_dist=5e5,
                                min_het_dist=1e5, 
                                gamma_logr=100,
                                length_adjacent=5e4,
                                gccorrectprefix=gccorrectprefix, 
                                repliccorrectprefix=repliccorrectprefix,
                                min_base_qual=min_base_qual,
                                min_map_qual=min_map_qual, 
                                allelecounter_exe=allelecounter_exe,
                                min_normal_depth=min_normal_depth,
                                skip_allele_counting=skip_allele_counting[sampleidx])
        } else if (analysis == "germline"){
          
          prepare_wgs_germline(chrom_names=chrom_names,
                               chrom_coord=chrom_coord_file,
                               germlinebam=sample_data_file,
                               germlinename=samplename,
                               g1000lociprefix=g1000prefix,
                               g1000allelesprefix=g1000allelesprefix,
                               gamma_ivd=1e5,
                               kmin_ivd=50,
                               centromere_noise_seg_size=1e6,
                               centromere_dist=5e5,
                               min_het_dist=2e3,
                               gamma_logr=100,
                               length_adjacent=5e4,
                               gccorrectprefix=gccorrectprefix,
                               repliccorrectprefix=repliccorrectprefix,
                               min_base_qual=min_base_qual,
                               min_map_qual=min_map_qual,
                               allelecounter_exe=allelecounter_exe,
                               min_normal_depth=min_normal_depth,
                               skip_allele_counting=skip_allele_counting[sampleidx])
        }
        
        
        # Kill the threads
        parallel::stopCluster(clp)
        
      } else if (data_type=="snp6" | data_type=="SNP6") {
        
        prepare_snp6(tumour_cel_file=sample_data_file[sampleidx],
                     normal_cel_file=normal_data_file,
                     tumourname=samplename[sampleidx],
                     chrom_names=chrom_names,
                     snp6_reference_info_file=snp6_reference_info_file,
                     apt.probeset.genotype.exe=apt.probeset.genotype.exe,
                     apt.probeset.summarize.exe=apt.probeset.summarize.exe,
                     norm.geno.clust.exe=norm.geno.clust.exe,
                     birdseed_report_file=birdseed_report_file,
                     genomebuild=genomebuild)
        
      } else {
        print("Unknown data type provided, please provide wgs or snp6")
        q(save="no", status=1)
      }
      preprocessing_done <- preprocessing_done + 1
    }
    
    if (data_type=="snp6" | data_type=="SNP6") {
      # Infer what the gender is - WGS requires it to be specified
      gender = infer_gender_birdseed(birdseed_report_file)
      ismale = gender == "male"
    }
    
    
    if (!skip_phasing[sampleidx]) {
      
      # if external phasing data is provided (as a vcf), split into chromosomes for use in haplotype reconstruction
      if (!is.na(externalhaplotypefile) && file.exists(externalhaplotypefile)) {
        externalhaplotypeprefix <- paste0(normalname, "_external_haplotypes_chr")
        
        # if these files exist already, no need to split again
        if (any(!file.exists(paste0(externalhaplotypeprefix, 1:length(chrom_names), ".vcf")))) {
          
          print(paste0("Splitting external phasing data from ", externalhaplotypefile))
          split_input_haplotypes(chrom_names = chrom_names,
                                 externalhaplotypefile = externalhaplotypefile,
                                 outprefix = externalhaplotypeprefix)
        } else {
          print("No need to split, external haplotype files per chromosome found")
        }
      } else {
        externalhaplotypeprefix <- NA
      }
      
      # Setup for parallel computing
      clp_setup = setup_cluster_for_debugging(nthreads=nthreads, libs=NULL, enable_debug=debug_parallel)
      clp = clp_setup$cluster
      
      # Reconstruct haplotypes
      # mclapply(1:length(chrom_names), function(chrom) {
      if (analysis=="germline"){
        foreach::foreach (i=1:length(chrom_names), .verbose=debug_parallel, .errorhandling="stop") %dopar% {
          chrom = chrom_names[i]
          print(chrom)
          
          run_haplotyping_germline(chrom=chrom,
                                   germlinename=samplename,
                                   normalname=normalname,
                                   ismale=ismale,
                                   imputeinfofile=imputeinfofile,
                                   problemloci=problemloci,
                                   impute_exe=impute_exe,
                                   min_normal_depth=min_normal_depth,
                                   chrom_names=chrom_names, 
                                   externalhaplotypeprefix = NA,
                                   use_previous_imputation=F,
                                   snp6_reference_info_file=NA, 
                                   heterozygousFilter=NA,
                                   usebeagle=usebeagle,
                                   beaglejar=beaglejar,
                                   beagleref=gsub("CHROMNAME",chrom,beagleref.template),
                                   beagleplink=gsub("CHROMNAME",chrom,beagleplink.template),
                                   beaglemaxmem=beaglemaxmem,
                                   beaglenthreads=beaglenthreads,
                                   beaglewindow=beaglewindow,
                                   beagleoverlap=beagleoverlap,
                                   region.size=impute_region_size,
                                   test_mode=test_mode,
                                   seed=seed)      
        }
      } else {
        foreach::foreach (i=1:length(chrom_names), .verbose=debug_parallel, .errorhandling="stop") %dopar% {
          chrom = chrom_names[i]
          print(chrom)      
          run_haplotyping(chrom=chrom,
                          tumourname=samplename[sampleidx],
                          normalname=normalname,
                          ismale=ismale,
                          imputeinfofile=imputeinfofile,
                          problemloci=problemloci,
                          impute_exe=impute_exe,
                          min_normal_depth=min_normal_depth,
                          chrom_names=chrom_names,
                          snp6_reference_info_file=snp6_reference_info_file,
                          heterozygousFilter=heterozygousFilter,
                          usebeagle=usebeagle,
                          beaglejar=beaglejar,
                          beagleref=gsub("CHROMNAME", chrom, beagleref.template),
                          beagleplink=gsub("CHROMNAME", chrom, beagleplink.template),
                          beaglemaxmem=beaglemaxmem,
                          beaglenthreads=beaglenthreads,
                          beaglewindow=beaglewindow,
                          beagleoverlap=beagleoverlap,
                          externalhaplotypeprefix=externalhaplotypeprefix,
                          use_previous_imputation=(sampleidx > 1),
                          region.size=impute_region_size,
                          test_mode=test_mode,
                          seed=seed)
        }
      }
      
      # Kill the threads as from here its all single core
      parallel::stopCluster(clp)
      
      # Combine all the BAF output into a single file
      combine.baf.files(inputfile.prefix=paste(samplename[sampleidx], "_chr", sep=""),
                        inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt",
                        outputfile=paste(samplename[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                        chr_names=chrom_names)
      phasing_done <- phasing_done + 1
    }

    print(samplename[sampleidx])
    
    if (!skip_segmentation) {
    print('SEGMENTING BAF HERE')

    # Segment the phased and haplotyped BAF data
    segment.baf.phased(samplename=samplename[sampleidx],
                       inputfile=paste(samplename[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                       outputfile=paste(samplename[sampleidx], ".BAFsegmented.txt", sep=""),
                       prior_breakpoints_file=prior_breakpoints_file,
                       gamma=segmentation_gamma,
                       phasegamma=phasing_gamma,
                       kmin=segmentation_kmin,
                       phasekmin=phasing_kmin,
                       calc_seg_baf_option=calc_seg_baf_option)
              single_seg_done <- single_seg_done + 1
    
    if (nsamples > 1 ) {
      # Write the Battenberg phasing information to disk as a vcf
      write_battenberg_phasing(tumourname = samplename[sampleidx],
                               SNPfiles = paste0(samplename[sampleidx], "_alleleFrequencies_chr", chrom_names, ".txt"),
                               imputedHaplotypeFiles = paste0(samplename[sampleidx], "_impute_output_chr", chrom_names, "_allHaplotypeInfo.txt"),
                               bafsegmented_file = paste0(samplename[sampleidx], ".BAFsegmented.txt"),
                               outprefix = paste0(samplename[sampleidx], "_Battenberg_phased_chr"),
                               chrom_names = chrom_names,
                               include_homozygous = F)
    }
    } else {
      single_seg_done <- single_seg_done + 1
    }
    
  }

  if (preprocessing_target > 0 && preprocessing_done == preprocessing_target) {
    emit_stage_event(
      stage = "preprocessing",
      status = "complete",
      stage_start = preprocessing_stage_start,
      extra_info = paste0("samples_completed=", preprocessing_done, "/", preprocessing_target)
    )
  } else if (preprocessing_target == 0) {
    emit_stage_event(
      stage = "preprocessing",
      status = "skipped",
      stage_start = preprocessing_stage_start,
      extra_info = "samples_completed=0/0"
    )
  }

  if (phasing_target > 0 && phasing_done == phasing_target) {
    emit_stage_event(
      stage = "phasing",
      status = "complete",
      stage_start = phasing_stage_start,
      extra_info = paste0("samples_completed=", phasing_done, "/", phasing_target)
    )
  } else if (phasing_target == 0) {
    emit_stage_event(
      stage = "phasing",
      status = "skipped",
      stage_start = phasing_stage_start,
      extra_info = "samples_completed=0/0"
    )
  }

  if (single_seg_done == single_seg_target) {
    emit_stage_event(
      stage = "single_sample_segmentation",
      status = if (skip_segmentation) "skipped" else "complete",
      stage_start = single_seg_stage_start,
      extra_info = paste0("samples_completed=", single_seg_done, "/", single_seg_target)
    )
  }
  
  # if this is a multisample run, combine the battenberg phasing outputs, incorporate it and resegment
  multisample_stage_start <- Sys.time()
  if (nsamples > 1 && !skip_segmentation) {
    print("Constructing multisample phasing")
    multisamplehaplotypeprefix <- paste0(normalname, "_multisample_haplotypes_chr")
    
    
      # Setup for parallel computing
      clp_setup = setup_cluster_for_debugging(nthreads=nthreads, libs=libs, enable_debug=debug_parallel)
      clp = clp_setup$cluster
    
    print(chrom_names)
    # Reconstruct haplotypes
    .libPaths()
    foreach::foreach (i=1:length(chrom_names), .verbose=debug_parallel, .errorhandling="stop") %dopar% {
      .libPaths()
      chrom = chrom_names[i]
      print(chrom)
      
      get_multisample_phasing(chrom = chrom,
                              bbphasingprefixes = paste0(samplename, "_Battenberg_phased_chr"),
                              maxlag = multisample_maxlag,
                              relative_weight_balanced = multisample_relative_weight_balanced,
                              outprefix = multisamplehaplotypeprefix)
    }
    
    # continue over all samples to incorporate the multisample phasing
    for (sampleidx in 1:nsamples) {
      
      # rename the original files without multisample phasing info
      MutBAFfiles <- paste0(samplename[sampleidx], "_chr", chrom_names, "_heterozygousMutBAFs_haplotyped.txt")
      heterozygousdatafiles <- paste0(samplename[sampleidx], "_chr", chrom_names, "_heterozygousData.png")
      raffiles <- paste0(samplename[sampleidx], "_RAFseg_chr", chrom_names, ".png")
      segfiles <- paste0(samplename[sampleidx], "_segment_chr", chrom_names, ".png")
      haplotypedandbafsegmentedfiles <- paste0(samplename[sampleidx], c("_heterozygousMutBAFs_haplotyped.txt", ".BAFsegmented.txt"))
      
      file.copy(from = MutBAFfiles, to = gsub(pattern = ".txt$", replacement = "_noMulti.txt", x = MutBAFfiles), overwrite = T)
      file.copy(from = heterozygousdatafiles, to = gsub(pattern = ".png$", replacement = "_noMulti.png", x = heterozygousdatafiles), overwrite = T)
      file.copy(from = raffiles, to = gsub(pattern = ".png$", replacement = "_noMulti.png", x = raffiles), overwrite = T)
      file.copy(from = segfiles, to = gsub(pattern = ".png$", replacement = "_noMulti.png", x = segfiles), overwrite = T)
      file.copy(from = haplotypedandbafsegmentedfiles, to = gsub(pattern = ".txt$", replacement = "_noMulti.txt", x = haplotypedandbafsegmentedfiles), overwrite = T)
      # done renaming, next sections will overwrite orignals
      
      
      foreach::foreach (i=1:length(chrom_names), .verbose=debug_parallel, .errorhandling="stop") %dopar% {
        chrom = chrom_names[i]
        print(chrom)
        
        input_known_haplotypes(chrom = chrom,
                               chrom_names = chrom_names,
                               imputedHaplotypeFile = paste0(samplename[sampleidx], "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
                               externalHaplotypeFile = paste0(multisamplehaplotypeprefix, chrom, ".vcf"),
                               oldfilesuffix = "_noMulti.txt")
        
        GetChromosomeBAFs(chrom=chrom,
                          SNP_file=paste0(samplename[sampleidx], "_alleleFrequencies_chr", chrom, ".txt"),
                          haplotypeFile=paste0(samplename[sampleidx], "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
                          samplename=samplename[sampleidx],
                          outfile=paste0(samplename[sampleidx], "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
                          chr_names=chrom_names,
                          minCounts=min_normal_depth)

	# Plot what we have until this point
        plot.haplotype.data(haplotyped.baf.file=paste0(samplename[sampleidx], "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
                            imageFileName=paste0(samplename[sampleidx],"_chr",chrom,"_heterozygousData.png"),
                            samplename=samplename[sampleidx],
                            chrom=chrom,
                            chr_names=chrom_names)
      }
      
    }
    
    # Kill the threads as from here its single core
    parallel::stopCluster(clp)
    
    for (sampleidx in 1:nsamples) {
      
      # Combine all the BAF output into a single file
      combine.baf.files(inputfile.prefix=paste0(samplename[sampleidx], "_chr"),
                        inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt",
                        outputfile=paste0(samplename[sampleidx], "_heterozygousMutBAFs_haplotyped.txt"), 
                        chr_names=chrom_names)
      
    }
    # Segment the phased and haplotyped BAF data
    segment.baf.phased.multisample(samplename=samplename,
                                   inputfile=paste(samplename, "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                                   outputfile=paste(samplename, ".BAFsegmented.txt", sep=""),
                                   prior_breakpoints_file=prior_breakpoints_file,
                                   gamma=segmentation_gamma_multisample,
                                   calc_seg_baf_option=calc_seg_baf_option,
                                   GENOMEBUILD=genomebuild)

    emit_stage_event(
      stage = "multisample_rephasing_segmentation",
      status = "complete",
      stage_start = multisample_stage_start,
      extra_info = paste0("samples_completed=", nsamples, "/", nsamples)
    )
  } else {
    emit_stage_event(
      stage = "multisample_rephasing_segmentation",
      status = "skipped",
      stage_start = multisample_stage_start,
      extra_info = "samples_completed=0/0"
    )
    
  }
  
  # Setup for parallel computing
  final_fit_stage_start <- Sys.time()
  clp_setup = setup_cluster_for_debugging(nthreads=min(nthreads, nsamples), libs=libs, enable_debug=debug_parallel)
  clp = clp_setup$cluster
  # for (sampleidx in 1:nsamples) {
  foreach::foreach (sampleidx=1:nsamples, .verbose=debug_parallel, .errorhandling="stop") %dopar% {
    tryCatch({
    print(paste0("Fitting final copy number and calling subclones for sample ", samplename[sampleidx]))
    
    if (data_type=="wgs" | data_type=="WGS") {
      logr_file = paste(samplename[sampleidx], "_mutantLogR_gcCorrected.tab", sep="")
      if (analysis=="paired") {
        allelecounts_file = paste(samplename[sampleidx], "_alleleCounts.tab", sep="")
      } else {
      # Not produced by a number of analysis and is required for some plots. Setting to NULL  makes the pipeline not attempt to create these plots
        allelecounts_file = NULL
      }
    }
    
    # Fit a clonal copy number profile
    fit.copy.number(samplename=samplename[sampleidx],
                    outputfile.prefix=paste(samplename[sampleidx], "_", sep=""),
                    inputfile.baf.segmented=paste(samplename[sampleidx], ".BAFsegmented.txt", sep=""),
                    inputfile.baf=paste(samplename[sampleidx],"_mutantBAF.tab", sep=""),
                    inputfile.logr=logr_file,
                    dist_choice=clonality_dist_metric,
                    ascat_dist_choice=ascat_dist_metric,
                    min.ploidy=min_ploidy,
                    max.ploidy=max_ploidy,
                    min.rho=min_rho,
                    max.rho=max_rho,
                    min.goodness=min_goodness,
                    uninformative_BAF_threshold=uninformative_BAF_threshold,
                    gamma_param=platform_gamma,
                    use_preset_rho_psi=F,
                    preset_rho=NA,
                    preset_psi=NA,
                    read_depth=30,
                    analysis=analysis,
                    nthreads=nthreads,
                    enhanced_grid_search=enhanced_grid_search,
                    PURPLE_purity_path=purple_path[sampleidx],
                    External_WES_purity_path = WES_solutions)
    
    # KT: need to iterate over callSubclones for each soltion that we find in fit.copy.numberd
    # first read in file with all soltions
    all_solutions <- read.table(paste0(samplename[sampleidx], "_all_solutions_rho_psi.txt"), head = T, sep = "\t")
    
    for(solution in 1:max(all_solutions$solution)){
      
      print(paste0("calling subclonal copy number for solution ", solution, " of ", max(all_solutions$solution)))
      
      rho <- all_solutions[all_solutions$solution == solution, "rho"]
      psi <- all_solutions[all_solutions$solution == solution, "psi"]
      solution_type <- all_solutions[all_solutions$solution == solution, "solution_type"]
      
      rho.psi.file.solution <- paste0(samplename[sampleidx], "_", solution_type, "_psi", psi, "_rho", rho, "_runclonalASCAT_rho_and_psi.txt")
      output.file.solution <- paste0(samplename[sampleidx], "_", solution_type, "_psi", psi, "_rho", rho, "_subclones.txt")
      output.figures.prefix.solution <- paste0(samplename[sampleidx], "_", solution_type, "_psi", psi, "_rho", rho, "_subclones_chr")
      output.gw.figures.prefix.solution  <- paste0(samplename[sampleidx], "_", solution_type, "_psi", psi, "_rho", rho, "_BattenbergProfile")
      masking_output_file <- paste0(samplename[sampleidx], "_", solution_type, "_psi", psi, "_rho", rho, "_segment_masking_details.txt")
      
      # Go over all segments, determine which segements are a mixture of two states and fit a second CN state
      callSubclones(sample.name=samplename[sampleidx],
                    baf.segmented.file=paste(samplename[sampleidx], ".BAFsegmented.txt", sep=""),
                    logr.file=logr_file,
                    rho.psi.file=rho.psi.file.solution,
                    output.file=output.file.solution,
                    output.figures.prefix=output.figures.prefix.solution,
                    output.gw.figures.prefix=output.gw.figures.prefix.solution,
                    masking_output_file=masking_output_file,
                    prior_breakpoints_file=prior_breakpoints_file,
                    chr_names=chrom_names, 
                    gamma=platform_gamma, 
                    segmentation.gamma=NA, 
                    siglevel=0.05, 
                    maxdist=0.01, 
                    max_allowed_state=max_allowed_state, 
                    cn_upper_limit=cn_upper_limit, 
                    noperms=noperms,
                    cn_confidence_level=cn_confidence_level,
                    seed=seed,
                    calc_seg_baf_option=calc_seg_baf_option,
                    RHO=rho,
                    PSI=psi,
                    solution_type=solution_type)
      
      # If patient is male, get copy number status of ChrX based only on logR segmentation (due to hemizygosity of SNPs)
      # Only do this when X chromosome is included
      if (ismale & "X" %in% chrom_names){
        callChrXsubclones(tumourname=samplename[sampleidx],
                          X_gamma=1000,
                          X_kmin=100,
                          genomebuild=genomebuild,
                          AR=TRUE,
                          prior_breakpoints_file=prior_breakpoints_file,
			                    chrom_names=chrom_names,
                          data_type=data_type,
                          RHO = rho,
                          PSI = psi,
                          solution_type = solution_type,
                          seed=seed)
      }
      
      # Make some post-hoc plots
      make_posthoc_plots(samplename=samplename[sampleidx],
                         logr_file=logr_file,
                         subclones_file=paste(samplename[sampleidx], "_", solution_type, "_psi", psi, "_rho", rho, "_subclones.txt", sep=""),
                         rho_psi_file=paste(samplename[sampleidx], "_", solution_type, "_psi", psi, "_rho", rho, "_runclonalASCAT_rho_and_psi.txt", sep=""),
                         bafsegmented_file=paste(samplename[sampleidx], ".BAFsegmented.txt", sep=""),
                         logrsegmented_file=paste(samplename[sampleidx], ".logRsegmented.txt", sep=""),
                         allelecounts_file=allelecounts_file,
                         RHO = rho,
                         PSI = psi,
                         solution_type = solution_type)
      
      # Save refit suggestions for a future rerun
      cnfit_to_refit_suggestions(samplename=samplename[sampleidx],
                                 subclones_file=paste(samplename[sampleidx], "_", solution_type, "_psi", psi, "_rho", rho, "_subclones.txt", sep=""),
                                 rho_psi_file=paste(samplename[sampleidx], "_", solution_type, "_psi", psi, "_rho", rho, "_runclonalASCAT_rho_and_psi.txt", sep=""),
                                 gamma_param=platform_gamma,
                                 RHO = rho,
                                 PSI = psi,
                                 solution_type = solution_type)

    }
    }, error = function(e) {
      msg <- paste0(
        "\n========== WORKER ERROR ==========\n",
        "Sample: ", samplename[sampleidx], " (sampleidx=", sampleidx, ")\n",
        "Error: ", conditionMessage(e), "\n",
        "Call: ", deparse(conditionCall(e)), "\n",
        "Traceback:\n",
        paste(capture.output(traceback(4)), collapse="\n"), "\n",
        "==================================\n"
      )
      cat(msg, file=stderr())
      cat(msg)
      stop(paste0("sample ", samplename[sampleidx], ": ", conditionMessage(e)))
    })
  }
  parallel::stopCluster(clp)

  emit_stage_event(
    stage = "final_fit_and_subclones",
    status = "complete",
    stage_start = final_fit_stage_start,
    extra_info = paste0("samples_completed=", nsamples, "/", nsamples)
  )

  #KT: compress files 
  compression_stage_start <- Sys.time()

  # If organizing output, move files into subdirectories first
  if (organize_output) {
    organize_output_files()
  }

  # Compress text, tab, and log files (search subdirectories if organized)
  gzip_files <- function(pattern) {
    files <- list.files(pattern = pattern, recursive = organize_output)
    n <- 0L
    for (f in files) {
      gz_name <- paste0(f, ".gz")
      if (!file.exists(gz_name)) {
        R.utils::gzip(f, gz_name)
        n <- n + 1L
      }
    }
    return(list(total = length(files), compressed = n))
  }

  txt_stats <- gzip_files("\\.txt$")
  tab_stats <- gzip_files("\\.tab$")
  log_stats <- gzip_files("\\.log$")

  emit_stage_event(
    stage = "compression",
    status = "complete",
    stage_start = compression_stage_start,
    extra_info = paste0("txt_files=", txt_stats$total,
                        " tab_files=", tab_stats$total,
                        " log_files=", log_stats$total,
                        if (organize_output) " organized=TRUE" else "")
  )

  emit_stage_event(
    stage = "pipeline",
    status = "complete",
    stage_start = pipeline_start_time,
    extra_info = paste0("samples_completed=", nsamples, "/", nsamples)
  )
}

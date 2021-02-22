#' Function to read NPX data into long format
#'
#' Modifications have been made to allow the reading of our specific dataset and
#' add GeneIDs in addition to the olink assay names. 
#'
#' @param filename Path to file NPX Manager output file.
#' @param sample_manifest Per individual phenotypic data
#' @param pheno Per sample phenotypic data
#' @param skip_mod How many lines to skip in the excel file
#' @param panel Which panel is currently being read
#' @param tab Which excel tab the data is being read from
#' @param this_exp Whether the data is plasma or serum
#' 
#' @return A tibble in long format.
#' @keywords NPX
#' @export
#' @examples \donttest{read_NPX("~/NPX data.xlsx")}
#' @import dplyr stringr tidyr biomaRt

read_NPX <- function(filename, this_exp="plasma", num_samples=31){

  raw_dat <-  data.frame(readxl::read_excel(filename, col_names = F,.name_repair="minimal"))
  
  plate_info_index <- which(grepl("Plate ID", raw_dat[2,]))
  plate_info <- raw_dat[,plate_info_index]
  raw_dat <- raw_dat[,-plate_info_index]
  plate_info <- plate_info[c(1, 5),]
  colnames(plate_info) <- plate_info[1,]
  plate_info <- plate_info[-1,]
  
  qc_info_index <- which(grepl("QC Warning", raw_dat[2,]))
  qc_info <- raw_dat[,c(1, qc_info_index)]
  raw_dat <- raw_dat[,-qc_info_index]
  qc_info <- qc_info[-c(2:4),]
  colnames(qc_info) <- qc_info[1,]
  qc_info <- qc_info[-1,]
  qc_info <- qc_info[1:num_samples,]
  rownames(qc_info) <- qc_info[,1]
  qc_info <- qc_info[,-1]
  qc_info$SampleID <- rownames(qc_info)
  
  colnames(raw_dat) <- raw_dat[4,]
  assay_info_index <- c(1:4, (num_samples+5):(num_samples+7))
  assay_info <- data.frame(t(raw_dat[assay_info_index,]))
  colnames(assay_info) <- assay_info[1,]
  assay_info <- assay_info[-1,]
  raw_dat <- raw_dat[-assay_info_index,]
  
  gene_ids <- uniprot_to_gene(data.frame(UniProt=assay_info$`Uniprot ID`, target=assay_info$Assay))
  gene_ids <- dplyr::select(gene_ids, UniProt, target, hgnc_symbol, prot.on.multiple.panel)
  assay_info <- dplyr::left_join(assay_info, gene_ids, by=c("Uniprot ID"="UniProt", "Assay"="target"))
  colnames(assay_info) <- c("Panel", "Assay", "Uniprot", "OlinkID", "LOD", "MissingFreq", "Normalisation", "GeneID", "ProtOnMultiplePanels")
  
  stopifnot(nrow(assay_info) == 460)
  
  rownames(raw_dat) <- raw_dat[,1]
  
  raw_dat <- tidyr::pivot_longer(raw_dat, cols = starts_with("OID"))
  colnames(raw_dat) <- c("SampleID", "OlinkID", "NPX")
  stopifnot(nrow(raw_dat) == num_samples * 460)
  
  raw_dat <- dplyr::left_join(raw_dat, assay_info)
  stopifnot(nrow(raw_dat) == num_samples * 460)
  
  long_panel_name_shortener <- function(long_name) {
    if (grepl("Cardiometabolic", long_name)) {
      return("CM")
    } else if (grepl("Cardiovascular III", long_name)) {
      return("CVD3")
    } else if (grepl("Cardiovascular II", long_name)) {
      return("CVD2")
    } else if (grepl("Immune Response", long_name)) {
      return("IR")
    } else if (grepl("Inflammation", long_name)) {
      return("Inf")
    } else {
      return(long_name)
    }
  }
  
  raw_dat$panel_name <- sapply(raw_dat$Panel, long_panel_name_shortener)
  colnames(qc_info) <- sapply(colnames(qc_info), long_panel_name_shortener)
  
  qc_info <- tidyr::pivot_longer(qc_info, c("CM", "CVD2", "CVD3", "IR", "Inf"))
  colnames(qc_info) <- c("SampleID", "panel_name", "QC_Warning")
  
  raw_dat <- dplyr::left_join(raw_dat, qc_info)
  stopifnot(nrow(raw_dat) == num_samples * 460)
  
  raw_dat$matrix <- this_exp
  
  return(raw_dat)
}

#' Map protein IDs to gene IDs. 
#' 
#' @param df meta_dat dataframe
#' 
#' @return Gene IDs

uniprot_to_gene <- function(df) {
  df$uniprot <- df$UniProt
  
  # data input error by Olink 'o' instead of 'O'
  df$target <- gsub("IL-2oRA", "IL-20RA", df$target)
  
  # turns out TWEAK labelled with an out of date or inferior UP id
  if (TRUE %in% grepl('TWEAK', df$target)) {
    df$uniprot[grep('TWEAK', df$target)] <- "O43508"
  }
  
  # no uniprot id for NT-proBNP - there is for proBNP
  if ("NT-proBNP" %in% df$target) {
    df[df$target == "NT-proBNP",]$uniprot <- "P16860"
  }
  # Clean up 1. identify bad entries: Olink have made some errors and Excel import causes some problems
  
  # clean whitespace
  df$uniprot <- gsub('\\\r\\\n', ";", df$uniprot, ignore.case = F)
  df$uniprot <- gsub(', |,', ";", df$uniprot, ignore.case = F)
  df$uniprot <- gsub("[[:space:]]", "", df$uniprot)
  
  # Clean up 2. '-' represents isoform notation eg O43521-2
  
  df$uniprot.isoform <- NA
  
  df$uniprot.isoform[grep('-', df$uniprot)] <- grep('-', df$uniprot, v=T)
  df$uniprot <- gsub("-[0-9]$", "", df$uniprot)

  # Special circumstances 2:two ids for protein complex eg IL12A-IL12B
  # uniprot ids sep by ';'
  # df[grep(";", df$uniprot), ]
  
  df$multiple.proteins <- FALSE
  df$multiple.proteins[grep(";", df$uniprot)] <- TRUE
  
  df$protein.1 <- df$uniprot  
  df$protein.2 <- NA
  
  df$protein.2[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot, pattern = ";[A-Z0-9]+")[which(df$multiple.proteins==T)]
  df$protein.2 <- gsub("^;", "", df$protein.2)
  
  df$protein.1[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot, pattern = "^[A-Z0-9]+;")[which(df$multiple.proteins==T)]
  df$protein.1 <- gsub(";$", "", df$protein.1)
  
  # where there are 2 uniprot ids (eg protein complex) the uniprot ids are not always in consistent order
  # lets make them in consistent alphabetical order
  df$uniprot.ordered <- NA
  df$uniprot.ordered[!df$multiple.proteins] <- df$uniprot[!df$multiple.proteins]
  
  alphabetize.up <- function(x) {
    if( !inherits(x, what='data.frame')){
      stop('argument must be a data.frame')
    }
    y <- paste(sort(c( x$protein.1, x$protein.2)),collapse=';')
    y
  }
  
  inds <- which(df$multiple.proteins)
  
  for (i in inds){
    df$uniprot.ordered[i] <- alphabetize.up(df[i,]) 
  }
  
  #annoying that p1 and p2 are arbitrary: now we've ordered things let's start over on this front
  df$protein.1 <- df$protein.2 <-NULL
  
  # now repeat the exercise for p1 and p2 using the alphabetized concatenation
  
  df$protein.1 <- df$uniprot.ordered  
  df$protein.2 <- NA
  
  df$protein.2[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot.ordered, pattern = ";[A-Z0-9]+")[which(df$multiple.proteins==T)]
  df$protein.2 <- gsub("^;", "", df$protein.2)
  
  df$protein.1[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot.ordered, pattern = "^[A-Z0-9]+;")[which(df$multiple.proteins==T)]
  df$protein.1 <- gsub(";$", "", df$protein.1)
  
  # col to identify dup proteins and which panels
  
  dup.prots <- union( which( duplicated(df$uniprot.ordered)), which( duplicated(df$uniprot.ordered, fromLast = T)) )
  df$prot.on.multiple.panel <- FALSE
  df$prot.on.multiple.panel[dup.prots] <- TRUE
  
  df$panels.with.prot <- NA
  
  
  tmp.list <- split( df[dup.prots,], f=df$uniprot.ordered[dup.prots] )
  
  mylist <- lapply(tmp.list, FUN = function(x) paste( as.character(x$panel), collapse=";" ) )

  for (i in dup.prots){
    uprot <- df$uniprot.ordered[i]
    df[i, "panels.with.prot"] <- mylist[[uprot]]
  }
  
  #--------------------- Gene symbol annotation ---------------------#
  
  # matching to gene symbols: do for p1 and p2
  
  #ensembl <- biomaRt::useMart(biomart="ensembl",
  #                   dataset="hsapiens_gene_ensembl",
  #                   host='http://jul2018.archive.ensembl.org')

  #filters <- biomaRt::listFilters(ensembl)

  #x <- biomaRt::getBM(attributes = c('uniprotswissprot', 'hgnc_symbol', 'entrezgene', 'chromosome_name'),
  #             filters = 'uniprotswissprot',
  #             values = df$protein.1,
  #             mart = ensembl)

  # some UP ids not found by BioMart: turns out we have outdated IDs
  #df[which(!df$protein.1 %in% x$uniprotswissprot),]
  
  #--------------------- Try an archived version of Ensembl ---------------------#
  
  # find urls for old ensembl versions
  biomaRt::listEnsemblArchives()
  
  # hg19/GRCh37
  ensembl.hg19 <- biomaRt::useMart(biomart= "ENSEMBL_MART_ENSEMBL",
                          dataset="hsapiens_gene_ensembl")
  
  # note attribute names differ in the older release
  gene.pos <- biomaRt::getBM(attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene',
                                   'chromosome_name', 'start_position', 'end_position'), 
                    filters = 'uniprotswissprot', 
                    values = unique(df$protein.1), 
                    mart = ensembl.hg19)
  
  # P0DOY2 / IGLC2 is present within db but swissprot id is NA
  missing_to_add <- biomaRt::getBM(
    attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene', 
                   'chromosome_name', 'start_position', 'end_position'), 
    filters = 'hgnc_symbol', 
    values = "IGLC2", 
    mart = ensembl.hg19)
  
  missing_to_add[missing_to_add$hgnc_symbol == "IGLC2",]$uniprotswissprot <- "P0DOY2"
  
  # same for hgnc_symbol
  missing_to_add <- rbind(missing_to_add, data.frame(biomaRt::getBM(
    attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene', 
                   'chromosome_name', 'start_position', 'end_position'), 
    filters = 'hgnc_symbol', 
    values = "PECAM1", 
    mart = ensembl.hg19)))
  
  missing_to_add[missing_to_add$hgnc_symbol == "PECAM1",]$uniprotswissprot <- "P16284"
  
  # same for CDSN
  missing_to_add <- rbind(missing_to_add, data.frame(biomaRt::getBM(
    attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene', 
                   'chromosome_name', 'start_position', 'end_position'), 
    filters = 'hgnc_symbol', 
    values = "CDSN", 
    mart = ensembl.hg19)))
  
  missing_to_add[missing_to_add$hgnc_symbol == "CDSN",]$uniprotswissprot <- "Q15517"

  gene.pos <- rbind(gene.pos, missing_to_add)

  # there are some duplicated genes
  dup.ind <- union( which(duplicated(gene.pos$hgnc_symbol)),
                    which(duplicated(gene.pos$hgnc_symbol, fromLast = T))
  )
  
  # strange chr names
  
  strange.ind <- which(!gene.pos$chromosome_name %in% c(1:22, 'X', 'Y'))
  
  to.cut <- intersect(dup.ind, strange.ind)
  
  gene.pos2 <- gene.pos[-to.cut,]

  #-------------------------------------------------------------------------------#
  # some proteins map to multiple genes
  # map_ids <- vector("integer", nrow(df)) 
  # for (i in 1:nrow(df)) {
  #   map_ids[i] <- which(gene.pos2$uniprotswissprot == df$protein.1[i])[1]
  # }
  # gene.pos2 <- gene.pos2[map_ids,]
  
  gene.pos2 <- gene.pos2[!is.na(gene.pos2$hgnc_symbol),]
  gene.pos2 <- gene.pos2[!duplicated(gene.pos2$uniprotswissprot),]
  
  
  df2 <- left_join(df, gene.pos2, by=c("protein.1" = "uniprotswissprot"))
  df2 <- df2[!duplicated(df2$target),]
  
  # is blank
  df2$hgnc_symbol[df2$UniProt == "P50225"] <- "SULT1A1"
  df2$hgnc_symbol[df2$UniProt == "Q8NHL6"] <- "LILRB1"
  df2$hgnc_symbol[df2$UniProt == "P01137"] <- "TGFB1"
  
  return(df2)
}

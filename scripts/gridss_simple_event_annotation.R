#!/usr/bin/env Rscript
# Annotates a raw gridss VCF (breakend-only) with INFO/SIMPLE_TYPE and
# INFO/SVLEN so dicast can read its DEL/DUP/INS events.
#
# Adapted from PapenfussLab/gridss example/simple-event-annotation.R: the
# classification logic is unmodified; SVLEN is declared in the header (the
# 2018 script relied on VariantAnnotation silently accepting an undeclared
# INFO write-back, which current versions reject) and the file names come
# from the command line.
#
# Usage:  Rscript gridss_simple_event_annotation.R gridss_raw.vcf gridss_simple.vcf [genome]
# Needs:  R with VariantAnnotation, StructuralVariantAnnotation, stringr
#         (conda: bioconductor-variantannotation bioconductor-structuralvariantannotation r-stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript gridss_simple_event_annotation.R <gridss_raw.vcf> <out.vcf> [genome=hg38]")
}
in_vcf <- args[1]
out_vcf <- args[2]
genome <- if (length(args) >= 3) args[3] else "hg38"

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(StructuralVariantAnnotation)
  library(stringr)
})

simpleEventType <- function(gr) {
  pgr <- partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomal
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))
}

vcf <- readVcf(in_vcf, genome)
info(header(vcf)) <- unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
  row.names = c("SIMPLE_TYPE", "SVLEN"),
  Number = c("1", "1"),
  Type = c("String", "Integer"),
  Description = c("Simple event type annotation based purely on breakend position and orientation.",
                  "Simple event length (from StructuralVariantAnnotation breakpointRanges)."))), "DataFrame"))
gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf)$SVLEN <- NA_integer_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen
writeVcf(vcf, out_vcf)

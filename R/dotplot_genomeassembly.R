# README
# This script reads and filters delta alignment files, and visualizes the filtered alignments using various ggplot2 plots.
# The script performs the following steps:
# 1. Installs and loads the necessary packages.
# 2. Defines functions to read and filter delta files.
# 3. Reads the delta file and filters contigs based on alignment size.
# 4. Creates basic and diagonal synteny plots.
# 5. Plots percent identity and coverage.

# Author: Bloodmark

# Input details:
# - Input delta file: "/home/bloodmark/Desktop/nucmer/chr10/out.delta"

# Output details:
# - Various plots visualizing the filtered alignments.

# Script functioning:
# 1. Installs and loads the necessary packages (GenomicRanges, dplyr, magrittr, knitr, ggplot2, tidyr).
# 2. Defines functions `readDelta` and `filterMum` to read and filter delta files.
# 3. Reads the delta file and filters contigs based on alignment size.
# 4. Creates basic and diagonal synteny plots.
# 5. Plots percent identity and coverage using ggplot2.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input delta file is in the specified path and properly formatted.

# Code starts below:

# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("GenomicRanges")
install.packages("dplyr")
install.packages("magrittr")
install.packages("knitr")
install.packages("ggplot2")
install.packages("tidyr")

# Load required libraries
library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)

# Function to read delta file
readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

# Function to filter contigs with poor alignments based on size of alignment
filterMum <- function(df, minl=1000, flanks=1e4){
  coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
    summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
    ungroup %>% arrange(desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
  merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
    mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

# Read delta file
mumgp = readDelta("/home/bloodmark/Desktop/nucmer/chr10/out.delta")

# Display the first few rows of the data
mumgp %>% head %>% kable

# Filter contigs with poor alignments based on size of alignment
mumgp.filt = filterMum(mumgp, minl=1e2)
mumgp.filt %>% head %>% kable

# Basic plot
ggplot(mumgp.filt, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) + geom_segment() +
  geom_point(alpha=0.1) + facet_grid(qid~., scales='free', space='free', switch='y') +
  theme_bw() + theme(strip.text.y=element_text(angle=180, size=10),
                     legend.position=c(.99,.01), legend.justification=c(1,0),
                     strip.background=element_blank(),
                     axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Pastel1')

# Make diagonal
diagMum <- function(df){
  # Find best qid order
  rid.o = df %>% group_by(qid, rid) %>% summarize(base=sum(abs(qe-qs)),
                                                  rs=weighted.mean(rs, abs(qe-qs))) %>%
    ungroup %>% arrange(desc(base)) %>% group_by(qid) %>% do(head(., 1)) %>%
    ungroup %>% arrange(desc(rid), desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid)))
  # Find best qid strand
  major.strand = df %>% group_by(qid) %>%
    summarize(major.strand=ifelse(sum(sign(qe-qs)*abs(qe-qs))>0, '+', '-'),
              maxQ=max(c(qe, qs)))
  merge(df, major.strand) %>% mutate(qs=ifelse(major.strand=='-', maxQ-qs, qs),
                                     qe=ifelse(major.strand=='-', maxQ-qe, qe),
                                     qid=factor(qid, levels=levels(rid.o$qid)))
}

mumgp.filt.diag = diagMum(mumgp.filt)

ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment() + geom_point(alpha=0.1) + theme_bw() + 
  facet_grid(qid~., scales='free', space='free', switch='y') +
  theme(strip.text.y=element_text(angle=180, size=15), strip.background=element_blank(),
        legend.position=c(.99,.01), legend.justification=c(1,0),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set2')

ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment() + geom_point(alpha=.5) + theme_bw() + 
  facet_grid(qid~rid, scales='free', space='free', switch='y') +
  theme(strip.text.y=element_text(angle=180, size=10), strip.background=element_blank(),
        legend.position=c(.99,.01), legend.justification=c(1,0),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set2')

# Plot Percent identity and coverage
mumgp %<>% mutate(similarity=1-error/abs(qe-qs))
mumgp.filt %<>% mutate(similarity=1-error/abs(qe-qs))

ggplot(mumgp, aes(x=rs, xend=re, y=similarity, yend=similarity)) + geom_segment() +
  theme_bw() + xlab('reference sequence') + ylab('similarity') + ggtitle('All contigs') +
  ylim(0,1)

ggplot(mumgp.filt, aes(x=rs, xend=re, y=similarity, yend=similarity)) + geom_segment() +
  theme_bw() + xlab('reference sequence') + ylab('similarity') +
  ggtitle('At least 10 Kbp aligned') + ylim(0,1)

maxSimilarityDisjoin <- function(df){
  ref.ir = GRanges('X', IRanges(df$rs, df$re), similarity=df$similarity)
  # Efficient clean up of low similarity within high similarity
  step = 1
  while(step>0){
    largealign = ref.ir[head(order(rank(-ref.ir$similarity), rank(-width(ref.ir))),step*1000)]
    ol = findOverlaps(ref.ir, largealign, type='within') %>% as.data.frame %>%
      mutate(simW=ref.ir$similarity[queryHits],
             simL=largealign$similarity[subjectHits]) %>% filter(simW<simL)
    if(length(largealign) == length(ref.ir)){
      step = 0
    } else {
      step = step + 1
    }
    ref.ir = ref.ir[-ol$queryHits]
  }
  # Disjoin and annotate with the max similarity
  ref.dj = disjoin(c(ref.ir, GRanges('X', IRanges(min(df$rs), max(df$rs)), similarity=0)))
  ol = findOverlaps(ref.ir, ref.dj) %>% as.data.frame %>%
    mutate(similarity=ref.ir$similarity[queryHits]) %>%
    group_by(subjectHits) %>% summarize(similarity=max(similarity))
  ref.dj$similarity = 0
  ref.dj$similarity[ol$subjectHits] = ol$similarity
  as.data.frame(ref.dj)
}

mumgp.sim = maxSimilarityDisjoin(mumgp)

mumgp.sim %>% select(similarity, start, end) %>% gather(end, pos, 2:3) %>%
  ggplot() + geom_line(aes(x=pos, y=similarity), alpha=.5, color='red') + theme_bw() +
  xlab('reference sequence') + ylab('similarity') + ggtitle('All contigs') + ylim(0,1) +
  geom_segment(aes(x=rs, xend=re, y=similarity, yend=similarity), data=mumgp)

ggplot(mumgp.sim) + geom_segment(aes(x=start, xend=end, yend=similarity, y=similarity),
                                 color='red', size=2) +
  theme_bw() + xlab('reference sequence') + ylab('similarity') + ylim(0,1) +
  geom_segment(aes(x=rs, xend=re, y=similarity, yend=similarity), data=mumgp)

mumgp.filt.m %>% filter(similarity==0) %>%
  ggplot(aes(x=start, xend=end, y=filter, yend=filter)) + geom_segment(size=10) +
  theme_bw() + xlab('reference sequence') + ylab('filter') +
  scale_colour_brewer(palette='Set1') + ggtitle('Reference regions not covered')


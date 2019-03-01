---
title: 'Badread: simulation of error-prone long reads'
tags:
  - long-read sequencing
  - oxford nanopore
  - ont
  - pacific biosciences
  - pacbio
authors:
  - name: Ryan R Wick
    orcid: 0000-0001-8349-0778
    affiliation: 1
affiliations:
 - name: Department of Infectious Diseases, Central Clinical School, Monash University, Melbourne, Victoria 3004, Australia
   index: 1
date: 1 March 2019
bibliography: paper.bib
---


# Background

DNA sequencing platforms aim to measure the sequence of nucleotides (A, C, G and T) in a sample of DNA. Sequencers made by Illumina have been the dominant technology for much of the past decade, but their platforms generate fragments of sequence ('reads') that are relatively small (~100–300 nucleotides in length). In contrast, Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio) produce 'long-read' sequencers that can generate sequence fragments with tens of thousands of nucleotides or more [@Eisenstein2017]. Long reads from these platforms can be very beneficial for genome assembly and other bioinformatic analyses [@Phillippy2017;@Koren2017]. ONT and PacBio sequencers achieve their long read lengths because they detect nucleotides in individual molecules of DNA, a.k.a. single-molecule sequencing [@Heather2016]. However, the stochastic nature of measuring at the single-molecule scale means that ONT and PacBio reads are 'noisy' – they contain a significant amount of errors.

Since sequencing reads from ONT and PacBio platforms are qualitatively different from Illumina reads (long and noisy vs short and accurate), they often require novel methods of analysis. The last few years have seen much research in this space, and one useful technique for evaluating new methods is read simulation: generating fake sequencing reads from a reference nucleotide sequence [@Huang2012]. This approach has some key advantages over using real sequencing data: it can be faster, more affordable and allow for a greater number of tests. Additionally, when using simulated reads, the reference nucleotide sequence provides a confident ground truth which may not be otherwise available.



# Summary

Here we introduce Badread, a software tool for _in silico_ simulation of long reads. Its primary aim is to generate simulated read sets for the purpose of evaluating tools or methods that take long reads as input. Badread differs from existing tools (e.g. PBSIM [@Ono2013], LongISLND [@Mu2016] and NanoSim [@Yang2017]) in two key ways. First, it can simulate types of read errors that other tools cannot. While other long read simulation tools focus on modelling read length and sequencing errors, Badread can additionally include chimeras (when a single read which consists of two or more non-contiguous sequences), adapters (additional sequences from the library preparation at the start or end of a read), glitches (localised regions of low accuracy) and junk reads (low-complexity repetitive sequences).

The second way Badread differs from existing tools is that it prioritises control over realism. Using read length as an example, other long read simulation tools may sample read lengths from a real read set, so their simulated reads follow a realistic distribution. Badread instead uses a gamma distribution for read lengths where the user specifies the mean and standard deviation – less realistic but highly tuneable. Users can therefore generate many read sets which quantitatively vary, e.g. mean lengths of 1000, 2000, 3000, etc. Other characteristics of the read set (read accuracy, chimera rate, glitch rate, etc.) can be similarly tuned in Badread, allowing users to systematically evaluate how they affect the performance of a tool or method.



# Availability

Badread is open-source and available via the GPLv3 license at [github.com/rrwick/Badread](https://github.com/rrwick/Badread).



# References

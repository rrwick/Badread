# Comparison to real reads and other simulators

* [Simulators tested](#simulators-tested)
* [Real read set](#real-read-set)
* [Training commands](#training-commands)
* [Simulation commands](#simulation-commands)
* [Results](#results)
    * [Read lengths](#read-lengths)
    * [Read identities](#read-identities)
    * [Identities across reads](#identities-across-reads)
    * [Qscore distribution](#qscore-distribution)
    * [Read identity vs mean qscore](#read-identity-vs-mean-qscore)
    * [Consensus accuracy](#consensus-accuracy)
* [Conclusions](#conclusions)




## Simulators tested

This comparison includes the following five long-read simulators:
* [Badread](https://github.com/rrwick/Badread) v0.1.0
* [NanoSim](https://github.com/bcgsc/nanosim) v2.1.0
* [PBSIM](https://github.com/pfaucon/PBSIM-PacBio-Simulator) v1.0.3
* [SiLiCO](https://github.com/ethanagbaker/SiLiCO) v1.0.1
* [LongISLND](https://github.com/bioinform/longislnd) v0.9.5

I was also interested in [loresim2](https://github.com/gt1/loresim2), but couldn't get it to compile.




## Real read set

I used a Nanopore MinION run on [this sample](https://www.ncbi.nlm.nih.gov/biosample/SAMEA3356997) as my real read set to compare against. It's a _Klebsiella pneumoniae_ genome about 5.5 Mbp in size. There are 115,042 reads which have a total size of 1.1 Gbp (200x depth) and an N50 read length of 21.1 kbp. The completed genome (a hybrid Illumina+Nanopore assembly) is `ref.fasta` in the commands below.

This is a slightly older sample from early 2017, so I redid the basecalling with a more recent version of Albacore (v2.3.1). I chose this sample because it was not part of a barcoded run – only the one genome was sequenced in its run. This means we should see the full spectrum of read quality, unlike for a barcoded sample where the demultiplexing process removes many of the lower quality reads.




## Training commands

Some read simulation tools can be trained for a more realistic model. Wherever possible, I used the real reads to train the simulators.

#### Badread:
```
minimap2 -t 32 -c -x map-ont ref.fasta real_reads.fastq.gz > alignments.paf
badread error_model --reference ref.fasta --reads real_reads.fastq.gz --alignment alignments.paf > badread_errors
badread qscore_model --reference ref.fasta --reads real_reads.fastq.gz --alignment alignments.paf > badread_qscores
```

#### NanoSim:
```
seqtk seq -A real_reads.fastq.gz > real_reads.fasta
read_analysis.py -i real_reads.fasta -r ref.fasta
```

#### LongISLND:
```
gunzip -c real_reads.fastq.gz > real_reads.fastq
minimap2 -a -x map-ont -t 16 ref.fasta real_reads.fastq.gz | samtools view -hu - | samtools sort -@ 4 -o real_reads.fastq.bam -
samtools index real_reads.fastq.bam
samtools faidx ref.fasta
sample.py --input_suffix fastq.bam --read_type fastq --model_dir longislnd_model --reference ref.fasta
```

PBSIM is also uses a real read set, but in the simulation command (see below) not in a separate training step.



## Simulation commands

While the sample genome has two plasmids, I chose to use only the chromosome for a simulation reference. This kept things simpler and allowed me to use NanoSim in [circular mode](https://github.com/bcgsc/NanoSim#2-simulation-stage):
```
seqtk seq ref.fasta | head -n 2 > ref_chromosome.fasta
```

I mostly used each tool's default parameters when simulating reads. The specific commands used are below:

#### Badread
``` 
badread simulate --quantity 100x --reference ref_chromosome.fasta --error_model badread_errors --qscore_model badread_qscores | gzip > badread_reads.fastq.gz
```

#### NanoSim
```
simulator.py circular -r ref_chromosome.fasta -c training -n 55000
gzip -c simulated_reads.fasta > nanosim_reads.fasta.gz
rm simulated_reads.fasta simulated.log simulated_error_profile
```

#### PBSIM
```
gunzip -c real_reads.fastq.gz > real_reads.fastq
pbsim --data-type CLR --depth 100 --sample-fastq real_reads.fastq --length-max 1000000 --accuracy-min 0 ref_chromosome.fasta
gzip -c sd_0001.fastq > pbsim_reads.fastq.gz
rm real_reads.fastq sd_0001.ref sd_0001.maf sd_0001.fastq
```

#### SiLiCO
```
echo ">chr1" > silico_ref.fasta
tail -n 1 ref_chromosome.fasta >> silico_ref.fasta
SiLiCO.py -i silico_ref.fasta -o silico_out -c 100 --nanopore --fasta
gzip -c silico_out/simulated_read_positions_trial_1.fa > silico_reads.fasta.gz
rm -r silico_out silico_ref.fasta
```

#### LongISLND
```
simulate.py --movie_id ONT --read_type fastq --model_dir longislnd_model --fasta ref_chromosome.fasta --coverage 100
cat out/*.fq | gzip > longislnd_reads.fastq.gz
rm -r out
```




## Results

I've analysed the read sets to qualitatively compare them against each other and against the real read set (in the top left of each plot). For each analysis, a good simulated read set is one that closely mimics the nature of the real read set.


### Read lengths

<p align="center"><img src="images/read_lengths.png" alt="Read lengths" width="85%"></p>

These plots show the distribution of read lengths in histogram form. The values are weighted by the read length, so this is actually a distribution of bases in reads (e.g. a 40 kbp read contributes 10x as much as a 4 kbp read). I find this to be more useful than the simpler distribution of read counts, as longer reads really are more sequencing than short reads. The read N50 is the median of these distributions.

NanoSim, PBSIM and LongISLND all sample from the real read length distribution, so their read lengths match the real reads very well. Badread generates its read lengths from a gamma distribution, so it's not quite as realistic. SiLiCO's reads are all close to 10 kbp in size.



### Read identities

<p align="center"><img src="images/read_identities.png" alt="Read identities" width="85%"></p>

Read identity was determined by aligning each read back to the reference genome from which it was made. Both NanoSim and LongISLND produce narrow identity distributions – i.e. all reads are similarly accurate. PBSIM's distribution is wider but still lacks many reads below about 75% identity. Only Badread covers the entire identity spectrum, with plenty of reads down to 70% identity and below. It's mean identity is actually worse than the real read set, but after all it is called <em>Bad</em>read. SiLiCO does not add any errors to its reads, which seriously limits its usefulness as a long-read simulator.



### Identities across reads

These plots may take a bit more explaining, so here are two examples:

<p align="center"><img src="images/window_identities_examples.png" alt="Identities across reads" width="55%"></p>

These show the read identity over a 250 bp sliding window across 10 kbp regions of a read. Both of the above examples have an overall identity of 87.5%, but the read on the left is more even – every one of its 250 bp windows stays within 80–92.5%. The read on the right is more variable with higher highs and lower lows.

With those examples in mind, here are plots showing 10 reads from each set. Again, these reads all have an average identity of 87.5% (I chose that value because all simulators, excluding SiLiCO, produced reads in that range):

<p align="center"><img src="images/window_identities.png" alt="Identities across reads" width="85%"></p>

The real reads have a moderate amount of variance from their mean identity. NanoSim looks pretty good here, largely matching the pattern of the real reads. Both PBSIM and LongISLND make reads that are a bit too even.

The real reads contain a few 'glitches', junk regions which case the window identity drastically fall. Badread is the only simulator which also includes glitches. It actually includes far more glitches than the real reads, but again, that's why it's called <em>Bad</em>read.



### Qscore distribution

<p align="center"><img src="images/qscores.png" alt="QScores" width="85%"></p>

NanoSim and SiLiCO output reads in FASTA format (which lack qscores) and so are excluded from this plot and the following plot. The other simulators all model their qscores using the real reads, and each achieves a very realistic looking distribution.



### Read identity vs mean qscore

<p align="center"><img src="images/qscore_identity.png" alt="Identity vs mean qscore" width="85%"></p>

Each point in these plots represents one read. Mean qscore should be positively correlated with read identity, i.e. qscores are informative about how good a read actually is. While this is true for the real reads, they seem to have a bimodal distribution with some reads appearing higher on the plot (more accurate than their qscore suggests) – I'm not sure what caused this behaviour.

PBSIM and Badread do best here, but neither quite matches the real read set. PBSIM's mean qscores tend to be too low, while Badread's mean qscore don't span the full range. Both PBSIM and Badread have unrealistically tight correlations. LongISLND's narrow identity distribution makes it difficult to assess, but the correlation looks weak.



### Consensus accuracy

<p align="center"><img src="images/consensus_accuracy.png" alt="Consensus accuracy" width="85%"></p>

These plots show the assembly accuracy (specifically, a [Rebaler](https://github.com/rrwick/Rebaler) assembly which uses Racon to polish) as a function of read depth. Note that the real reads never achieve a particular high accuracy due to systematic basecalling errors (mainly getting homopolymer lengths wrong and botching locations with modified bases).

Read simulators which approach 100% assembly accuracy (NanoSim and PBSIM) are not modelling these systematic errors very well. LongISLND's assemblies have the lowest accuracy of the simulated reads, suggesting that it is modelling the read errors. Badread is middle-of-the-pack in this analysis, approaching 100% accuracy but at a higher depth.



## Conclusions

None of the simulators excel in every category – i.e. each produces reads which could be distinguished from real reads. That being said, PBSIM is a very good all-rounder, despite being the oldest tool tested (developed way back in 2012), and was also very fast to run. NanoSim is best at varying identity across the length of a read. LongISLND seems to be best at modelling systematic basecalling errors. Badread (true to its name) is best at modelling many read problems, like glitches and low identities.

I didn't make Badread with the intention of simulating reads as realistically as possible. Rather, Badread strives to give users control over various long-read faults, so they can torture test tools which take long reads as input. That being said, a few tweaks to its parameters (like tightening up the read identity distribution) could make it pretty good at imitating real Nanopore reads.

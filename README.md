<p align="center"><picture><source srcset="images/logo-dark.png" media="(prefers-color-scheme: dark)"><img src="images/logo.png" alt="Badread" width="75%"></picture></p>

Badread is a long-read simulator tool that makes – you guessed it – bad reads! It can imitate many kinds of problems one might encounter in real long-read sets: chimeras, low-quality regions, systematic basecalling errors and more.

Badread does not try to be best at imitating real reads (though it's not too bad, see [this comparison between Badread and other long-read simulators](comparison)). Rather, it was intended to give users _control_ over the quality of its simulated reads. I made Badread for the purpose of testing tools which take long reads as input. With it, one can increase the rate of different types of read problems, to see what effect it has.

Badread is published in the [Journal of Open Source Software](http://joss.theoj.org). If you use it in your research, please cite this manuscript:<br>
> Wick RR. Badread: simulation of error-prone long reads. _Journal of Open Source Software_. 2019;4(36):1316. doi:[10.21105/joss.01316](https://doi.org/10.21105/joss.01316).

[![Build Status](https://travis-ci.com/rrwick/Badread.svg?token=71gNPkycVbFoEsJC4qcj&branch=master)](https://travis-ci.com/rrwick/Badread) [![Coverage Status](https://coveralls.io/repos/github/rrwick/Badread/badge.svg?branch=master&service=github)](https://coveralls.io/github/rrwick/Badread?branch=master) [![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![status](http://joss.theoj.org/papers/a9ae38312391668801949ddc59e10cb1/status.svg)](http://joss.theoj.org/papers/a9ae38312391668801949ddc59e10cb1) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2622997.svg)](https://doi.org/10.5281/zenodo.2622997)


## Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Quick usage](#quick-usage)
* [Method](#method)
* [Detailed usage](#detailed-usage)
  * [Command line](#command-line)
  * [Reference FASTA](#reference-fasta)
  * [Fragment lengths](#fragment-lengths)
  * [Read identities](#read-identities)
  * [Error model](#error-model)
  * [QScore model](#qscore-model)
  * [Adapters](#adapters)
  * [Junk and random reads](#junk-and-random-reads)
  * [Chimeras](#chimeras)
  * [Small plasmid bias](#small-plasmid-bias)
  * [Glitches](#glitches)
* [Contributing](#contributing)
* [License](#license)




## Requirements

Badread runs on MacOS and Linux. It may not work natively on Windows (I haven't tried) but can be run using the [Windows subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10). It requires [Python](https://www.python.org/) 3.6 or later.

To install Badread you'll need [pip](https://pypi.org/project/pip/) and [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git). It also uses a few Python packages ([Edlib](https://github.com/Martinsos/edlib/tree/master/bindings/python), [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/) and [Matplotlib](https://matplotlib.org/)) but these should be taken care of by the installation process.




## Installation

### Install from source

You can install Badread using pip, either from a local copy:
```bash
git clone https://github.com/rrwick/Badread.git
pip3 install ./Badread
badread --help
```

Or directly from GitHub:
```bash
pip3 install git+https://github.com/rrwick/Badread.git
badread --help
```

If these installation commands aren't working for you (e.g. an error message like `Command 'pip3' not found` or `command 'gcc' failed with exit status 1`) then check out the [installation issues page on the wiki](https://github.com/rrwick/Badread/wiki/Installation-issues).


### Run without installation

Badread can also be run directly from its repository by using the `badread-runner.py` script, no installation required:

```bash
git clone https://github.com/rrwick/Badread.git
Badread/badread-runner.py -h
```

If you run Badread this way, it's up to you to make sure that all [necessary Python packages](#requirements) are installed.




## Quick usage

If you need a reference genome to try out Badread, you can download [this file](https://bit.ly/2Ejwr6V) which is an assembly of the [_Klebsiella pneumoniae_ SGH10](https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN06112188) genome – a nasty hypervirulent strain ([read more about it here](https://www.nature.com/articles/s41467-018-05114-7)).

Badread's default settings correspond to Oxford Nanopore R10.4.1 reads of mediocre quality:
```bash
badread simulate --reference ref.fasta --quantity 50x \
    | gzip > reads.fastq.gz
```

To simulate older Oxford Nanopore reads (R9.4.1, worse basecalling):
```bash
badread simulate --reference ref.fasta --quantity 50x \
    --error_model nanopore2020 --qscore_model nanopore2020 --identity 90,98,5 \
    | gzip > reads.fastq.gz
```

To simulate PacBio HiFi reads:
```bash
badread simulate --reference ref.fasta --quantity 50x \
    --error_model pacbio2021 --qscore_model pacbio2021 --identity 30,3 \
    | gzip > reads.fastq.gz
```

Very bad reads:
```bash
badread simulate --reference ref.fasta --quantity 50x --glitches 1000,100,100 \
    --junk_reads 5 --random_reads 5 --chimeras 10 --identity 80,90,6 --length 4000,2000 \
    | gzip > reads.fastq.gz
```

Pretty good reads:
```bash
badread simulate --reference ref.fasta --quantity 50x --glitches 10000,10,10 \
    --junk_reads 0.1 --random_reads 0.1 --chimeras 0.1 --identity 20,3 \
    | gzip > reads.fastq.gz
```

Very good reads:
```bash
badread simulate --reference ref.fasta --quantity 50x --error_model random \
    --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
    --identity 30,3 --length 40000,20000 --start_adapter_seq "" --end_adapter_seq "" \
    | gzip > reads.fastq.gz
```




## Method

Badread simulates reads by mimicking the process of real sequencing: breaking the DNA into fragments, adding adapters and then reading the fragments into nucleotide sequences.

Here is an overview of how Badread makes each of its reads:

1. Use the [fragment length distribution](#fragment-lengths) to choose a length for the read.

2. Choose a type of fragment:
    * Most will be fragments of sequence from the [reference FASTA](#reference-fasta). These are equally likely to come from either strand, and can loop around circular references. If there are multiple reference sequences with different depths, then the likelihood of the fragment coming from each sequence is proportional to that sequence's depth.
    * Depending on the settings, some fragments may also be [junk or random sequence](#junk-and-random-reads).

3. Add adapter sequences to the start and end of the fragment, based on the [adapter settings](#adapters).

4. As determined by the [chimera rate](#chimeras), there is a chance that Badread will make another fragment and concatenate it onto the current fragment (possibly with adapter sequences in between, possibly not).

5. Add glitches to the fragment, based on the [glitch settings](#glitches).

6. Choose a percent identity for the read using the [read identity distribution](#read-identities).

7. 'Sequence' the fragment by adding errors until it has the target percent identity.
    * Errors are chosen using the [error model](#error-model) and are added at random positions in the read.
    * This step performs periodic alignments between the original fragment and the error-added sequence, so Badread can track the read's actual identity. This allow it to be precise (if Badread is aiming for a 91.5% identity read, it will be very close to 91.5% identity) but slow. If you find that Badread is too slow, check out [the wiki page on running it in parallel](https://github.com/rrwick/Badread/wiki/Running-in-parallel).

8. Generate quality scores for each base using the [qscore model](#error-model).

9. Output the read and quality in FASTQ format.




## Detailed usage

### Command line

```
usage: badread simulate --reference REFERENCE --quantity QUANTITY [--length LENGTH]
                        [--identity IDENTITY] [--error_model ERROR_MODEL]
                        [--qscore_model QSCORE_MODEL] [--seed SEED] [--start_adapter START_ADAPTER]
                        [--end_adapter END_ADAPTER] [--start_adapter_seq START_ADAPTER_SEQ]
                        [--end_adapter_seq END_ADAPTER_SEQ] [--junk_reads JUNK_READS]
                        [--random_reads RANDOM_READS] [--chimeras CHIMERAS] [--glitches GLITCHES]
                        [--small_plasmid_bias] [-h] [--version]

Generate fake long reads

Required arguments:
  --reference REFERENCE           Reference FASTA file (can be gzipped)
  --quantity QUANTITY             Either an absolute value (e.g. 250M) or a relative depth (e.g. 25x)

Simulation parameters:
  Length and identity and error distributions

  --length LENGTH                 Fragment length distribution (mean and stdev, default: 15000,13000)
  --identity IDENTITY             Sequencing identity distribution (mean,max,stdev for beta
                                  distribution or mean,stdev for normal qscore distribution, default:
                                  95,99,2.5)
  --error_model ERROR_MODEL       Can be "nanopore2018", "nanopore2020", "nanopore2023", "pacbio2016",
                                  "pacbio2021", "random" or a model filename (default: nanopore2023)
  --qscore_model QSCORE_MODEL     Can be "nanopore2018", "nanopore2020", "nanopore2023", "pacbio2016",
                                  "pacbio2021", "random", "ideal" or a model filename (default:
                                  nanopore2023)
  --seed SEED                     Random number generator seed for deterministic output (default:
                                  different output each time)

Adapters:
  Controls adapter sequences on the start and end of reads

  --start_adapter START_ADAPTER   Adapter parameters for read starts (rate and amount, default: 90,60)
  --end_adapter END_ADAPTER       Adapter parameters for read ends (rate and amount, default: 50,20)
  --start_adapter_seq START_ADAPTER_SEQ
                                  Adapter sequence for read starts (default:
                                  AATGTACTTCGTTCAGTTACGTATTGCT)
  --end_adapter_seq END_ADAPTER_SEQ
                                  Adapter sequence for read ends (default: GCAATACGTAACTGAACGAAGT)

Problems:
  Ways reads can go wrong

  --junk_reads JUNK_READS         This percentage of reads will be low-complexity junk (default: 1)
  --random_reads RANDOM_READS     This percentage of reads will be random sequence (default: 1)
  --chimeras CHIMERAS             Percentage at which separate fragments join together (default: 1)
  --glitches GLITCHES             Read glitch parameters (rate, size and skip, default: 10000,25,25)
  --small_plasmid_bias            If set, then small circular plasmids are lost when the fragment
                                  length is too high (default: small plasmids are included regardless
                                  of fragment length)

Other:
  -h, --help                      Show this help message and exit
  --version                       Show program's version number and exit
```



### Reference FASTA

The reference genome must be given as a FASTA file (either gzipped or not) using the `--reference` argument.

Each sequence's depth can be specified in the FASTA header, e.g. using `depth=1.1` or `depth=15`. Badread will use this to determine the relative abundance of each sequence. This can be useful for both bacterial genomes (where plasmids may be higher depth than the chromosome) and eukaryote genomes (where chloroplast/mitochondrial genomes may be higher depth than the rest of the genome).

Circular sequences are indicated by including `circular=true` in the FASTA header. This allows reads to loop past the end and back to the start of the sequence.

For a couple of examples, check out [the reference FASTA page on the wiki](https://github.com/rrwick/Badread/wiki/Example-reference-FASTAs).



### Fragment lengths

Badread generates fragment lengths from a [gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution). While a gamma distribution is usually parameterised with shape and scale (_k_ and _θ_) or shape and rate (_α_ and _β_), I don't find these particularly intuitive. So Badread instead defines fragment lengths using mean and standard deviation.

There are two ways to think about fragment lengths: the distribution of the fragment lengths and the distribution of _bases_ in the fragments. The latter distribution is higher because larger fragments contribute more bases. The read N50 is the median of the base distribution – half the bases will be in reads shorter than this and half in longer reads.

Badread's default is `--length 15000,13000` (mean=15000, stdev=13000) which corresponds to a decent Nanopore run (N50=22.6 kbp). To see the equations and interactively explore how different parameters affect the distributions, check out [this Desmos plot](https://www.desmos.com/calculator/z2a3yqssie).

Note that these parameters control the length of the _fragments_, not the final _reads_. These differ because: adapters are added to fragments, glitches can lengthen/shorten fragments, adding read errors can change the length (especially if the error model is biased towards insertions or deletions) and chimeras are made by concatenating multiple fragments together.



### Read identities

Badread can generate read identities in two alternative ways.

The first way uses a [beta distribution](https://en.wikipedia.org/wiki/Beta_distribution) to sample read identities. Like with fragment lengths, Badread defines the distribution with mean and standard deviation instead of using the more formal (and less intuitive) shape parameters (_α_ and _β_). In addition, Badread scales the distribution down using a maximum value, for a total of three parameters. To use this method, give three comma-delimited values (identity mean, max, stdev), e.g. `--identity 95,99,2.5`.

The second way uses a [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution) to sample read [qscores](https://en.wikipedia.org/wiki/Phred_quality_score). To use this method, give two comma-delimited values (qscore mean and stdev), e.g. `--identity 20,4`. This approach is better suited to high accuracy reads like ONT duplex or PacBio HiFi.

Badread's default is `--identity 95,99,2.5` which corresponds to an okay (but not great) R10.4.1 Nanopore sequencing run. To see the equations and interactively explore how different parameters affect the distribution, check out these Desmos plots:
* [three-parameter identity beta distribution](https://www.desmos.com/calculator/u9fivqmisa)
* [two-parameter qscore normal distribution (qscore on x-axis)](https://www.desmos.com/calculator/kpganphxj0)
* [two-parameter qscore normal distribution (identity on x-axis)](https://www.desmos.com/calculator/kg9s0dvbup)

For detail on how Badread defines identity, check out [this page on the wiki](https://github.com/rrwick/Badread/wiki/Definition-of-identity).



### Error model

The possible values for the `--error_model` argument are:
* `nanopore2023`: a model trained on ONT R10.4.1 reads from 2023 (the default)
* `nanopore2020`: a model trained on ONT R9.4.1 reads from 2020
* `nanopore2018`: a model trained on ONT R9.4/R9.4.1 reads from 2018
* `pacbio2016`: a model trained on PacBio RS II reads from 2016
* `pacbio2021`: a model trained on PacBio Sequel IIe HiFi reads from 2021
* `random`: a random error model with 1/3 chance each of insertion, deletion and substitution
* a file path for a trained model

For more information on how error models work, see [this page on the wiki](https://github.com/rrwick/Badread/wiki/Error-models). For instructions on building your own error model, see [this page](https://github.com/rrwick/Badread/wiki/Generating-error-and-qscore-models).



### QScore model

The possible values for the `--qscore_model` argument are:
* `nanopore2023`: a model trained on ONT R10.4.1 reads from 2023 (the default)
* `nanopore2020`: a model trained on ONT R9.4.1 reads from 2020
* `nanopore2018`: a model trained on ONT R9.4/R9.4.1 reads from 2018
* `pacbio2016`: a model trained on PacBio RS II reads from 2016
* `pacbio2021`: a model trained on PacBio Sequel IIe HiFi reads from 2021
* `random`: a model where qscores are meaningless and give no indication of read/base quality
* `ideal`: a model where scores are unrealistically informative about read/base quality
* a file path for a trained model

For more information on how qscore models work, see [this page on the wiki](https://github.com/rrwick/Badread/wiki/QScore-models). For instructions on building your own qscore model, see [this page](https://github.com/rrwick/Badread/wiki/Generating-error-and-qscore-models).



### Adapters

Adapter sequences are controlled with the `--start_adapter_seq` and `--end_adapter_seq` options. The default adapters are those for the Nanopore ligation adapters. To see what those sequences are and some alternatives, check out the [adapter sequences page on the wiki](https://github.com/rrwick/Badread/wiki/Adapter-sequences). If you supply numbers for the adapter sequences (e.g. `--start_adapter_seq 20`, then Badread will make a random sequence of that length to be the adapter.

How much adapter is added to the start/end of a read is controlled by two parameters: rate and amount. These are set using the `--start_adapter` and `--end_adapter` options, with each taking two comma-delimited numbers (rate first and amount second). Rate is the percent chance that the adapter will appear at all. E.g. a start-adapter rate of 90 means that 10% of reads will have no adapter at their start and 90% of reads will have some adapter at their start. Think of it like a Bernoulli distribution.

Amount controls how much of the adapter, on average, appears on the read. E.g. a start-adapter amount of 60 means that when an adapter is on the start of the read, it has an expected length of 60% its full length. Start-adapters are truncated at the start and end-adapters are truncated at the end. The amount of adapter for each read is controlled by a beta distribution, and you can interactively explore different values using [this Desmos plot](https://www.desmos.com/calculator/qzza86553k).

To turn off adapters entirely, set the sequences to nothing:<br>
`--start_adapter_seq "" --end_adapter_seq ""`




### Junk and random reads

Badread can add two types of completely wrong reads to the output: junk and random. Junk reads are low-complexity sequence and simulate something wrong with the sequencer. A junk fragment might look like `AATAATAATAATAATAATAAT` (and so on). Random reads are made of random sequence (25% chance of each base at each position). By default, Badread includes 1% of each type, so 2% of the total reads will be junk or random.



### Chimeras

Chimeric reads can occur in real datasets for two possible reasons: 1) fragments of DNA were actually ligated together before sequencing (probably more common in library preps that use ligase), and 2) two or more reads were sequenced in quick succession such that the sequencing software didn't recognise them as separate (an _in silico_ chimera). Chimeras can occur on both PacBio and Nanopore sequencing platforms.

As an example, imagine you used `--chimeras 2` to set chimeras to 2%. After making a sequence fragment, Badread then has a 2% chance of making another fragment and concatenating on to the first. It then has a 2% chance of concatenating on yet another fragment, and so on. This means that about 2% of the reads will be chimeras of two or more fragments, (2%)<sup>2</sup> will be chimeras of three or more, (2%)<sup>3</sup> will be chimeras of four or more, and so on. If you are using start/end adapters, they will sometimes (but not always) be added between the fragments.



### Small plasmid bias

Small circular plasmids can be underrepresented in long read sequencing – a topic addressed in [Completing bacterial genome assemblies with multiplex MinION sequencing](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000132). In short, it is necessary to avoid excessive DNA shearing in order to achieve long read lengths. But (at least for ligation-based preps) an unsheared circular plasmid can't be sequenced because it has no blunt ends for adapter ligation. So plasmids on the low end of the read-length distribution will be sequenced at lower-than-expected depth.

Badread simulates this effect if you use the `--small_plasmid_bias` option. When turned on, Badread tosses out fragments that land in a circular reference which is smaller than the fragment size. The degree of bias is therefore strongly dependent on the fragment length distribution (specifically how much of the distribution is less than the plasmid's size). Note that this only affects _circular_ sequences – linear sequences are unaffected.



### Glitches

Glitches are points in the read where the sequence is briefly messed up. They are controlled by three parameters:
* rate: how often glitches occur
* size: how much random sequence is added to the read
* skip: how much read sequence is lost

These are specified with the `--glitches` option by giving all three parameters in a comma-delimited list (no spaces). E.g. `--glitches 5000,100,100`. Each of these parameters is a mean for a [geometric random variable](https://en.wikipedia.org/wiki/Geometric_distribution). E.g. a glitch rate of 1000 doesn't mean glitches evenly occur at 1000 bp intervals, it means glitches are _on average_ 1000 bp apart. Turn glitches off entirely with `--glitches 0,0,0`.

Take a look at the [glitches page on the wiki](https://github.com/rrwick/Badread/wiki/Glitches) to see some dotplots which illustrate the concept.



## Contributing

If you are interested in contributing to Badread, please take a look at the [contribution guidelines](CONTRIBUTING.md).



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)

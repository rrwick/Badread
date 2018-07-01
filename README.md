<p align="center"><img src="images/logo.png" alt="Badread" width="75%"></p>

Badread is a long read simulator tool that makes – you guessed it – bad reads! It can imitate many kinds of problems one might encounter in real read sets: chimeric reads, low-quality regions, systematic basecalling errors and more.

Badread is pretty good at generating realistic simulated long read sets, but its focus is more on providing users with _control_ over the simulated reads. I made Badread for the purpose of testing long read assemblers. With it, one can increase the rate of different types of read problems, to see how they affect assembly quality.




## Table of contents

  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Quick usage](#quick-usage)
  * [Method](#method)
  * [Detailed usage](#detailed-usage)
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
  * [Other tools](#other-tools)
     * [Generating an error model](#generating-an-error-model)
     * [Visualising error rates in reads](#visualising-error-rates-in-reads)
  * [License](#license)




## Requirements

Badread runs on MacOS and Linux. Its only dependencies are some Python packages ([Edlib](https://github.com/Martinsos/edlib/tree/master/bindings/python), [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/) and [Matplotlib](https://matplotlib.org/)) but these should be taken care of by pip when installing Badread.




## Installation

### Install from source

Running the `setup.py` script will install a `badread` executable:

```bash
git clone https://github.com/rrwick/Badread.git
cd Badread
python3 setup.py install
badread -h
```

* If the `python3 setup.py install` command complains about permissions, you may need to run it with `sudo`.
* Install just for your user: `python3 setup.py install --user`
    * If you get a strange 'can't combine user with prefix' error, read [this](http://stackoverflow.com/questions/4495120).
* Install to a specific location: `python3 setup.py install --prefix=$HOME/.local`
* Install with pip (local copy): `pip3 install path/to/Badread`
* Install with pip (from GitHub): `pip3 install git+https://github.com/rrwick/Badread.git`


### Run without installation

Badread can be run directly from its repository by using the `badread-runner.py` script, no installation required:

```bash
git clone https://github.com/rrwick/Badread.git
Badread/badread-runner.py -h
```

If you run Badread this way, it's up to you to make sure that all [necessary Python packages](#requirements) are installed.




## Quick usage

__Generate fake Nanopore reads at 50x depth:__
```
badread simulate --reference ref.fasta --quantity 50x --error_model nanopore --qscore_model nanopore | gzip > reads.fastq.gz
```

__Generate 200 Mbp of fake PacBio reads:__
```
badread simulate --reference ref.fasta --quantity 200M --error_model pacbio --qscore_model pacbio | gzip > reads.fastq.gz
```

__Generate very bad reads:__
```
badread simulate --reference ref.fasta --quantity 50x --error_model nanopore --qscore_model nanopore --glitches 1000,100,100 --junk_reads 5 --random_reads 5 --chimeras 10 --identity 75,90,8 | gzip > reads.fastq.gz
```

__Generate very nice reads:__
```
badread simulate --reference ref.fasta --quantity 50x --error_model random --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 --identity 95,100,4 --start_adapter 0,0 --end_adapter 0,0 | gzip > reads.fastq.gz
```



## Method

Badread simulates reads by roughly following the process of sequencing real DNA: breaking the DNA into fragments, adding adapters and then reading the fragments into nucleotide sequences.

Here is an overview of all the steps in Badread's read-generation:
* Choose a length for a sequence fragment using the [fragment length distribution](#fragment-lengths).
* Choose a type of fragment:
  * Most will be fragments of sequence from the [reference FASTA](#reference-fasta). These are equally likely to come from either strand, and can loop around circular references.
  * Depending on the settings, some fragments may also be [junk or random sequence](#junk-and-random-reads).
* Add adapter sequences to the start and end of the fragment, based on the [adapter settings](#adapters). 
* As determined by the [chimera setting](#chimeras), there is a chance that Badread will make another fragment and concatenate it onto the current fragment (possibly with adapter sequences in between, possibly not).
* Add glitches to the fragment, based on the [glitch settings](#glitches).
* Choose a percent identity for the read using the [read identity distribution](#read-identities).
* 'Sequence' the fragment by adding errors to the sequence until it has the target percent identity.
  * Errors are added at random positions, leading to a somewhat variable identity across the span of the read.
  * Errors are chosen using the [error model](#error-model).
* Generate quality scores for each base using the [qscore model](#error-model).
* Output the finished read, and repeat until the total volume of reads reaches the target amount.


## Detailed usage

### Reference FASTA

The reference genome must be given in FASTA format using the `--reference` argument.

Each sequence's depth can be specified in the FASTA header using `depth=1.1` or `depth=15`, etc. Badread will use this to determine the relative abundance of each sequence. This can be useful for both bacterial genomes (where plasmids may be higher depth than the chromosome) and eukaryote genomes (where chloroplast/mitochondrial genomes may be higher depth than the rest of the genome).

Circular sequences are indicated by including `circular=true` in the FASTA header. This allows reads to loop past the end and back to the start of the sequence.

An example bacterial genome reference FASTA might look like this:
```
>chromosome depth=1.0 circular=true
GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGCA
...
>plasmid depth=1.3 circular=true
ATGACGAGCGAAAATAACAGCTTACTTCTGAACCTTCAGGAAGTTGATAAGACAACCGGCGAAGTTGTTA
...
```

A eukaryote genome reference FASTA might look like this:
```
>chromosome_I depth=1.0
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACA
...
>chromosome_II depth=1.0
AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGATGTTCAACCAAAAGCTACTT
...
>chromosome_III depth=1.0
CCCACACACCACACCCACACCACACCCACACACCACACACACCACACCCACACACCCACACCACACCACA
...
...
...
>mitochondria depth=100.0 circular=true
TTCATAATTAATTTTTTATATATATATTATATTATAATATTAATTTATATTATAAAAATAATATTTATTA
...
```



### Fragment lengths

Badread generates fragment lengths from a [gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution). Instead of providing the gamma distributions shape and rate (which are not very intuitive paramaters), Badread defines the distribution using the mean and standard deviation.

Note that these parameters control the length of the _fragments_, not the final _reads_. These differ because: adapters are added to fragments, glitches can lengthen/shorten fragments, adding read errors can change the length (especially if the error model is biased towards insertions or deletions) and chimeras are made by concatenating multiple fragments together.

There are two ways to think about fragment lengths: the distribution of the fragment lengths and the distribution of the amount of bases in the fragments. The latter distribution is higher because larger fragments contribute more bases. The read N50 is the median of the base (red) distribution – half the bases will be in reads shorter than this and half in reads longer.

<table>
    <tr>
        <td>
            <img align="right" src="images/default_frag_lengths.png" alt="Default length distribution" width="400">
            Badread's default is <code>--length 15000,13000</code> (mean=15000, stdev=13000) which corresponds to a decent Nanopore run (N50=22.6 kbp). The fragment length distribution is in blue, while the base distribution is in red.
        </td>
    </tr>
</table>

You can interactively explore different values using [this Desmos plot](https://www.desmos.com/calculator/xrkqgzt4o5).



### Read identities

#### Definition of identity

Badread defines identity the same way as BLAST does: the number of matching bases over the length of the alignment. Take this example of a 24 bp read which originated from a 24 bp fragment of DNA. The read has 3 errors: one deletion, one substitution and one insertion. This read's identity is 22 / 25 = 88%. Note that the denominator is the not the length of the _read_ but rather the length of the _alignment_.

```
             Read:  ACGAC-CAGCAGTCGCGACTAGCTT
                    ||||| |||||| || |||||||||
Original sequence:  ACGACTCAGCAGACG-GACTAGCTT
```

Since DNA has only a 4 letter alphabet, two completely random sequences can typically align with >50% identity. As an example, here are two random sequences aligned to each other which match in 32 places over 59 alignment positions, giving an identity of 54%:
```
AAT-CGGCGCGTCCCGCGTTTCGGAAATTGA-C-ACTCTGACG-GTT---AGCACAG--
| | ||| | | |  || ||  ||   | || | | | ||| | |||   || | ||  
ATTACGG-GAG-C--GC-TTA-GGC--T-GAACTATTATGATGCGTTGCGAGAAAAGGA
```

This means that any read with less than about 60% identity is difficult to distinguish from random sequence.


#### Badread identity distributions

Badread generates read identities from a [beta distribution](https://en.wikipedia.org/wiki/Beta_distribution). There are three parameters: `mean,max,stdev`. Max sets the upper end of the distribution. Stdev controls the shape: smaller values create a tighter distribution around the mean, while larger values make a broader distribution.

<table>
    <tr>
        <td>
            <img align="right" src="images/default_identities.png" alt="Default identity distribution" width="400">
            Badread's default is <code>--identity 85,95,5</code> which corresponds to a decent Nanopore sequencing run.
        </td>
    </tr>
</table>

You can interactively explore different values using [this Desmos plot](https://www.desmos.com/calculator/q7qw6rq2lb).


### Error model

The possible values for the `--error_model` argument are:
* `nanopore`: a model trained on real Nanopore reads (the default)
* `pacbio`: a model trained on real PacBio reads
* `random`: a random error model with 1/3 chance each of insertion, deletion and substitution
* a filepath for a trained model

For more information, see the [Error models page on the wiki](https://github.com/rrwick/Badread/wiki/Error-models).


### QScore model

The possible values for the `--qscore_model` argument are:
* `nanopore`: a model trained on real Nanopore reads (the default)
* `pacbio`: a model trained on real PacBio reads
* `random`: a model where qscores are meaningless, give no indication of read or base quality
* `ideal`: a model where scores are unrealistically informative about read and base quality

For more information, see the [QScore models page on the wiki](https://github.com/rrwick/Badread/wiki/QScore-models).



### Adapters

Adapter sequences are controlled with the `--start_adapter_seq` and `--end_adapter_seq` options. The default adapters are those for the Nanopore ligation adapters. To see what those sequences are and some alternatives, check out the [Adapter sequences page on the wiki](https://github.com/rrwick/Badread/wiki/Adapter-sequences).

How much adapter is added to the start/end of a read is controlled by two parameters: rate and amount. Rate is the percent chance that the adapter will appear at all. E.g. a start-adapter rate of 90 means that 10% of reads will have no adapter at their start and 90% of reads will have some adapter at their start. Think of it like a Bernoulli distribution.

Amount controls how much of the adapter, on average, appears on the read. E.g. a start-adapter rate of 60 means that when an adapter is on the start of the read, it has an expected length of 60% its full length. Start-adapters are truncated at the start and end-adapters are truncated at the end. The actual distribution of amount is controlled by a beta distribution, and you can interactively explore different values using [this Desmos plot](https://www.desmos.com/calculator/qzza86553k).

To turn off adapters entirely, you can either set the sequences to nothing:<br>
`--start_adapter_seq "" --end_adapter_seq ""`<br>
Or set the rate/amount to 0:<br>
`--start_adapter 0,0 --end_adapter 0,0`



### Junk and random reads

Badread can add two types of completely wrong reads to the output: junk and random. Junk reads are low-complexity sequence and simulate something wrong with the sequencer. A junk fragment might look like `AATAATAATAATAATAATAATAATAAT` and so on. Random reads are made of random sequence and can serve to simulate external contamination. By default, Badread includes 1% of each type, so 2% of the total reads will be junk/random.



### Chimeras

Chimeric reads can occur in real datasets for two possible reasons: 1) fragments of DNA were actually ligated together before sequence (probably more common library preps that use ligase) and 2) more than one read was sequenced in quick succession such that the software didn't recognise them as separate (an in-silico chimera). They can occur on both PacBio and Nanopore sequencing platforms.

The `--chimeras` option controls the rate of fragment-joining in Badread. It takes a percentage, e.g. `--chimeras 5` means that 5% of reads will be made of multiple separate fragments.



### Small plasmid bias

Small circular plasmids can be underrepresented in long read sequencing (a topic addressed in [Completing bacterial genome assemblies with multiplex MinION sequencing](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000132)). The most straightforward explanation is that circular DNA is unavailable for sequencing because it has no free ends on which adapters can ligate.

Badread simulates this effect if you use the `--small_plasmid_bias` option. When turned on, Badread tosses out any fragment that lands in a circular reference which is smaller than the fragment size. The degree of bias is therefore strongly dependent on the fragment length distribution (specifically how much of the distribution is less than the plasmid's size). Note that this only affects _circular_ sequences – linear sequences are unaffected.



### Glitches

Glitches are points in the read where the sequence is briefly messed up. They are controlled by three parameters:
* rate: how often glitches occur
* size: how much random sequence is added to the read
* skip: how much read sequence is lost

These are specified with the `--glitches` option by giving all three parameters in a comma-delimited list (no spaces). E.g. `--glitches 5000,100,100`. Each of these parameters is a mean for a [geometric random variable](https://en.wikipedia.org/wiki/Geometric_distribution). E.g. a glitch rate of 1000 doesn't mean glitches occur at 1000 bp intervals, it means glitches are _on average_ 1000 bp apart. Turn glitches off with `--glitches 0,0,0`.

To better understand glitches, take a look at the [Glitches page on the wiki](https://github.com/rrwick/Badread/wiki/Glitches) to see some dotplots which illustrate the concept.



## Other tools

Badread comes with some other functionality, in addition to its main read-generation.


### Generating an error model

Badread comes with two k-mer-based errors: one that I built with Oxford Nanopore reads (MinION, R9.4 flowcell) and one that I built with PacBio reads (PacBio RS II, CLR). If you'd like to build your own model, keep reading!

Requirements:
* Long reads at least a Gbp would be good.
* A high quality reference FASTA. Ideally this is a high quality FASTA of the same genome as the reads came from. Illumina-polished genomes are probably best.
* [minimap2](https://github.com/lh3/minimap2) (my favourite long read aligner)

First, you must align your long reads to your reference. Make sure to use minimap2's `-c` option so it includes the CIGAR string in the output:
```
minimap2 -c -x map-ont reference.fasta.gz reads.fastq.gz | gzip > alignments.paf.gz
```

Now give the files to `badread model` to build the model:
```
badread model --reference reference.fasta.gz --reads reads.fastq.gz --alignment alignments.paf.gz > new_model
```

Options:
* `--k_size`: Smaller values result in smaller model files. Larger values give bigger files but more more specific error modelling. Larger k-mers will require more reads/alignments to build. Default: 7.
* `--max_alignments`: Limits the model-building to only using this many alignments. Useful if you have a lot of reads/alignments and are running out of RAM. Default: use all alignments.
* `--max_alt`: How many alternatives to save for each k-mer. Larger values give more specific error models but larger model files. Default: 25.



### Visualising error rates in reads

This tool was mainly for me to use when developing Badread, but I left it in case anyone else finds in useful. Use it to see a plot for each read's identity over a sliding window:
```
minimap2 -c -x map-ont reference.fasta.gz reads.fastq.gz | gzip > alignments.paf.gz
badread plot --reference reference.fasta.gz --reads reads.fastq.gz --alignment alignments.paf.gz
```




## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)

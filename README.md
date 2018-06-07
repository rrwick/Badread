<p align="center"><img src="images/logo.png" alt="Deepbinner" width="75%"></p>

Badread is a long read simulator tool that makes – you guessed it – bad reads! It can imitate many kinds of problems one might encounter in real read sets: chimeric reads, low-quality regions, systematic basecalling errors and more. You can use Badread to generate a realistic read set (with a normal amount of problems) or turn up the problems for a truly terrible read set. I made it for the purposes of testing long read assemblers. With it, one can see how different types of read problems affect assembly quality.




## Table of contents

  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Quick usage](#quick-usage)
  * [Simulating reads](#simulating-reads)
  * [Generating an error model](#generating-an-error-model)
  * [Generating an error model](#generating-an-error-model-1)
  * [Visualising error rates in reads](#visualising-error-rates-in-reads)
  * [License](#license)



## Requirements




## Installation




## Quick usage




## Simulating reads

### Fragment lengths

desmos.com/calculator/rddlqip1if

### Read identities

desmos.com/calculator/t03zr2thap

### Error model

### Adapters

### Junk and random reads

### Chimeras

### Small plasmid bias

### Glitches

Glitches are points in the read where the sequence is briefly messed up. They are controlled by three parameters:
* glitch rate: how often glitches occur
* glitch size: how much random sequence is added to the read
* glitch skip: how much read sequence is lost

These are specified with the `--glitches` option by giving all three parameters in a comma-delimited list (no spaces). E.g. `--glitches 5000,100,100`. Each of these parameters is a mean for a [geometric random variable](https://en.wikipedia.org/wiki/Geometric_distribution). E.g. a glitch rate of 1000 doesn't mean glitches occur at 1000 bp intervals, it means glitches are _on average_ 1000 bp apart.

To better understand glitches, let's look at some dotplots comparing a sequence before glitches are added (x axis) and after (y axis):

| 5000,100,100 | 1000,100,100 | 200,100,100 |
| --- | --- | --- |
| <img src="images/glitch_dotplots/5000-100-100.png" width="200" /> | <img src="images/glitch_dotplots/1000-100-100.png" width="200" /> | <img src="images/glitch_dotplots/200-100-100.png" width="200" /> |

The three dotplots above show 10 kbp reads where the glitch rate is varied. You can see that the smaller the value, the smaller the space between glitches.


| 1000,100,0 | 1000,0,100 | 1000,100,100 |
| --- | --- | --- |
| <img src="images/glitch_dotplots/1000-100-0.png" width="200" /> | <img src="images/glitch_dotplots/1000-0-100.png" width="200" /> | <img src="images/glitch_dotplots/1000-100-100.png" width="200" /> |


These dotplots show the effect of the size and skip parameters. The first has a size of 100 and a skip of 0. This means that each glitch adds random sequence but loses no real sequence, resulting in a longer glitched sequence. The second is the opposite, a size of 0 and a skip of 100, where the glitched sequence is shorter than the original. The third has size and skip both equal to 100, so the glitched sequence should on average be the same length as the original.

The examples shown above are particularly glitchy to illustrate the concept. Badread's default value is `5000,50,50` for a modest amount of glitches. You can use `0,0,0` to turn glitches off entirely:

| 5000,50,50 | 0,0,0 |
| --- | --- |
| <img src="images/glitch_dotplots/5000-50-50.png" width="200" /> | <img src="images/glitch_dotplots/0-0-0.png" width="200" /> |


## Generating an error model




## Generating an error model




## Visualising error rates in reads




## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)

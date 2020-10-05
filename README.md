# **A**lignment-free **p**hylogenetic **p**lacement algorithm based on **SPA**ced-word **M**atches

_**A**lignment-free **p**hylogenetic **p**lacement algorithm based on **SPA**ced-word **M**atches_ (_App-SpaM_) is a software for performing _phylogenetic placement_. Phylogenetic placement is the task of placing (usually short) query sequences of unknown taxonomic origin into an existing phylogeny of reference sequences. The input normally consists of three files:
- A `.fasta` file with reference sequences.
- A `.newick` file containing the reference phylogeny. The phylogeny comprises the exact same references as the `.fasta` file of references. 
- A `.fasta` file with query sequences that are to be placed in the reference phylogeny.

_App-SpaM_ will place each of the query sequences into the reference phylogeny at a phylogenetically appropriate position. The placement is based on the concept of _Spaced Word Matches_ (_FSWM_[^FSWM]). Depending on the chosen placement heuristic, _App-SpaM_ uses either the number of spaced word matches between query and references, or the estimated number of nucleotide substitutions per sequence position between query and references. The output is a _JPlace_[^JPlace] file containing all query placements.

We are planning to have _App-SpaM_ also available in _PEWO_, the _**P**lacement **E**valuation **WO**rkflows_[^PEWO], a tool developed to rapidly test and compare tools for phylogenetic placement.

When using _App-SpaM_, please cite:
_In preparation_

## Installation
### Prerequisites
If not already installed on your system, install Git and CMake; on Ubuntu e.g.:
```bash
sudo apt-get install git
sudo apt-get install cmake
```

At the moment, _App-SpaM_'s parallelization is performed via OpenMP. Most modern compilers are supporting OpenMP, but it is advisable to update your compiler to the newest version. If you experience problems during the compilation do not hesitate to contact us.

### Installing _App-SpaM_
On the command line, download the newest version of _App-SpaM_ (alternatively download and unpack the `.zip`-Archive from this Github page):
```bash
git clone https://github.com/matthiasblanke/_App-SpaM_
```
Navigate into the _App-SpaM_ directory, create a build folder, and navigate into it:
```bash
cd path/to/appspam
mkdir build
cd build
```
Build the program with CMake:
```bash
cmake ..
make
```

## Running _App-SpaM_
From within the build directory you can just run _App-SpaM_ like so:
```bash
./appspam -h
```
If you want to run it from anywhere, add it to your path:
```bash
export PATH=$PATH:~/path/to/appspam
```
If you want this change to be permanent add this line to your `~/.profile` or `~/.bash_profile`.

### Placing query reads
When running _App-SpaM_, at least you need to specify the three input file parameters: The reference sequences (`-s`), the reference phylogeny (`-t`) and the query sequences (`-q`):
```
./appspam -s path/to/references.fasta -t path/to/referencetree.nwk -q path/to/query.fasta
```
The paths can be either absolut paths, or relative to your current working directory. All other parameters will be set to default values. All output files will be placed in your current working directory. You can specify the output location and file name with the flag `-o`, e.g.:
```
./appspam -s references.fasta -t referencetree.nwk -q query.fasta -o path/to/output.jplace
```
If other output files are produced (see below) they will be placed in the same folder as the _JPlace_ file.

### Using unassembled references
Will be updated soon.

### Parameters
There are several other parameters that can influence the accuracy, speed, and output of _App-SpaM_:

**Performance**
| Parameter | Full name | Default | Info |
| -------- | -------- | -------- | -------- |
| `-w`     | `--weight`     | `12`     | Weight of pattern (number of _match positions_ (1s)). Higher weight generally leads to faster computation, but on small datasets it may result in too few spaced words, resulting in low accuracy. |
| `-d`     | `--dontCare`     | `32`     | Number of _don't care positions_ in pattern (number of 0s). |
| `-g`     | `--mode`     | `LCACOUNT`     | Assignment mode determines how a placement position is chosen from the calculated reference-query distances. For more information see paper. Possible values are: `MINDIST`,`SPAMCOUNT`,`LCADIST`,`LCACOUNT`, `APPLES`...|
| `-u`     | `--unassembled`     |     | Enables support for unassembled references, see below. |
|      | `--delimiter`     | `"-"`     | Specifies delimiter in reference names when unassembled mode is executed. All reads from the same reference should have this delimiter in their name. They are then regarded as one reference sequence. |
| `-h`     | `--help`     |     | Show help and exit. |

**Further**
| Parameter | Full name | Default | Info |
| -------- | -------- | -------- | -------- |
| `-b`     | `--read_block_size`     | `100000`     | Weight of pattern (number of match positions). Higher weight generally leads to faster computation, but on small datasets it may result in too few spaced words, resulting in low accuracy. |
| `-v`     | `--verbose`     |       | Outputs additional information about the current run on the standard output. |
|      | `--threads`       | `1`     | Specify number of threads to use. |
|      | `--write-histogram`     |    | Write a histogram of all spaced word matches to file `histogram.txt`. |
|      | `--write-scoring`       |     | Write file with all pairwise distances between references and queries to file `scoring_table.txt`. |
|      | `--threshold`     | `0`     | Specifies filtering threshold of spaced word filtering procedure. |

### Further help
Write to matthias.blanke@biologie.uni-goettingen.de

## Publication
### Associated publication
follows
### Data sets used for comparison
follows
### Additional materials
follows

## Links
[^FSWM]: C.-A. Leimeister, S. Sohrabi-Jahromi, B. Morgenstern (2017)
    Fast and Accurate Phylogeny Reconstruction using Filtered Spaced-Word Matches
    Bioinfomatics 33, 971-979, DOI 10.1093/bioinformatics/btw776

    https://github.com/burkhard-morgenstern/FSWM

[^JPlace]:  Matsen FA, Hoffman NG, Gallagher A, Stamatakis A (2012)
    A Format for Phylogenetic Placements. 
    PLoS ONE 7(2): e31009. https://doi.org/10.1371/journal.pone.0031009

[^PEWO]: Benjamin Linard, Nikolai Romashchenko, Fabio Pardi, Eric Rivals
    PEWO: a collection of workflows to benchmark phylogenetic placement
    Bioinformatics, btaa657, https://doi.org/10.1093/bioinformatics/btaa657

    https://github.com/phylo42/PEWO

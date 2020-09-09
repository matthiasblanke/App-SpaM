# APP-SpaM

APP-SpaM is a program for performing **A**lignment-free **P**hylogenetic **P**lacement based on **SPA**ced-word **M**atches. The input to the program consists of three files: 
1. A fasta file with reference sequences.
2. A newick file containing the reference phylogeny of the reference sequences. 
3. A fasta file of query sequences. 

APP-SpaM will place each of the query sequences in the reference phylogeny at a position that it proposes to be the closest to the correct phylogenetic position with respect to the given references. 
The output is a JPlace file (link) with all query placements.

APP-SpaM is also available in PEWO, the **P**lacement **E**valuation **WO**rkflows, a tool developed to rapidly test and compare tools for phylogenetic placement (link).

When using, please cite:
link

## Setup
1. Download the newest APP-SpaM version:
`git clone https://github.com/matthiasblanke/APP-SpaM`
2. Navigate into the local APP-SpaM directory, create a build folder and go into it:
`cd path/to/appspam`
`mkdir build`
`cd build`
3. Build the program with CMake:
`cmake ..`
`make`

## Running the program
From within the build directory you can just run APP-SpaM like so:
```bash
./appspam -h
```
If you want to run it from anywhere, add it to your PATH:
```bash
export PATH=$PATH:~/path/to/appspam
```
If you want this change to be permanent add this line to your `~/.profile` or `~/.bash_profile`.

### Standard run
When running APP-SpaM you need to at least specify the three input parameters, namely the reference sequences (`-s`), the reference phylogeny (`-t`) and the query sequences (`-q`):
```
./appspam -s path/to/references -t path/to/referencetree -q path/to/query
```
The paths can be either absolut paths, or relative to your current working directory. All other parameters will be set to default values. All output files will be placed in your current working directory. You can specify the output location and file name with the flag `-o`, e.g.:
```
./appspam -s path/to/references.fasta -t path/to/referencetree.nwk -q path/to/query.fasta -o path/to/output.myjplace
```
will use the specified file as output JPlace. If other output files are produced (look below) they will be placed in the same folder as the JPlace file.

### Parameter choices
There are several other parameters that can influence the performance or output of APPSPAM.

**Performance**
| Parameter | Full name | Default | Info |
| -------- | -------- | -------- | -------- |
| `-w`     | `--weight`     | `12`     | Weight of pattern (number of match positions or 1s). Higher weight generally leads to faster computation, but on small datasets it may result in too few spaced words, resulting in low accuracy. |
| `-d`     | `--dontCare`     | `32`     | Number of don't care positions in pattern (number of 0s). |
| `-g`     | `--assignment_mode`     | `LCACOUNT`     | Assignment mode determines how a placement position is chosen from the calculated reference-query distances. For more information see paper. Possible values are: `LCACOUNT`,`BESTCOUNT`,`SAC`,`APPLES`,...|
| `-u`     | `--unassembled`     | `False`     | If this flag is set, the references can be unassembled, see below. |
| `-h`     | `--help`     | `False`     | Show help and exit. |

**Further**
| Parameter | Full name | Default | Info |
| -------- | -------- | -------- | -------- |
| `-b`     | `--read_block_size`     | `100000`     | Weight of pattern (number of match positions). Higher weight generally leads to faster computation, but on small datasets it may result in too few spaced words, resulting in low accuracy. |
| `-v`     | `--verbose`     | `False`     | Outputs additional information about the current run on the standard output. Additional information include: Number of patterns used. Patterns used. Number of spaced words extracted. Number of spaced word matches. |
|      | `--threads`       | `1`     | Specify number of threads to use. |
|      | `--histogram`     | `False`     | Writes a histogram of all spaced word matches to file `histogram.txt`. |
|      | `--scoring`       | `False`     | Writes file with all pairwise distances between references and queries to file `scoring_table.txt`. |
|      | `--threshold`     | `0`     | Specifies filtering threshold of spaced word filtering procedure. |
|      | `--delimiter`     | `"-"`     | Specifies delimiter in reference names when unassembled mode is executed. All reads from the same reference should have this delimiter in their name. They are then regarded as one reference sequence. |

### Further help
Write to matthias.blanke@biologie.uni-goettingen.de

## Publication
### Associated publication
follows
### Data sets used for comparison
follows
### Additional materials
follows


# K-LEB

## K-mer Layers Estimation using Bloom filters

## Bird-eye view
Reads k-mers from files and inserts them in a BloomFilter per file. Once every input file has been processed, Bloom Filters are aggregated to get a count of k-mer presence per colors.
Outputs a CSV that can be used to generate a plot of k-mers abundances by color inside your data.

## Compilation

```sh
git clone   https://github.com/TimRouze/UpSet
cargo build
```

## Usage

```sh
cargo r -r input_fof.txt
```
### Options

- ```-o``` or ```--output```: Output file (Default = Out.csv) 
- ```-t``` or ```--threads```: number of threads (Default = 1)
- ```-m``` or ```--memory```: Memory (in Gb) allocated to each Bloom Filters (1 BF per file + 2 ABF) (Default = 4Gb/BF)
- ```-h``` or ```--hashes```: number of Hashes used in Bloom Filters (Default = 1)
- ```-M``` or ```--modimizer```: Modimizer modulo value (Default = 1)
- ```-h``` or ```--help```: Display help

## Acknowledgements
Many tanks to Igor Martayan who helped me discovering programming in rust and allowed me to re-use lots of code from this repo: [BRRR](https://www.github.com/imartayan/BRRR).

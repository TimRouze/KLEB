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
- ```-m``` or ```--memory```: Memory (in Gb) allocated to each Bloom Filters (2 BF per file) (Default = 4Gb)
- ```-h``` or ```--hashes```: number of Hashes used in Bloom Filters (Default = 1)
- ```-M``` or ```--modimizer```: Modimizer modulo value (Default = 1)

## Acknowledgements
Many tanks to Igor Martayan who helped me discovering programming in rust and allowed me to re-use lots of code from this repo: [BRRR](github.com/imartayan/BRRR).

# UpSet

## UpSet plot generator for k-mer membership in multiple genomes

## Bird-eye view
Reads k-mers from files and inserts them in a BloomFilter. Once every input file has been processed, Bloom Filters are aggregated to get a count of k-mer presence in input files.
Outputs a CSV that can be used to generate an UpSet plot of k-mers abundances inside your data.

## Compilation

```sh
git clone   https://github.com/TimRouze/UpSet
cargo r -r fof.txt
```

## Acknowledgements
Many tanks to Igor Martayan from whom I took a lot of code from his repo: github.com/imartayan/BRRR.

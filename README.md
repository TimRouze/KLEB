# K-LEB

## K-mer Layers Estimation using Bloom filters

## Bird-eye view
Reads k-mers from files and inserts them in a BloomFilter per file. Once every input file has been processed, Bloom Filters are aggregated to get a count of k-mer presence per colors.
Outputs a CSV that can be used to generate a plot of k-mers abundances by color inside your data.

## Compilation

```sh
git clone   https://github.com/TimRouze/UpSet
cargo r -r fof.txt
```

## Acknowledgements
Many tanks to Igor Martayan from whom I took a lot of code from his repo: [BRRR](github.com/imartayan/BRRR).

### Introduction
Locality-sensitive bucketing (LSB) functions generalize the widely-used
locality-sensitive hashing (LSH) methods that are designed to be
able to recognize similar but not necessarily identical sequences.
Such functionalities are useful in many bioinformatics applications
including homology detection, overlap graph construction, and
phylogenetic tree reconstruction, especially when the sequencing
data has a high error rate.

A $(d_1, d_2)$-sensitive bucketing function sends a sequence into
multiple buckets such that any two sequences
within an edit distance of $d_1$ are guaranteed to share at least one bucket,
and any two sequences with an edit distance at least $d_2$
are guaranteed to be mapped into disjoint buckets.
Here we provide implementation of an optimal $(1, 2)$-sensitive bucketing
function that assigns to each length-$n$ sequence $n$ buckets.
If the user is interested in a small subset of all length-$n$ sequences,
we also provide an efficient function that computes the buckets in the
above $(1,2)-sensitive bucketing function for a specific sequence.

A subset of fixed-length sequences is said to be $(1, 1)$-guaranteed
if they can be used to label all the buckets such that the resulting
bucketing function is $(1, 3)$-sensitive.
Here we provide implementation of both the construction of a 
minimum $(1,1)$-guaranteed subset and the efficient linear time 
membership query of such a set.
### Installation
```
git clone https://github.com/Shao-Group/lsbucketing.git
cd lsbucketing
make
```
### Usage
- To generate buckets for all length $n$ sequences, run
`./assignBuckets.out n` where `n` is the length of the sequences.
The results are written in a file named `buckets-n.txt`.
This program also verifies the correctness of 
the efficient algorithm that generates
buckets for a specific sequence without a global counter.
This algorithm is provided as the function 
`assignBuckets` in `src/assignBuckets.c`.

- To generate a $(1,1)$-guaranteed subset, run
`./genSampleD1.out n` where `n` is the length of the sequences.
The results are written in a file named `n01.sample`.
The first line of the file is the number of length-$n$ sequences
in this subset, which equals to $4^{n-1}$ for the default alphabet
{A, C, G, T}.
The remaining lines each contains one sequence.

- `./LSB-statistics.out` can be used to reproduce the 
experimental results in [our paper](https://arxiv.org/abs/2206.03097).
It takes three parameters:
  - `n`: integer, the length of the sequences.
  - `r`: integer, the radius of the neighborhood of a sequence.
  - `w|s`: char, the option `w` uses all the length-$n$ sequences as the 
  bucketing set; the option `s` uses a $(1,1)$-guaranteed subset.
  The program utilizes the efficient membership query function
  `isInSampleD1` defined in `lib/util.h` so it is not needed to 
  explicitly generate the $(1,1)$-guaranteed subset.

  Results for three LSB functions are given in the paper, the corresponding
  parameters are
  ```
  ./LSB-statistics.out 20 1 w
  ./LSB-statistics.out 20 1 s
  ./LSB-statistics.out 20 2 s
  ```
  The program writes to standard output which can be redirected
  ```
  ./LSB-statistics.out 20 1 w > output.txt &
  ```

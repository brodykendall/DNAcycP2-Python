DNAcycP Python package 
================

**Maintainer**: Ji-Ping Wang, \<<jzwang@northwestern.edu>\>; Brody Kendall \<<curtiskendall2025@u.northwestern.edu>\>; Keren Li, \<<keren.li@northwestern.edu>\>

**License**: GPLv3

**Cite DNAcycP package**:

TODO: Update citation when applicable

Li, K., Carroll, M., Vafabakhsh, R., Wang, X.A. and Wang, J.-P., DNAcycP: A Deep Learning Tool for DNA Cyclizability Prediction, *Nucleic Acids Research*, 2021

## What is DNAcycP?

**DNAcycP**, short for **DNA** **cyc**lizablity **P**rediction, is a Python package for accurate prediction of DNA intrinsic cyclizablity score. It was built upon a deep learning architecture with a hybrid of Inception and Residual network structure and an LSTM layer. The original DNAcycP was trained based on loop-seq data from Basu et al 2021 (see below). An updated version (DNAcycP2) was trained based on smoothed predictions of this loop-seq data. The predicted score (for either DNAcycP or DNAcycP2), termed **C-score** achieves high accuracy compared to the experimentally measured cyclizablity score by loop-seq assay.

## Available format of DNAcycP

TODO: update reference to R package
TODO: update web server

DNAcycP is available in three formats: A web server available at http://DNAcycP.stats.northwestern.edu for real-time prediction and visualization of C-score up to 20K bp, a standalone Python package available for free download from https://github.com/jipingw/DNAcycP, and an R package (coming soon).


## Architecture of DNAcycP

The core of DNAcycP is a deep learning architecture mixed with an Inception-ResNet structure and an LSTM layer (IR+LSTM, Fig 1b) that processes the sequence and its reverse complement separately, the results from which are averaged and detrended to reach the predicted intrinsic score. (Fig 1a).

IR+LSTM starts with a convolutional layer for dimension reduction such that the encoded sequence space is reduced from 2D to 1D. The output is fed into an inception module that contains two parallel branches, each having two sequentially connected convolutional layers with branch-specific kernels to capture sequence features of different scale. The first branch has kernel dimension 3x1 for both layers and the second has kernel dimension 11x1 and 21x1 sequentially. The output of the inception module is combined by concatenation and added back to the input of the inception module to form a short circuit or residual network. Finally, the IR+LSTM concludes with a dense layer to predict output with linear activation function. 

![A diagram of DNAcycP.](./figures/Figure1.png)

## DNAcycP required packages

*The recommended python version is 3.11*

* `numpy==1.26.1`
* `pandas==2.1.2`
* `tensorflow==2.14.0`
* `keras==2.14.0`
* `bio==1.7.1`
* `docopt==0.6.2`


## Installation

**DNAcycP** Python package requires specific versions of dependencies. We recommend to install and run **DNAcycP** in a virtual environment. For example, suppose the downloaded DNAcycP package is unpacked as a folder `dnacycp-main`. We can install DNAcycP in a virtual environment as below:

```bash
cd dnacycp-main
python3 -m venv env
source env/bin/activate test
pip install -e .
```

Run `dnacycp-cli ` to see whether it is installed properly.

*Note: You may need to deactivate then re-activate the virtual environment prior to this step (see below)*

```bash
dnacycp-cli 
```

Once done with DNAcycP for prediction, you can close the virtual environment by using:
```bash
deactivate
```

Once the virtual environment is deactivated, you need to re-activate it before you run another session of prediciotn as follows:
```bash
cd dnacycp-main
source env/bin/activate test
```

## Usage

DNAcycP supports the input sequence in two formats: FASTA format (with sequence name line beginning with “>”) or plain TXT format. Unlike in the web server version where only one sequence is allowed in input for prediction, the Python package allows multiple sequences in the same input file. In particular for the TXT format, each line (can be of different length) in the file is regarded as one input sequence for prediction. 

The main funciton in DNAcycP is `dnacycp-cli`, which using one of the following lines:
```bash
dnacycp-cli -f -s <inputfile> <basename> [-L <chunk_length>] [-n <num_cores>]
dnacycp-cli -f <inputfile> <basename> [-L <chunk_length>] [-n <num_cores>]
dnacycp-cli -t -s <inputfile> <basename>
dnacycp-cli -t <inputfile> <basename>
```

where 
  * `-f/-t`: indicates the input file name in FASTA or TXT format respectively; either one must be specified.
  * `-s`: (optional) indicates the updated model trained on smoothed data (DNAcycP2) should be used. If `-s` is omitted, the model trained on the original data (DNAcycP) will be used.
  * `<inputfile>`: is the name of the intput file;
  * `<basename>`: is the name base for the output file.
  * `-L <chunk_length>`: is the length of sequence that a given core will be predicting on at any given time (default 100000; only applicable with `-f`)
  * `-n <num_cores>`: is the number of cores to be used in parallel (default 1; only applicable with `-f`)

### Example 1:

```bash
dnacycp-cli -f -s ./data/raw/ex1.fasta ./data/raw/ex1_smooth -L 1000 -n 2
dnacycp-cli -f ./data/raw/ex1.fasta ./data/raw/ex1_original -L 1000 -n 2
```

The `-f` option specifies that the input file named "ex1.fasta" is in fasta format. 

The `-s` option specifies that the DNAcycP2 model should be used for prediction, while omitting the `-s` argument specifies that the original DNAcycP model should be used for prediction.

The `-L` option specifies that prediction will occur on sequences of length 1000 for a given core

The `-n` option specifies that prediction will occur on 2 cores in parallel

The `./data/raw/ex1.fasta` is the sequence file path and name, and `./data/raw/ex1` specifies the output file will be saved in the directory `./data/raw` with file name initialized with `ex1`.
For example, `ex1.fasta` contains two sequences with IDs "1" and "2" respectively.
The output files containing DNAcycP2 predictions will be named as "ex1_smooth_cycle_1.txt" and "ex1_smooth_cycle_2.txt" for the first and second sequences respectively, while the output files containing DNAcycP predictions will be named as "ex1_original_cycle_1.txt" and "ex1_original_cycle_2.txt".

 Each output file contains three columns: `position`, `C_score_norm`, `C_score_unnorm`. The `C_score_norm` is the predicted C-score from the model trained based on the standardized loop-seq score (in the case of DNAcycP) or the standardized smoothed intrinsic cyclizability estimate (in the case of DNAcycP2) of the tiling library of Basu et al 2021 (i.e. 0 mean unit variance). When predictions are made using the original DNAcycP, the `C_score_unnorm` is the predicted C-score recovered to the original scale of loop-seq score in the tiling library data from Basu et el 2021. When predictions are made using the updated DNAcycP2 (`-s`), the `C_score_unnorm` is the predicted C-score recovered to the scale of standardized raw cyclizability scores of the tiling library data. The standardized loop-seq score provides two advantages. As loop-seq may be subject to a library-specific constant, standardized C-score is defined with a unified baseline as yeast genome (i.e. 0 mean in yeast genome). Secondly, the C-score provides statisitcal significance indicator, i.e. a C-score of 1.96 indicates 97.5% in the distribution.


### Example 2:

```bash
dnacycp-cli -t -s ./data/raw/ex2.txt ./data/raw/ex2_smooth
dnacycp-cli -t ./data/raw/ex2.txt ./data/raw/ex2_original
```
With `-t` option, the input file is regarded as in TXT format, each line representing a sequence without sequence name line that begins with ">".

The `-s` option again specifies that the DNAcycP2 model should be used for prediction, while omitting the `-s` argument specifies that the original DNAcycP model should be used for prediction.

The predicted C-scores will be saved into two files, one with `_unnorm.txt` and the other with `_norm.txt` for unnormalized and normalized C-score, with C-scores in each line corresponding to the sequence in the input file in the same order.

For any input sequence, DNAcycP predicts the C-score for every 50 bp. Regardless of the input sequence format the first C-score in the output file corresponds to the sequence from position 1-50, second for 2-51 and so forth.

### Run prediction within Python interactive session

```python
from dnacycp import cycle_fasta, cycle_txt

# Smooth prediction using DNAcycP2:
cycle_fasta("data/raw/ex1.fasta","ex1_smooth", smooth=True, chunk_size=1000, num_threads=2)
cycle_txt("data/raw/ex2.txt","ex2_smooth",smooth=True)

# Original prediction using DNAcycP:
cycle_fasta("data/raw/ex1.fasta","ex1_original",smooth=False, chunk_size=1000, num_threads=2)
cycle_txt("data/raw/ex2.txt","ex2_original",smooth=False)
```


## Other References

* Basu, A., Bobrovnikov, D.G., Qureshi, Z., Kayikcioglu, T., Ngo, T.T.M., Ranjan, A., Eustermann, S., Cieza, B., Morgan, M.T., Hejna, M. et al. (2021) Measuring DNA mechanics on the genome scale. Nature, 589, 462-467.



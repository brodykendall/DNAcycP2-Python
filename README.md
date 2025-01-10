DNAcycP2 Python package 
================

**Maintainer**: Ji-Ping Wang, \<<jzwang@northwestern.edu>\>; Brody Kendall \<<curtiskendall2025@u.northwestern.edu>\>; Keren Li, \<<keren.li@northwestern.edu>\>

**License**: GPLv3

**Cite DNAcycP2 package**:

TODO: Update citation when applicable

## What is DNAcycP2?

**DNAcycP2**, short for **DNA** **cyc**lizablity **P**rediction v**2**, is a Python package for accurate, unbiased prediction of DNA intrinsic cyclizablity score. It was built upon a deep learning architecture with a hybrid of Inception and Residual network structure and an LSTM layer. DNAcycP2 is an updated version of DNAcycP, released by Li et al 2021 (see below). DNAcycP was trained based on loop-seq data from Basu et al 2021 (see below), while DNAcycP2 was trained based on smoothed predictions of this loop-seq data. The predicted score (for either DNAcycP or DNAcycP2), termed **C-score** achieves high accuracy compared to the experimentally measured cyclizablity score by loop-seq assay.

## Key differences between DNAcycP2 and DNAcycP

Following the release of DNAcycP, it was discovered that the training data contained residual measurement bias, leading to biased predictions. To correct this bias in the data, we employed a data augmentation + periodic smoothing approach to generate new, unbiased estimates of intrinsic DNA cyclizability for each sequence in the original training dataset. We then trained a new model on the unbiased data with architecture identical to that of DNAcycP, named DNAcycP2. More details on this process can be found in the following paper: (CITATION).

TODO: fill in citation above

Previously, the measurement bias introduced by the location of the biotin tether was not adequately accounted for. By employing data augmentation and smoothing with a moving average approach over the length of 1 full helical repeat at 1bp resolution in the genome, we can remove this bias while still maintaining high resolution, accurate estimates of intrinsic cyclizability.

![Visualization of difference between DNAcycP2 and DNAcycP.](./figures/Figure7.png)

## Available formats of DNAcycP2 and DNAcycP

DNAcycP2 is available in three formats: A web server available at http://DNAcycP.stats.northwestern.edu for real-time prediction and visualization of C-score up to 20K bp, a standalone Python package avilable for free download from https://github.com/jipingw/DNAcycP2, and an R package available for free download from https://github.com/jipingw/dnacycp2-R.

TODO: update web server - possible selection on server of which model to use?

DNAcycP is still available in its two original formats: A web server available at http://DNAcycP.stats.northwestern.edu for real-time prediction and visualization of C-score up to 20K bp, and a standalone Python package available for free download from https://github.com/jipingw/DNAcycP

## Architecture of DNAcycP2

The core of DNAcycP2 is a deep learning architecture mixed with an Inception-ResNet structure and an LSTM layer (IR+LSTM, Fig 1b) that processes the sequence and its reverse complement separately, the results from which are averaged and detrended to reach the predicted intrinsic score. (Fig 1a).

IR+LSTM starts with a convolutional layer for dimension reduction such that the encoded sequence space is reduced from 2D to 1D. The output is fed into an inception module that contains two parallel branches, each having two sequentially connected convolutional layers with branch-specific kernels to capture sequence features of different scale. The first branch has kernel dimension 3x1 for both layers and the second has kernel dimension 11x1 and 21x1 sequentially. The output of the inception module is combined by concatenation and added back to the input of the inception module to form a short circuit or residual network. Finally, the IR+LSTM concludes with a dense layer to predict output with linear activation function. 

![A diagram of DNAcycP2.](./figures/Figure1.png)

## DNAcycP2 required packages

*The recommended python version is 3.11*

* `numpy==1.26.1`
* `pandas==2.1.2`
* `tensorflow==2.14.0`
* `keras==2.14.0`
* `bio==1.7.1`
* `docopt==0.6.2`


## Installation

**DNAcycP2** Python package requires specific versions of dependencies. We recommend to install and run **DNAcycP2** in a virtual environment. For example, suppose the downloaded DNAcycP2 package is unpacked as a folder `dnacycp2-main`. We can install DNAcycP2 in a virtual environment as below:

```bash
cd dnacycp2-main
python3 -m venv env
source env/bin/activate test
pip install -e .
```

Run `dnacycp2-cli ` to see whether it is installed properly.

*Note: You may need to deactivate then re-activate the virtual environment prior to this step (see below)*

```bash
dnacycp2-cli 
```

Once done with DNAcycP2 for prediction, you can close the virtual environment by using:
```bash
deactivate
```

Once the virtual environment is deactivated, you need to re-activate it before you run another session of prediction as follows:
```bash
cd dnacycp2-main
source env/bin/activate test
```

## Usage

DNAcycP2 supports the input sequence in two formats: FASTA format (with sequence name line beginning with “>”) or plain TXT format. Unlike in the web server version where only one sequence is allowed in input for prediction, the Python package allows multiple sequences in the same input file. In particular for the TXT format, each line (can be of different length) in the file is regarded as one input sequence for prediction, however the computation is significantly more efficient when every sequence has length exactly 50bp. 

The main function in DNAcycP2 is `dnacycp2-cli`, which is called through one of the following lines:
```bash
dnacycp2-cli -f -s <inputfile> <basename> [-L <chunk_length>] [-n <num_cores>]
dnacycp2-cli -f <inputfile> <basename> [-L <chunk_length>] [-n <num_cores>]
dnacycp2-cli -t -s <inputfile> <basename>
dnacycp2-cli -t <inputfile> <basename>
```

where 
  * `-f/-t`: indicates the input file name in FASTA or TXT format respectively; either one must be specified.
  * `-s`: (optional) indicates the updated model trained on smoothed data (**DNAcycP2**) should be used. **If `-s` is omitted, the model trained on the original, biased data (DNAcycP) will be used**.
  * `<inputfile>`: is the name of the intput file;
  * `<basename>`: is the name base for the output file.
  * `-L <chunk_length>`: is the length of sequence that a given core will be predicting on at any given time (default 100000; only applicable with `-f`)
  * `-n <num_cores>`: is the number of cores to be used in parallel (default 1 - no parallelization; only applicable with `-f`)

The `-f` setting (FASTA format) is designed for larger files, so it has added parallelization capability. To utilize this capability, specify the number of cores to be greater than 1 using the `n_cores` argument (default 1 - no parallelization). You can also specify the length of the sequence that each core will predict on at a given time using the `chunk_length` argument (default 100000).

For reference, on a personal computer (16 Gb RAM, M1 chip with 8-core CPU), prediction at full parallelization directly on the yeast genome FASTA file completes in 12 minutes, and on the hg38 human genome Chromosome I FASTA file in just over 4 hours. In our experience, selection of parallelization parameters (-L and -n) has little affect when making predictions on a personal computer, but if using the package on a high-performance compute cluster, prediction time should decrease as the number of cores increases. If you do run into memory issues, we suggest first reducing -L

### Example 1:

```bash
dnacycp2-cli -f -s ./data/raw/ex1.fasta ./data/raw/ex1_smooth -L 1000 -n 2
dnacycp2-cli -f ./data/raw/ex1.fasta ./data/raw/ex1_original -L 1000 -n 2
```

The `-f` option specifies that the input file named "ex1.fasta" is in fasta format. 

The `-s` option specifies that the DNAcycP2 model should be used for prediction, while omitting the `-s` argument specifies that the original DNAcycP model should be used for prediction.

The `-L` option specifies that prediction will occur on sequences of length 1000 for a given core

The `-n` option specifies that prediction will occur on 2 cores in parallel

The `./data/raw/ex1.fasta` is the sequence file path and name, and `./data/raw/ex1` specifies the output file will be saved in the directory `./data/raw` with file name initialized with `ex1`.
For example, `ex1.fasta` contains two sequences with IDs "1" and "2" respectively.
The output files containing DNAcycP2 predictions will be named as "ex1_smooth_cycle_1.txt" and "ex1_smooth_cycle_2.txt" for the first and second sequences respectively, while the output files containing DNAcycP predictions will be named as "ex1_original_cycle_1.txt" and "ex1_original_cycle_2.txt".

Each output file contains three columns. The first columns is always `position`. When `-s` is specified, the next columns are `C0S_score_norm`, `C0S_score_unnorm`, and when `-s` is not specified, the next columns are `C0_score_norm` and `C0_score_unnorm`. The predicted C-score for either model is the normalized output (`C0S_score_norm` and `C0_score_norm`), the predictions from the model trained based on the standardized loop-seq score (in the case of DNAcycP) or the standardized smoothed intrinsic cyclizability estimate (in the case of DNAcycP2) of the Tiling library of Basu et al 2021 (i.e. 0 mean unit variance). When predictions are made using the original DNAcycP, the `C0_score_unnorm` is the predicted C-score recovered to the original scale of loop-seq score in the Tiling library data from Basu et el 2021. When predictions are made using the updated DNAcycP2 (`-s`), the `C0S_score_unnorm` is the predicted C-score recovered to the scale of standardized raw cyclizability scores of the Tiling library data. The standardized scores provide two advantages. As loop-seq may be subject to a library-specific constant, standardized C-score is defined with a unified baseline as yeast genome (i.e. 0 mean in yeast genome). Secondly, the C-score provides statisitcal significance indicator, i.e. a C-score of 1.96 indicates 97.5% in the distribution.


### Example 2:

```bash
dnacycp2-cli -t -s ./data/raw/ex2.txt ./data/raw/ex2_smooth
dnacycp2-cli -t ./data/raw/ex2.txt ./data/raw/ex2_original
```
With `-t` option, the input file is regarded as in TXT format, each line representing a sequence without sequence name line that begins with ">".

The `-s` option again specifies that the DNAcycP2 model should be used for prediction, while omitting the `-s` argument specifies that the original DNAcycP model should be used for prediction.

The predicted C-scores will be saved into two files, one with `_C0S_unnorm.txt` and the other with `_C0S_norm.txt` in the DNAcycP2 (`-s`) case, or `_C0_unnorm.txt` and `_C0_norm.txt` in the DNAcycP case. C-scores in each line correspond to the sequence in the input file in the same order.

For any input sequence, DNAcycP2 predicts the C-score for every 50 bp. Regardless of the input sequence format the first C-score in the output file corresponds to the sequence from position 1-50, second for 2-51 and so forth.

### Run prediction within Python interactive session

```python
from dnacycp2 import cycle_fasta, cycle_txt

# Smooth prediction using DNAcycP2:
cycle_fasta("data/raw/ex1.fasta","ex1_smooth", smooth=True, chunk_size=1000, num_threads=2)
cycle_txt("data/raw/ex2.txt","ex2_smooth",smooth=True)

# Original prediction using DNAcycP:
cycle_fasta("data/raw/ex1.fasta","ex1_original",smooth=False, chunk_size=1000, num_threads=2)
cycle_txt("data/raw/ex2.txt","ex2_original",smooth=False)
```


## Other References

* Li, K., Carroll, M., Vafabakhsh, R., Wang, X.A. and Wang, J.-P., DNAcycP: A Deep Learning Tool for DNA Cyclizability Prediction, *Nucleic Acids Research*, 2021

* Basu, A., Bobrovnikov, D.G., Qureshi, Z., Kayikcioglu, T., Ngo, T.T.M., Ranjan, A., Eustermann, S., Cieza, B., Morgan, M.T., Hejna, M. et al. (2021) Measuring DNA mechanics on the genome scale. Nature, 589, 462-467.



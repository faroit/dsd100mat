# DSD100 Matlab

MATLAB scripts to parse and process the demixing secrets dataset (DSD100) as
part of the [Signal Separation Evaluation Campaign (SISEC)](https://sisec.inria.fr/).
This scripts are intended to perform the full evaluation of the separation
quality of your estimates on the DSD100.

## Usage

We provide two functions

Function Name  | Description
------------- | -------------
`DSD100_separate_and_eval.m`  | Parse the DSD100 and generates estimates with a user provided function. Evaluates with BSS_EVAL and saves results.
`DSD100_eval_only.m`  | Only evaluates existing estimates (also multiple methods) with BSS_EVAL and saves results. Good in combination with the [DSD100 python package](https://github.com/faroit/dsd100py)

### Separate and Evaluate

The function `DSD100_separate_and_eval.m` should be used along with the
Demixing Secret Dataset 100 (DSD100) for the purpose of music source separation.
Before you start please set the `dataset_folder` to point to the DSD100 root folder.

The separation function should be named `myfunction.m` placed in the
root folder, and have the following syntax:

```matlab
[bass, drums, other, vocals, accompaniment] = myfunction(mixture, fs)
```
where `mixture` is a matrix of size `[#samples, #channels]` corresponding
to the mixture, `fs` is the corresponding sampling frequency in Hz, and
`bass`, `drums`, `other`, `vocals` and `accompaniment` are matrices of
same size as the mixture corresponding to the estimates, i.e., the bass,
the drums, the other instruments, the vocals, and the full accompaniment
(i.e., bass+drums+other), respectively. If one or more sources are not
meant to be estimated, the function should return an empty matrix (i.e.,
`[]`). Any other parameter of the algorithm should be defined internally.

The evaluation function should then be called simply as follows:
```matlab
DSD100_separate_and_eval.m
```
The function loops over all the 100 songs of the MSD100 data set, and,
for each song, for both the "Dev" and "Test" subsets, performs source
separation on the mixture `mixture.wav` from the folder "Mixtures" using
the separation function `myfunction.m` and saves the estimates as
`bass.wav,` `drums.wav,` `other.wav,` and `vocals.wav` (if estimated) to
the folder `Estimates.` The function then measures the separation
performance using the BSS Eval toolbox 3.0 (included in this function)
and the sources from the folder "Sources," and saves the results (i.e.,
SDR, ISR, SIR, and SAR) in the file "results.mat," including the song
name and the processing time, along with the estimates to the folder
"Estimates". The function also saves the results for all the songs in a
single file "result.mat" to the estimates folder, along with this function.

### Evaluate only

If you already have generated the estimates before (e.g. by using the [DSD100 python package](https://github.com/faroit/dsd100py)
you can run `DSD100_eval_only.m` separately. Before you start please set the
`dataset_folder` to point to the DSD100 folder.

The variable ```base_estimates_directory``` stand for the root folder in
which the script should find subdirectories containing the results of the
methods you want to evaluate. each of these subdirectories must contain
the exact same file structure than the DSD dataset, as produced by the
DSD100_separate_and_eval_parallel.m script or the [DSD100 python package](https://github.com/faroit/dsd100py).
The matching is case sensitive. There is the possibility of not including
all sources, and also to include the "accompaniment" source, defined as
the sum of all sources except vocals.

The evaluation function should then be called simply as follows:
```matlab
DSD100_eval_only.m
```

The function loops over all the 100 songs of the MSD100 data set, and,
for each song, for both the "Dev" and "Test" subsets, performs evaluation
using the BSS Eval toolbox 3.0 (included in this function) and saves the
results (i.e. SDR, ISR, SIR, and SAR) in the file "results.mat," including
the song name. The function also saves the results for all the songs in a
single file "resultX.mat" to the estimates folder, along with this function.

A first evaluation is performed for the 4 sources vocals/bass/drums
and other, and a second is performed for accompaniment.

### References

We would like to thank [Emmanuel Vincent](http://www.loria.fr/~evincent/) for giving us the permission to
use the [BSS Eval toolbox 3.0](http://bass-db.gforge.inria.fr/bss_eval/)

If you use this script, please reference the following paper
```latex
@inproceedings{SiSEC2015,
    TITLE = {{The 2015 Signal Separation Evaluation Campaign}},
    AUTHOR = {N. Ono and Z. Rafii and D. Kitamura and N. Ito and A. Liutkus},
    BOOKTITLE = {{International Conference on Latent Variable Analysis and Signal Separation  (LVA/ICA)}},
    ADDRESS = {Liberec, France},
    SERIES = {Latent Variable Analysis and Signal Separation},
    VOLUME = {9237},
    PAGES = {387-395},
    YEAR = {2015},
    MONTH = Aug,
}
```

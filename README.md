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
`DSD100_only_eval.m`  | Only evaluates existing estimates folder with BSS_EVAL and saves results. Good in combination with the [DSD100 python package](https://github.com/faroit/dsd100py)

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
you can run `DSD100_only_eval.m` separately. You therefore should place the file in a folder,
along with your estimates, containing your results to evaluate. Before you start please set the `dataset_folder` to point to the DSD100  folder and the `estimates_folder` to point to your estimates.

The `estimates_folder` should have exactly the same structure as the
DSD100/Sources folder, the matching is case sensitive. In each directory,
there should be the separated sources estimates whose quality is to be
evaluated. There is the possibility of not including all sources, and
also to include the "accompaniment" source, defined as the sum of all
sources except vocals.

The evaluation function should then be called simply as follows:
```matlab
DSD100_only_eval.m
```

The function loops over all the 100 songs of the MSD100 data set, and,
for each song, for both the "Dev" and "Test" subsets, performs evaluation
using the BSS Eval toolbox 3.0 (included in this function) and saves the
results (i.e. SDR, ISR, SIR, and SAR) in the file "results.mat," including
the song name. The function also saves the results for all the songs in a
single file "resultX.mat" to the estimates folder, along with this function.

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

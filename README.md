# oku2023dmn
Code for a manuscipt of the default mode network in common marmoset (2023)

## Requirements: Software
* MATLAB R2019b or later
* [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn)

Please download the [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn) and "Add Path" in the MATLAB before using this code.

## Installation
1. Download this code and [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn) zip files.
2. Extract zip files under your working directory <work_path>.
3. Run the MATLAB software, and "Add Path" extracted directories (i.e. <work_path>/vardnn-main and <work_path>/oku2023dmn-main).
4. Move to <work_path>/oku2023dmn-main directory and run the following demos.

## Demo Codes
<b>Demo</b><br>
The first demo shows the calculation of GLM analysis for auditory task-fMRI data in common marmoset.<br>
Pre-processed NIfTI files should be downloaded from [zendo](https://zendo/) and extracted under 'data' directory before running this code.
~~~
>> marmoAudGLMindividual
loading : data/s34waM3_1.nii.gz
apply mask atlas...
apply highpass filter (0.0078125 Hz) : tfMRI and design matrix...
process GLM with Tukey-Taper(8) estimation ...
done t=5.3043sec
P-value=0.05, T-value=1.6525
Tmax of GLM6 marmoAuCube1s34waM3_1CTukey8 : audio tmax=9.2829, tcnt=8111, mrv=1.6545
...
~~~

After calculation of GLM analysis for individual sessions, mixed-effects model (2nd analysis) could be applied.
~~~
>> marmoAudGLMmixed
process GLM with Tukey-Taper(8) estimation ...
done t=0.5378sec
P-value=0.001, T-value=4.1437
Tmax of GLM6marmoAudD 2nd-mix-Tukey8full : audio tmax=50.7377, tcnt=1759, mrv=0.0036622
~~~

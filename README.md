# SLPred: a multi-view subcellular localization prediction tool for human proteins
* SLPred is a multi-view subcellular localization prediction tool for human proteins.
* The tool consists of nine independently developed model for the proteins which have annotation with nine subcellular locations: **Cytoplasm, Nucleus, Cell Membrane, Mitochondrion, Extra cellular, Endoplasmic reticulum, Golgi apparatus, Lysosome and Peroxisome.** 
* SLPred exploits the features of nineteen different protein descriptors from the publicly available tools: POSSUM, iFeature and SPMAP.
* Support Vector Machine (SVM) is used to construct probabilistic prediction models, which produces probabilistic scores indicating the localization probability for a query protein sequence. 
* A weighted score is calculated based on the obtained probabilistic scores from seven feature-based probabilistic prediction models (SVMs) by employing weighted mean voting.
* Binary prediction is given by applying thresholding on the weighted score.

* The following figure shows the proposed method
![alt text](https://github.com/gozsari/SLPred/blob/master/images/model_architecture.png)

## Installation

SLPred is a command-line prediction tool written in Python 3.7.1. SLPred was developed and tested in Ubuntu 20.04 LTS. Please run the below commands to install requirements. Dependencies are available in requirements.txt file.

```
conda create -n slpred_env python=3.7
conda activate slpred_env
conda install -r requirements.txt
```

## How to run SLPred to obtain the predictions 

* Clone the Git Repository
* Download Trust dataset (the dataset we created for training)
* Download saved models (these are the pre-trained models)
* **cd SLPred/ncbi-blast** (navigate ncbi-blast folder which is aldready inside SLPred)
* **chmod 777 psiblast** (give necessary permissions to use psiblast command inside the code)
* **cd ..** (navigate back to SLPred folder)
* Put the fasta file (that you want to take predictions) under the folder **fasta_files**. Fasta file may contain any number of sequences. 
* Run **run_SLPred.py** script as shown below 
## Fasta file format

```
>sp|A9WZ33|14KL_BRUSI
MNSFRKTCAGALALIFGATSIVPTVAAPMNMDRPAINQNVIQARAHYRPQNYNRGHRPGY
WHGHRGYRHYRHGYRRHNDGWWYPLAAFGAGAIIGGAISQPRPVYRAPAGSPHVQWCYSR
YKSYRASDNTFQPYNGPRKQCRSPYSR
>sp|C0JAT6|A1H3_LOXHI Dermonecrotic toxin LhSicTox-alphaIA1iii (Fragment) OS=Loxosceles hirsuta OX=571525 PE=2 SV=1
WIMGHMVNAIGQIDEFVNLGANSIETDVSFDSSANPEYTYHGIPCDCGRNCKKWENFNDF
```
## Explanation of Parameters
* **--file**: this is the file name of the fasta file. For example if fasta file name is **input.fasta**, this argument must be just **input**

### To run SLPred the command is as follows:
```
python run_SLPred.py --file input 
```



Yeast-MetaTwin
======

<p align="center">
  <img  src="doc/logo.png" width = "600">
</p>


## Introduction
We designed a systematic workflow for mining underground metabolism, which combines rule-based retrosynthesis approach with deep learning-based enzyme annotation approach. Using this workflow, we constructed **Yeast-MetaTwin**, the first genome-scale metabolic model that systematically integrates underground networks, Yeast-MetaTwin encompasses 84% of the predicted metabolic enzymes and 92% of the metabolome in yeast.

## Usage
  - Download the Yeast-MetaTwin package
  
         git clone https://github.com/LiLabTsinghua/Yeast-MetaTwin.git
  
  - Create and activate enviroment
  
         conda create -n  Yeast_MT python=3.7.16
         conda activate Yeast_MT

  - Download required Python package
         
         conda install ipykernel
         pip install biopython
         pip install cobra
         pip install fair-esm
         pip install gurobipy  
         pip install matplotlib
         pip install numpy
         pip install pandas
         pip install plotly    
         pip install pubchempy
         pip install rdchiral
         pip install rdkit  
         pip install rdkit-pypi
         pip install rxnmapper
         pip install scikit-learn
         pip install seaborn        
         pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu113
          

## Citation
Please cite the paper [To be updated] ()

Contact
-------

-   Feiran Li ([@feiranl](https://github.com/feiranl)), Tsinghua Shenzhen International Graduate School, Tsinghua University, Shenzhen, China
-   Ke Wu ([@wuke](https://github.com/wuke0714)), Tsinghua Shenzhen International Graduate School, Tsinghua University, Shenzhen, China


Last update: 2024-07-15

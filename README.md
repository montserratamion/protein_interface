


## Overview

In this project I implemented python modules to analyse the interface of SARS-CoV-2 spike protein (prot_A) which plays a key role in the receptor recognition to the host receptor angiotensin-converting enzyme 2 (ACE2, prot_B). This preliminary study looks for the binding of prot_A with different proteins to diminish binding with Antibody P5A3C8 (prot_H and prot_L) and maintain or increase binding with ACE2 (prot_B). Studying their interfaces will provide insights that could inform future interface redesign efforts for prot_A.

Usage is in the Jupyter notebok : example.ipynb


<img src="./images/prot_AB_chain1.png" width= "80%">



## 1. Setup

Python libraries used:

- MDAnalysis

The set up can be easily done with conda enviroments using the YAML file.

```
conda env create -f analysis_env.yaml
```

## 2. Python modules overview and approach


 ### Algorithm and approach
<img src="./images/approach.png" width= "100%">


### python module

```mermaid
flowchart TB
%%% Colors %%%
linkStyle default stroke-width:3px,stroke:black
classDef lightblue fill:#d0e1f3,stroke:#292d3f,stroke-width:2px,color:#052003
classDef red fill:#f2d8dd,stroke:#292d3f,stroke-width:2px,color:#052003
classDef green fill:#aabc7d,stroke:#292d3f,stroke-width:2px,color:#052003
classDef green2 fill:#5ac1c0,stroke:#292d3f,stroke-width:2px,color:#052003
classDef orange fill:#d08c23,stroke:#292d3f,stroke-width:2px,color:#052003


%%% main branches of the tool %%%
    %%% link 0
    I[(PDB file)]:::lightblue ---> T(my_tools):::green

    T(my_tools) ---> E(MyPDBTool Class):::green2 & X(ProteinAProteinB Class):::green2
    E(MyPDBTool Class)---> ML2(describe):::red & ML3(select_chain):::red & ML4(dimensions):::red


    X(ProteinAProteinB Class)---> ML5(get_com_residue_residue_distance):::red & ML6(get_com_protein_residue_distance):::red & ML7(search_protein_A_protein_B_interface):::red

```



## 3. Outputs

PDB outputs

<img src="./images/result.png" width= "100%">

## 4. Usage
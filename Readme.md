# Landscape and transition path approach for distributed working memory

This repository provides the necessary resources and instructions to reproduce the results presented in the paper titled "[Controlling brain dynamics: landscape and transition path for working memory](
https://doi.org/10.48550/arXiv.2209.05002)".

## Paper Abstract

Understanding the underlying dynamical mechanisms of the brain and controlling it is a crucial issue in brain science. The energy landscape and transition path approach provides a possible route to address these challenges. Here, taking working memory as an example, we quantified its landscape based on a large-scale macaque model. The working memory function is governed by the change of landscape and brain-wide state switching in response to the task demands. The kinetic transition path reveals that information flow follows the direction of hierarchical structure. Importantly, we propose a landscape control approach to manipulate brain state transition by modulating external stimulation or inter-areal connectivity, demonstrating the crucial roles of associative areas, especially prefrontal and parietal cortical areas in working memory performance. Our findings provide new insights into the dynamical mechanism of cognitive function, and the landscape control approach helps to develop therapeutic strategies for brain disorders.

## Requirements

All simulations are implemented at MATLAB R2020b.

## Dataset

The dataset used in this study is available at [https://core-nets.org](https://core-nets.org). But we also provided the preprocessed data in `model/`, including `spine_count.mat`, `SLN_matrix.mat` and `FLN_matrix.mat`.



## Instructions

Follow these steps to reproduce the results:

1. Clone this repository to your local machine:

```bash
git clone [repository_url]
```

2. Download the dataset from [https://core-nets.org](https://core-nets.org), place it in the specified directory and preprocessing (Optional).

3. Run `model/parameters.m` to set the model parameters.

4. Run `simulation/SimulationIPulse.m`. You will get the simulation result with different initial conditions under specified parameters, which is used in the construction of potential landscape.

5. Run `main.m` to reproduce the Fig.1 G and Fig. 3D in our paper.

## Citation  

If you use this code or replicate the results of the paper, please cite the original paper:  


## Contact  

If you have any questions or need further assistance, please feel free to contact yeyeleijun@gmail.com.
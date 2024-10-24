# Beyond Diagonal RIS: Passive Maximum Ratio Transmission and Interference Nulling Enabler

This repository contains the code and data associated with the paper:

H. Yahya, H. Li, M. Nerini, B. Clerckx and M. Debbah, "[Beyond Diagonal RIS: Passive Maximum Ratio Transmission and Interference Nulling Enabler](https://arxiv.org/abs/2408.09887)," Archive, 2024.

## Table of Contents
- [Introduction](#introduction)
- [Code Overview](#code-overview)
- [Routines](#routines)
- [Results](#results)
- [Citation](#citation)

## Introduction
This repository provides the implementation of the proposed designs and benchmarks presented in the paper. The goal is to demonstrate the effectiveness of passive MRT and interference nulling using BD-RIS.

## Code Overview
The provided code includes the following:
- **func_MRT_GC.m**: A MATLAB script that designs RIS scattering matrix to achieve passive MRT.
- **func_Null_GC.m**: A MATLAB script that designs RIS scattering matrix to achieve passive interference nulling.
- **func_MaxF.m**: A MATLAB script that designs RIS scattering matrix based on the Max-F benchmark.
- **func_Joint.m**: A MATLAB script that designs RIS scattering matrix based on the ([Joint Benchmark](https://ieeexplore.ieee.org/abstract/document/10158988)).
- **func_Prec_UP.m**: A MATLAB script that designs BS precoding matrix based on uniform power allocation.
- **func_Prec_RM.m**: A MATLAB script that designs BS precoding matrix based on RM allocation.
- **func_Prec_WF.m**: A MATLAB script that designs BS precoding matrix based on water-filling power allocation.
- **func_Prec_ZF.m**: A MATLAB script that designs BS precoding matrix based on zero-forcing.
- **func_Prec_MRT.m**: A MATLAB script that designs BS precoding matrix based on MRT.
- **alg_symuni.m**: A MATLAB script that achieves ([Symmetric Unitary Projection](https://github.com/YijieLinaMao/BD-RIS-low-complexity/tree/main)).
- **kr.m**: A MATLAB script that computes ([Khatri-Rao product](https://uk.mathworks.com/matlabcentral/fileexchange/28872-khatri-rao-product)).

### Dependencies
- MATLAB

## Routines
The routine files are used to complete the Monte-Carlo simulations. The raw results are stored in ResultsSaved folder.

## Results
The results are plotted and stored in Figures_pub folder.

## Citation
If you find this work useful in your research, please consider citing our paper:

@article{yahya2024beyond,
  title={Beyond Diagonal RIS: Passive Maximum Ratio Transmission and Interference Nulling Enabler},
  author={Yahya, Hamad and Li, Hongyu and Nerini, Matteo and Clerckx, Bruno and Debbah, Merouane},
  journal={arXiv preprint arXiv:2408.09887},
  year={2024}
}

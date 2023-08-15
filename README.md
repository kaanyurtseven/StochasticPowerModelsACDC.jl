# StochasticPowerModelsACDC.jl

StochasticPowerModelsACDC.jl is a Julia/JuMP/PowerModels package which contains an open-source code for a mathematical chance-constrained framework for solving the stochastic optimal power flow and stochastic optimal transmission switching problems for AC/DC grids. This framework takes into account HVDC controllability and uncertainties brought about by load and renewable energy generation, by using an intrusive general polynomial chaos approach.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Installation

The latest stable release of PowerModelACDC can be installed using the Julia package manager with

```
] add https://github.com/kaanyurtseven/StochasticPowerModelsACDC.jl.git
```
### Core Problem Specification
* Stochastic Optimal Power Flow for AC and AC/DC grids
* Stochastic Optimal Transmission Switching for AC and AC/DC grids

### Core Network Formulation
The code currently supports only the IVR-formulated optimal power flow and optimal transmission switching problem.

## Contributing

We welcome contributions from the community. If you would like to contribute to this project, please fork the repository and submit a pull request.

## Acknowledgements

The primary developer is Kaan Yurtseven ([@kaanyurtseven](https://github.com/kaanyurtseven)), with support from the following contributor(s):
* Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala)), KU Leuven

## Citing StochasticPowerModelsACDC

If you find StochasticPowerModelsACDC useful in your work, we kindly request that you cite the following publication:

Yurtseven, Kaan; Koirala, Arpan; Ergun, Hakan; Van Hertem, Dirk (2023): Stochastic Optimal Power Flow For Hybrid ACDC Grids Using General Polynomial Chaos. TechRxiv. Preprint. https://doi.org/10.36227/techrxiv.22606840.v1

## License

This code is provided under a BSD license.

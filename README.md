# opt-distillation

Implementation of optimization problems, test scenarios and solution strategies described in thesis "Optimization of Distillation Systems" by Juliane Ritter.

## Abstract

A chemical plant uses chemical processes in order to transform feedstock materials into products.
Chemical processes can be divided roughly into two categories: Transformations use chemical reactions or biological processes to convert available compounds into new ones, while separations exploit differences in the chemical or physical properties of the constituents of a mixture to produce two or more mixtures with different compositions.

In Germany the chemical industry accounts for about 10% of the overall primary energy demand, as can be seen in Figure 1.1. It is estimated that over 40% of the energy demand of a chemical plant is due to separation processes [1]. Therefore improving the efficiency of separation processes is an important step in reducing the overall energy demand.

In this work, we focus exclusively on the process of distillation, which is used for around 95% of all separations tasks in the chemical industry [1]. Large-scale industrial distillation is typically performed as a continuous distillation in a distillation column. Often, in order to separate a mixture of multiple components into products of desired purity, a system of distillation columns is needed. The number
of potential structures grows rapidly in the number of components present.

The goal of this work is to examine methods for finding the optimal distillation system for a given separation task. To achieve this, we follow the following steps:

In Part I we lay the foundations: We explain the thermodynamic principles behind the method of distillation and derive equations for modelling a distillation column. Further, we introduce some shortcut methods for modelling a distillation column based on suitable simplifications. Additionally, we introduce some basic definitions and results of mathematical optimization and examine different types of optimization problems and their properties and typical algorithms used in solvers.

In Part II we review literature that has treated the problem of finding the optimal distillation system. We first compare different definitions of a feasible distillation system and give results that justify restricting the feasible set to a specific type. We then examine formulations for describing this search space mathematically. Further, we briefly describe and analyse different algorithms that have been developed to search for the optimal distillation system.

In Part III we describe the two approaches we have developed and tested in this work: The first approach combines the rigorous models for distillation, introduced in Part I, with the superstructure formulation for the feasible set of distillation systems, discussed in Part II. The second approach is based on the Underwood shortcut model for distillation, introduced in Part I, and a matrix-based formulation for the feasible set of distillation systems, examined in Part II.
Firstly, we give a detailed derivation of all variables and constraints arising from both approaches. Then, we discuss the technical framework in which they were implemented and solved numerically. Finally, we compare and discuss the results obtained in both cases for multiple test scenarios and using a variety of solvers.

We finish this thesis with an outlook on further problems that arise from this work and may be examined in the future.

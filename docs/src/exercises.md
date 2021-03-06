# Exercises

## 1

Add support for equality constraints to one or more of the algorithms in `AnimatedOptimization.jl`.

## 2

Take your favoriate optimization problem (preferrably one that is not too easy, but also does not take more than a few minutes to solve), write code to evaluate the objective function and constraints, and then experiment with using various solvers. Either 

  1. Add your problem as a test to `AnimatedOptimization.jl` or
   
  2. Benchmark the performance of various solvers (Optim, Ipopt, NLopt, Knitro, JuliaSmoothOptimizers) and/or a variety of solver options.    

## 3
  
Adapt one or more of the algorithms in `AnimatedOptimization.jl` to use the API and infrastructure of one of the popular Julia optimization packages 

The popular Julia optimization packages are: 
  - [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl/)      
  - [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) which uses [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl)        
  - [JuliaSmoothOptimizers](https://juliasmoothoptimizers.github.io/#home) which uses [NLPModels.jl](https://github.com/JuliaSmoothOptimizers/NLPModels.jl)         
There is no sequential quadratic programming algorithm in any of these packages, so that would be a good one to work on.

## 4

Implement an optimization algorithm not included in `AnimatedOptimization.jl`. Of the algorithms mentioned, but not implemented, SLSQP, Augmented Lagrangian, and CMA-ES seem the most interesting and potentially useful (as far as I know, there are no quality Julia implementations for any of them, except perhaps CMA-ES).
   

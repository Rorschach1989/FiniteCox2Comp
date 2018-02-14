# FiniteCox2Comp
Finite Mixture-of-Experts Proportional Hazard Regression Model with 2 mixture components

Implementations of algorithms in [Subgroup Analysis with Time-to-Event Data Under a Logistic-Cox Mixture Model](http://onlinelibrary.wiley.com/doi/10.1111/sjos.12213/abstract)

## Toy Simulations
2 scripts for small-scale simulation of EM-Test and NPMLE estimation algorithms are available as [simulation_Testing.R](simulation_Testing.R)
and [simulation_Estimation.R](simulation_Estimation.R), see comments therein for instructions

- In EM-Test configurations, one may make inference based on asymptotic distributions of the test statistics under some situations
but the paper pointed out this is inaccurate in small-sample cases, bootstrap testing algorithms are safter choices, which results in heavy computational burden

- In NPMLE estimations of Finite Mixture-of-Experts Proportional Hazard Regression Models, estimation of the asymptotic variance of the estimator
is possible via inverting the Fisher information matrix of corresponding to the observed likelihood function, which is usually O(n) by O(n) in size, where n is the sample size. While this memory complexity and computational complexity of inverting the matrix
could be reduced greatly via careful decomposition of the matrix and applying some matrix computation tricks, they're currently not implemented here. For moderate sample-sizes (up to ~1000), brute-force implementation works not that bad. Bootstrap estimation is provided as well. 

## Notes
- 2 components finite-mixture seems not "expressive enough" for a finite-mixture use case, general estimation and influencial schemes are trivial generalizations of the approach used here
- Cox model is inadequate or inappropriate under phenomenons like crossing hazards, we're gonna wrap a much more general interface concerning
semiparametric linear transformation models or even frailty models into a fully-functioned R package. While I'm not sure about the release date :(

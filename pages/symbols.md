# Guide to Common Symbols Used in Notation

## λ
Lambda. Typical label used for the instantaneous rate in a Poisson process. For us, it will usually mean the instantaneous birth rate of lineages, per time-unit, per lineage. This 'birth' rate is also commonly known as the 'per-capita' origination rate, speciation rate, branching rate, cladogenesis rate, etc, depending on the context. 

The 'per time-unit, per lineage' aspect means the rates are relative to time and the number of lineages extant at that time (this is what is meant when rates are referred to sometimes as 'per capita'). For example, if there is a group with a hundred lineages, and another group containing only a single lineage, but they had identical birth rates and were diversifiying over the same time interval, the first group is a *hundred* times more likely to give arise to additional lineages, relative to the group with the single lineage. This variable is sometimes denoted as 'p' in paleontological literature.

## μ
Mu. Typically used as label for the instantaneous death rate (extinction rate) per time-unit, per lineage. This variable is sometimes denoted as 'q' in paleontological literature.
	
## ψ
Psi. Typically used in Fossilized Birth Death models for the instantaneous 'serial sampling' rate per time-unit, per lineage. Serial sampling means the sampling process over time - i.e. the chance of a lineage leaving behind remains that get preserved and survive to the modern day, so that we observe them. The serial sampling rate is sometimes denoted as 'r' in paleontological literature, and is related to (but not the same as) the per-time-interval sampling probability, often denoted 'R' (e.g. Foote & Raup, 1996). 

Sometimes, Psi is also used as a symbol for the dated tree we are trying to infer (sometimes this is instead denoted with Tau, below). Occassionally, authors will use Psi for both the serial sampling rate and the dated tree in the same manual (This is typical for RevBayes documentation).
	
## ρ
Rho. Typically used for probability of sampling a species at time=zero (e.g. the modern).
	
## τ
Tau. Often used as a symbol representing the dated phylogeny with branch lengths that we are trying to infer. (Some literature uses Psi to denote the tree.)
	
## Γ
Gamma. Typically used as an abbreviation for the Gamma distribution.
	
## β
Beta. Typically used as an abbreviation for the Beta distribution.
	
## σ
Sigma (lower-case). Symbol used for various parameters, such as the standard deviation of a normal distribution, the log-standard deviation of a lognormal distribution, or the rate of trait change in a non-directional model of trait change (essentially, the expected variance across a time-unit), etc. 
	
## α
Alpha.  In models of Markov processes, such as the model used to describe character change, alpha is used to denote the instantaneous rate of any particular transition between states (Lewis, 2001). In models of continuous trait evolution, alpha is often used to denote the strength of attractors.
	
## r 
Lower-case r, from the Latin Alphabet. Small r is often used as the notation for the probability of serial sampling events also causing an immediate death event in that lineage. In the fossil record, this would that preservation was also tied to the extinction of that lineage. When r=1, there can be no sampled ancestors. Most implementations of the the Fossilized Birth Death model treat r=0, and it is more commonly referenced in models of disease transmission.
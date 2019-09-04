This repository will serve as a place for developing and storing materials for the Revbayes-based short course in tip-dating phylogenetics. for GSA 2019, presented by The Paleontological Society.

This hands-on computational workshop is organized by Sandra Carlson and Peter Wagner (@PeterJWagner3), instructors April Wright (@wrightaprilm), Laura Soul (@laurasoul), David "Davey" F. Wright (@daveyfwright), and David Bapst (@dwbapst), with instructional material development from Rachel Warnock (@rachelwarnock). Additionally, this course will utilize a number of 'roving' teaching assistants (tentative list, as follows: William Gearty (@willgearty), Lucy Chang (@lucymchang), ...).

## Current Schedule:

- 08:00 AM - Module 1 - Introduction: Tree-Thinking
	- Sandy C., Pete W.
	- 15 min, Lecture.	
- 08:15 AM - Module 2 - Meet Our Data: Smith and Zamora, 2009
	- Davey W.
	- 10 min, Lecture.
	- Brief discussion about data set
	- Explaining what cinctans are and why echinoderms are cool
- 08:25 AM - Module 3 - Morphological Character Coding
	- Sandy C., Davey W., Pete W.
	- 30 min, Lecture & exercise
	- Mesquite exercise with S&Z 2009 data set.
	- Experiment - Examining character coding decisions
	- Possible consequences of different coding options
	- Sandy: I think the first hour could be reduced to half an hour, or even 20 minutes if that would help with the rest of the timing.	
	- From Pete:
		- A. Emphasis: Ways to identify particular conditions that we think are shared (or potentially shared) via common ancestry;
		- Z. But even if my characters are “good,” what if they do not all evolve in the same way all of the time?
- 09:00 AM - Module 4 - Intro to RevBayes, Graphical Modeling
	- April W, ?
	- 45 min, Lecture & exercise
	- Plate notation?
	- From Pete:
		- A. Emphasis: We want a “flow chart” detailing all of the variables that affect how character states are distributed in the fossil record.
		- B. Show RevBayes symbols for different types of variables
		- Z. But how do I get a tree without assuming that all of these other parameters have particular conditions?
- 09:45 AM - BREAK (15 min)	
- 10:00 AM - Module 5 - Tripartite Model 1: Morphological Character Change Models 
	- April W, ?
	- Lecture & Hands-On with RevBayes/RevNotebook.	
	- Likelihood concepts, Bayesian Methods, MCMC
	- Simple Mk, filtering
	- Teasers of how to relax some assumptions of the model:
		- transition rate symmetry?
		- across-site/char heterogeneity?
	- April: Learners first contact with RB; stuff will go wrong. Also a part where there is likely to be significant discussion, since no one really believes the Mk model. 
		- DWB: We've already seen this in comments in our registration survey. There will definitely be non-believers about Mk.
	- Push up right before lunch so people can yell at April afterwards
	- April: I think I'll have them work as a table: one person do an MCMC with the strict Mk model; others can each pick one thing to change. Maybe choose a different distribution for among-site rate variation, or allow asymmetrical transitions. They can compare results during/after lunch.
	- April: I'm going to turn my morph models talk in to an R Markdown document so that I can interweave activities with lecture.
	- April's materials from previous workshop: 		https://github.com/danlwarren/Evolution-2019-Phylogenetic-Methods-Workshop/tree/master/3%20-%20April%20Wright%20-%20Teaching%20and%20resesarch%20with%20RevBayes
	- From Pete:
		- MCMC
			- A. Emphasis: Way to search complicated hypothesis space and and marginalize effects of unknown rate parameters.
			- Z. But how do I actually do the math at any step?
		- MK Models (Tripartite Model 1a)
			- A. Emphasis: Probability of deriving a matrix given a particular phylogeny (cladistic relationships + divergence times)
			- B.  MCMC allows us to “integrate” over different general rates.
			- Z. But what if rates vary among characters?
		- Gamma & lognormal variation & hyperpriors (Tripartite Model 1b)
			- A. Emphasis: Assume distribution of rates and thus vary rates without specifying which characters are “good” or “bad”.
			- B. Introduces  hyperpriors for varying lognormal or Gamma distributions 
			- Z. But I keep seeing “time” in these equations: what if these rates are not consistent over time?	
- 12:00 PM - LUNCH (1.5 hour)
- 01:30 PM - Module 6 - Tripartite Model 2: Clock Models for Character Change
	- Who? Davey W.?
	- 45 min, Hands-On with RevBayes/RevNotebook.
	- April: light, hands-on ... this is just sorta boring
	- Pete: tease looking at different rates across morph partitions?
	- From Pete:
		- A. Emphasis: similar to the “set rate” vs. “distributed rate” issue and with similar solutions.
		- Z. But what if my clade spans “exceptional” events (e.g., major radiations or major global shifts affecting lots of clades), or what if my characters are “cheating” in some other way beyond what among-branch & among-character variation models?	
- 02:00 PM - Module 7 - Tripartite Model 3: Fossilized Birth Death
	- Dave B., Davey W.?, Laura S.?
	- 1:30 min, Lecture & Hands-On with RevBayes/RevNotebook
	- Fossilized Birth Death Models 
		- Stadler, 2010; Heath et al., 2014; Stadler et al., 2018
		- FBD, Range-FBD
	- Anc-Desc Relationships, Sampled Ancestors in Tip-Dating
		- Sampled Ancestor Moves in MCMC
		- Gavryushkina et al. 2014; 2016
	- Try to simplify, do not get lost in covering historical perspective
		- April: Simulate under a pure-birth process and ask whose data looks like this? Simulate under the FBD and show that this process allows sampled ancestors? 	
	- Unifies morph, clocks, tree priors
	- Goes from start to end of an analytical process
	- From Pete:
		- BD & FBD (BDS) models for putting prior probabilities on tree topologies offer independent check on characters (Tripartite Model 3).
		- A. Emphasis: Probabilities getting fossil record given phylogeny based on likelihood of unsampled ancestors AND probability of no sampled sister taxa (aka, Raup’s Revenge).
		- B. Our old friend hyperpriors used to create distributions for these rates.
		- Z. Wow, that’s cool: but I’m now using everything and the kitchen sink; can I do anything with the results?	
- 03:30 PM - BREAK (15 min)	
- 03:45 PM - Module 8 - Worked Example with Phylo Comp Methods
	- Laura S., Davey W.
	- 1:15 min, Lecture & Hands-on with RevBayes
	- From Pete:
		- A. How to test macroevolutionary ideas while admitting that we can never know the “true” tree
		- Z. But what about X, Y, Z, etc.?
	- Laura: Being able to show something (even if its simple) that one can do whether you made your own tree, generated a supertree, or made most of a tree but want to include 10 taxa you can't code, would be the most broadly useful to the (clearly very mixed experience) group we have coming. From an earlier chat with Davey, he pointed out that we could also instead take the approach of doing a light touch overview of what is now possible in the paleo-PCM world and highlighting the R packages that are available to do it. We could also stick to modelling rates of discrete character change on the tree if switching to continuous characters feels too 'from scratch'.
	- Peter: As Laura noted, those approaches rely heavily on divergence times, which means that you have to rehash (and reinforce!) the FBD & tip-dating material from earlier parts; many PCM approaches integrate over trees, too, which gets back to the MCMC approach.  And both the FBD/tip-dating and methods for assessing how characters evolve over trees get back to the morphological evolution modules.  There could be a lot of “As April / Laura / David / Davey showed us earlier, we do XXX here” parts that could really tie the whole thing together.  And, of course, it would also get back to something you certainly will raise from the start: why we want decent estimations of phylogeny in the first place.
	- Sandy: 
		- Cooper, Natalie, Gavin H. Thomas, and Richard G. FitzJohn. "Shedding light on the ‘dark side’of phylogenetic comparative methods." Methods in ecology and evolution 7.6 (2016): 693-699.
		- https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12533
	- Laura: Yup this article highlights important points, there are a couple of others that highlight similar things (e.g. Uyeda). It will be important for us to emphasise that throwing loads of pcms at a single topology in R is not a great idea, but that integrating over tree and trait uncertainly for carefully chosen methods that are appropriate for the question is. 
	- I’d like Davey and I to be able to emphasise their power and necessity, but also the importance of choosing your approach carefully. 
	- April: I think you two have exactly the right approach. It's really easy to get bogged down in the PCM debates. Should we be using BAMM? Does Ohrenstein-Uhlenbeck even work? Those are important questions! But it's probably important that learners walk away with some roadmap to even basic points like which are process-based, and which are more neutral to the underlying evolutionary model? 
		- If you don't want to do a single point estimate tree, how can you avoid that (maybe treeplyr in R for working with tree vectors, co-estimation of  model + tree in RevBayes) see here for an example: 
		- https://revbayes.github.io/tutorials/chromo/ 
	- Peter: What might be done with the cinctan dataset?  Smith & Zamora separated “feeding” characters (or at least food groove characters) from others, and there certainly have been PCM studies contrasting how one class of characters evolves vs. another.  We probably could gather basic size data, too: e.g., dimensions of the theca or whatever you call the “body” of these things.  The benefit of that is that there might be a prize for the 1 millionth tree-based test of body size evolution!  
- End-of-Day Summary 
	- Who? Peter W.? Laura S.? Davey W.?
	- When? Unknown duration? 1.25 hrs is gonna be very tight for PCMs as it is...
	- Conclusions, Inspirations
	- What Paleontologists Can Contribute to Systematics
	
## Readings

Attendees will be asked to read over three overviews of the field beforehand:

Holder, M., and P. O. Lewis. 2003. Phylogeny estimation: traditional and Bayesian approaches. Nat Rev Genet 4(4):275-284.

Nee, S. 2006. Birth-death models in macroevolutin. Annual Reviews in Ecology, Evolution and Systematics 37:1-17.

Wright, A. M. 2019. A Systematist’s Guide to Estimating Bayesian Phylogenies From Morphological Data. Insect Systematics and Diversity 3(3).
	Link: https://academic.oup.com/isd/article/3/3/2/5519658
	Preprint Here: https://paleorxiv.org/jupva

Wright, D. F. 2017. Bayesian estimation of fossil phylogenies and the evolution of early to middle Paleozoic crinoids (Echinodermata). Journal of Paleontology 91(4):799-814.



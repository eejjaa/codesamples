Selected recent code samples by Emerson Arehart, PhD candidate at the University of Utah.


emerson.arehart@utah.edu

Power Analysis and Experimental Design for Salivary Proteins (RS_pwr_analysis.R) << work in progress >>

Our collaborators are designing an experiment in which they maintain a group of rats on a 'control' diet for several days, then switch them to an experimental diet. Each day, a saliva sample is taken from each rat and assessed for the volume of each size class of protein. We will be building a functional model of the processes which underpin the change in salival protein profile associated with each diet, but for now (in this code) I use data from previous experiments to build generalized additive models of the change in volume (for each protein size class) over time, and use a likelihood ratio test to compute the appropriate p-value. By re-sampling the individuals from the data set at different sample sizes, I can determine the minimum sample size needed to obtain significant results for all bands.

Simulation of a Population of Wild-Harvested Mushrooms Under Habitat Change (loyo_example.R) << work in progress >>

This code models the population dynamics of Boletus loyo, an archetypal wild edible mushroom endemic to the Patagonian forests of Chile, where it is an obligate mycorrhizal symbiont of Nothofagus trees. Our model incorporates the effects of spatial distribution and genetic diversity to model the growth and spread (or decline) of a population of B. loyo, under various harvesting strategies and habitat change scenarios. By parameterizing the model with data from the field, we hope to determine which harvesting strategies are most sustainable under different climate change scenarios. 

Glomerulus - Kenyon Distance Metrics (gk_distance.R) << work in progress >>

In the fruit fly olfactory system, odors activate different olfactory receptor neurons, which activate corresponding glomeruli in the brain. As such, a given odor is 'represented' by the combination of glomeruli that it activates. However, the glomeruli are in turn connected to Kenyon Cells (KCs) in the Mushroom Body of the brain, a region responsible for associative learning. An odor is thus also 'represented' by a combination of KCs, and that representation depends on the wiring of these glomeruli-KC connections. This code generates a wiring pattern, uses it to translate odors from glomerulus to KC representations, and then applies metrics to the representations - with a goal of understanding the effect that different wirings have on attributes of the KC representations.

Handicap Principle Simulation (handicap.R)

This project simulates a population of asexually-reproducing organisms which use (and inherit) different strategies for parent-offspring signaling. Some offspring accurately signal their need for food to the parent, while some exaggerate; some parents match the food request of their offspring exactly, others underfeed. This simulation serves as a sandbox for ideas related to honest signaling theory and the Handicap Principle, as proposed by Zahavi, Godfray, and others.

Heterogenous Neuron Decision Making (heterogenous_neurons.R)

This code simulates a population of neurons with different 'strengths' of connectivity to each other as they respond to a stimulus (and to input from each other) to reach a collective decision. The code generates plots for various metrics of the decisions at different levels of background noise, heterogeneity in connections between neurons, and overall strength of connectivity.

Cooperation and 'Cheating' in Yeast (yeast_coop.R)

This project, my first effort in the R language, simulates a population of yeast as they multiply (by budding) and consume resources that diffuse through a 2D environment. Some of the yeast 'cooperate' by expending energy to convert sucrose into glucose, which is the food source for the yeast; others 'cheat' or 'free-load' by consuming the glucose produced by other yeast, but not expending energy to convert any sucrose. The simulation can be run with different intial populations of cooperators and cheaters, demonstrating how the population dynamics play out under different blends of strategies.

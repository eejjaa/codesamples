Selected recent code samples by Emerson Arehart, PhD candidate at the University of Utah. 
emerson.arehart@utah.edu




Glomerulus - Kenyon Distance Metrics (gk_distance.R) << work in progress >>

In the fruit fly olfactory system, odors activate different olfactory receptor neurons, which activate corresponding glomeruli in the brain. As such, a given odor is 'represented' by the combination of glomeruli that it activates. However, the glomeruli are in turn connected to Kenyon Cells (KCs) in the Mushroom Body of the brain, a region responsible for associative learning. An odor is thus also 'represented' by a combination of KCs, and that representation depends on the wiring of these glomeruli-KC connections. This code generates a wiring pattern, uses it to translate odors from glomerulus to KC representations, and then applies metrics to the representations - with a goal of understanding the effect that different wirings have on attributes of the KC representations.

Handicap Principle Simulation (handicap.R)

This project simulates a population of asexually-reproducing organisms which use (and inherit) different strategies for parent-offspring signaling. Some offspring accurately signal their need for food to the parent, while some exaggerate; some parents match the food request of their offspring exactly, others underfeed. This simulation serves as a sandbox for ideas related to honest signaling theory and the Handicap Principle, as proposed by Zahavi, Godfray, and others.

Heterogenous Neuron Decision Making (heterogenous_neurons.R)

This code simulates a population of neurons with different 'strengths' of connectivity to each other as they respond to a stimulus (and to input from each other) to reach a collective decision. The code generates plots for various metrics of the decisions at different levels of background noise, heterogeneity in connections between neurons, and overall strength of connectivity.

Cooperation and 'Cheating' in Yeast (yeast_coop.R)

This project, my first effort in the R language, simulates a population of yeast as they multiply (by budding) and consume resources that diffuse through a 2D environment. Some of the yeast 'cooperate' by expending energy to convert sucrose into glucose, which is the food source for the yeast; others 'cheat' or 'free-load' by consuming the glucose produced by other yeast, but not expending energy to convert any sucrose. The simulation can be run with different intial populations of cooperators and cheaters, demonstrating how the population dynamics play out under different blends of strategies.
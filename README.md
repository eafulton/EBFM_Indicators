This quarto (R) code and supporting documentation describes how to generate some files useful for considering the state and function of marine ecosystems, which can be useful in taking an ecosystem oriented approach to fisheries management. It assumes you are obtaining your “data” from Ecopath with Ecosim (EwE) and processing it with R quarto code EBFM_indicator_calculation.qmd

There are 2 stages to this process.
1.	Generating the base information – typically done via EwE (as described here), though it is possible to generate it from data if you have sufficiently long and comprehensive data time series (this is not described in this instance)
2.	Calculating the hub species and ecosystem indicators (using the R code in EBFM_indicator_calculation.qmd)

Full explanation of all the steps is provided in the documentation file - Calculating Hub index and ETI from Ecopath with Ecosim.docx
A short worked example for SE Australia is provided so you can see what the set-up looks for a model and some example files. You will need to execute it on your machine to create a Rendered html file of output, but example png files are in the output folder so you can see what they should look like as you create your own.

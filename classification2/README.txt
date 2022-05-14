Files in this directory:

simulation.root : copy of events generated from the S and B models
dataBigSig.root : data sample with largest S/B ratio
dataSmSig.root  : data sample with small S/B ratio
dataVSmSig.root : data sample with very small S/B ratio

VarAnalyzer.C : Example code to analyze S,B variables
	        Creates simulation_vars.root
Usage in ROOT:
		root[] .L VarAnalyzer.C+
		root[] VarAnalyzer()
		// plots S,B variables, correlations, calculates and plots
		// s/sqrt(s+b) for different cut values

		root[] VarAnalyzer(0.015)
		// Run VarAnalyzer with a scale factor and the signal will
		// be scaled.  This allows finding optimal cuts with smaller
		// S/B ratios
		
		
Selection.C : Example code to apply selections to S,B variables
	        Calculates signal significance after applying seletions
		to the S and B samples and creates simulation_sel.root
Usage in ROOT:
		root[] .L Selection.C+
		root[] Selection()
		root[] Selection(0.015)
		// as above, performs calculations w/ signal sample weighted


SelectD.C : Example code to apply selections to Data variables
		creates <datafileName>_sel.root containing the data
		distributions following the selections.
		You cause the distributions in simulation_sel.root and
		<datafileName>_sel.root to perform templated fits to measure
		the signal.
Usage in ROOT:
		root[] .L Selection.C+
		root[] Selection()
		root[] Selection(0.015)
		// as above, performs calculations w/ signal sample weighted

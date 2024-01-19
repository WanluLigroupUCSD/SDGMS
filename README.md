# SDGMS
Modified Basin Hopping algorithm for fast global minimum optimization

Usage:
-b basis name
-m method name (DFT,HF,etc.)
-l previous file name
-f output file
-n taskName (can be anything you want, no spaces)
-comp composition (a composition like B3 or H2O)
-pr radial criteria percent (this is the percentage of similarity allowed to determine if two structures are identical. should be lower than percent change i)
-pi percent change i (the lowest percent Xnew deviates from S)
-pf percent change f (the highest percent Xnew deviates from S)
-charge charge
-s state
-o optimization cycles (for gaussian optimization. recommended is 20-30)
-t time program can run, an float number of hours
-stepsC coordination steps
-stepsL local minima steps
-xyz x y z box that the molecule is bound to
-ABC cell shape. As there is only currently support for Gaussian this does not do much, but in the future its needed for CP2K.

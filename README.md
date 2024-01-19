# SDGMS
Modified Basin Hopping algorithm for fast global minimum optimization

Usage:
<br>-b basis name
<br>-m method name (DFT,HF,etc.)
<br>-l previous file name
<br>-f output file
<br>-n taskName (can be anything you want, no spaces)
<br>-comp composition (a composition like B3 or H2O)
<br>-pr radial criteria percent (this is the percentage of similarity allowed to determine if two structures are identical. should be lower than percent change i)
<br>-pi percent change i (the lowest percent Xnew deviates from S)
<br>-pf percent change f (the highest percent Xnew deviates from S)
<br>-charge charge
<br>-s state
<br>-o optimization cycles (for gaussian optimization. recommended is 20-30)
<br>-t time program can run, an float number of hours
<br>-stepsC coordination steps
<br>-stepsL local minima steps
<br>-xyz x y z box that the molecule is bound to
<br>-ABC cell shape. As there is only currently support for Gaussian this does not do much, but in the future its needed for CP2K.

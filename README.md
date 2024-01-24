# SDGMS
Modified Basin Hopping algorithm for fast global minimum optimization

Usage:
<br>-b basis name
<br>-m method name (PBE,HF,etc.)
<br>-d calculator used for energy calculations (g for gaussian, b for basin hopping formula. No other support at this moment)
<br>-l previous file name
<br>-f output file
<br>-n taskName (can be anything you want, no spaces)
<br>-comp composition (a composition like B3 or H2O)
<br>-pr radial criteria percent (percentage of deviation from known covalent radii of elements in structure allowed. Recommended at least 40%)
<br>-pi percent change i (the lowest percent Xnew deviates from S based on xyz)
<br>-pf percent change f (the highest percent Xnew deviates from S based on xyz)
<br>-charge charge
<br>-s state
<br>-o optimization cycles (for gaussian optimization. recommended is 20-30)
<br>-t time program can run, an float number of hours
<br>-stepsC coordination steps
<br>-stepsL local minima steps
<br>-stepsH steps before switching from hartree fock to method of choice
<br>-stepsS steps for maximum scf in optimization
<br>-xyz x y z is the maximum x, y ,z values that a 3d coordinate can take during a displacement.
<br>-ABC cell shape. As there is only currently support for Gaussian this does not do much, but in the future its needed for CP2K.

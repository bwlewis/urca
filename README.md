# urca
Pfaff/Stigler Unit Root and Cointegration Analysis R package

## Warning

This is a work in progress...don't install this package (yet)! Use the one from
CRAN. Unless that is, you want to help, then please have at it.

## Intention

We're revising the package to enhance some of the numeric computations for
stability and efficiency along the lines suggested by Doornik and O'Brien in
their superb paper http://www.doornik.com/research/CointNum2_final.pdf.

Our variation will use a generalized SVD. Because Pfaff's package is so
comprehensive at handling trends, constants, etc., there are quite a few
details to work through but the basic idea looks solid.

We're starting with the `ca.jo` function first. Eventually, almost everywhere
you find the function `solve` in the source code will include new routines.

We intend to make the new code paths optional to make it easy to compare with
the standard Johansen procedure.



## Status

No code yet! But my hand-written notes are online. Code will follow...

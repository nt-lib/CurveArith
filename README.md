# CurveArith

A Magma package for computing class groups and gonalities of curves over finite fields.

## Requirements

This packages requires Magma version v2.29-1 or higher.

After installation, it is recommended to run the automated tests, as described below, in order to see whether CurveArith works correctly with your version of Magma.

Some older versions are only partially supported. The end user code has been verified to work in Magma versions as old as v2.28-3. If you want to reproduce the timing data then Magma v2.28-19 or higher is needed due to the usage of [Fork](https://magma.maths.usyd.edu.au/magma/handbook/text/42#426).

## Usage

To load the package into Magma, run
```
AttachSpec("<path to installation>/CurveArith.spec");
```
Alternatively, you can add the above to your `.magmarc` file to load the package automatically upon Magma startup.

This package provides three intrinsics:
- `CAClassGroup(C) : Crv[FldFin] -> GrpAb`  
    `CAClassGroup(F) : FldFun -> GrpAb`

    Computes the divisor class group of the curve C or the function field F.

    Takes three optional parameters:
    - `BaseDivisor: DivCrvElt/DivFunElt`  
        The base divisor defining the Riemann-Roch space in which to look for class group relations. If this is too small, the computation may not terminate. If omitted, an appropriate divisor is chosen heuristically.
    - `FactorBasisDegree: RngIntElt`  
        The degree bound for the places in the factor basis. If this is too small, the computation may return an incorrect result. If omitted, the smallest degree giving a provably correct result is used.
    - `MaximumTime: RngReSubElt`  
        The maximum time (in seconds) after which the calculation is aborted. If omitted, no time limit is used.

- `CAGonality(C) : Crv[FldFin] -> RngIntElt`  
    `CAGonality(F) : FldFun -> RngIntElt`

    Computes the gonality of the curve C or the function field F.

    Takes three optional parameters:
    - `Bound: RngIntElt`  
        An upper bound for the degrees of functions to look for. If provided, the return value will be `Bound + 1` if the gonality exceeds `Bound`. Otherwise, or if omitted, the return value equals the gonality.
    - `Al: MonStgElt`  
        The method to use for computing Riemann-Roch dimensions. The default, `"LinAlg"`, uses a linear algebra-based procedure described in [[1]](#citing). Can be set to `"Hess"` instead to use Magma's `Dimension` intrinsic instead, which may be faster on very small examples.
    - `MaximumTime: RngReSubElt`  
        The maximum time (in seconds) after which the calculation is aborted. If omitted, no time limit is used.

- `CAHasFunctionOfDegreeAtMost(C, d) : Crv[FldFin], RngIntElt -> BoolElt`  
    `CAHasFunctionOfDegreeAtMost(F, d) : FldFun, RngIntElt -> BoolElt`

    Returns whether there is a function on the curve C, or in the function field F, of degree at most d.

    Takes four optional parameters:
    - `Al: MonStgElt`  
        The method to use for computing Riemann-Roch dimensions. The default, `"LinAlg"`, uses a linear algebra-based procedure described in [[1]](#citing). Can be set to `"Hess"` instead to use Magma's `Dimension` intrinsic instead, which may be faster on very small examples.
    - `MaximumTime: RngReSubElt`  
        The maximum time (in seconds) after which the calculation is aborted. If omitted, no time limit is used.
    - `TimingData : BoolElt`  
        If set, returns an additional value containing detailed timing data for the calculation.
    - `StopAfterFirst : BoolElt`  
        If set, iterates through all possible divisors of degree d, even after a function of degree at most d has been found. Mainly used for timing scripts.

## Automated testing

The code in this repository can be tested automatically. If you have `make` installed this can be done as follows:
```shell
make test
```

The tests can also be run manually:
```shell
cd tests
magma -n src/test_all.m
```
At the moment it is only possible to run the tests from the `tests` folder. Running them from any other folder will result in an error. This is due to the way `load` works in Magma.

## Timing

The timing data used in [[1]](#citing) can be found in `timings/output`. The timing data can be reproduced by running the following from the root directory:
```shell
cd timings
magma -n src/time_all.m
```
At the moment it is only possible to run the timings from the `timings` folder. Running them from any other folder will result in an error. This is due to the way `load` works in Magma.

## Citing

If this repository has been useful in your research, please cite the following preprint.

[1] Maarten Derickx and Kenji Terao, *Computing class groups and gonalities of algebraic curves over finite fields*. 2026. arXiv: [2602.17417](https://arxiv.org/abs/2602.17417)

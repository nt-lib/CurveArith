# CurveArith
Computes gonalities and class groups of curves over finite fields using Magma

## Requirements

This package requires magma v2.28-3 or higher. If you want to reproduce the timing data then magma v2.28-19 or higher is needed due to the usage of [Fork](https://magma.maths.usyd.edu.au/magma/handbook/text/42#426) .

## Usage

Start magma from the root directory of this repository and run:

```
AttachSpec("CurveArith.spec");
```

in order to load this magma package.

## Automated Testing

The code in this repository can be tested automatically. If you have make installed this can be done as follows:

```shell
make test
```

The tests can also be run manually:

```shell
cd tests
magma -n src/test_all.m
```

At the momemnt it is only possible to run the tests from the tests folder. Running them from any other folder will result in an error. The reason for this is how `load` works in magma.

## Timing information

The timing information can be reproduced as follows.

```shell
cd timings
magma -n src/time_all.m
```

At the momemnt it is only possible to run the timings from the timings folder. Running them from any other folder will result in an error. The reason for this is how `load` works in magma.




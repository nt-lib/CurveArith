# CurveArith
Computes gonalities and class groups of curves over finite fields using Magma

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




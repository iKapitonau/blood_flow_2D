# Blood flow modelling
2D blood flow numerical modelling program in a vessel without defects.

## Requirements

* gcc
* make
* libconfig
* openmp

## Compilation
Release target is default:
```Bash
make
```
Parallel version (with OpenMP):
```Bash
make parallel
```
Debug version:
```Bash
make debug
```

## Execution
```Bash
./blood_flow_2d [config_filename]
```
Output values in every layer of the calculation grid will be written to `text` file.

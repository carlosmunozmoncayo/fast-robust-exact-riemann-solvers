#!/bin/bash
rm modules/*
for elemento in $(ls hstar*)
do
    foo=${elemento#"hstar_"}
    foo=${foo%".f90"}
    python3 -m numpy.f2py -c --fcompiler=gnu95 $elemento -m  $foo
    mv *$foo*.so modules
done

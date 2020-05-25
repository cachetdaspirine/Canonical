#!/bin/bash

SimNum=0

rm -rf "Res/Sim$SimNum"
mkdir "Res/Sim$SimNum"

cp Parameter.py "Res/Sim$SimNum/Parameter.py"

sbatch CanonicalAnnealing.pbs $SimNum

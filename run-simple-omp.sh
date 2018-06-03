#!/bin/bash

## Specifies the interpreting shell for the job.
#$ -S /bin/bash

## Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -V

## Specifies the parallel environment
#$ -pe smp 4

## Execute the job from the current working directory.
#$ -cwd 

## The  name  of  the  job.
#$ -N omp3a32th

##send an email when the job ends
#$ -m e

##email addrees notification
#$ -M email@alumnes.udl.cat

##Passes an environment variable to the job
#$ -v  OMP_NUM_THREADS=32

## The folders to save the standard and error outputs.
#$ -o $HOME
#$ -e $HOME

## In this line you have to write the command that will execute your application.
./ompcompiled test_trees3.txt $HOME/omp3a32th




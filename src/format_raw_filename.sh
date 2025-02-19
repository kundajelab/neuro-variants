#!/bin/bash

indir=$1
echo $indir

rename 's/ /_/g' $indir/*

rename "s/^'//" $indir/*

rename "s/'$//" $indir/*

rename 's/\?//g' $indir/*


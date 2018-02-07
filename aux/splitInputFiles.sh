#!/bin/bash
##num=0;for i in `seq 1 20 230`; do num=$((num + 1)); echo "Writing Split File VidFilesToProcessSplit$num.txt";sed -n "$i,$((i+procsize-1)) p" startVidsToProcess.txt > "startVidFilesToProcessSplit$num.txt"; done

INPUTFILE=$1
PARPROC=$2
lines=`wc -l $INPUTFILE | cut -f1 -d' '`
#procsize=$lines

procsize=`echo $lines $PARPROC | awk '{print int(($1+1)/$2)}'`
echo "Splitting for $PARPROC , each $procsize"

num=0
##Break It Down To Processing Files

for i in `seq 1 $procsize $lines`; do num=$((num + 1)); echo "Writing Split File $INPUTFILE Split $num";sed -n "$i,$((i+procsize-1)) p" $INPUTFILE > "$INPUTFILE"Split"$num"; done

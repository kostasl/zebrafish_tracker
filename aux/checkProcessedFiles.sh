#!/bin/bash
echo "Kostas Lagogiannis 2017 - Check Track Videos Progress Utility"
echo "Finds Which videos at param dir do not have tracked .csv file in the current directory . ie finds which videos have not been processed yet"


DATDIR=$1
VIDDIR=$2
PARPROC=$3

echo $VIDDIR
echo $DATDIR

rm datafileslist.txt
rm vidfilelist.txt
find $DATDIR -maxdepth 3 -name "*.csv" | sed 's,^[^/]*/,,' | sed s/\.[^\_]*$// | sed s/\.[^\_]*$// | sed 's/.*\///' |sed '/^$/d' | sort | uniq > datafileslist.txt
find $VIDDIR -maxdepth 3 > vidfilesfullpath.txt
cat vidfilesfullpath.txt | sed 's,^[^/]*/,,' | sed s/\.[^\.]*$// | sed 's/.*\///' |sed '/^$/d' | sort | uniq > vidfilelist.txt
diff --changed-group-format='%>' --unchanged-group-format='' datafileslist.txt vidfilelist.txt > unprocessedfiles.txt


echo "There are :"
wc -l unprocessedfiles.txt
echo " unprocessed video files in $VIDDIR" 
echo "Here are the list of video file that have not been tracked/processed yet:"
echo "--------------------------"

grep -f unprocessedfiles.txt vidfilesfullpath.txt > VidFilesToProcess.txt

lines=`wc -l VidFilesToProcess.txt | cut -f1 -d' '`
#procsize=$lines

procsize=`echo $lines $PARPROC | awk '{print int(1+($1)/$2)}'`
echo "Splitting for $PARPROC , each $procsize"

num=0
##Break It Down To Processing Files

for i in `seq 1 $procsize $lines`; do num=$((num + 1)); echo "Writing Split File VidFilesToProcessSplit$num.txt";sed -n "$i,$((i+procsize-1)) p" VidFilesToProcess.txt > "VidFilesToProcessSplit$num.txt"; done

#more unprocessedfiles.txt


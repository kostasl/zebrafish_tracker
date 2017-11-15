#!/bin/bash
echo "Kostas Lagogiannis 2017 - Check Track Videos Progress Utility"
echo "Finds Which videos at param dir do not have tracked .csv file in the current directory . ie finds which videos have not been processed yet"


DATDIR=$1
VIDDIR=$2

echo $VIDDIR
echo $DATDIR

rm datafileslist.txt
rm vidfilelist.txt
find $DATDIR -maxdepth 1 | sed 's,^[^/]*/,,' | sed s/\.[^\_]*$// | sed s/\.[^\_]*$// | sed 's/.*\///' | sort | uniq > datafileslist.txt
find $VIDDIR -maxdepth 1 | sed 's,^[^/]*/,,' | sed s/\.[^\.]*$// | sed 's/.*\///' | sort | uniq > vidfilelist.txt
diff --changed-group-format='%>' --unchanged-group-format='' datafileslist.txt vidfilelist.txt > unprocessedfiles.txt

echo "There are :"
wc -l unprocessedfiles.txt
echo " unprocessed video files in $VIDDIR" 
echo "Here are the list of video file that have not been tracked/processed yet:"
echo "--------------------------"

more unprocessedfiles.txt


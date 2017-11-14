#!/bin/bash
echo "Kostas Lagogiannis 2017 - Check Track Videos Progress Utility"
echo "Finds Which videos at param dir do not have tracked .csv file in the current directory . ie finds which videos have not been processed yet"


DATDIR=$1
VIDDIR=$2

echo $VIDDIR
echo $DATDIR

rm datafileslist.txt
rm vidfilelist.txt
find $DATDIR -maxdepth 1 | sed 's,^[^/]*/,,' | sed s/\.[^\_]*$// | sed s/\.[^\_]*$// | sort | uniq > datafileslist.txt
find $VIDDIR -maxdepth 1 | sed 's,^[^/]*/,,' | sed s/\.[^\.]*$// | sed 's/.*\///' | sort | uniq > vidfilelist.txt
diff datafileslist.txt vidfilelist.txt --changed-group-format='%>'

echo "**file listed are the video file names in source video directory for which a datafile is missing.** "
echo "--------------------"

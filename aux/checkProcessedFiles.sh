#!/bin/bash

rm datafileslist.txt
rm vidfilelist.txt
find . -maxdepth 1 | sed 's,^[^/]*/,,' | sed s/\.[^\_]*$// | sed s/\.[^\_]*$// | sort | uniq -u > datafileslist.txt
find /media/kostasl/extStore/ExpData/zebrafish_preycapturesetup/AutoSet_02-11-17 -maxdepth 1 | sed 's,^[^/]*/,,' | sed s/\.[^\.]*$// | sed 's/.*\///' | sort | uniq -u > vidfilelist.txt
diff datafileslist.txt vidfilelist.txt

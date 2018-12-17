#!/usr/bin/env python3

import os
import subprocess
import sys

##requires mediainfo
directory = sys.argv[1] # "/media/LinuxDat/expDataKostas/AutoSet420fps_07-12-17LRes/"

# list the files in the directory
files_tosort = list()

for (dirpath, dirnames, filenames) in os.walk(directory):
    files_tosort += [os.path.join(dirpath, file) for file in filenames]


filedata = []
for file in files_tosort:
    # combine filepath and file, take care of the whitespaces 
    
    filepath = file; command = "mediainfo "+"'"+filepath+"'"
#    print(filepath)
    # get the file's data
    data = subprocess.check_output(["/bin/bash", "-c", command]).decode("utf-8")
    # extract the bitrate from the output
    if data.find("Bit rate") > 0:
	    bitrate = [line[line.find(":")+2:].replace(" ", "") \
		       for line in data.splitlines() if "Bit rate" in line][0]

	    bitdepth = [line[line.find(":")+2:].replace(" ", "") \
		       for line in data.splitlines() if "Bits/(Pixel*Frame)" in line][0]

	    # add the found bitrate+filename to he list
	    filedata.append((float(bitdepth),bitrate, file))
    else:
	    print(data)

# sort the list by the bitrate
filedata.sort(key=lambda item: item[0])
# print out
print("bitdepth \t  bitrate \t filename")
for item in filedata:
    print(str(item[0])+"\t" + item[1] + "\t"+item[2])
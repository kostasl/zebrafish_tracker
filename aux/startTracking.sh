#!/bin/bash
for i in `seq 1 1 12`;
do num=$((num + 1));
echo "start $i";
echo "./zebraprey_track --invideolist=startVidsToProcess.txtSplit$i --outputdir=/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackedInitVids/ --logtofile=/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackedInitVids/logs/log$i.txt "
`./zebraprey_track --invideolist=startVidsToProcess.txtSplit$i --outputdir=/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackedInitVids/ --logtofile=/mnt/4E9CF34B9CF32BD9/kostasl/Dropbox/Calculations/zebrafishtrackerData/TrackedInitVids/logs/log$i.txt` &
done

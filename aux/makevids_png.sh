echo ".BATCH Convert to MP4 from Pgm. 2018 kostasl"
echo "Make HD Vids from Image Sequence Dirs.."
echo "example : ./makeVids.sh InputDir fps OutputDir"
fps=$2
outdir=$3
for dir in `find $1 -type d`
do
  #test -d "$dir" || continue
	files=($dir/*.png)
	if [ ${#files[@]} -gt 10 ]; then
#	 echo $dir;
	 echo "Found PNG Image FILES..in $dir";
#	 Make Video
	  filename=${dir//[\/]/_}
	  filename=${filename//[.]/}
	  echo $filename

	firstFrame=`ls $dir -1v | head -n 1 | sed 's/[^0-9]*//g'`
	echo "First Frame $firstFrame"
	   ffmpeg -framerate   $fps -start_number $firstFrame -i $dir/%05d.png -c:v libx264 -crf 18  -crf_max 35 $outdir/$filename.mp4
#	   ffmpeg -framerate   $fps  -pattern_type glob -i '*.png' -c:v libx264 -crf 18 -start_number $firstFrame -crf_max 35 $outdir/$filename.mp4
#	   ffmpeg -framerate $fps $(printf -- "-i %s ") '*.png' -c:v libx264 -crf 18 -crf_max 35 $outdir/$filename.mp4

	fi


# Do something with $dir...
done


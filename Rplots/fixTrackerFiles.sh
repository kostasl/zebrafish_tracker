

grep -rl   '		 templateScore' *.csv | xargs sed -i 's/		 templateScore/	 templateScore/g'
##Remove  Rotifer COunt only llines - starting with tab
sed -i '/^\t/d' AutoSet420fps_*.csv

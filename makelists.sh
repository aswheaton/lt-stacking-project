for year in "2012" "2013" "2014" "2015" "2016"
do
for band in "R" "G" "U"
do
ls | grep ^year | grep band. &>> fits_lists.txt
done
done

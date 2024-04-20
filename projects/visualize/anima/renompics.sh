#!/bin/bash/

# ejemplo:
# sh renompics.sh primavera

cd $1/
mkdir num
jj=0
for ifile in *.png;
do
	jj=$(($jj+1))
	cp $ifile num/$jj.png
	rm $ifile
	mv num/$jj.png .
done
cd ..

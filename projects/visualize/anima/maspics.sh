#!/bin/bash/

# ejemplo:
# sh maspics.sh 29 monstruos

for ii in $(seq 1 $1);
do
	jj=$(($ii+$1))
	cp $2/$ii.png $2/$jj.png
done

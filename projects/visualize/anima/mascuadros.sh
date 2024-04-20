#!/bin/bash/

# ejemplo:
# sh mascuadros.sh 11 25 corredor

for ii in $(seq $1 $2);
do
	jj=$(($ii+$1))
	cp $3/$ii.png $3/$jja.png
done

#!/bin/bash
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------
printf "\n"
printf "\n           do some ✏️ on some 📁\n\n"
printf "                 ■⎔■⎔■\n"
printf "                 ⎔■⎔■⎔\n"
printf "                 ◁▷◁◁◎\n"
printf "\n"
# ------------------------------------------------------------------------------
# get address to access from 1ˢᵗ & 2ⁿᵈ line in paths.txt
dirread=$(sed '1q;d' paths.txt)
dirsave=$(sed '2q;d' paths.txt)
# ------------------------------------------------------------------------------
# loop through files (or directories) in the chosen directory address,
# and write the address of those files (or directories) in
# the 3ʳᵈ line of paths.txt
# each line written in the 3ʳᵈ line of  paths.txt will be pushed to the
# 4ᵗʰ line of paths.txt...
# so at the end there will be a record of what happened.
# ------------------------------------------------------------------------------
for dire in $dirread*; do
  printf "\n"
  printf "       ■⎔■⎔■⎔■⎔■⎔⎔■⎔■⎔⎔■⎔■⎔⎔■⎔■⎔⎔■⎔■⎔⎔■⎔■⎔⎔■⎔■⎔\n"
  printf "       ⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■⎔■\n"
  printf "       ◁▷◁◁◎◁▷◁◁◎◁▷◁◁◎◁▷◁◁◎◁▷◁◁◎◁▷◁◁◎◁▷◁◁◎◁▷◁◁◎\n"
  printf "\n"
  # this skips txt ending dire in dirread
  if [ "${dire: -3}" = "txt" ]; then
    echo "skip the" ${dire: -3}
    continue
  fi
  # write dire address in the third line of paths.txt
  echo $(basename $dire) | sed -i '2r /dev/stdin' paths.txt
  # read & print the third line of paths.txt
  dirread_=$(sed '3q;d' paths.txt)
  echo "  👉 gonna do some work on ⟶  " $dirread_
  # here you can call other stuff to work on the address stored in
  # the third line of paths.txt
  # 👷👷👷
done
# ------------------------------------------------------------------------------

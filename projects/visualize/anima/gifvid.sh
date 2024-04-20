#!/bin/bash/

gifski -o paloma.gif --fps 6 *.png

ffmpeg -i paloma.gif -b:v 1M paloma.mp4

# ffmpeg no sirve para instagram :(
https://ezgif.com/gif-to-mp4

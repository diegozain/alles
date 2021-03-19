function rgb=rainbow2_cb(n)
% diego domenzain
% spring 2021 @ CSM
% ------------------------------------------------------------------------------
% to get even better colormaps, concatenate them!
% the factor n will move the white towards the bottom if n>1.
% cmap=rainbow2(8);colormap(cmap)
rgb=rainbow(0,[0.9412,0.8941,0.2588; 0.0471,0.4824,0.8627; 0.865,0.865,0.865],64);
rgb_=rainbow(0,[0.865,0.865,0.865; 0.8627,0.1961,0.1255; 0,0,0],fix(64*n));
rgb=[rgb;rgb_]; 
end
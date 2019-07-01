%function that implement data augmentation technique to randomize the
%boundary points of the bone mask to generate more medium for kwave
function [mask,ind] = rnd_bone(image) %the initial index for first image

%use data augmentation:
%randomize the boundary points
r_max =  max(ind, [], 'all');
pts = (ind ~=0);
rl = normrnd(0,10);
ind(pts) = ind(pts) + rl;
%for pts with larger value, the variance can be larger
ph = (ind > round(r_max/2));
rh = normrnd(10,30);
ind(ph) = ind(ph) + rh;
%truncate the outlier points:
ind(ind <0) = 0;

plot(ind)
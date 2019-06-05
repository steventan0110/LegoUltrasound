load('mask.mat','M')
[image, ind] = resize(M, 1000, 1000);
rnd_bone(image,10, ind)
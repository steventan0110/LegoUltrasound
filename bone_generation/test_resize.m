load('mask1.mat', 'M')

image = resize(M, 1000, 1000); %use 1000 pixels as Lx and Ly
imshow(image);
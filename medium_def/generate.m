%the main function to generate n bones given the data I cropped
%the only function called by the kwave simulation file when defining the
%meidum

function [image, inds] = generate(Lx, Ly) %generate num of bones

load('mask.mat', 'M')
%resize original image
image = resize(M, Lx, Ly);
%shift the image data by an random amount:
amount = randi([0,Ly],1);
image = circshift(image, [0,amount]);
%50% chance flip the image horizontally 
if (rand > 0.5)
    image = fliplr(image);
end
%randomize the boundaries
[image, inds] = rnd_bone(image);
%dimension of images/inds would add one more dimension of size num


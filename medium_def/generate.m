%the main function to generate n bones given the data I cropped
%the only function called by the kwave simulation file when defining the
%meidum

function [image, inds] = generate(Lx, Ly) %generate num of bones

load('mask.mat', 'M')
[r,c] = size(M);
%resize original image (randomize part incorporated in it)
[image] = resize(M, Lx, Ly);
imshow(M);
%shift the image data by an random amount:
amount = randi([0,round((Ly-c)/2)],1); 
image = circshift(image, [0,amount]);
%50% chance flip the image horizontally 
if (rand > 0.5)
    image = fliplr(image);
end

%retrive the index for label
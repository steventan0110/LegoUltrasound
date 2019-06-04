%the main function to generate n bones given the data I cropped
%the only function called by the kwave simulation file
function [images, inds] = generate(num, Lx, Ly) %generate num of bones

load('mask.mat', 'M')
%resize original image
[image, ind] = resize(M, Lx, Ly); %the index will be used as label for supervised learning
% imshow(image) for test purpose
% hold on
% plot(ind)



%randomize the boundaries
[images, inds] = rnd_bone(image, num);
%dimension of images/inds would add one more dimension of size num


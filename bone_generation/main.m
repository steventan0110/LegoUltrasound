img = imread('C:\Users\Steve\Documents\LegoUltrasound\bone_generation\img1.png');
I = rgb2gray(img);
%crop the useful region out:
I = I(170:630,250:650);

Bmask = img_convert(I);

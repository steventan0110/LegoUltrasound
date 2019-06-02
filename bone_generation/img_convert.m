% this script convert an image into medium usable in kwave simulation
% process
function Bmask = img_convert(img)


%imadjust to enhance the contrast

% pout_imadjust = imadjust(img);
% pout_histeq = histeq(img);
pout_adapthisteq = adapthisteq(img);
% montage({img,pout_imadjust,pout_histeq,pout_adapthisteq},'Size',[1 4])
% title("Original Image and Enhanced Images using imadjust, histeq, and adapthisteq")

%from comparison, the adapthisteq technique is the most useful one
%Now find the binary mask
image = pout_adapthisteq;
imshow(image) 
h1 = drawline('SelectedColor','yellow');
h1.Selected = true;
h2 = drawline('SelectedColor','yellow');
h2.Selected = true;
h3 = drawline('SelectedColor','yellow');
h3.Selected = true;
h4 = drawline('SelectedColor','yellow');
h4.Selected = true;
h5 = drawline('SelectedColor','yellow');
h5.Selected = true;
h6 = drawline('SelectedColor','yellow');
h6.Selected = true;
h7 = drawline('SelectedColor','yellow');
h7.Selected = true;
h8 = drawline('SelectedColor','yellow');
h8.Selected = true;
h9 = drawline('SelectedColor','yellow');
h9.Selected = true;
h10 = drawline('SelectedColor','yellow');
h10.Selected = true;
h11 = drawline('SelectedColor','yellow');
h11.Selected = true;
M = createMask(h1) | createMask(h2)|createMask(h3)|createMask(h4)|createMask(h5)...
    |createMask(h6)|createMask(h7)|createMask(h8)|createMask(h9)|createMask(h10)...
    |createMask(h11);

save('mask1.mat', 'M')

Bmask = M;

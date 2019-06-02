%% Converting image into mask M
img = imread('/Users/weitingtan/Documents/LegoUltrasound/bone_generation/img1.png');
I = rgb2gray(img);
% crop the useful region out:
I = I(170:630,250:650);

Bmask = img_convert(I);


%% Fill in the region under the mask
load('mask1.mat', 'M')
[r,c] = size(M); 
%find the upper 1:
[~, ind] = max(M, [],1);
[~, min_ind] = min(M, [],1);
for i = 1:c
    if ind(i) == min_ind(i)
        ind(i) = r;
    end
end
ind = r - ind;
plot(1:c, ind)


%extrapolate the max point and from an area
img = M;
%find the first value 1


% for i = 1:r
%     for j = 1:c
%         if j < ind(j)
%             img(i,j) =1;
%         else
%             img(i,j) =0;
%         end
%     end 
% end
% imshow(img)
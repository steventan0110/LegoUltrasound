%Resize the region into correct medium dimension 
function [image, new_ind] = resize(img, Lx, Ly)
% Lx and Ly are given in kwave main function to contain the input image 
image = ones(Lx,Ly);
[r,c] = size(img);
%position the bone in the middle of the image
%check dimension:
if Lx < r || Ly < c
    disp('dimension too small')
    return
end 

%find start and position in lateral axis:
mid = Lx/2;
start_x = round(mid - c/2);
end_x = round(mid + c/2);

%fill in the new image:
ind = get_ind(img); 


ind = r - ind;
for i = 1:r  
    for j = start_x:end_x
        index = round(j- start_x) + 1;
        if index > c %prevent oversize due to round
            index = c;
        end 
        if i < ind(index)
            image(i,j) = 0;
        end
    end
end
image = flipud(image);

%adjust index into correct dimensions

new_ind = zeros(1,Lx);
for i = start_x: end_x
    index = round(i- start_x) + 1;
    if index > c %prevent oversize due to round
        index = c;
    end  
    new_ind(i) = ind(index);
end



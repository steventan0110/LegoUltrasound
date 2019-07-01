function ind = get_ind(M)
[r,c] = size(M); 
%find the upper 1:
[~, ind] = max(M, [],1);
[~, min_ind] = min(M, [],1);
for i = 1:c
    if ind(i) == min_ind(i)
        ind(i) = r;
    end
end
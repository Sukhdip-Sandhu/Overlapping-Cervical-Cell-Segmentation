function [bool] = connected_in_clump(dims, clump, nroi)
%CONNECTED_IN_CLUMP Summary of this function goes here
%   Detailed explanation goes here
img1 = zeros(dims);
img2 = zeros(dims);

img1(clump) = 1;
img2(nroi) = 1;

compare = img1 & img2;
compare = any(any(compare));

bool = compare;
end


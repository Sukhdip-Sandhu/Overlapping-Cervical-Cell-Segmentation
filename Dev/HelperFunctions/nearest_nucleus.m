function [nearest_nucleus] = nearest_nucleus(idx, superpixel, nucleus, dims)

tmp = zeros(dims);
tmp(idx{superpixel})= 255;
rp = regionprops(tmp, 'Centroid');
superpixel_centroid = rp(255).Centroid;

distance_list = [];
for i = 1 : size(nucleus, 1)
    nucleus_centroid = nucleus(i).Centroid;
    distance = euclidian_distance(superpixel_centroid, nucleus_centroid);
    distance_list = [distance_list, distance];
end

dl = distance_list;
dl = sort(dl);
nn = find(distance_list == dl(1));
nearest_nucleus = nn;
end


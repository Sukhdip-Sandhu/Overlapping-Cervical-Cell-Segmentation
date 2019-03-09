function [dsc] = display_superpixel_clusters_from_list(image, idx, list)
display_image = zeros(size(image));
cat = [];
for i = 1 : size(list,2)
    val = idx{list(i)};
    cat = vertcat(cat, val);
end
display_image(cat) = 255;
% figure, imshow(display_image)
dsc = display_image;
end
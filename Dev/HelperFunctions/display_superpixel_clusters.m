function [dsc] = display_superpixel_clusters(image, label)
display_image = zeros(size(image));
display_image(label) = 255;
% figure, imshow(display_image)
dsc = display_image;
end
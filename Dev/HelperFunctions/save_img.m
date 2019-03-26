function save_img(im, fn)
save_path = '.\SavedImages\';
save_path = strcat(save_path, fn, '.tif');
imwrite(im, save_path);
end


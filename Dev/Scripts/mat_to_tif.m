%{
Simply script that MATLAB mat data as individual tif files
%}

path_to_data = '../Dataset/Training/Train45Test90/isbi_train.mat';
path_to_write = '../Dataset/Training/Train45Test90/IndividualImages/';

mat = load(path_to_data);
mat_data = mat.ISBI_Train;

% for each 'cell', save the image as a tif
for i=1 : size(mat_data, 1)
    C = mat_data{i};
    save_image = strcat(path_to_write, string(i), '.tif');
    imwrite(C, save_image);
end



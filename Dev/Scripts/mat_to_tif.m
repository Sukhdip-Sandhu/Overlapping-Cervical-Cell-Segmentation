%{
Simply script that MATLAB mat data as individual tif files
%}

% path_to_data = '../../../Dataset/Training/Train45Test90/isbi_train.mat';
path_to_data = '../../../Dataset/Synthetic/trainset_GT.mat';
path_to_write = 'C:\Users\Sukhdip\Desktop\435Project\Dataset\Synthetic\trainset_GT\train_cytoplasm\';

mat = load(path_to_data);
mat_data = mat.train_Cytoplasm;

% for each 'cell', save the image as a tif
for i=1 : size(mat_data, 1)
    C = mat_data{i};
    for j = 1 : size(C,1)
        CC = C{j};
        save_image = strcat(path_to_write, string(i), string(j), '.tif');
        imwrite(CC, save_image);

    end
end



% =========================================================================
%% Housekeeping and Initialization
% =========================================================================
tic; clear; close all hidden; clc;
addpath(genpath('Dataset'));
% path to .mat containing images
path_to_mat_file = 'trainset.mat';
% getting images
images = load(path_to_mat_file);
images = images.trainset;
% number of images
num_images = numel(images);


% =========================================================================
%% MAIN 
% =========================================================================
for image_iter = 1 : num_images
    % CURRENT IMAGE
    im = images{image_iter};
    [curr_image, h, w] = read_image(im);
    [preprocessed_image] = preprocess_image(curr_image);
    [mean_superpixels, label_matrix, idx, region_boundaries] = oversegment_image(preprocessed_image, h);
    
    % image histogram
    img_hist = imhist(mean_superpixels);
    
    % calculate triangle value (Zack GW, Rogers WE, Latt SA (1977))
    actual_threshold_value = triangle_th_actual(img_hist, 256);
    
    % calculate triangle value (David (2018))
    david_threshold_value = triangle_th_david(img_hist, 256);
    
    % calculate difference between the two
    diff = actual_threshold_value - david_threshold_value;
    
%     % UNCOMMENT THIS IF YOU THINK ITS READY!
%     if diff ~= 0, fprintf('*** ERROR IN ITERATION %f ***\n', image_iter); return; end
    
    fprintf('Actual: %f\t David: %f\t Difference: %f\n', actual_threshold_value, david_threshold_value, diff);
    
end


% =========================================================================
%% Functions
% =========================================================================
%% Return Image from Path
function [image, height, width] = read_image(image)
    % image border relected
    image(1,:) = image(2,:);
    image(:,1) = image(:,2);
    image(end,:) = image(end-1,:);
    image(:,end) = image(:,end-1);
    [height, width] = size(image);
end

%% Returns the Preprocessed Image
function [preprocessed_image] = preprocess_image(image)
    % apply 2D median filter across image
    preprocessed_image = medfilt2(image);
end

%% Superpixel Function
function [mean_superpixels, label_matrix, idx, region_boundaries] = oversegment_image(image, height)
    % generate superpixels
    % label_matrix = labels of each superpixel region
    [label_matrix, num_labels] = superpixels(image, height*2);
    
    % superpixel regions
    region_boundaries = boundarymask(label_matrix);

    % show the superpixel overlay
    % figure, imshow(imoverlay(image, region_boundaries, 'K'),'InitialMagnification', 100)

    % mean of superpixels initialization
    mean_superpixels = zeros(size(image),'like',image);

    % linear indices
    idx = label2idx(label_matrix);

    % calculate means for superpixel
    for label = 1:num_labels
        val = idx{label};
        mean_superpixels(val) = mean(image(val));
    end
end
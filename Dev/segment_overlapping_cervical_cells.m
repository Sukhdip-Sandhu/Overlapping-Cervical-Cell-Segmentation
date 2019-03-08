%{
    Main Script
     Dev Branch
%}
%% Housekeeping and Setup
tic
clear all %#ok<CLALL>
close all
clc
addpath('EvalCode', 'HelperFunctions', 'Scripts', 'Dataset/Training/Train45Test90/IndividualImages/')


%% Hardcoded Example Variables
training_image_path = 'synth.tif';

%% Image Preprocessing
image = imread(training_image_path);
figure, imshow(image);
% apply 2D median filter across image
preprocessed_image = medfilt2(image);
% preprocessed_image = imgaussfilt(preprocessed_image);


%% Cell Mass Detection
% generate superpixels
[label_matrix, num_labels] = superpixels(preprocessed_image, size(preprocessed_image,1)*2);
region_boundaries = boundarymask(label_matrix);

% show the superpixel overlay
% figure, imshow(imoverlay(preprocessed_image, region_boundaries, 'K'),'InitialMagnification', 100)

% mean of superpixels initialization
mean_superpixels = zeros(size(preprocessed_image),'like',preprocessed_image);

% linear indices
idx = label2idx(label_matrix);

% calculate means for superpixel
for label = 1:num_labels
    val = idx{label};
    mean_superpixels(val) = mean(preprocessed_image(val));
end

figure, imshow(imoverlay(mean_superpixels, region_boundaries, 'K'),'InitialMagnification', 100)

% image histogram
img_hist = imhist(mean_superpixels);

% apply triangle thresholding (Zack GW, Rogers WE, Latt SA (1977))
triangle_theshold_image = triangle_th(img_hist, 256);

% ROI
binarized_image = ~imbinarize(mean_superpixels, triangle_theshold_image);
binarized_image = double(binarized_image);
figure, imshow(binarized_image);
cell_mass = label_matrix .* binarized_image;

roi = unique(cell_mass);
roi = roi(2:end); % remove unique '0'
num_rois = size(roi,1);

%% 2 Nucleus Extraction
% distance metric for determining local window of each superpixel
d = 0.1 * size(preprocessed_image,1);
% d = 7;
% Key: superpixel label. Value: X,Y coordinate of superpixel centroid
superpixel_centroid_dictionary = containers.Map('KeyType','double', 'ValueType', 'any');
hold on
% For each label, find its centroid (X,Y) and store into dictionary
for label = 1:num_rois
    tmp = zeros(size(preprocessed_image));
    val = idx{roi(label)};
    tmp(val) = 255;
    props = regionprops(tmp);
    superpixel_centroid_dictionary(label) = props(255).Centroid;
%     plot(props(255).Centroid(1), props(255).Centroid(2), 'X', 'MarkerEdgeColor', 'r') % used for visualization...
end

% Key: superpixel label. Value: list of local window candidates
window_dictionary = containers.Map('KeyType','double','ValueType', 'any');

% compare every superpixel with every other superpixel, and if centroids
% are within distance d, add it as part of local window
for label = 1 : num_rois
    count = 1;
    superpixel_local_window = [];
    for j = 1 : num_rois
%         if label == j
%             continue
%         end
        if (euclidian_distance(superpixel_centroid_dictionary(label), superpixel_centroid_dictionary(j)) < d)
            superpixel_local_window(count) = j;
            count = count + 1;
        end
    end
    window_dictionary(label) = superpixel_local_window;
end

% Key: superpixel label. Value: local window
local_windows_dictionary = containers.Map('KeyType','double','ValueType', 'any');

% for each superpixel label, calculate its local window
for label = 1 : num_rois
    local_window_concat = [];
    n = window_dictionary(label);
    for j = 1 : size(n, 2)
        val = idx{roi(n(j))};
        local_window_concat = vertcat(local_window_concat, val);
    end
    local_windows_dictionary(label) = local_window_concat; 
end

% paper paramater
p = 2;
q = 10;
k = 0.25;
s = 255;

count = 0;
D = zeros(size(mean_superpixels));

% find nuclei candidates
nuclei_candidates = [];
nuclei_candidates_idx = [];
for label = 1 : num_rois
    mew = mean(mean_superpixels(local_windows_dictionary(label)));
    sigma = std(double(mean_superpixels(local_windows_dictionary(label))));
    T = mew * (1 + p*exp(-q*mew) + k *((sigma/s) - 1));
    
    if mean(mean_superpixels(idx{roi(label)})) < T
        nuclei_candidates = vertcat(nuclei_candidates, idx{roi(label)});
        nuclei_candidates_idx = [nuclei_candidates_idx, roi(label)];
        count = count + 1;
    end 
end

%% Nucleus Fine Tuning
% 1 break superpixels - Not complete
nc = display_superpixel_clusters(image, nuclei_candidates);
ncl = logical(nc);
cc = bwconncomp(ncl);
rp = regionprops(ncl);

% 2 reject tiny superpixels - Not complete
for i = 1 : numel(rp)
    candidate_nuclei_area = rp(i).Area;
    if candidate_nuclei_area < 30
        % disregard it
    end
end

% 3 reject low circularity - Not complete
perm = regionprops(cc, 'perimeter');
for i = 1 : numel(rp)
    A = rp(i).Area;
    P = perm(i);
    circularity = (4*pi*A)/(P.Perimeter^2);
    if circularity < 0.5
        % disregard it
    end
end

%% Cytoplasm Segmentation



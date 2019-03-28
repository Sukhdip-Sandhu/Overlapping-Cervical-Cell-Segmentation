%{
    Main Script
	
%}
%% Housekeeping and Setup
tic
clear all %#ok<CLALL>
close all
clc
addpath('EvalCode', 'HelperFunctions', 'Scripts', 'Dataset/Training/Train45Test90/IndividualImages/')


%% Hardcoded Variables
training_image_path = '9.tif';

%% Image Preprocessing
image = imread(training_image_path);
[w,h] = size(image);
image = image(2:w-1, 2:h-1);
% figure, imshow(image);
% save_img(image, '1.original_image')
% apply 2D median filter across image
preprocessed_image = medfilt2(image);
% preprocessed_image = imgaussfilt(preprocessed_image);
% preprocessed_image = imsharpen(preprocessed_image);


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

% figure, imshow(imoverlay(mean_superpixels, region_boundaries, 'K'),'InitialMagnification', 100)
% save_img(imoverlay(mean_superpixels, region_boundaries, 'K'), '3.region_boundaries')
% save_img(mean_superpixels, '4.mean_superpixels')
% image histogram
img_hist = imhist(mean_superpixels);

% apply triangle thresholding (Zack GW, Rogers WE, Latt SA (1977))
triangle_theshold_image = triangle_th(img_hist, 256);

% ROI
binarized_image = ~imbinarize(mean_superpixels, triangle_theshold_image);
binarized_image = double(binarized_image);
% figure, imshow(binarized_image);
binarized_image = bwareaopen(binarized_image,175); 
% save_img(imoverlay(mean_superpixels, edge(binarized_image), 'r'), '5.cell_mass');
% figure, imshow(binarized_image);
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
% hold on
% For each label, find its centroid (X,Y) and store into dictionary
for label = 1:num_rois
    tmp = zeros(size(preprocessed_image));
    val = idx{roi(label)};
    tmp(val) = 255;
    rp = regionprops(tmp);
    superpixel_centroid_dictionary(label) = rp(255).Centroid;
%     plot(rp(255).Centroid(1), rp(255).Centroid(2), 'X', 'MarkerEdgeColor', 'r') % used for visualization...
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
% 1 Breaking Superpixel Clusters
nc = display_superpixel_clusters(image, nuclei_candidates);
ncl = logical(nc);
cc = bwconncomp(ncl);

num_comp = 2;

nuclei_candidates_idx_new = [];
for i = 1 : cc.NumObjects
    connected_components = [];
    clump = cc.PixelIdxList{i};
    for j = 1 : numel(nuclei_candidates_idx)
        nuclei_roi = idx(nuclei_candidates_idx(j));
        nuclei_roi = nuclei_roi{1};
        if connected_in_clump(size(image) ,clump, nuclei_roi)
            connected_components = [connected_components, nuclei_candidates_idx(j)];
        end
    end
    if size(connected_components, 2) > num_comp
       nuclei_candidates_idx_new = [nuclei_candidates_idx_new, remove_connected_clump(connected_components, idx, mean_superpixels, num_comp)];
    else
       nuclei_candidates_idx_new = [nuclei_candidates_idx_new, connected_components];
    end    
end

%% 2 Reject Tiny Superpixels
% (if this fails, look at step 3 indices to remove alg)
nc = display_superpixel_clusters_from_list(image, idx, nuclei_candidates_idx_new);
ncl = logical(nc);
cc = bwconncomp(ncl);
rp = regionprops(ncl);

for i = 1 : numel(nuclei_candidates_idx_new)
    tmp = zeros(size(image));
    tmp(idx{nuclei_candidates_idx_new(i)}) = 255;
    rp = regionprops(tmp);
    if rp(255).Area < 30
        nuclei_candidates_idx_new(i) = [];
    end
end

%% 3 Reject Low Circularity
nc = display_superpixel_clusters_from_list(image, idx, nuclei_candidates_idx_new);
ncl = logical(nc);
cc = bwconncomp(ncl);
rp = regionprops(ncl);

A = regionprops(cc, 'area');
P = regionprops(cc, 'perimeter');
indices_to_remove = [];
for i = 1 : numel(A)
    circularity = (4*pi*A(i).Area)/(P(i).Perimeter^2);
    if circularity < 0.75
        for j = 1 : size(nuclei_candidates_idx_new, 2)
            nuclei_roi = idx(nuclei_candidates_idx_new(j));
            nuclei_roi = nuclei_roi{1};
            if connected_in_clump(size(image) , cc.PixelIdxList{i}, nuclei_roi)
                indices_to_remove = [indices_to_remove, j];
            end
        end
    end
end

for i = size(indices_to_remove,2):-1:1
    nuclei_candidates_idx_new(indices_to_remove(i)) = [];
end

nuclei = nuclei_candidates_idx_new;
nc = display_superpixel_clusters_from_list(image, idx, nuclei);
% figure, imshow(nc);
% save_img(nc, '6.nucleus_extraction')
%% Nucleus_Superpixels (Map To) Region
ncl = logical(nc);
cc = bwconncomp(ncl);
rp = regionprops(ncl);
nucleus_superpixels_to_region_dictionary = containers.Map('KeyType','double', 'ValueType', 'any');

for i = 1 : numel(nuclei)
    nuclei_roi = idx(nuclei(i));
    nuclei_roi = nuclei_roi{1};
    for j = 1 : cc.NumObjects
        if connected_in_clump(size(image) , cc.PixelIdxList{j}, nuclei_roi)
            nucleus_superpixels_to_region_dictionary(j) = nuclei(i);
        end
    end
end


%% Cytoplasm Segmentation
%% Superpixel Partitioning
nearest_nucleus_dictionary = containers.Map('KeyType','double', 'ValueType', 'any');
for i = 1 : num_rois
    superpixel = roi(i);
    nn = nearest_nucleus(idx, superpixel, rp, size(image));
    try
        nearest_nucleus_dictionary(nn) = [nearest_nucleus_dictionary(nn), superpixel];
    catch
        nearest_nucleus_dictionary(nn) = superpixel;
    end
end

for i = 1 : cc.NumObjects
    nc = display_superpixel_clusters_from_list(image, idx, nearest_nucleus_dictionary(i));
%     figure, imshow(nc);
    bw = edge(nc);
    edge_color = edge_LUT(i);
    if i == 1
        segmentation_overlay = imoverlay(image, bw, edge_color);
    else
        segmentation_overlay = imoverlay(segmentation_overlay, bw, edge_color);
    end
end

figure, imshow(segmentation_overlay)
% save_img(segmentation_overlay, '7.segmentation_overlay')

%% Cellwise Contour Refinement
%% END
toc


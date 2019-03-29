%% Housekeeping and Initialization
tic; clear; close all; clc;
addpath('EvaluationCode', 'HelperFunctions', genpath('Dataset'), 'Results')

%% MAIN
tic
path_to_mat_file = 'testset.mat';
images = load(path_to_mat_file);
images = images.testset;
num_images = numel(images);
master_cell = cell([num_images,1]);

waitbar_handle = waitbar(0, 'Segmenting Cytoplasms... Please Wait :)');

for image_iter = 1 : num_images
    %% GET CURRENT IMAGE
    [curr_image, h, w] = read_image(images{image_iter});
    
    %% CELL MASS DETECTION
    [preprocessed_image] = preprocess_image(curr_image);
    [mean_superpixels, label_matrix, idx, region_boundaries] = oversegment_image(preprocessed_image, h);
    [cell_mass, roi, num_roi] = triangle_threshold(mean_superpixels, label_matrix);
    [superpixel_area_LUT] = label_to_area_dictionary(cell_mass, roi, num_roi);

    %% NUCLEUS EXTRACTION
    [candidate_labels, nuclei_candidates] = initial_nuclei_candidates(superpixel_area_LUT, roi, num_roi, idx, w, mean_superpixels);
    [broken_cluster_candidates, nuclei_candidates] = breaking_superpixel_clusters(superpixel_area_LUT, candidate_labels, nuclei_candidates, mean_superpixels, idx);
    [reject_tiny_candidates, new_nuclei_candidates] = reject_tiny_superpixels(superpixel_area_LUT, broken_cluster_candidates);
    [nuclei_labels, nuclei] = reject_low_circularity(superpixel_area_LUT, reject_tiny_candidates, new_nuclei_candidates, idx);

    %% CYTOPLASM SEGMENTATION
    [initial_segmentation] = superpixel_partitioning(superpixel_area_LUT, nuclei, roi, num_roi);
    [cellwise_contour_segmentation] = cellwise_contour_refinement(superpixel_area_LUT, initial_segmentation, curr_image, label_matrix);

    [initial_boundaries] = draw_boundary_initial(superpixel_area_LUT, initial_segmentation, curr_image);
    [cellwise_contour_boundaries] = draw_boundary_refinement(cellwise_contour_segmentation, curr_image);

    %% SAVE
    [master_cell] = save_to_mat(cellwise_contour_segmentation, master_cell, image_iter);
    waitbar(image_iter/num_images);
end
% safely close the waitbar
close(waitbar_handle);

% save segmentation results
save(strcat('.\Results\Segmentations_', path_to_mat_file(1:end-4)), 'master_cell');

toc


%% Functions
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

%% Triangle Threshold Function
function [cell_mass, roi, num_roi] = triangle_threshold(mean_superpixels, label_matrix)
    % image histogram
    img_hist = imhist(mean_superpixels);
    
    % calculate triangle value (Zack GW, Rogers WE, Latt SA (1977))
    triangle_theshold_value = triangle_th(img_hist, 256);
    
    % binarize image, convert to double image
    binary_image = double(~imbinarize(mean_superpixels, triangle_theshold_value));
    
    % remove all components with area < 250 pixels
    binary_image = bwareaopen(binary_image, 250); 
    
    % keep the labels are are part of the binary cell image
    cell_mass = label_matrix .* binary_image;

    % array of the unique labels
    roi = nonzeros(unique(cell_mass));
    % number of unique labels
    num_roi = size(roi,1); 
end

%% Look Up Table
function [dict] = label_to_area_dictionary(cell_mass, roi, num_roi)
    % dictionary {KEY = ROI}, {VALUE = Area}
    dict = containers.Map('KeyType','double','ValueType', 'any');
    
    % for each roi (label)
    for label = 1 : num_roi
        % area = region corresponding to label
        area = double(cell_mass == roi(label))*255;
        % add area value to label key
        dict(roi(label)) = area;
    end
end

%% 
function [candidate_labels, nuclei_candidates] = initial_nuclei_candidates(superpixel_area_LUT, roi, num_roi, idx, width, mean_superpixels)
    [p, q, k, s, d] = paper_parameters(width);
    superpixel_centroids = determine_superpixel_centroids(superpixel_area_LUT, roi, num_roi);
    local_windows_dictionary = determine_local_window(superpixel_centroids, roi, num_roi, idx, d);
    candidate_labels = find_nuclei_candidates(p, q, k, s, roi, num_roi, local_windows_dictionary, idx, mean_superpixels);
    nuclei_candidates = image_from_labels(superpixel_area_LUT, candidate_labels);
end

%% Return Parameters Specified in Research Paper
function [p, q, k, s, d] = paper_parameters(width)
    % weighting variables given in paper
    p = 2;
    q = 10;
    k = 0.25;
    s = 255;
    % distance between centroids as part of comparison window
    d = 0.1 * width;
end

%% Construct Dictionary of Label to Centroid
function [centroid_dict] = determine_superpixel_centroids(superpixel_area_LUT, roi, num_roi)
    % dictionary [key = label] [value = centroid]
    centroid_dict = containers.Map('KeyType','double', 'ValueType', 'any');
    
    for label = 1:num_roi
        % hold on
        rp = regionprops(superpixel_area_LUT(roi(label)));
        centroid_dict(label) = rp(255).Centroid;
        % plot(rp(255).Centroid(1), rp(255).Centroid(2), 'X', 'MarkerEdgeColor', 'r') % used for visualization...
    end
end

%% Construct Dictionary 
function [local_windows_dictionary] = determine_local_window(superpixel_centroids, roi, num_roi, idx, d)
    % dictionary [key = index] [value = list of neighbours]
    window_dictionary = containers.Map('KeyType','double','ValueType', 'any');

    % compare every superpixel with every other superpixel, and if centroids
    % are within distance d, add it as part of local window
    for label1 = 1 : num_roi
        count = 1;
        superpixel_local_window = [];
        for label2 = 1 : num_roi
            if (euclidian_distance(superpixel_centroids(label1), superpixel_centroids(label2)) < d)
                superpixel_local_window(count) = roi(label2);
                count = count + 1;
            end
        end
        % [key = index] [value = list of neighbours]
        window_dictionary(label1) = superpixel_local_window;
    end
    
    % Key: superpixel label. Value: local window
    local_windows_dictionary = containers.Map('KeyType','double','ValueType', 'any');

    % for each superpixel label, calculate its local window
    for label = 1 : num_roi
        
        local_window_concat = [];
        neighbours = window_dictionary(label);
        
        % for each
        for i = 1 : size(neighbours, 2)
            val = idx{neighbours(i)};
            local_window_concat = vertcat(local_window_concat, val);
        end
        
        local_windows_dictionary(label) = local_window_concat;
    
    end
end

%% Function Returns Nuclei Candidates
function [candidate_labels] = find_nuclei_candidates(p, q, k, s, roi, num_roi, local_windows_dictionary, idx, mean_superpixels)
    % thresholding algorithm based off paper

    % nuclei candidates
    candidate_labels = [];

    % for each superpixel
    for label = 1 : num_roi
        % calculate mean and sigma of its determined neighbours
        mew = mean(mean_superpixels(local_windows_dictionary(label)));
        sigma = std(double(mean_superpixels(local_windows_dictionary(label))));
        % determine threshold based off superpixels determined neighbours
        T = mew * (1 + p*exp(-q*mew) + k *((sigma/s) - 1));

        % if below threshold (signifyings it a nucleus candidate)
        if mean(mean_superpixels(idx{roi(label)})) < T
            % track label as a candidate
            candidate_labels = [candidate_labels, roi(label)];
        end 
    end
end

%% Function Breaks Superpixel Clusters of Clumps > 2
function [broken_cluster_candidates, new_nuclei_candidates] = breaking_superpixel_clusters(superpixel_area_LUT, candidate_labels, nuclei_candidates, mean_superpixels, idx)
    % specified clump size ...
    CLUMP_SIZE = 2;    
    
    % convert candidates to logical for next step
    nuclei_candidates_logical = logical(nuclei_candidates);
    
    % connected components
    cc = bwconncomp(nuclei_candidates_logical);
    
    % parameter for connected_in_clump function
    region_size = size(nuclei_candidates);
    
    % candidates after breaking clusters
    broken_cluster_candidates = [];
    
    % for each clump
    for i = 1 : cc.NumObjects
        % clumped nuclei list
        clumped_nuclei = [];
        % ares corresponding to clump
        clump = cc.PixelIdxList{i};
        % for each candidate nuclei
        for j = 1 : numel(candidate_labels)
            % area corresponding to candidate nuclei
            nuclei_roi = idx{candidate_labels(j)};
            if connected_in_clump(region_size, clump, nuclei_roi)
                clumped_nuclei = [clumped_nuclei, candidate_labels(j)];
            end
        end
        
        % if the nuclei clumped greater than alloted clumpsize
        if numel(clumped_nuclei) > CLUMP_SIZE
           broken_cluster_candidates = [broken_cluster_candidates, remove_connected_clump(clumped_nuclei, idx, mean_superpixels, CLUMP_SIZE)];
        else
           broken_cluster_candidates = [broken_cluster_candidates, clumped_nuclei];
        end    
    end
    new_nuclei_candidates = image_from_labels(superpixel_area_LUT, broken_cluster_candidates);
end

%% Function Reject Tiny Superpixels
function [reject_tiny_candidates, new_nuclei_candidates] = reject_tiny_superpixels(superpixel_area_LUT, broken_cluster_candidates)
    % specified area size ...
    AREA_SIZE = 30;
    
    % for each nuclei superpixel
    for i = 1 : numel(broken_cluster_candidates)
        superpixel = superpixel_area_LUT(broken_cluster_candidates(i));
        rp = regionprops(superpixel);
        % if area of superpixel less than AREA_SIZE
        if rp(255).Area < AREA_SIZE
            % flag it as -1
            broken_cluster_candidates(i) = 0;
        end
    end
    % return all candidates not flagged with -1 (>0)
    reject_tiny_candidates = broken_cluster_candidates(broken_cluster_candidates > 0);
    new_nuclei_candidates = image_from_labels(superpixel_area_LUT, reject_tiny_candidates);
end

function [nuclei_labels, nuclei] = reject_low_circularity(superpixel_area_LUT, reject_tiny_candidates, new_nuclei_candidates, idx)
    % specified circularity...
     CIRCULARITY_VALUE = 0.75;
    
    % convert candidates to logical for next step
    nuclei_candidates_logical = logical(new_nuclei_candidates);
    
    % connected components
    cc = bwconncomp(nuclei_candidates_logical);
    
    % parameter for connected_in_clump function
    region_size = size(new_nuclei_candidates);
    
    % calculate area and perimeter of connected components
    area = regionprops(cc, 'area');
    perimeter = regionprops(cc, 'perimeter');
    
    % for each candidate nuclie
    for i = 1 : numel(area)
        % calculate its circularity
        circularity = (4*pi*area(i).Area)/(perimeter(i).Perimeter^2);
        % if under CIRCULARITY_VALUE
        if circularity < CIRCULARITY_VALUE
            % remove superpixels corresponding to nuclei
            for j = 1 : numel(reject_tiny_candidates)
                if(reject_tiny_candidates(j) == 0)
                    continue
                end
                nuclei_roi = idx{reject_tiny_candidates(j)};
                % if nuclie part of nuclei clump to remove
                if connected_in_clump(region_size, cc.PixelIdxList{i}, nuclei_roi)
                    % flag index with -1
                    reject_tiny_candidates(j) = 0;
                end
            end
        end
    end

    nuclei_labels = reject_tiny_candidates(reject_tiny_candidates > 0);
    nuclei = image_from_labels(superpixel_area_LUT, nuclei_labels);
    
end

%% Function For Nearest Nucleus
function [nearest_cytoplasm, cc] = superpixel_partitioning(superpixel_area_LUT, nuclei, roi, num_roi)
    
    % conneced componens and region properties 
    nuclei_logical = logical(nuclei);
    cc = bwconncomp(nuclei_logical);
    rp = regionprops(nuclei_logical);
    
    % dictionary [key = nucleus label] [value = list of nearest superpixels]
    nearest_cytoplasm = containers.Map('KeyType','double', 'ValueType', 'any');
    
    for i = 1 : num_roi
        % for each superpixel
        superpixel = roi(i);
        % find its nearest nucleus
        nn = nearest_nucleus(superpixel_area_LUT, superpixel, rp);
        % add list of superpixels nearest to nuclues
        try
            nearest_cytoplasm(nn) = [nearest_cytoplasm(nn), superpixel];
        catch
            nearest_cytoplasm(nn) = superpixel;
        end
    end
end


%% Function Returns the Nearest Nucleus (bad data structure right now...)
function [nn] = nearest_nucleus(superpixel_area_LUT, superpixel, nuclei)
    % find the centroid of the given superpixel
    current_superpixel = superpixel_area_LUT(superpixel);
    rp = regionprops(current_superpixel, 'Centroid');
    superpixel_centroid = rp(255).Centroid;
    
    % initialize a list of distances between superpixel and nuclie
    distance_list = [];
    
    % for each nuclei
    for i = 1 : numel(nuclei)
        % calculate its centroid
        nucleus_centroid = nuclei(i).Centroid;
        % calculate distances between superpixel and nuclei centroid
        distance = euclidian_distance(superpixel_centroid, nucleus_centroid);
        distance_list = [distance_list, distance];
    end
    
    % return the index of the nearest nucleus
    sorted_distance_list = sort(distance_list);
    nn = find(distance_list == sorted_distance_list(1));
end

%% Function Returns the Initial Boundary Image
function [overlayed_boundary] = draw_boundary_initial(superpixel_area_LUT, initial_segmentation, image)
    % for each nucleus
    for i = 1 : initial_segmentation.Count
        % determine superpixel_partitioning
        initial_cell_cytoplasm = image_from_labels(superpixel_area_LUT, initial_segmentation(i));
        % determine boundary of partitioning
        cell_boundary = edge(initial_cell_cytoplasm);
        % assign color to edge
        edge_color = edge_LUT(i);
        if i == 1
            overlayed_boundary = imoverlay(image, cell_boundary, edge_color);
        else
            overlayed_boundary = imoverlay(overlayed_boundary, cell_boundary, edge_color);
        end
    end
end

%% Function Returns the Boundary Image
function [overlayed_boundary] = draw_boundary_refinement(initial_segmentation, image)
    % for each nucleus
    for i = 1 : initial_segmentation.Count
        % determine superpixel_partitioning
        initial_cell_cytoplasm = initial_segmentation(i);
        % determine boundary of partitioning
        cell_boundary = edge(initial_cell_cytoplasm);
        % assign color to edge
        edge_color = edge_LUT(i);
        if i == 1
            overlayed_boundary = imoverlay(image, cell_boundary, edge_color);
        else
            overlayed_boundary = imoverlay(overlayed_boundary, cell_boundary, edge_color);
        end
    end
end


%% Function Performs Cellwise Contour Refinement
function [contour_refined_cytoplasm] = cellwise_contour_refinement(superpixel_area_LUT, initial_segmentation, image, L)

    % dictionary [key = nucleus label] [value = list of nearest superpixels]
    contour_refined_cytoplasm = containers.Map('KeyType','double', 'ValueType', 'any');
    
    for i = 1 : initial_segmentation.Count
        % determine superpixel_partitioning
        initial_cell_cytoplasm = image_from_labels(superpixel_area_LUT, initial_segmentation(i));
        
        % morpological structuring elements
        SE_dilate = strel('square', 25);
        SE_erode = strel('square', 25);
        
        % specify what is forsure foreground and forsure background based
        % off morpological operations
        foremask = imerode(initial_cell_cytoplasm, SE_erode);
        backmask = ~imdilate(initial_cell_cytoplasm, SE_dilate);
        
        % grabcut segmentation :)
        BW = grabcut(image, L, logical(~backmask), logical(foremask), logical(backmask));
        contour_refined_cytoplasm(i) = BW;
     end
end

%% Function Stores Segmention to Mask
function [master_cell] = save_to_mat(segmentation, master_cell, image_iter)
    master_cell{image_iter} = cell([segmentation.Count, 1]);
    for i = 1 : segmentation.Count
        master_cell{image_iter}{i} = segmentation(i);
    end
end

%% HELPER FUNCTIONS
%% Output an image based off labels
function [image] = image_from_labels(superpixel_area_LUT, labels)
    image = zeros(size(superpixel_area_LUT(labels(1))));
    for i = 1 : numel(labels)
        image = image + superpixel_area_LUT(labels(i));
    end
end

%% Function that returns Euclidian Distance between two (X,Y) coordinates.
function [dist] = euclidian_distance(p1, p2)
    diff = p2 - p1;
    dist = sqrt(diff * diff');
end

%% Function Returns True if Candidate Connected to Clump
function [bool] = connected_in_clump(size, cluser, nuclei)
    % initial two regions for each corresponding object
    region1 = zeros(size);
    region2 = zeros(size);
    % set area equal to 1
    region1(cluser) = 1;
    region2(nuclei) = 1;
    % if there is any overlap, return true (connected)
    compare = region1 & region2;
    compare = any(any(compare));
    bool = compare;
end

%% Function Iteratively Removes Clusted Superpixels By Highest Intensity
function [updated_labels] = remove_connected_clump(clustered_superpixels, idx, mean_superpixels, clump_size)
    % initialize intensity list 
    intensity_list = zeros(1, numel(clustered_superpixels));

    % calculate the intensities of clustered superpixels
    for i = 1 : numel(clustered_superpixels)
        intensity_list(i) = mean_superpixels(idx{clustered_superpixels(i)}(1));
    end

    % sort list in accending order
    sorted_list = sort(intensity_list);

    % return the first 'clump_size' number of elements
    indices = sorted_list(1:clump_size);

    % create list of labels to be returned
    updated_labels = zeros(1, numel(indices));

    for i = 1 : numel(indices)
        % index of clustered_superpixels corresponding to low intensity
        index = intensity_list == indices(i);

        % superpixel label 
        superpixel_label = clustered_superpixels(index);

        % if regions have the same intensity, return the first one (arbitrary)
        if size(superpixel_label,2) ~= 1
            superpixel_label = superpixel_label(1);
        end

        updated_labels(i) = superpixel_label;
    end
end

%% Function Assigns Color to Cell Boundary
function [color] = edge_LUT(i)
    num_colors = 3;
    r = mod(i,num_colors);
    if r == 0, color = 'r'; return, end
    if r == 1, color = 'g'; return, end
    if r == 2, color = 'b'; return, end
end
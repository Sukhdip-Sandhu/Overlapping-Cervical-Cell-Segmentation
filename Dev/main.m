%{
    Main Dev Script
%}
%% Housekeeping
clear all %#ok<CLALL>
close all
clc
%% Hardcoded Example Variables
training_image_path = './Dataset/Training/Train45Test90/Individual_Images/single.tif';
%% Image Preprocessing
image = imread(training_image_path);
% apply 2D median filter across image
preprocessed_image = medfilt2(image);
%% 1 : CELL MASS DETECTION
% generate superpixels
[L,N] = superpixels(preprocessed_image, size(preprocessed_image,1));
BW = boundarymask(L);

% show the superpixel overlay
imshow(imoverlay(preprocessed_image, BW, 'K'),'InitialMagnification', 100)

% mean of superpixels 
B = zeros(size(preprocessed_image),'like',preprocessed_image);
idx = label2idx(L);

for labelVal = 1:N
    val = idx{labelVal};
    B(val) = mean(preprocessed_image(val));
%     imshow(B)
end

% figure, imshow(B);  
img_hist = imhist(B);
output = triangle_th(img_hist, 256);
bi = ~imbinarize(B, output);
% figure, imshow(bi)
bi_uint8 = uint8(bi);
masked = B .* bi_uint8;
% imshow(masked)

%% 2 : NUCLEUS EXTRACTION
d = 0.1 * size(preprocessed_image,1);

dict = containers.Map('KeyType','double','ValueType', 'any');
hold
for labelVal = 1:N
    C = zeros(size(preprocessed_image),'like',preprocessed_image);
    val = idx{labelVal};
    C(val) = 255;
    props = regionprops(C);
    plot(props(255).Centroid(1), props(255).Centroid(2), 'X', 'MarkerEdgeColor', 'r')
    dict(labelVal) = props(255).Centroid;
end

window_dict = containers.Map('KeyType','double','ValueType', 'any');
for i = 1 : N
    count = 1;
    local_window = [];
    for j = 1 : N
        if (euclidian_distance(dict(i), dict(j)) < d)
%             sprintf('%d - %d', i, j)
            local_window(count) = j;
            count = count + 1;
        end
    end
    window_dict(i) = local_window;
end

local_windows = containers.Map('KeyType','double','ValueType', 'any');
for i = 1 : N
    concat_list = [];
    n = window_dict(i);
    for j = 1 : size(n,2)
        val = idx{n(j)};
        concat_list = vertcat(concat_list, val);
    end
    local_windows(i) = concat_list; 
end

p = 2;
q = 10;
k = 0.25;
s = 1;
count = 1;
D = zeros(size(preprocessed_image),'like',preprocessed_image);
for i = 1 : N
    mew = mean(preprocessed_image(local_windows(i)));
    sigma = std(double(preprocessed_image(local_windows(i))));
    T = mew * ( 1 + (p^(-q*mew)) + k *((sigma/s) - 1));
    if mean(preprocessed_image(idx{labelVal})) < T
        sprintf('Nucleus Candidate! - %d', count)
        count = count + 1;
    end
end


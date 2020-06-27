%%%%%%%%%%%%%%%%
%%%WORMFINDER%%%
%%%%%%%%%%%%%%%%

%This script is designed to detect worms of different stages (stitched brightfield
%images of individual wells in 96 well plates) and return their body area and length. 
function wf = wormF(mypath)
%specify directory with images
path_to_images=mypath;
files = dir(strcat(path_to_images, '/*.png'));

%make output directories if none exist
if ~exist(strcat(path_to_images, '/results'),'dir') mkdir(strcat(path_to_images, '/results')); end
if ~exist(strcat(path_to_images, '/results_txt'),'dir') mkdir(strcat(path_to_images, '/results_txt')); end
dir2 = strcat(path_to_images, '/results');
dir3 = strcat(path_to_images, '/results_txt');


for file = files'   
message = sprintf('Analyzing image %s ...', file.name)
full_path = fullfile(path_to_images, file.name);
image_w=imread(full_path);
image_gray=image_w;

%mask the image based on V hsv component
minValue=0;
maxValue=0.35;
v = mat2gray(image_w);
valueMask = v > minValue & v < maxValue;

%close holes in the mask
se = strel('disk',50);
valueMask = imdilate(valueMask, se);
se = strel('disk',50);
valueMask = imdilate(valueMask, se);

%pick the largest object in the mask
valueMask = imerode(valueMask, se);
valueMask = imerode(valueMask, se);
valueMask = imerode(valueMask, se);
valueMask=imcomplement(valueMask);
valueMask = bwareafilt(valueMask, 1, 'largest');
valueMask=imfill(valueMask, 'holes');

%apply mask to the gray image
image_m=immultiply(image_gray, valueMask);

%segment worms by texture filter
E = entropyfilt(image_m);
Eim = mat2gray(E);

%apply threshold
bw = im2bw(Eim,0.7);
bwm = bwareaopen(bw, 100);
bwm(imcomplement(valueMask))=255;
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWsdil = imdilate(bwm, [se90 se0]);
bw1=BWsdil-bwareaopen(BWsdil, 100000);

BWdfill = imfill(bw1, 'holes');
seD = strel('diamond',1);
BWfinal = imerode(BWdfill,[se90 se0]);
BWfinal = imerode(BWfinal,[se90 se0]);
BWfinal=im2bw(BWfinal);
BWfinal=bwareaopen(BWfinal, 500);

%remove boundary objects
BWfinal=imclearborder(BWfinal);

skeletons=bwmorph(BWfinal, 'thin', Inf);

skeletons_pruned = bwmorph(skeletons, 'spur', 15);
%skeletons_pruned = bwmorph(skeletons_pruned, 'spur');

[labels, num_labels] = bwlabel(skeletons_pruned);
full_path1 = fullfile(dir2, sprintf('%s_skeletons.jpg', file.name));
imwrite(skeletons_pruned, full_path1,'jpg', 'Quality', 100);

%hold on;
j=[];
for i = 1:num_labels
    %// Find the ith object mask
    mask = labels == i;
    %// Find if worm has branching points
    br=bwmorph(mask,'branchpoints');
    br = sum(br(:));
     %// If area is less than a threshold
    %// don't process this object
    if br > 0
        continue;
    end
    %cut-off for minimal skeleton length
    if bwarea(mask)<20
        continue;
    end
        j=[j i];
end

%remove branching worms
[lab, num_lab] = bwlabel(BWfinal);

%straight_skels = ismember(labels, j);
straight = ismember(lab, j);

%Label each blob
[labeledImage, num_image] = bwlabel(straight);     
skelExport = regionprops(bwmorph(labeledImage, 'thin', Inf), 'Area');

%assign color to each blob
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle');

%get all the blob properties
blobMeasurements = regionprops(labeledImage, 'all');
numberOfBlobs = size(blobMeasurements, 1);
blobExport = regionprops(labeledImage, 'Area', 'Perimeter', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength');
skelExport = regionprops(bwmorph(labeledImage, 'thin', Inf), 'Area');

%filter objects by shape (eccentricity) and size
allBlobEcc = [blobMeasurements.Eccentricity];
allBlobAreas = [blobMeasurements.Area];
allowableEccIndexes = (allBlobEcc >= 0.8) ;

%allowableAreaIndexes = (allBlobAreas > 2000) & (allBlobAreas <5000);
med = median(allBlobAreas);
md=mad(allBlobAreas, 1);
med_length = median([skelExport.Area]);
allowableAreaIndexes = (allBlobAreas > med/10)& (allBlobAreas <= med*10);
%allowableAreaIndexes = (allBlobAreas >= med/2) & (allBlobAreas <= med*2);

%get actual indexes
keeperIndexes = find(allowableEccIndexes & allowableAreaIndexes);

%extract blobs that meet our criteri
keeperBlobsImage = ismember(labeledImage, keeperIndexes);

% re-label selected blobs
labeledDimeImage = bwlabel(keeperBlobsImage, 8);

wormArea = regionprops(labeledDimeImage,'area');
wormCent = regionprops(labeledDimeImage,'centroid');

%label rows in "results.txt"
rowNum = strsplit(num2str([1:numel(wormCent)]));

coloredDimeLabels = label2rgb (labeledDimeImage, 'hsv', 'white', 'shuffle'); % pseudo random color labels
imshow(coloredDimeLabels, []);

%add numbered labels to picture
hold on
for k = 1:numel(wormCent)
    c = wormCent(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 20, 'Color', 'k');
end
 
hFrame = getframe(gca);

%hold off
full_path1 = fullfile(dir2, sprintf('%s_segmentation2.jpg', file.name));
imwrite(hFrame.cdata, full_path1,'jpg', 'Quality', 100); %to print labels on pictures
hold off

%write the outline of identified objects (draw perimeter)
BWoutline = imdilate(bwperim(labeledDimeImage), strel('disk', 1));
Segout = image_m;
Segout(BWoutline) = 0;

%write segmentation image
full_path1 = fullfile(dir2, sprintf('%s_segmentation.jpg', file.name));
imwrite(Segout, full_path1,'jpg', 'Quality', 100);

%export
blobExport1 = blobExport(keeperIndexes);
skelExport1 = skelExport(keeperIndexes);
[blobExport1(:).Length] = deal(skelExport1.Area);

%%%%
%OUTPUT
full_path2 = fullfile(dir3, sprintf('%s.txt', file.name));
if ~isempty(blobExport1) 
writetable(struct2table(blobExport1, 'RowNames', rowNum),full_path2,'Delimiter','\t','WriteRowNames',true);
end

end
message='All done!'
end
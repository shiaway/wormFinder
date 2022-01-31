%%%%%%%%%%%%%%%%
%%%WORMFINDER%%%
%%%%%%%%%%%%%%%%
%Version 2, updated January 2022
%This script is designed to detect worms of different stages using 2X brightfield images
%of individual wells in 96 well plates and return their body area and length. Worms need 
%to be paralyzed and washed into an empty 96 well plate.

function wf = wormFind(mypath)
%specify directory with images
path_to_images=mypath;
files = dir(strcat(path_to_images, '/*.JPG'));

%make output directories if none exist
if ~exist(strcat(path_to_images, '/results'),'dir') mkdir(strcat(path_to_images, '/results')); end
if ~exist(strcat(path_to_images, '/results_txt'),'dir') mkdir(strcat(path_to_images, '/results_txt')); end
dir2 = strcat(path_to_images, '/results');
dir3 = strcat(path_to_images, '/results_txt');


for file = files'   
message = sprintf('Analyzing image %s ...', file.name)
full_path = fullfile(path_to_images, file.name);
image_w=imread(full_path);
% Adjust data to span data range.
image_gray = rgb2gray(image_w);
image_gray = imadjust(image_gray);
%% create a mask
% Threshold image - adaptive threshold
BW = imbinarize(image_gray, 'adaptive', 'Sensitivity', 0.7);

% Invert mask
BW = imcomplement(BW);

% Fill holes
BW = imfill(BW, 'holes');

% Open mask with line
length = 3.000000;
angle = 0.000000;
se = strel('line', length, angle);
BW = imopen(BW, se);

% Create masked image.
maskedImage = image_gray;
maskedImage(~BW) = 0;
%% create mask
valueMask = BW;
% filter non-worm componenets or multiple worms
valueMask = bwareafilt(valueMask, [50 4500]);
%% clean up worms
% remove boundary objects
BWfinal=imclearborder(valueMask);

%Label each blob
[labeledImage, num_image] = bwlabel(BWfinal);
%assign color to each blob
coloredLabels = label2rgb(labeledImage, 'hsv', 'k', 'shuffle');

%get all the blob properties
blobExport = regionprops('struct',labeledImage, 'Area', 'Perimeter', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength');
% find the longest skeleton
skeleton = false(size(BWfinal,1),size(BWfinal,2));
for i = 1:num_image
    %// Find the ith object mask
    mask = labeledImage == i;
    try 
        [maskOut,thinnedImg] = longestConstrainedPath(mask);
    catch 
        % shrinked to a single point
        maskOut = bwmorph(mask, 'thin', Inf);
    end
    bodyLen = bwarea(maskOut);
    bodyArea = bwarea(mask);

    skeleton(maskOut) = 1;
    blobExport(i).Length = bodyLen;
end

full_path1 = fullfile(dir2, sprintf('skeletons_%s.jpg', file.name));
imwrite(skeleton, full_path1,'jpg', 'Quality', 100);


wormCent = regionprops(labeledImage, 'Centroid');
%label rows in "results.txt"
rowNum = strsplit(num2str([1:numel(wormCent)]));

imshow(coloredLabels, []);
%add numbered labels to picture
hold on
for k = 1:numel(wormCent)
    c = wormCent(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 20, 'Color', 'w');
end
 
hFrame = getframe(gca);

%hold off
full_path1 = fullfile(dir2, sprintf('segmentation_%s.jpg', file.name));
imwrite(hFrame.cdata, full_path1,'jpg', 'Quality', 100); %to print labels on pictures
hold off

%write the outline of identified objects (draw perimeter)
BWoutline = imdilate(bwperim(labeledImage), strel('disk', 1));
% apply mask to the gray image
Segout = im2double(cat(3, image_gray, image_gray, image_gray));
col = [BWoutline.*1,BWoutline.*0,BWoutline.*0];
Segout(cat(3,BWoutline,BWoutline,BWoutline)) = [ones(sum(sum(BWoutline)),1);zeros(sum(sum(BWoutline)),1);zeros(sum(sum(BWoutline)),1)];
%write segmentation image
full_path1 = fullfile(dir2, sprintf('segmentation2_%s.jpg', file.name));
imwrite(Segout, full_path1,'jpg', 'Quality', 100);

%export
blobExport1 = blobExport;


%%%%
%OUTPUT
full_path2 = fullfile(dir3, sprintf('%s.txt', file.name));
if ~isempty(blobExport1) 
writetable(struct2table(blobExport1, 'RowNames', rowNum),full_path2,'Delimiter','\t','WriteRowNames',true);
end

end
message='All done!'
end

%%
function [bwOut,thinnedImg] = longestConstrainedPath(bwin,varargin)
% BWOUT = LONGESTCONSTRAINEDPATH(BW)
% BWOUT = LONGESTCONSTRAINEDPATH(BW,'thinOpt',thinOption)
%
% Calculates the longest continuous path in the infinitely thinned bw
% image, following the calculation of the bwdistgeodesic transform.
% Robustly ignores spurs. Note that only a single path is detected; if
% there are multiple paths the same length, only one will be returned.
%
% INPUTS:
% bwin:    2D binary input image.
%
% PV PAIRS:
% 'thinOpt': {'Thin','Skel'}. Thinning option.
%            'Thin' uses infinite thinning (bwmorph(bwin, 'thin', Inf);
%            'Skel' uses infinite skeletonization (bwmorph(bwin,
%            'skeleton', Inf); DEFAULT: 'Thin'.
% 'geodesicMethod':
%            'Method' used in call to bwdistgeodesic. One of:
%            {'cityblock', 'chessboard','quasi-euclidean'}. DEFAULT:
%            'quasi-euclidean'.
%
% OUTPUT:
% bwOut:   2D binary image showing the longest calculated path.
%
% thinnedImg: 2D thinned image (using infinite 'thinning' or
%          'skeletonization'  in bwmorph).
%
% % Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 
% See also: bwmorph bwdistgeodesic

% Copyright 2015 The MathWorks, Inc.

validateattributes(bwin,{'numeric' 'logical'},{'real' 'nonsparse' '2d'}, ...
	mfilename, '', 1);
[thinOpt, geodesicMethod] = parseInputs(varargin{:});
% 8-connected only:
M = size(bwin, 1);
neighborOffsets = [-1, M, 1, -M, M + 1, M-1, -M + 1, -M - 1]; %8-connected
thinOpt = lower(thinOpt);

switch thinOpt
	case 'skel'
		thinnedImg = bwmorph(bwin, 'skeleton', Inf);
	case 'thin'
		%This is the default; seems to perform better than skel
		thinnedImg = bwmorph(bwin, 'thin', Inf);
end
endpoints = find(bwmorph(thinnedImg, 'endpoints'));
if numel(endpoints)==2
	bwOut = thinnedImg;
	return
end
mask = false(size(thinnedImg));
mask(bwmorph(thinnedImg, 'endpoints')) = true; %endpoints mask
bwdg = bwdistgeodesic(thinnedImg,mask,geodesicMethod);
bwdg(bwdg==0)= NaN;

%Now the maximum position value of tmp must be on the longest path, and we
%can pare the graph by keeping only its largest two neighbors
bwOut = false(size(thinnedImg));
startPoint = find(bwdg==max(bwdg(:)));% SEED?
startPoint = startPoint(1); %In case there are multiple paths...
bwdg(startPoint) = NaN;
bwOut(startPoint)= true;
neighbors = bsxfun(@plus,startPoint,neighborOffsets);
[sortedNeighbors,inds] = sort(bwdg(neighbors));
bothNeighbors = neighbors(inds(1:2));
for ii = 2:-1:1
	activePixel = bothNeighbors(ii);
	while ~isempty(activePixel)
		bwOut(activePixel)= true;
		bwdg(activePixel) = NaN;
		
		neighbors = bsxfun(@plus,activePixel,neighborOffsets);
		activePixel = neighbors(bwdg(neighbors)==max(bwdg(neighbors)));
		if ~isempty(activePixel)% What to do for dupes?
			activePixel = activePixel(1);
		end
	end
end

	function [thinOpt,geodesicMethod] = parseInputs(varargin)
		% Setup parser with defaults
		parser = inputParser;
		parser.CaseSensitive = false;
		parser.FunctionName  = 'longestConstrainedPath';
		parser.addParameter('thinOpt','Thin');
		parser.addParameter('geodesicMethod','quasi-euclidean');
		% Parse input
		parser.parse(varargin{:});
		% Assign outputs
		r = parser.Results;
		[thinOpt,geodesicMethod] = deal(r.thinOpt,r.geodesicMethod);
	end

end
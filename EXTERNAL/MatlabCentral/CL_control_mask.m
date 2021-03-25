
infold = '/home/cchunharas/Desktop/objects_control';
outfold = '/home/cchunharas/Desktop/objects_control_mask';


imgs = dir(infold);
imgs(1:2) = [];

%%

black_name = '/home/cchunharas/Downloads/blankb.jpg';
black = imread(black_name);
black(444,:,:) = [];
black   = rgb2gray(black);


mask_name = '/home/cchunharas/Downloads/blank.jpg';
mask = imread(mask_name);
mask(444,:,:) = [];
% mask   = rgb2gray(mask);


black(ismember(black,(35:255))) = 255;
black(ismember(black,(1:34))) = 0;
imtool(black)
windex = ismember(black,255);

windex = double(windex);
windex = repmat(windex,[1 1 3]);
windex(:,:,3) = windex(:,:,1);
mask = double(mask);

mask(~windex) = 129;
mask(~~windex) = 255;

% set to rbg grey 
mask2 = 129*ones(478,498,3);
mask2(7:449,100:399,:) = mask;

for m = 1 : numel(imgs);
    img_name = sprintf('%s/%s',infold,imgs(m).name);
    img = imread(img_name);
    img = double(img);
    final_pic = mask2;
    final_pic(ismember(mask2,255)) = img(ismember(mask,255));
    final_pic = uint8(final_pic);
    % imtool(final_pic);
    imwrite(final_pic,[outfold '/' imgs(m).name],'Quality',100);
    
close all force

end


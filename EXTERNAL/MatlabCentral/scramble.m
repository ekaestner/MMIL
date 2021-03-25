%read and rescale (0-1) image
imgs = dir('/home/cchunharas/Desktop/objects_mask');
imgs(1:62) = [];
final_scramble = cell(1,numel(imgs));
final_R = zeros(1,numel(imgs));
final_m = zeros(1,numel(imgs));
for m = 2:numel(imgs);
img_name = sprintf('/home/cchunharas/Desktop/objects_mask/%s',imgs(m).name);
img1 = imread(img_name);
% img = img1(7:449,100:399,:);
% imtool(img)

close all force
bestR = .2;
for i = 1:5000
 [ImScrambled RandomPhase] = randomize_image_phase(img, .1, [],[]);
 R = corr(ImScrambled(:), im2double(img(:)));
 if R > bestR
   bestR = R
   best_img = ImScrambled;
 end
 if R > .6;
   break;
 end
end
final_scramble{m} = best_img;
final_R(m) = bestR;
final_m(m) = m;
% figure()
% fig = imtool(best_img);
imwrite(best_img,sprintf('/home/cchunharas/Desktop/objects_control/%s',imgs(m).name),'Quality',100);
end

% % % 
% foo.m = final_m;
% foo.R = final_R;
% foo.scramble = final_scramble;
% vstimCL.scramble = final_scramble;
% vstimCL.R = final_R;
% vstimCL.m = final_m;
% save('/home/cchunharas/CL/vstimCL.m','vstimCL','-v7.3');
% load('/space/mdeh4/1/halgdev/projects/cchunharas/note/CL/vstimCL.m');
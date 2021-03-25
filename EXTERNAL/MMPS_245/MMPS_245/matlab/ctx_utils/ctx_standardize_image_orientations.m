function ctx_vol_out = ctx_standardize_image_orientations(ctx_vol_in)

slice_normal = ctx_vol_in.Mvxl2lph(1:3,3);
slice_normal = slice_normal/norm(slice_normal);
sag_corr = dot(slice_normal,[1; 0; 0]);   % we want sagittal slice order to be R->L
cor_corr = dot(slice_normal,[0; -1; 0]);  % we want coronal slice order to go from P->A
hor_corr = dot(slice_normal,[0; 0; 1]);   % we want axial slice order to be I->S

if abs(sag_corr)>=0.5
  reverse_slice_flag = sag_corr<0;
elseif abs(cor_corr)>=0.5
  reverse_slice_flag = cor_corr<0;
else
  reverse_slice_flag = hor_corr<0;
end

ctx_vol_out = ctx_vol_in;

if reverse_slice_flag
  ctx_vol_out.imgs = ctx_vol_out.imgs(:,:,end:-1:1);
  ctx_vol_out.Mvxl2lph = ctx_vol_out.Mvxl2lph*[[1 0 0 0]' [0 1 0 0]' [0 0 -1 0]' [0 0 0 1]'];
  % JCR: I think the following line could be redone as 
  %   ctx_vol_out.Mvxl2lph(:,4) = ctx_vol_in.Mvxl2lph* [0 0 size(ctx_vol_out.imgs,3)+1 1]'; 
  ctx_vol_out.Mvxl2lph = [[1 0 0 0]' [0 1 0 0]' [0 0 1 0]' ([0 0 0 1]'+ctx_vol_in.Mvxl2lph*[1 1 1 1]'-ctx_vol_out.Mvxl2lph*[1 1 size(ctx_vol_out.imgs,3) 1]')]*ctx_vol_out.Mvxl2lph;
end

% The following 2 lines were commented out the evening of 12/15/05. They were uncommented at noon on 12/20/05.
ctx_vol_out.imgs = permute(ctx_vol_out.imgs,[2 1 3]); % Transpose images in-plane -- should check Mvxl2lph matrix
ctx_vol_out.Mvxl2lph = ctx_vol_out.Mvxl2lph(:,[2 1 3 4],:);

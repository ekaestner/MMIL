%% 2D CNN Epilepsy vs controls
% Leo Bonilha
% April 2020
% ref: https://www.mathworks.com/help/deeplearning/ref/connectlayers.html
%
% P -> prj_dir
% P_controls -> con_dir
% P_patients -> epd_dir
% regions = sub_reg, tot_reg
%
%

%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; clc

prj_dir = '/home/ekaestne/PROJECTS/DATA/FSL_segmented/FSL_segmented';

tissue_type_folder = 'Smoothed_8/gray'; % 'brain_'; % 'Smoothed_8/gray';
tissue_type_prefix = 'smoothed8'; % ''; % 'smoothed8'; % Change for tissue type

chosen_color       = @gray; % @jet; % will also be used to save 2D images below

% %%%%
con_dir = [ prj_dir '/' 'Controls' '/' tissue_type_folder '/' ];

epd_dir = [ prj_dir '/' 'Patients' '/' 'left' '/' tissue_type_folder  ];

resPath = [ prj_dir '/' 'results_GM_smoothed_with_Atlas_masked' ];

Atlas =  [ prj_dir '/' 'aal.nii' ]; % 'C:\Users\bonil\Dropbox\Atlases\aal.nii';

% %%%%

if exist(resPath); rmdir(resPath,'s'); end; mkdir(resPath)

%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data %%%%
fprintf('Getting image data from %s\n',con_dir)
[Con_matrix, con_names] = extract_image_data(con_dir, tissue_type_prefix);

fprintf('Getting image data from %s\n',epd_dir)
[Pat_matrix, pat_names] = extract_image_data(epd_dir, tissue_type_prefix);

% Mask Data %%%%%
regions = [29 31 33 35  37  39  41  47  55  77  81  85  89]; % regions = 1:90;
   
fprintf('Masking volumes from %s\n',con_dir)
[masked_Con_matrix] = mask_volumes(con_dir,Atlas,Con_matrix,[regions regions+1]);

fprintf('Masking volumes from %s\n',epd_dir)
[masked_Pat_matrix] = mask_volumes(epd_dir,Atlas,Pat_matrix,[regions regions+1]);

%% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
county = 1;

for jj = 95 %15:10:145 % loop through coronal planes posterior to anterior
    
    countx = 1;
    
    for tt = 12 % 10:30:60 %loop through windows from left to right
        
        countz = 1;
        
        for zz = 12 % 25:50:75 %loop through windows from ventral to dorsal

            %% Define the 2D slices            
            siz = 140;
            range_cor = tt:tt+siz; range_vert = zz:zz+siz; %range_vert = range_vert(1:numel(range_cor));
            
            %% Get the 2D slice from T1 data            
            fprintf('Getting the 2D slice from T1 data')            
            Coronal_controls = extract_slice(Con_matrix,jj,range_cor ,range_vert);  % masked_Con_matrix         
            Coronal_patients = extract_slice(Pat_matrix,jj,range_cor ,range_vert);  % masked_Pat_matrix         
            fprintf('----Done\n')
            
            %% place in conventional radiology display            
            fprintf('Placing in conventional radiology view')            
            Coronal_controls = fliplr(flipud(permute(Coronal_controls,[2,1,3])));            
            Coronal_patients = fliplr(flipud(permute(Coronal_patients,[2,1,3])));            
            fprintf('----Done\n')
            
            %% Display images to check
            fprintf('Displaying images to check')
            show_slices(Coronal_controls, 20, resPath,['Slice_y_' num2str(jj) '_and_x_' num2str(tt) '_and_z_' num2str(zz)], chosen_color);
            fprintf('----Done\n')
            
            %% Organizing data and splitting the images into train, validation and test groups (with smoting)
            fprintf('Organizing data and splitting the images into train, validation and test groups (with smoting)')
            
            tic
            X = cat(3, Coronal_controls,Coronal_patients);
            y = [ones(1,size(Coronal_controls,3)) zeros(1,size(Coronal_patients,3))];
            
            num_fold = 10;
            
            [xtrain, xval, xtest, ...
                ytrain, yval, ytest, ...
                xtrain_smoted, xval_smoted, ...
                ytrain_smoted, yval_smoted] = split(X, y, num_fold);
            fprintf('----Done\n')
            toc
            
            %% run network            
            fprintf('Running the models')
            
            for i = 1:numel(xtrain_smoted)
                
                [lgraph, options] = make_net_not_images(siz, xval(i).x, categorical(yval(i).y));
                
                [net] = trainNetwork(xtrain(i).x, categorical(ytrain(i).y), lgraph, options);
                
                ypred = classify(net,xtest);
                
                accuracy(i) = sum(ytest == ypred')/numel(ytest);             
                                
                % print activation weights
                
                fprintf('Saving activation weights\n')                
                countx = countx + 1; county = county + 1; countz = countz + 1;                
                close all force
                
                list_layers = {'conv_1','conv_2','conv_3'};      
                
                print_activation_weights(net, xtest(:,:,1), list_layers, resPath,['Slice_y_' num2str(jj) '_and_x_' num2str(tt) '_and_z_' num2str(zz)])
            
            end
           
             % Report mean accuracy    
             fprintf('Mean accuracy with %d fold cross-validation = %0.4f (SD = %0.4f)\n',numel(xtrain_smoted), mean(accuracy), std(accuracy))
            
             %% Run random models
             fprintf('Running random models')
             
             for i = 1:numel(xtrain_smoted)                 
                 [lgraph, options] = make_net_not_images(siz, xval(i).x,shuffle_labels(categorical(yval(i).y)));
                 
                 [net] = trainNetwork(xtrain(i).x,shuffle_labels(categorical(ytrain(i).y)), lgraph, options);
                 
                 ypred = classify(net,xtest);
                 
                 random_accuracy(i) = sum(ytest == ypred')/numel(ytest);
             end
             
             % Report mean accuracy
             fprintf('Mean accuracy with %d fold cross-validation = %0.4f (SD = %0.4f)\n',numel(xtrain_smoted), mean(random_accuracy), std(random_accuracy))
             
             %% Report probabilities
             [p, higher_than] = compare_real_random_dist(mean(accuracy), random_accuracy);
             save_image_real_vs_random_dist(resPath,random_accuracy,mean(accuracy), jj, zz, tt,higher_than,p)
         
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% extract_image_data %%%%%
function [FourD, names] = extract_image_data(folder, pattern)

co = dir(folder);

tic
count = 1;
for i = 1:numel(co)    
    if strfind(co(i).name,pattern)
        
        S = load_nii(fullfile(folder,co(i).name));
        
        FourD(:,:,:,count) = S.img;
        
        names{count} = co(i).name;
        
        count = count + 1;
        
        if rem(i,25) == 0            
            disp(sprintf('%d out of %d',i,numel(co)))            
        end
    end
end
toc

end
%% mask_volumes %%%%%
function [masked_fourD] = mask_volumes(folder, atlas, fourD, atlas_indices)

co = dir(folder);

sample_image = fullfile(folder,co(10).name);

I{2} = atlas; I{1} = sample_image;

flags.interp = 0;

spm_reslice(I,flags)

delete(fullfile(folder,['r' co(10).name]));
delete(fullfile(folder,['mean' co(10).name]));

[a,b,c] = fileparts(atlas);

ratlas = fullfile(a,['r' b c]);

R = load_nii(ratlas);

Rmask = zeros(size(R.img));

tempmask = zeros(size(R.img));

for i = 1:size(atlas_indices,2)    
    tempmask(R.img==(atlas_indices(i))) = 1;    
    Rmask = Rmask + tempmask;    
end

Rmask = double(Rmask>0);

for i = 1:size(fourD,4)    
    masked_fourD(:,:,:,i) = fourD(:,:,:,i).*Rmask;    
end

end
%% extract_slice %%%%%
function [threeDslice] = extract_slice(matrix,slice, rangex,rangez)
threeDslice = squeeze(matrix(rangex,slice,rangez,:));
end
%% show_slices %%%%%
function show_slices(group_slices, num_plots, path, name, chosen_color)

rows = 4;

columns = round(num_plots/rows);

f = figure;

set(f,'color',[1 1 1 ])

for i = 1:rows*columns    
    subplot(rows,columns,i)    
    imagesc(group_slices(:,:,i))    
    colormap(chosen_color(1000))    
end

print(gcf,fullfile(path,[name '_sample']),'-dpng','-r300');close all

end
%% split %%%%%
function [xtrain, xval, xtest, ...    
ytrain, yval, ytest, ...
xtrain_smoted, xval_smoted,...
ytrain_smoted, yval_smoted] = split(X, y, num_fold)

% define the imbalance between groups
cat1 = sum(y);
cat2 = size(y,2) - sum(y);

if cat1<cat2    
    smallercat = 1;    
    perc_to_add_in_smote = cat1 / size(y,2);    
else    
    smallercat = 0;    
    perc_to_add_in_smote = cat2 / size(y,2);    
end

% remove testing sample from the data
c = cvpartition(y,'Kfold',num_fold);
xtest = X(:,:,test(c,1));
ytest = categorical(y(test(c,1)));
xtest = add_fourdim_to_ThreeD(xtest);

% split the remaining data into training and validation
xmodel = X(:,:,training(c,1));
ymodel = categorical(y(training(c,1)));
xmodel = add_fourdim_to_ThreeD(xmodel);
cmodel = cvpartition(ymodel,'Kfold',num_fold);

%
for i = 1:cmodel.NumTestSets
    
    xtrain_temp = xmodel(:,:,training(cmodel,i));    
    ytrain(i).y = categorical(ymodel(training(cmodel,i)));    
    xtrain(i).x = add_fourdim_to_ThreeD(xtrain_temp);    
    xval_temp = xmodel(:,:,test(cmodel,1));    
    yval(i).y = categorical(ymodel(test(cmodel,1)));    
    xval(i).x = add_fourdim_to_ThreeD(xval_temp);
    
    %  to smote the data, define the category with fewer cases    
    [xtrain_sm, ytrain_sm] = make_smoted_data(from_matrix_to_vector_group(xtrain(i).x),ytrain(i).y, smallercat);
    [xval_sm, yval_sm] = make_smoted_data(from_matrix_to_vector_group( xval(i).x),yval(i).y, smallercat);

    xtrain_smoted(i).x = xtrain_sm;
    ytrain_smoted(i).y = ytrain_sm;
    
    xval_smoted(i).x = xval_sm;
    yval_smoted(i).y = yval_sm;

end

end
%% add_fourdim_to_ThreeD %%%%%
function FourD = add_fourdim_to_ThreeD(X)
FourD = reshape(X,[size(X,1), size(X,2), 1, size(X,3)]);
end
%% make_smoted_data %%%%%
function [xsmoted, ysmoted] = make_smoted_data(V_xtrain, ytrain, smallercat)

%%
X_full_smote_smallercat = mySMOTE(V_xtrain(ytrain == categorical(smallercat),:), 200, 5);
X_full_smote = mySMOTE(V_xtrain(ytrain ~= categorical(smallercat),:), 100, 5);

original_size = size(V_xtrain(ytrain == categorical(smallercat)),2);
num_to_add = size(V_xtrain(ytrain ~= categorical(smallercat)),2);

xsmoted = [V_xtrain(ytrain ~= categorical(smallercat),:) ; X_full_smote_smallercat(1: num_to_add,:)];

ysmoted = [ytrain categorical(repmat(smallercat,1,num_to_add-original_size))];

end
%% mySMOTE %%%%%
function X_smote = mySMOTE(X, N, k)
% from https://github.com/kedarps/MATLAB-SMOTE
% from Kedar Prabhudesai
% mySMOTE  Synthetic Minority Oversampling Technique. A technique to
% generate synthetic samples as given in: https://www.jair.org/media/953/live-953-2037-jair.pdf
%   Usage:
%   X_smote = mySMOTE(X, N, k)
%   Inputs:
%   X: Original dataset
%   N: Percentage of data-augmentation intended, Typically, N > 100, if N < 100, then N is set to 100.
%   k: number of nearest neighbors to consider while performing
%   augmentation
%   Outputs:
%   X_smote: augmented dataset containing original data as well.
%   See also datasample, randsample

T = size(X, 1);

if N < 100    
    N = 100;    
end

N = ceil(N / 100);

X_smote = X;

for i = 1:T    
    y = X(i,:);    
    % find k-nearest samples    
    [idx, ~] = knnsearch(X,y,'k',k);    
    % retain only N out of k nearest samples    
    idx = datasample(idx, N);    
    x_nearest = X(idx,:);    
    x_syn = bsxfun(@plus, bsxfun(@times, bsxfun(@minus,x_nearest,y), rand(N,1)), y);    
    X_smote = cat(1, X_smote, x_syn);    
end

end
%% from_matrix_to_vector_group %%%%
function V = from_matrix_to_vector_group(X)

for i = 1: size(X,4)    
    V(i,:) = from_matrix_to_vector_ind(X(:,:,1,i));    
end

end
%% from_matrix_to_vector_group %%%%
function Vind = from_matrix_to_vector_ind(rgb)
Vind = reshape(rgb,[1,size(rgb,1)*size(rgb,2)*size(rgb,3)*size(rgb,4)]);
end
%% make_net_not_images %%%%%
function [lgraph, options] = make_net_not_images(siz, xval, yval)

layers = [ imageInputLayer([siz+1 siz+1 1],'Name','input')
           convolution2dLayer(10,16,'Padding','same','Name','conv_1')
           batchNormalizationLayer('Name','BN_1')
           reluLayer('Name','relu_1')
           convolution2dLayer(3,32,'Padding','same','Stride',2,'Name','conv_2')
           batchNormalizationLayer('Name','BN_2')
           reluLayer('Name','relu_2')
           convolution2dLayer(3,32,'Padding','same','Name','conv_3')
           batchNormalizationLayer('Name','BN_3')
           reluLayer('Name','relu_3')
           additionLayer(2,'Name','add')
           averagePooling2dLayer(2,'Stride',2,'Name','avpool')
           fullyConnectedLayer(2,'Name','fc')
           softmaxLayer('Name','softmax')
           classificationLayer('Name','classOutput')];

lgraph = layerGraph(layers);
skipConv = convolution2dLayer(1,32,'Stride',2,'Name','skipConv');
lgraph = addLayers(lgraph,skipConv);
lgraph = connectLayers(lgraph,'relu_1','skipConv');
lgraph = connectLayers(lgraph,'skipConv','add/in2');

options = trainingOptions("sgdm",...    
                          'MaxEpochs',60, ...    
                          "ExecutionEnvironment","auto",...    
                          "InitialLearnRate",0.001,...    
                          "Shuffle","every-epoch",...    
                          "ValidationFrequency",30,...    
                          "Plots","none",...    
                          'ValidationData',{xval,yval});
end
%% shuffle_labels %%%%%
function shuffled_y = shuffle_labels(y)

idx = randperm(numel(y));
shuffled_y = y(idx);

end
%% print_activation_weights %%%%%
function print_activation_weights(net, im, list_layers,path,name)

for i = 1 : numel(list_layers)
    
    act1 = activations(net,im,list_layers{i});
    % The activations are returned as a 3-D array, with the third dimension indexing the channel on the conv1 layer. To show these activations using the imtile function, reshape the array to 4-D. The third dimension in the input to imtile represents the image color. Set the third dimension to have size 1 because the activations do not have color. The fourth dimension indexes the channel.
    sz = size(act1);
    act1 = reshape(act1,[sz(1) sz(2) 1 sz(3)]);
    
    % Now you can show the activations. Each activation can take any value, so normalize the output using mat2gray. All activations are scaled so that the minimum activation is 0 and the maximum is 1. Display the 64 images on an 8-by-8 grid, one for each channel in the layer.
    
    I = imtile(mat2gray(act1),'GridSize',[2 8]);
    
    imshow(I)
    
    print(gcf,fullfile(path,[name '_layer_weights_' list_layers{i}]),'-dpng','-r300')
    
    close all
    
    imgSize = size(im);
    imgSize = imgSize(1:2);
    [maxValue,maxValueIndex] = max(max(max(act1)));
    
    act1chMax = act1(:,:,:,maxValueIndex);
    act1chMax = mat2gray(act1chMax);
    act1chMax = imresize(act1chMax,imgSize);
    
    I = imtile({im,act1chMax});
    
    imshow(I)
    
    print(gcf,fullfile(path,[name '_max_weight_' list_layers{i}]),'-dpng','-r300')
    
    close all
    
end

end
%% compare_real_random_dist %%%%%
function [p, higher_than] = compare_real_random_dist(real_accuracy, rand_dist)

W = find(mean(real_accuracy) > rand_dist);
higher_than= 100* numel(W) / numel(rand_dist);
p = 1- (numel(W) / numel(rand_dist));

end
%% save_image_real_vs_random_dist %%%%%
function save_image_real_vs_random_dist(path,rand_accuracy,real_accuracy,jj,tt,zz, higher_than, p)

f= figure;
set(f,'color',[1 1 1])
h = histogram(rand_accuracy,10);
h.FaceAlpha = 0.4;
xl = xline(real_accuracy,'k:');
xl.LineWidth = 5;
title(sprintf('Real accuracy = %0.1f%%\nGreater than %0.1f%% of random accuracy (p = %0.3f)',100*real_accuracy,higher_than, p))
print(gcf,fullfile(path,['Compared_with_random_Slice_y_' num2str(jj) '_and_x_' num2str(tt) '_and_z_' num2str(zz)]),'-dpng','-r300')

close all

end





%% 2D CNN Epilepsy vs controls

% April 2020

% ref: https://www.mathworks.com/help/deeplearning/ref/connectlayers.html

clearvars; clc

P = '/home/ekaestne/PROJECTS/DATA/FSL_segmented/FSL_segmented/';

%% Change below for tissue type

tissue_type_prefix = 'smoothed8';

%% Need to define chosen color map

chosen_color = @gray; % will also be used to save 2D images below

%% Change below for source and saving results

P_controls = '/home/ekaestne/PROJECTS/DATA/FSL_segmented/FSL_segmented/Controls/Smoothed_8/gray';

P_patients = '/home/ekaestne/PROJECTS/DATA/FSL_segmented/FSL_segmented/Patients/left/Smoothed_8/gray';

resPath = '/home/ekaestne/PROJECTS/DATA/FSL_segmented/FSL_segmented/new_results';

if exist(resPath)
    
    cd(P)
    
    rmdir(resPath,'s')
    
end

mkdir(resPath)

% Atlas = '/home/ekaestne/PROJECTS/DATA/FSL_segmented/FSL_segmented/aal.nii';

%% Get T1 data

fprintf('Getting image data from %s',P_controls)

[Con_matrix, con_names] = extract_image_data(P_controls, tissue_type_prefix);

fprintf('Getting image data from %s',P_patients)

[Pat_matrix, pat_names] = extract_image_data(P_patients, tissue_type_prefix);

%% Mask volumes

regions = [29 31 33 35  37  39  41  47  55  77  81  85  89];

regions = 1:90;

masking = 0;

if masking
    
    fprintf('Masking volumes from %s',P_controls)
    
    [masked_Con_matrix] = mask_volumes(P_controls,Atlas,Con_matrix,[regions regions+1]);
    
    fprintf('Masking volumes from %s',P_patients)
    
    [masked_Pat_matrix] = mask_volumes(P_patients,Atlas,Pat_matrix,[regions regions+1]);
    
else
    
    masked_Con_matrix = Con_matrix;
    
    masked_Pat_matrix = Pat_matrix;
    
end

%%

county = 1;

for jj = 95 %15:10:145 % loop through coronal planes posterior to anterior
    
    countx = 1;
    
    for tt = 12 % 10:30:60 %loop through windows from left to right
        
        countz = 1;
        
        for zz = 12 % 25:50:75 %loop through windows from ventral to dorsal
            
            for iSB = 1:1000 % - THIS IS THE ACCURACY FOR LOOP
                
                %% Define the 2D slices
                
                siz = 140;
                
                range_cor = tt:tt+siz; range_vert = zz:zz+siz; %range_vert = range_vert(1:numel(range_cor));
                
                %% Get the 2D slice from T1 data
                
                fprintf('Getting the 2D slice from T1 data')
                
                Coronal_controls = extract_slice(masked_Con_matrix,jj,range_cor ,range_vert);
                
                Coronal_patients = extract_slice(masked_Pat_matrix,jj,range_cor ,range_vert);
                
                fprintf('----Done')
                
                %% place in conventional radiology display
                
                fprintf('Placing in conventional radiology view')
                
                Coronal_controls = fliplr(flipud(permute(Coronal_controls,[2,1,3])));
                
                Coronal_patients = fliplr(flipud(permute(Coronal_patients,[2,1,3])));
                
                fprintf('----Done')
                
                %% Display images to check
                
                fprintf('Displaying images to check')
                
%                 show_slices(Coronal_controls, 20, resPath,...
%                     ['Slice_y_' num2str(jj) '_and_x_' num2str(tt) '_and_z_' num2str(zz)], chosen_color);
                
                fprintf('----Done')
                
                
                
                %% Organizing data and splitting the images into train, validation and test groups (with smoting)
                
                fprintf('Organizing data and splitting the images into train, validation and test groups (with smoting)')
                
                tic
                
                X = cat(3, Coronal_controls,Coronal_patients);
                
                y = [ones(1,size(Coronal_controls,3)) zeros(1,size(Coronal_patients,3))];
                
                %
                
                num_fold = 10;
                
                [xtrain, xval, xtest, ...
                    ytrain, yval, ytest, ...
                    xtrain_smoted, xval_smoted,...
                    ytrain_smoted, yval_smoted] = split(X, y, num_fold);
                
                fprintf('----Done')
                
                toc
                
                %% run network
                
                
                
                fprintf('Running the models')
                
                f = {@make_simple_cnn, @make_dag_net, @make_alex_net, @make_resnet18, @make_googlenet};
                
                for ff = 2 %1:numel(f)
                    
                    fnames = {'Simple_CNN', 'DAGNet', 'AlexNet', 'ResNet18', 'GoogleNet'};
                    
                    for i = 1:numel(xtrain_smoted)
                        
                        [lgraph, options] = f{ff}(siz, xval(i).x, categorical(yval(i).y));
                        
                        [net] = trainNetwork(xtrain(i).x, categorical(ytrain(i).y), lgraph, options);
                        
                        ypred = classify(net,xtest);
                        
                        accuracy(i) = sum(ytest == ypred')/numel(ytest);
                        
                        % print activation weights
                        
                        fprintf('Saving activation weights')
                        
                        countx = countx + 1; county = county + 1; countz = countz + 1;
                        
                        close all force
                        
                        list_layers = {'conv_1','conv_2','conv_3'};
                        
%                         print_activation_weights(net, xtest(:,:,1), list_layers, resPath,...
%                             [fnames{ff} 'Slice_y_' num2str(jj) '_and_x_' num2str(tt) '_and_z_' num2str(zz)])
                        
                    end
                    
                    % Report mean accuracy
                    
                    fprintf('%s -- mean accuracy  with %d fold cross-validation = %0.4f (SD = %0.4f)',...
                        fnames{ff},numel(xtrain_smoted), mean(accuracy), std(accuracy))
                    
                    
                    
                    %% Run random models
                    
                    fprintf('Running random models')
                    
                    for i = 1:numel(xtrain_smoted)
                        
                        [lgraph, options] = f{ff}(siz, xval(i).x,...
                            shuffle_labels(categorical(yval(i).y)));
                        
                        [net] = trainNetwork(xtrain(i).x,...
                            shuffle_labels(categorical(ytrain(i).y)), lgraph, options);
                        
                        ypred = classify(net,xtest);
                        
                        random_accuracy(i) = sum(ytest == ypred')/numel(ytest);
                        
                    end
                    
                    % Report mean accuracy
                    
                    fprintf('Random %s Mean accuracy with %d fold cross-validation = %0.4f (SD = %0.4f)',...
                        fnames{ff}, numel(xtrain_smoted), mean(random_accuracy), std(random_accuracy))
                    
                    
                    
                    %% Report probabilities
                    
%                     [p, higher_than] = compare_real_random_dist(mean(accuracy), random_accuracy);
                    
%                     save_image_real_vs_random_dist(fnames{ff}, resPath,random_accuracy,mean(accuracy), jj, zz, tt,higher_than,p)
                    
                end
                
                hld_rel_acc_avg{iSB} = mean(accuracy);
                hld_rel_acc_std{iSB} = std(accuracy);
                hld_shf_acc_avg{iSB} = mean(random_accuracy);
                hld_shf_acc_std{iSB} = std(random_accuracy);
                
                clear random_accuracy accuracy
                
                cell2csv('/home/ekaestner/Downloads/leavitrunning.csv',[ hld_rel_acc_avg ; hld_rel_acc_std ; hld_shf_acc_avg ; hld_shf_acc_std ])
                
            end
        end
        
    end
    
end

%% Supporting functions

function shuffled_y = shuffle_labels(y)

idx = randperm(numel(y));

shuffled_y = y(idx);

end

function [xtrain, xval, xtest, ...
    ytrain, yval, ytest, ...
    xtrain_smoted, xval_smoted,...
    ytrain_smoted, yval_smoted] = split(X, y, num_fold)

%%

%% define the imbalance between groups

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

%%

for i = 1:cmodel.NumTestSets
    
    xtrain_temp = xmodel(:,:,training(cmodel,i));
    
    ytrain(i).y = categorical(ymodel(training(cmodel,i)));
    
    xtrain(i).x = add_fourdim_to_ThreeD(xtrain_temp);
    
    xval_temp = xmodel(:,:,test(cmodel,1));
    
    yval(i).y = categorical(ymodel(test(cmodel,1)));
    
    xval(i).x = add_fourdim_to_ThreeD(xval_temp);
    
    %  to smote the data, define the category with fewer cases
    
    [xtrain_sm, ytrain_sm] = make_smoted_data(from_matrix_to_vector_group(xtrain(i).x),...
        ytrain(i).y, smallercat);
    
    [xval_sm, yval_sm] = make_smoted_data(from_matrix_to_vector_group( xval(i).x),...
        yval(i).y, smallercat);
    
    
    
    xtrain_smoted(i).x = xtrain_sm;
    
    ytrain_smoted(i).y = ytrain_sm;
    
    
    
    xval_smoted(i).x = xval_sm;
    
    yval_smoted(i).y = yval_sm;
    
    
    
end

end

function [xsmoted, ysmoted] = make_smoted_data(V_xtrain, ytrain, smallercat)

%%

X_full_smote_smallercat = mySMOTE(V_xtrain(ytrain == categorical(smallercat),:), 200, 5);

X_full_smote = mySMOTE(V_xtrain(ytrain ~= categorical(smallercat),:), 100, 5);



original_size = size(V_xtrain(ytrain == categorical(smallercat)),2);

num_to_add = size(V_xtrain(ytrain ~= categorical(smallercat)),2);



xsmoted = [V_xtrain(ytrain ~= categorical(smallercat),:);...
    
X_full_smote_smallercat(1: num_to_add,:)];

ysmoted = [ytrain categorical(repmat(smallercat,1,num_to_add-original_size))];

end



function V = from_matrix_to_vector_group(X)

%%

for i = 1: size(X,4)
    
    V(i,:) = from_matrix_to_vector_ind(X(:,:,1,i));
    
end

end

function Vind = from_matrix_to_vector_ind(rgb)

Vind = reshape(rgb,[1,size(rgb,1)*size(rgb,2)*size(rgb,3)*size(rgb,4)]);

end

function FourD = add_fourdim_to_ThreeD(X)

FourD = reshape(X,[size(X,1), size(X,2), 1, size(X,3)]);

end

function [layers, options] = make_simple_cnn(siz, xval, yval)

% from https://www.mathworks.com/help/deeplearning/ug/create-simple-deep-learning-network-for-classification.html



layers = [
    
imageInputLayer([siz+1 siz+1 1])

convolution2dLayer(3,8,'Padding','same','Name','conv_1')

batchNormalizationLayer

reluLayer

maxPooling2dLayer(2,'Stride',2)

convolution2dLayer(3,16,'Padding','same','Name','conv_2')

batchNormalizationLayer

reluLayer

maxPooling2dLayer(2,'Stride',2)

convolution2dLayer(3,32,'Padding','same','Name','conv_3')

batchNormalizationLayer

reluLayer

fullyConnectedLayer(2)

softmaxLayer

classificationLayer];





options = trainingOptions("sgdm",...
    'MaxEpochs',60, ...
    "ExecutionEnvironment","auto",...
    "InitialLearnRate",0.001,...
    "Shuffle","every-epoch",...
    "ValidationFrequency",30,...
    "Plots","none",...
    'ValidationData',{xval,yval});

end

function [lgraph, options] = make_dag_net(siz, xval, yval)

layers = [
    
imageInputLayer([siz+1 siz+1 1],'Name','input')

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

function [lgraph, options] = make_alex_net(siz, xval, yval)

layers = [
    
imageInputLayer([siz+1 siz+1 1],"Name","imageinput")

convolution2dLayer([3 3],32,"Name","conv_1","Padding","same")

reluLayer("Name","relu1")

crossChannelNormalizationLayer(5,"Name","norm1","K",1)

maxPooling2dLayer([3 3],"Name","pool1","Stride",[2 2])

convolution2dLayer([3 3],32,"Name","conv_2","Padding","same")

reluLayer("Name","relu2")

crossChannelNormalizationLayer(5,"Name","norm2","K",1)

maxPooling2dLayer([3 3],"Name","pool2","Stride",[2 2])

convolution2dLayer([3 3],32,"Name","conv_3","Padding","same")

reluLayer("Name","relu3")

convolution2dLayer([3 3],32,"Name","conv_4","Padding","same")

reluLayer("Name","relu4")

convolution2dLayer([1 1],2,"Name","conv_5","BiasLearnRateFactor",10,"Padding","same","WeightLearnRateFactor",10)

reluLayer("Name","relu5")

maxPooling2dLayer([3 3],"Name","pool5","Stride",[2 2])

fullyConnectedLayer(10,"Name","fc_1")

reluLayer("Name","relu6")

dropoutLayer(0.5,"Name","drop6")

fullyConnectedLayer(10,"Name","fc_2")

reluLayer("Name","relu7")

dropoutLayer(0.5,"Name","drop7")

fullyConnectedLayer(2,"Name","fc_3")

softmaxLayer("Name","prob")

classificationLayer("Name","classoutput")];



lgraph = layerGraph(layers);

options = trainingOptions("sgdm",...
    'MaxEpochs',60, ...
    "ExecutionEnvironment","auto",...
    "InitialLearnRate",0.001,...
    "Shuffle","every-epoch",...
    "ValidationFrequency",30,...
    "Plots","none",...
    'ValidationData',{xval,yval});

end

function [lgraph, options] = make_resnet18(siz, xval, yval)

lgraph = layerGraph();

tempLayers = [
    
imageInputLayer([siz+1 siz+1 1],"Name","data","Normalization","zscore")

convolution2dLayer([7 7],64,"Name","conv_1","BiasLearnRateFactor",0,"Padding",[3 3 3 3],"Stride",[2 2])

batchNormalizationLayer("Name","bn_conv1")

reluLayer("Name","conv1_relu")

maxPooling2dLayer([3 3],"Name","pool1","Padding",[1 1 1 1],"Stride",[2 2])];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([3 3],64,"Name","conv_2","BiasLearnRateFactor",0,"Padding",[1 1 1 1]) %res2a_branch2a

batchNormalizationLayer("Name","bn2a_branch2a")

reluLayer("Name","conv_2_relu")

convolution2dLayer([3 3],64,"Name","conv_3","BiasLearnRateFactor",0,"Padding",[1 1 1 1]) %res2a_branch2b

batchNormalizationLayer("Name","bn2a_branch2b")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
additionLayer(2,"Name","res2a")

reluLayer("Name","res2a_relu")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([3 3],64,"Name","res2b_branch2a","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn2b_branch2a")

reluLayer("Name","res2b_branch2a_relu")

convolution2dLayer([3 3],64,"Name","res2b_branch2b","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn2b_branch2b")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
additionLayer(2,"Name","res2b")

reluLayer("Name","res2b_relu")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([3 3],128,"Name","res3a_branch2a","BiasLearnRateFactor",0,"Padding",[1 1 1 1],"Stride",[2 2])

batchNormalizationLayer("Name","bn3a_branch2a")

reluLayer("Name","res3a_branch2a_relu")

convolution2dLayer([3 3],128,"Name","res3a_branch2b","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn3a_branch2b")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],128,"Name","res3a_branch1","BiasLearnRateFactor",0,"Stride",[2 2])

batchNormalizationLayer("Name","bn3a_branch1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
additionLayer(2,"Name","res3a")

reluLayer("Name","res3a_relu")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([3 3],128,"Name","res3b_branch2a","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn3b_branch2a")

reluLayer("Name","res3b_branch2a_relu")

convolution2dLayer([3 3],128,"Name","res3b_branch2b","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn3b_branch2b")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
additionLayer(2,"Name","res3b")

reluLayer("Name","res3b_relu")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],256,"Name","res4a_branch1","BiasLearnRateFactor",0,"Stride",[2 2])

batchNormalizationLayer("Name","bn4a_branch1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([3 3],256,"Name","res4a_branch2a","BiasLearnRateFactor",0,"Padding",[1 1 1 1],"Stride",[2 2])

batchNormalizationLayer("Name","bn4a_branch2a")

reluLayer("Name","res4a_branch2a_relu")

convolution2dLayer([3 3],256,"Name","res4a_branch2b","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn4a_branch2b")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
additionLayer(2,"Name","res4a")

reluLayer("Name","res4a_relu")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([3 3],256,"Name","res4b_branch2a","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn4b_branch2a")

reluLayer("Name","res4b_branch2a_relu")

convolution2dLayer([3 3],256,"Name","res4b_branch2b","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn4b_branch2b")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
additionLayer(2,"Name","res4b")

reluLayer("Name","res4b_relu")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],512,"Name","res5a_branch1","BiasLearnRateFactor",0,"Stride",[2 2])

batchNormalizationLayer("Name","bn5a_branch1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([3 3],512,"Name","res5a_branch2a","BiasLearnRateFactor",0,"Padding",[1 1 1 1],"Stride",[2 2])

batchNormalizationLayer("Name","bn5a_branch2a")

reluLayer("Name","res5a_branch2a_relu")

convolution2dLayer([3 3],512,"Name","res5a_branch2b","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn5a_branch2b")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
additionLayer(2,"Name","res5a")

reluLayer("Name","res5a_relu")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([3 3],512,"Name","res5b_branch2a","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn5b_branch2a")

reluLayer("Name","res5b_branch2a_relu")

convolution2dLayer([3 3],512,"Name","res5b_branch2b","BiasLearnRateFactor",0,"Padding",[1 1 1 1])

batchNormalizationLayer("Name","bn5b_branch2b")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
additionLayer(2,"Name","res5b")

reluLayer("Name","res5b_relu")

globalAveragePooling2dLayer("Name","pool5")

fullyConnectedLayer(2,"Name","fc1000")

softmaxLayer("Name","prob")

classificationLayer("Name","ClassificationLayer_predictions")];

lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable

clear tempLayers;



lgraph = connectLayers(lgraph,"pool1","conv_2");

lgraph = connectLayers(lgraph,"pool1","res2a/in2");

lgraph = connectLayers(lgraph,"bn2a_branch2b","res2a/in1");

lgraph = connectLayers(lgraph,"res2a_relu","res2b_branch2a");

lgraph = connectLayers(lgraph,"res2a_relu","res2b/in2");

lgraph = connectLayers(lgraph,"bn2b_branch2b","res2b/in1");

lgraph = connectLayers(lgraph,"res2b_relu","res3a_branch2a");

lgraph = connectLayers(lgraph,"res2b_relu","res3a_branch1");

lgraph = connectLayers(lgraph,"bn3a_branch1","res3a/in2");

lgraph = connectLayers(lgraph,"bn3a_branch2b","res3a/in1");

lgraph = connectLayers(lgraph,"res3a_relu","res3b_branch2a");

lgraph = connectLayers(lgraph,"res3a_relu","res3b/in2");

lgraph = connectLayers(lgraph,"bn3b_branch2b","res3b/in1");

lgraph = connectLayers(lgraph,"res3b_relu","res4a_branch1");

lgraph = connectLayers(lgraph,"res3b_relu","res4a_branch2a");

lgraph = connectLayers(lgraph,"bn4a_branch2b","res4a/in1");

lgraph = connectLayers(lgraph,"bn4a_branch1","res4a/in2");

lgraph = connectLayers(lgraph,"res4a_relu","res4b_branch2a");

lgraph = connectLayers(lgraph,"res4a_relu","res4b/in2");

lgraph = connectLayers(lgraph,"bn4b_branch2b","res4b/in1");

lgraph = connectLayers(lgraph,"res4b_relu","res5a_branch1");

lgraph = connectLayers(lgraph,"res4b_relu","res5a_branch2a");

lgraph = connectLayers(lgraph,"bn5a_branch2b","res5a/in1");

lgraph = connectLayers(lgraph,"bn5a_branch1","res5a/in2");

lgraph = connectLayers(lgraph,"res5a_relu","res5b_branch2a");

lgraph = connectLayers(lgraph,"res5a_relu","res5b/in2");

lgraph = connectLayers(lgraph,"bn5b_branch2b","res5b/in1");



options = trainingOptions("sgdm",...
    'MaxEpochs',60, ...
    "ExecutionEnvironment","auto",...
    "InitialLearnRate",0.001,...
    "Shuffle","every-epoch",...
    "ValidationFrequency",30,...
    "Plots","none",...
    'ValidationData',{xval,yval});

end

function [lgraph, options] = make_googlenet(siz, xval, yval)

lgraph = layerGraph();

tempLayers = [
    
imageInputLayer([siz+1 siz+1 1],"Name","data")

convolution2dLayer([7 7],64,"Name","conv_1","BiasLearnRateFactor",2,"Padding",[3 3 3 3],"Stride",[2 2]) %conv1-7x7_s2

reluLayer("Name","conv1-relu_7x7")

maxPooling2dLayer([3 3],"Name","pool1-3x3_s2","Padding",[0 1 0 1],"Stride",[2 2])

crossChannelNormalizationLayer(5,"Name","pool1-norm1","K",1)

convolution2dLayer([1 1],64,"Name","conv_2","BiasLearnRateFactor",2) %conv2-3x3_reduce

reluLayer("Name","conv2-relu_3x3_reduce")

convolution2dLayer([3 3],192,"Name","conv_3","BiasLearnRateFactor",2,"Padding",[1 1 1 1]) %conv2-3x3

reluLayer("Name","conv2-relu_3x3")

crossChannelNormalizationLayer(5,"Name","conv2-norm2","K",1)

maxPooling2dLayer([3 3],"Name","pool2-3x3_s2","Padding",[0 1 0 1],"Stride",[2 2])];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],96,"Name","inception_3a-3x3_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_3a-relu_3x3_reduce")

convolution2dLayer([3 3],128,"Name","inception_3a-3x3","BiasLearnRateFactor",2,"Padding",[1 1 1 1])

reluLayer("Name","inception_3a-relu_3x3")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],16,"Name","inception_3a-5x5_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_3a-relu_5x5_reduce")

convolution2dLayer([5 5],32,"Name","inception_3a-5x5","BiasLearnRateFactor",2,"Padding",[2 2 2 2])

reluLayer("Name","inception_3a-relu_5x5")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
maxPooling2dLayer([3 3],"Name","inception_3a-pool","Padding",[1 1 1 1])

convolution2dLayer([1 1],32,"Name","inception_3a-pool_proj","BiasLearnRateFactor",2)

reluLayer("Name","inception_3a-relu_pool_proj")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],64,"Name","inception_3a-1x1","BiasLearnRateFactor",2)

reluLayer("Name","inception_3a-relu_1x1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = depthConcatenationLayer(4,"Name","inception_3a-output");

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
maxPooling2dLayer([3 3],"Name","inception_3b-pool","Padding",[1 1 1 1])

convolution2dLayer([1 1],64,"Name","inception_3b-pool_proj","BiasLearnRateFactor",2)

reluLayer("Name","inception_3b-relu_pool_proj")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],128,"Name","inception_3b-1x1","BiasLearnRateFactor",2)

reluLayer("Name","inception_3b-relu_1x1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],32,"Name","inception_3b-5x5_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_3b-relu_5x5_reduce")

convolution2dLayer([5 5],96,"Name","inception_3b-5x5","BiasLearnRateFactor",2,"Padding",[2 2 2 2])

reluLayer("Name","inception_3b-relu_5x5")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],128,"Name","inception_3b-3x3_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_3b-relu_3x3_reduce")

convolution2dLayer([3 3],192,"Name","inception_3b-3x3","BiasLearnRateFactor",2,"Padding",[1 1 1 1])

reluLayer("Name","inception_3b-relu_3x3")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
depthConcatenationLayer(4,"Name","inception_3b-output")

maxPooling2dLayer([3 3],"Name","pool3-3x3_s2","Padding",[0 1 0 1],"Stride",[2 2])];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],192,"Name","inception_4a-1x1","BiasLearnRateFactor",2)

reluLayer("Name","inception_4a-relu_1x1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],96,"Name","inception_4a-3x3_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4a-relu_3x3_reduce")

convolution2dLayer([3 3],208,"Name","inception_4a-3x3","BiasLearnRateFactor",2,"Padding",[1 1 1 1])

reluLayer("Name","inception_4a-relu_3x3")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],16,"Name","inception_4a-5x5_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4a-relu_5x5_reduce")

convolution2dLayer([5 5],48,"Name","inception_4a-5x5","BiasLearnRateFactor",2,"Padding",[2 2 2 2])

reluLayer("Name","inception_4a-relu_5x5")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
maxPooling2dLayer([3 3],"Name","inception_4a-pool","Padding",[1 1 1 1])

convolution2dLayer([1 1],64,"Name","inception_4a-pool_proj","BiasLearnRateFactor",2)

reluLayer("Name","inception_4a-relu_pool_proj")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = depthConcatenationLayer(4,"Name","inception_4a-output");

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],160,"Name","inception_4b-1x1","BiasLearnRateFactor",2)

reluLayer("Name","inception_4b-relu_1x1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],24,"Name","inception_4b-5x5_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4b-relu_5x5_reduce")

convolution2dLayer([5 5],64,"Name","inception_4b-5x5","BiasLearnRateFactor",2,"Padding",[2 2 2 2])

reluLayer("Name","inception_4b-relu_5x5")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],112,"Name","inception_4b-3x3_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4b-relu_3x3_reduce")

convolution2dLayer([3 3],224,"Name","inception_4b-3x3","BiasLearnRateFactor",2,"Padding",[1 1 1 1])

reluLayer("Name","inception_4b-relu_3x3")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
maxPooling2dLayer([3 3],"Name","inception_4b-pool","Padding",[1 1 1 1])

convolution2dLayer([1 1],64,"Name","inception_4b-pool_proj","BiasLearnRateFactor",2)

reluLayer("Name","inception_4b-relu_pool_proj")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = depthConcatenationLayer(4,"Name","inception_4b-output");

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
maxPooling2dLayer([3 3],"Name","inception_4c-pool","Padding",[1 1 1 1])

convolution2dLayer([1 1],64,"Name","inception_4c-pool_proj","BiasLearnRateFactor",2)

reluLayer("Name","inception_4c-relu_pool_proj")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],128,"Name","inception_4c-1x1","BiasLearnRateFactor",2)

reluLayer("Name","inception_4c-relu_1x1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],128,"Name","inception_4c-3x3_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4c-relu_3x3_reduce")

convolution2dLayer([3 3],256,"Name","inception_4c-3x3","BiasLearnRateFactor",2,"Padding",[1 1 1 1])

reluLayer("Name","inception_4c-relu_3x3")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],24,"Name","inception_4c-5x5_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4c-relu_5x5_reduce")

convolution2dLayer([5 5],64,"Name","inception_4c-5x5","BiasLearnRateFactor",2,"Padding",[2 2 2 2])

reluLayer("Name","inception_4c-relu_5x5")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = depthConcatenationLayer(4,"Name","inception_4c-output");

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
maxPooling2dLayer([3 3],"Name","inception_4d-pool","Padding",[1 1 1 1])

convolution2dLayer([1 1],64,"Name","inception_4d-pool_proj","BiasLearnRateFactor",2)

reluLayer("Name","inception_4d-relu_pool_proj")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],112,"Name","inception_4d-1x1","BiasLearnRateFactor",2)

reluLayer("Name","inception_4d-relu_1x1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],144,"Name","inception_4d-3x3_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4d-relu_3x3_reduce")

convolution2dLayer([3 3],288,"Name","inception_4d-3x3","BiasLearnRateFactor",2,"Padding",[1 1 1 1])

reluLayer("Name","inception_4d-relu_3x3")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],32,"Name","inception_4d-5x5_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4d-relu_5x5_reduce")

convolution2dLayer([5 5],64,"Name","inception_4d-5x5","BiasLearnRateFactor",2,"Padding",[2 2 2 2])

reluLayer("Name","inception_4d-relu_5x5")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = depthConcatenationLayer(4,"Name","inception_4d-output");

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],256,"Name","inception_4e-1x1","BiasLearnRateFactor",2)

reluLayer("Name","inception_4e-relu_1x1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
maxPooling2dLayer([3 3],"Name","inception_4e-pool","Padding",[1 1 1 1])

convolution2dLayer([1 1],128,"Name","inception_4e-pool_proj","BiasLearnRateFactor",2)

reluLayer("Name","inception_4e-relu_pool_proj")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],160,"Name","inception_4e-3x3_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4e-relu_3x3_reduce")

convolution2dLayer([3 3],320,"Name","inception_4e-3x3","BiasLearnRateFactor",2,"Padding",[1 1 1 1])

reluLayer("Name","inception_4e-relu_3x3")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],32,"Name","inception_4e-5x5_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_4e-relu_5x5_reduce")

convolution2dLayer([5 5],128,"Name","inception_4e-5x5","BiasLearnRateFactor",2,"Padding",[2 2 2 2])

reluLayer("Name","inception_4e-relu_5x5")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
depthConcatenationLayer(4,"Name","inception_4e-output")

maxPooling2dLayer([3 3],"Name","pool4-3x3_s2","Padding",[0 1 0 1],"Stride",[2 2])];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],32,"Name","inception_5a-5x5_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_5a-relu_5x5_reduce")

convolution2dLayer([5 5],128,"Name","inception_5a-5x5","BiasLearnRateFactor",2,"Padding",[2 2 2 2])

reluLayer("Name","inception_5a-relu_5x5")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],256,"Name","inception_5a-1x1","BiasLearnRateFactor",2)

reluLayer("Name","inception_5a-relu_1x1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
maxPooling2dLayer([3 3],"Name","inception_5a-pool","Padding",[1 1 1 1])

convolution2dLayer([1 1],128,"Name","inception_5a-pool_proj","BiasLearnRateFactor",2)

reluLayer("Name","inception_5a-relu_pool_proj")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],160,"Name","inception_5a-3x3_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_5a-relu_3x3_reduce")

convolution2dLayer([3 3],320,"Name","inception_5a-3x3","BiasLearnRateFactor",2,"Padding",[1 1 1 1])

reluLayer("Name","inception_5a-relu_3x3")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = depthConcatenationLayer(4,"Name","inception_5a-output");

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],48,"Name","inception_5b-5x5_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_5b-relu_5x5_reduce")

convolution2dLayer([5 5],128,"Name","inception_5b-5x5","BiasLearnRateFactor",2,"Padding",[2 2 2 2])

reluLayer("Name","inception_5b-relu_5x5")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],192,"Name","inception_5b-3x3_reduce","BiasLearnRateFactor",2)

reluLayer("Name","inception_5b-relu_3x3_reduce")

convolution2dLayer([3 3],384,"Name","inception_5b-3x3","BiasLearnRateFactor",2,"Padding",[1 1 1 1])

reluLayer("Name","inception_5b-relu_3x3")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
convolution2dLayer([1 1],384,"Name","inception_5b-1x1","BiasLearnRateFactor",2)

reluLayer("Name","inception_5b-relu_1x1")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
maxPooling2dLayer([3 3],"Name","inception_5b-pool","Padding",[1 1 1 1])

convolution2dLayer([1 1],128,"Name","inception_5b-pool_proj","BiasLearnRateFactor",2)

reluLayer("Name","inception_5b-relu_pool_proj")];

lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    
depthConcatenationLayer(4,"Name","inception_5b-output")

globalAveragePooling2dLayer("Name","pool5-7x7_s1")

dropoutLayer(0.4,"Name","pool5-drop_7x7_s1")

fullyConnectedLayer(2,"Name","loss3-classifier","BiasLearnRateFactor",2)

softmaxLayer("Name","prob")

classificationLayer("Name","output")];

lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable

clear tempLayers;

%     Connect Layer Branches

%     Connect all the branches of the network to create the network graph.

lgraph = connectLayers(lgraph,"pool2-3x3_s2","inception_3a-3x3_reduce");

lgraph = connectLayers(lgraph,"pool2-3x3_s2","inception_3a-5x5_reduce");

lgraph = connectLayers(lgraph,"pool2-3x3_s2","inception_3a-pool");

lgraph = connectLayers(lgraph,"pool2-3x3_s2","inception_3a-1x1");

lgraph = connectLayers(lgraph,"inception_3a-relu_pool_proj","inception_3a-output/in4");

lgraph = connectLayers(lgraph,"inception_3a-relu_3x3","inception_3a-output/in2");

lgraph = connectLayers(lgraph,"inception_3a-relu_5x5","inception_3a-output/in3");

lgraph = connectLayers(lgraph,"inception_3a-relu_1x1","inception_3a-output/in1");

lgraph = connectLayers(lgraph,"inception_3a-output","inception_3b-pool");

lgraph = connectLayers(lgraph,"inception_3a-output","inception_3b-1x1");

lgraph = connectLayers(lgraph,"inception_3a-output","inception_3b-5x5_reduce");

lgraph = connectLayers(lgraph,"inception_3a-output","inception_3b-3x3_reduce");

lgraph = connectLayers(lgraph,"inception_3b-relu_1x1","inception_3b-output/in1");

lgraph = connectLayers(lgraph,"inception_3b-relu_5x5","inception_3b-output/in3");

lgraph = connectLayers(lgraph,"inception_3b-relu_3x3","inception_3b-output/in2");

lgraph = connectLayers(lgraph,"inception_3b-relu_pool_proj","inception_3b-output/in4");

lgraph = connectLayers(lgraph,"pool3-3x3_s2","inception_4a-1x1");

lgraph = connectLayers(lgraph,"pool3-3x3_s2","inception_4a-3x3_reduce");

lgraph = connectLayers(lgraph,"pool3-3x3_s2","inception_4a-5x5_reduce");

lgraph = connectLayers(lgraph,"pool3-3x3_s2","inception_4a-pool");

lgraph = connectLayers(lgraph,"inception_4a-relu_5x5","inception_4a-output/in3");

lgraph = connectLayers(lgraph,"inception_4a-relu_pool_proj","inception_4a-output/in4");

lgraph = connectLayers(lgraph,"inception_4a-relu_1x1","inception_4a-output/in1");

lgraph = connectLayers(lgraph,"inception_4a-relu_3x3","inception_4a-output/in2");

lgraph = connectLayers(lgraph,"inception_4a-output","inception_4b-1x1");

lgraph = connectLayers(lgraph,"inception_4a-output","inception_4b-5x5_reduce");

lgraph = connectLayers(lgraph,"inception_4a-output","inception_4b-3x3_reduce");

lgraph = connectLayers(lgraph,"inception_4a-output","inception_4b-pool");

lgraph = connectLayers(lgraph,"inception_4b-relu_5x5","inception_4b-output/in3");

lgraph = connectLayers(lgraph,"inception_4b-relu_pool_proj","inception_4b-output/in4");

lgraph = connectLayers(lgraph,"inception_4b-relu_1x1","inception_4b-output/in1");

lgraph = connectLayers(lgraph,"inception_4b-relu_3x3","inception_4b-output/in2");

lgraph = connectLayers(lgraph,"inception_4b-output","inception_4c-pool");

lgraph = connectLayers(lgraph,"inception_4b-output","inception_4c-1x1");

lgraph = connectLayers(lgraph,"inception_4b-output","inception_4c-3x3_reduce");

lgraph = connectLayers(lgraph,"inception_4b-output","inception_4c-5x5_reduce");

lgraph = connectLayers(lgraph,"inception_4c-relu_1x1","inception_4c-output/in1");

lgraph = connectLayers(lgraph,"inception_4c-relu_pool_proj","inception_4c-output/in4");

lgraph = connectLayers(lgraph,"inception_4c-relu_3x3","inception_4c-output/in2");

lgraph = connectLayers(lgraph,"inception_4c-relu_5x5","inception_4c-output/in3");

lgraph = connectLayers(lgraph,"inception_4c-output","inception_4d-pool");

lgraph = connectLayers(lgraph,"inception_4c-output","inception_4d-1x1");

lgraph = connectLayers(lgraph,"inception_4c-output","inception_4d-3x3_reduce");

lgraph = connectLayers(lgraph,"inception_4c-output","inception_4d-5x5_reduce");

lgraph = connectLayers(lgraph,"inception_4d-relu_pool_proj","inception_4d-output/in4");

lgraph = connectLayers(lgraph,"inception_4d-relu_1x1","inception_4d-output/in1");

lgraph = connectLayers(lgraph,"inception_4d-relu_3x3","inception_4d-output/in2");

lgraph = connectLayers(lgraph,"inception_4d-relu_5x5","inception_4d-output/in3");

lgraph = connectLayers(lgraph,"inception_4d-output","inception_4e-1x1");

lgraph = connectLayers(lgraph,"inception_4d-output","inception_4e-pool");

lgraph = connectLayers(lgraph,"inception_4d-output","inception_4e-3x3_reduce");

lgraph = connectLayers(lgraph,"inception_4d-output","inception_4e-5x5_reduce");

lgraph = connectLayers(lgraph,"inception_4e-relu_1x1","inception_4e-output/in1");

lgraph = connectLayers(lgraph,"inception_4e-relu_pool_proj","inception_4e-output/in4");

lgraph = connectLayers(lgraph,"inception_4e-relu_5x5","inception_4e-output/in3");

lgraph = connectLayers(lgraph,"inception_4e-relu_3x3","inception_4e-output/in2");

lgraph = connectLayers(lgraph,"pool4-3x3_s2","inception_5a-5x5_reduce");

lgraph = connectLayers(lgraph,"pool4-3x3_s2","inception_5a-1x1");

lgraph = connectLayers(lgraph,"pool4-3x3_s2","inception_5a-pool");

lgraph = connectLayers(lgraph,"pool4-3x3_s2","inception_5a-3x3_reduce");

lgraph = connectLayers(lgraph,"inception_5a-relu_1x1","inception_5a-output/in1");

lgraph = connectLayers(lgraph,"inception_5a-relu_pool_proj","inception_5a-output/in4");

lgraph = connectLayers(lgraph,"inception_5a-relu_5x5","inception_5a-output/in3");

lgraph = connectLayers(lgraph,"inception_5a-relu_3x3","inception_5a-output/in2");

lgraph = connectLayers(lgraph,"inception_5a-output","inception_5b-5x5_reduce");

lgraph = connectLayers(lgraph,"inception_5a-output","inception_5b-3x3_reduce");

lgraph = connectLayers(lgraph,"inception_5a-output","inception_5b-1x1");

lgraph = connectLayers(lgraph,"inception_5a-output","inception_5b-pool");

lgraph = connectLayers(lgraph,"inception_5b-relu_5x5","inception_5b-output/in3");

lgraph = connectLayers(lgraph,"inception_5b-relu_3x3","inception_5b-output/in2");

lgraph = connectLayers(lgraph,"inception_5b-relu_1x1","inception_5b-output/in1");

lgraph = connectLayers(lgraph,"inception_5b-relu_pool_proj","inception_5b-output/in4");

options = trainingOptions("sgdm",...
    'MaxEpochs',60, ...
    "ExecutionEnvironment","auto",...
    "InitialLearnRate",0.001,...
    "Shuffle","every-epoch",...
    "ValidationFrequency",30,...
    "Plots","none",...
    'ValidationData',{xval,yval});

end

function print_activation_weights(net, im, list_layers,path,name)

for i = 1 : numel(list_layers)
    
    act1 = activations(net,im,list_layers{i});
    
    % The activations are returned as a 3-D array, with the third dimension indexing the channel on the conv1 layer. To show these activations using the imtile function, reshape the array to 4-D. The third dimension in the input to imtile represents the image color. Set the third dimension to have size 1 because the activations do not have color. The fourth dimension indexes the channel.
    
    sz = size(act1);
    
    act1 = reshape(act1,[sz(1) sz(2) 1 sz(3)]);
    
    % Now you can show the activations. Each activation can take any value, so normalize the output using mat2gray. All activations are scaled so that the minimum activation is 0 and the maximum is 1. Display the 64 images on an 8-by-8 grid, one for each channel in the layer.
    
    I = imtile(mat2gray(act1),'GridSize',[2 8]);
    
    imshow(I)
    
    print(gcf,...
        fullfile(path,[name '_layer_weights_' list_layers{i}]),...
        '-dpng','-r300')
    
    close all
    
    
    
    imgSize = size(im);
    
    imgSize = imgSize(1:2);
    
    [maxValue,maxValueIndex] = max(max(max(act1)));
    
    act1chMax = act1(:,:,:,maxValueIndex);
    
    act1chMax = mat2gray(act1chMax);
    
    act1chMax = imresize(act1chMax,imgSize);
    
    I = imtile({im,act1chMax});
    
    imshow(I)
    
    print(gcf,...
        fullfile(path,[name '_max_weight_' list_layers{i}]),...
        '-dpng','-r300')
    
    close all
    
end

end

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

print(gcf,...
    fullfile(path,[name '_sample']),...
    '-dpng','-r300')

close all



end

function [threeDslice] = extract_slice(matrix,slice, rangex,rangez)

threeDslice = squeeze(matrix(rangex,slice,rangez,:));

end

%

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

function [p, higher_than] = compare_real_random_dist(real_accuracy, rand_dist)

W = find(mean(real_accuracy) > rand_dist);

higher_than= 100* numel(W) / numel(rand_dist);

p = 1- (numel(W) / numel(rand_dist));

end

function save_image_real_vs_random_dist(fname, path,rand_accuracy,real_accuracy,jj,tt,zz, higher_than, p)

f= figure;

set(f,'color',[1 1 1])

h = histogram(rand_accuracy,10);

h.FaceAlpha = 0.4;

xl = xline(real_accuracy,'k:');

xl.LineWidth = 5;

title(sprintf('Real accuracy = %0.1f%%Greater than %0.1f%% of random accuracy (p = %0.3f)',100*real_accuracy,higher_than, p))

print(gcf,...
    fullfile(path,[fname 'Compared_with_random_Slice_y_' num2str(jj) '_and_x_' num2str(tt) '_and_z_' num2str(zz)]),...
    '-dpng','-r300')

close all

end



%% functions downloaded from the web

% from https://github.com/kedarps/MATLAB-SMOTE

% from Kedar Prabhudesai

function X_smote = mySMOTE(X, N, k)

% mySMOTE  Synthetic Minority Oversampling Technique. A technique to

% generate synthetic samples as given in: https://www.jair.org/media/953/live-953-2037-jair.pdf

%   Usage:

%   X_smote = mySMOTE(X, N, k)

%

%   Inputs:

%   X: Original dataset

%   N: Percentage of data-augmentation intended, Typically, N > 100, if N < 100, then N is set to 100.

%   k: number of nearest neighbors to consider while performing

%   augmentation

%

%   Outputs:

%   X_smote: augmented dataset containing original data as well.

%

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
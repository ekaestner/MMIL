function [Sclass] = CVclassify(XMat, class_labels, cls_options)

% CVclassify.m - Perform binary, fully cross-validated classification on N samples, each 
%                having P dimensons. Classification methods include linear discriminant
%                analysis (LDA) and quadratic discriminant analysis (QDA).  Feature
%                selection via stepwise discriminant analysis can be applied to training
%                set samples prior to classification of test samples. 
%                  Cross-validation methods include K-fold and Leave-N-Out CV.
%                  Output includes ROC curve data, feature selection probability, and 
%                classification scores for all test samples.
%
%  This function is very much so in "alpha" devel stage. Please report bugs and unexpected
%  results/features to Cooper Roddey: cooper.roddey@gmail.com
%
%
%  Usage: [Sclass] = CVclassify(XMat, class_labels, cls_options)
%
%  Inputs:
%       XMat         :  N x P data matrix.
%       class_labels :  vector (length=N) of class labels (0 or 1) for rows of XMat
%       cls_options  :  structure containing control parameters for classification
%
%   cls_options fields:
%              modeltype: functional form of classifier ['linear' or 'quadratic']
%                CIsamps: number of independent re-partitions of input data (only applies 
%                          when estim_type = 'Kfolds')
%             estim_type: ['LeaveMOut' or 'Kfolds']
%             classifier: ['DiscrimAnalysis', 'Logistic']    ('SVM' coming soon)
%                 kfolds: number of k-folds (only applies when estim_type = 'Kfolds')
%              leaveMout: number of test samples per classifier when estim_type = 'LeaveMOut' 
%    feature_select_flag: perform feature selection prior to classification [true or false]
%           f_thresh_del: feature selection feature-remove F-threshold (e.g. 2.71)
%           f_thresh_add: feature selection feature-add F-threshold (e.g. 3.84). If -1, then
%                          perform limited search for best F-thresholds.
%             tol_thresh: feature selection correlation tolerance threshold (e.g. 1.0e-03)
%             ROC_points: number of points to calculate on ROC curve
%
%  Outputs: 
%       Sclass  :   struct containing classification results.
%
%    Sclass fields:
%            features: struct w/ the following fields: 
%                      'list' = cell array where each element contains list of features
%                                 selected for a particular training subset of XMat.
%                      'prob' = probability of feature selection (length P vector)
%                 ROC: struct w/ the following fields:
%                      'cutvalues' = classification score cutoff values for generating diff
%                                    ROC curve points.
%                      'TPF' = true positive fraction values for ROC curve points
%                      'FPF' = false positive fraction values for ROC curve points
%                      'errate' = the error rate for each ROC curve point
%                      'AUC' = area under the ROC curve
%                      'AUC_std' = std deviation of ROC curve (Hanley and McNeil)
%                  CV: struct w/ the following fields:
%                      'indtest' = test sample indeces
%                      'testlabel' = labels of test samples
%                      'cls_scores' = test sample classification scores
%                      'cls_prob' = probability that test sample belongs to group 1.
%

nsubjects = size(XMat,1);
nfeatures = size(XMat,2);

classifier = cls_options.classifier;
if isequal(classifier,'Logistic') || isequal(classifier,'SVM')
   cls_options.modeltype = 'linear';
end

if cls_options.ROC_points <= 1,
   cls_options.ROC_points = 25;
end

if cls_options.feature_select_flag,
   rand('twister',sum(100*clock));
end

prior.prob = [0.5 0.5];
prior.group = unique(class_labels);
cls_options.prior = prior;

switch cls_options.estim_type
  case 'err632plus'

%    Sclass = Boot632plusStats(XMat,class_labels, classifier, cls_options);
     error('err632plus method not yet QC checked.');
  case 'BCV'

     error('BCV method not yet implemented.');

  case 'HoldOut'

     error('HoldOut method not yet implemented.');

  case 'LeaveMOut'

    flist = [1:size(XMat,2)];
    Sclass.feature.list = {flist};

    if cls_options.leaveMout < 1,

       error('cls_options.leaveMout must be positive.');

    elseif cls_options.leaveMout >= 1,

       testlabel = [];
       indtest = []; 
       cls_scores = [];
       cls_prob = [];
       cls_out = [];

       if cls_options.leaveMout == 1,
           niters = nsubjects;
       else
           error('LeaveMout > 1 case not yet implemented.');
           %niters = cls_options.CIsamps;
       end

       K = zeros(niters,1);

       switch(cls_options.feature_select_mode)
           case 'None'
               flist = [1:size(XMat,2)];
               Sclass.feature.list = {flist};
           case 'SPSS_stepwise'
               flist = getFeatureList(XMat,class_labels,cls_options);
               Sclass.feature.list = {flist};
       end

       for ii=1:niters,

          fprintf(1,'%d ', ii);

          if cls_options.leaveMout == 1,
              train = setdiff([1:nsubjects],ii);
              test = ii;
          else
              % not yet fully implemented
              % [train,test] = crossvalind('LeaveMOut',class_labels,cls_options.leaveMout);
          end

          if strcmp(cls_options.feature_select_mode,'Nested_stepwise')
             flist = getFeatureList(XMat(train,:),class_labels(train),cls_options);
             Sclass.feature.list(ii) = {flist};
          end

          if ~exist('flist'), error('FLIST does not exist'); end

          [classout, class_prob, coeff] = classifySamples(XMat(train,flist), XMat(test,flist), class_labels(train), cls_options);

          cls_out = [cls_out; classout];
          testlabel = [testlabel; class_labels(test)];
          indtest = [indtest; test];

          if strcmp(cls_options.classifier,'DiscrimAnalysis') || strcmp(cls_options.classifier,'Logistic'),
             cls_prob = [cls_prob; class_prob];
             cscores = discrim_scores(XMat(test,flist),coeff,cls_options.modeltype);
             cls_scores = [cls_scores; cscores];
             K(ii) = coeff.const;
          end
%          trainlabel = [trainlabel; class_labels(train)];
%          indtrain = [indtrain; train];

%          tscores = discrim_scores(XMat(train,flist),coeff(2,1),cls_options.modeltype);
%          train_scores = [train_scores; tscores];        

       end
       fprintf(1,'\n');

       Sclass.CV.indtest = indtest;
       Sclass.CV.testlabel = testlabel;
       Sclass.CV.cls_out = cls_out;

       if strcmp(cls_options.classifier,'DiscrimAnalysis') || strcmp(cls_options.classifier,'Logistic'),
          interval = (max(cls_scores)-min(cls_scores))/(cls_options.ROC_points-1);
          cutoffvalues = [min(cls_scores):interval:max(cls_scores)];
          Sclass.ROC = sample_ROC_curve(cls_scores,testlabel,cutoffvalues);

          Sclass.CV.cls_scores = 2*cls_scores/mean(abs(K));   % scales scores such that score=1 corresp. to stdev=1
          Sclass.CV.cls_prob = cls_prob;   % prob of label=0, when prior = [0.5 0.5]
       end
    end

  case 'Kfolds'

     Sclass.feature.list = [];

     for ii=1:cls_options.CIsamps

        fprintf(1,'CIsamp = %d\n', ii);

        [cv(ii), feature_list] = KfoldCVclassify(XMat, class_labels, cls_options);
        Sclass.feature.list = [Sclass.feature.list; feature_list];

        cls_scores = cv(ii).cls_scores;
        if ii==1,
           interval = (max(cls_scores)-min(cls_scores))/(cls_options.ROC_points-1);
           cutoffvalues = [min(cls_scores):interval:max(cls_scores)];
        end

        roc(ii) = sample_ROC_curve(cls_scores,cv(ii).testlabel,cutoffvalues);
     end

     npoints = length(cutoffvalues);
     errate = reshape([roc(:).errate],npoints,length(roc))';
     TPF = reshape([roc(:).TPF],npoints,length(roc))';
     FPF = reshape([roc(:).FPF],npoints,length(roc))';
     Sclass.ROC.errate = mean(errate,1);
     Sclass.ROC.errate_std = std(errate);
     Sclass.ROC.TPF = mean(TPF,1);
     Sclass.ROC.TPF_std = std(TPF);
     Sclass.ROC.FPF = mean(FPF,1);
     Sclass.ROC.FPF_std = std(FPF);
     Sclass.ROC.AUC = mean([roc(:).AUC]);
     Sclass.ROC.AUC_std = mean([roc(:).AUC_std]);

     Sclass.CV = cv;

  case 'Stratified_Kfolds'

     error('Stratified_Kfolds method not yet implemented.');

  otherwise
     error(sprintf('%s is not a valid CV method type.', cls_options.estim_type));
end

Sclass.feature.prob = zeros(size(XMat,2),1);
for jj=1:size(XMat,2)
    Sclass.feature.prob(jj) = length(find([Sclass.feature.list{:}] == jj))/length(Sclass.feature.list);
end

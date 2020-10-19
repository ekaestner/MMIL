function [classout, sel_features] = KfoldCVclassify(Xmat,class_labels,cls_options)
% KfoldCVclassify.m
%
%  [classout, sel_features] = KfoldCVclassify(Xmat,class_labels,cls_options)
%

sel_features = {[1:size(Xmat,2)]};  % feature list
flist = sel_features{1};
classifier = cls_options.classifier;

testlabel = []; indtest = [];
cls_prob = []; cls_scores = [];
K=zeros(cls_options.kfolds,1);

% make each fold contain approx. the same class proportions as in Xmat
indices = crossvalind('Kfold', class_labels, cls_options.kfolds);
for ii = 1:cls_options.kfolds

   test = (indices == ii); train = ~test;
   test = find(test); train = find(train);

   if cls_options.feature_select_flag   
      [flist] = getFeatureList(Xmat(train,:),class_labels(train),cls_options);
      sel_features(ii) = {flist};
   end

   % cross-validation
   [classout, pprob, coeff] = classifySamples(XMat(train,flist), XMat(test,flist), class_labels(train), cls_options);
   testlabel = [testlabel; class_labels(test)];
   indtest = [indtest; test];
   cls_prob = [cls_prob; pprob];
   K(ii) = coeff.const;

   cscores = discrim_scores(Xmat(test,flist),coeff,cls_options.modeltype);
   cls_scores = [cls_scores; cscores];
   
   % resubstitution
   % [classout, pprob, coeff] = classifySamples(XMat(train,flist), XMat(train,flist), class_labels(train), cls_options);
%   classout.train = [classout.train; class_labels(train)];
%   classout.indtrain = [classout.indtrain; train];

%   tscores = discrim_scores(Xmat(train,flist), coeff(2,1), cls_options.modeltype);
%   train_scores = [train_scores; tscores];

end

classout.indtest = indtest;
classout.testlabel = testlabel;
classout.cls_scores = cls_scores/mean(abs(K));
classout.cls_prob = cls_prob;
%classout.train_scores = train_scores;

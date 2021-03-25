function [addThresh,delThresh] = optimizeWilks(Xmat,cls_labels,cls_options)

%  [addThresh,delThresh] = optimizeWilks(Xmat,cls_labels,cls_options)
%
%    Find optimal add and rm F-threshhold vals for Wilks lambda-based feature selection
%    selection.
%
%

F_add_threshvals = [2.7 3.84 5.0]; % [3.84 5.0 7.0];         % [1.05 2.7 3.84 5.0];
F_del_threshvals = [2.0 2.71 3.2];  % [2.71 3.5 5.3];         % [0.5 2.0 2.71 3.2];

classifier = cls_options.classifier;

kfolds = 10;
indices = crossvalind('Kfold', cls_labels, kfolds);

prior.group = unique(cls_labels);
prior_gp1 = length(find(cls_labels == prior.group(1)))/length(cls_labels);
prior.prob = [prior_gp1; 1-prior_gp1];
cls_options.prior = prior;


validerr = zeros(length(F_add_threshvals),1);
for indf = 1:length(F_add_threshvals),

   cls_options.f_thresh_add = F_add_threshvals(indf);
   cls_options.f_thresh_del = F_del_threshvals(indf);

   predclass = []; testlabel = [];
   for kk=1:kfolds,
      test = (indices == kk); train = ~test;
      sfeats = selectfeatures(Xmat(train,:),cls_labels(train),cls_options);
      %classout = classifier(Xmat(test,sfeats),Xmat(train,sfeats),cls_labels(train),cls_options.type,cls_options.prior);
      [classout, pprob, coeff] = classifySamples(Xmat(train,sfeats), Xmat(test,sfeats), cls_labels(train), cls_options);
      predclass = [predclass; classout];
      testlabel = [testlabel; cls_labels(test)];
   end

   validerr(indf) = length(find(ne(predclass,testlabel)))/length(testlabel);
end

[minerr,ind] = min(validerr);

delThresh = F_del_threshvals(ind);
addThresh = F_add_threshvals(ind);

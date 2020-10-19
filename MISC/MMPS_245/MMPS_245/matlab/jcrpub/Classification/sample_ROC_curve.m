function [ROC] = sample_roc_curve(cls_scores,cls_label,cutvalues)

% CLS_LABEL is more likely to be zero for lower CLS_SCORES and to be one for higher CLS_SCORES

% If CLS_LABEL's are all the same value, then this function will not work.

npoints = length(cutvalues);

if npoints == 0,
   error('CUTVALUES cannot be empty.');
end

ROC.cutvalues = cutvalues;
ROC.TPF = zeros(npoints,1);
ROC.FPF = zeros(npoints,1);
for ii=1:npoints,
   classout = zeros(length(cls_scores),1);
   classout(find(cls_scores > cutvalues(ii))) = 1;
   ROC.errate(ii) = length(find(ne(classout,cls_label)))/length(cls_label);
   ROC.TPF(ii) = length(find(classout==1 & cls_label==1))/length(find(cls_label==1));
   ROC.FPF(ii) = length(find(classout==1 & cls_label~=1))/length(find(cls_label~=1));
end

ROC.AUC = area_under_curve(ROC.FPF,ROC.TPF);
A = ROC.AUC;
label_list = unique(cls_label);   % our two class labels (bigger label will be 2nd)
ng1 = length(find(cls_label == label_list(1)));
ng2 = length(cls_label)-ng1;

% get std dev for AUC using method of Hanley-McNeil. Radiology, 1982.
ROC.AUC_std = sqrt((A*(1-A) + (ng2-1)*(A/(2-A) - A*A) + (ng1-1)*(2*A*A/(1+A) - A*A))/(ng1*ng2));

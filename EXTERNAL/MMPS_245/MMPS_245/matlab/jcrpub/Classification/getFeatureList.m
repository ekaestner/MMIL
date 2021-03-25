function [flist] = getFeatureList(Xmat,class_labels,cls_options)

if cls_options.f_thresh_add <= 0,
     % get optimal add and rm F-threshhold vals for Wilks lambda-based feature selection
     % (actually, the rm F-threshhold val is fixed for now)
     [addF,delF] = optimizeWilks(Xmat,class_labels,cls_options);
     cls_options.f_thresh_add = addF;
     cls_options.f_thresh_del = delF;
     fprintf(1,'f_thresh_add = %f   f_thresh_del = %f \n', addF, delF);
end

% get 'best' features for this training set
flist = selectfeatures(Xmat,class_labels,cls_options);
%flist = selectfeatures2(Xmat,class_labels,cls_options)

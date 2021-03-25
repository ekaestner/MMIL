function writeCVclass2csv(R,SubjectID,featureName,basename)

filename = sprintf('%s_feature.csv',basename);
fid=fopen(filename,'w');
if fid == -1,
   error(sprintf('Cannot open %s.', filename)); 
end
fprintf(fid,'feature_name, feature_index, feature_selection_prob\n');
prob=R.feature.prob;
for ii=1:length(prob)
   fprintf(fid,'%s, %d, %.4f\n', featureName{ii}, ii, prob(ii));
end
fclose(fid);

filename = sprintf('%s_ROC.csv',basename);
fid=fopen(filename,'w');
if fid == -1,
   error(sprintf('Cannot open %s.', filename));
end
fprintf(fid,'false_pos_fract, true_pos_fract, error_rate, discrim_cutoff, ROC_AUC, ROC_AUC_stdev\n');
TPF=R.ROC.TPF;
FPF=R.ROC.FPF;
cutval=R.ROC.cutvalues;
errate=R.ROC.errate;
fprintf(fid,'%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n',FPF(1),TPF(1),errate(1),cutval(1),R.ROC.AUC,R.ROC.AUC_std);
for ii=2:length(TPF)
   fprintf(fid,'%.4f, %.4f, %.4f, %.4f\n',FPF(ii),TPF(ii),errate(ii),cutval(ii));
end
fclose(fid);

filename = sprintf('%s_scores.csv',basename);
fid=fopen(filename,'w');
if fid == -1,
   error(sprintf('Cannot open %s.', filename));
end
fprintf(fid,'SubjectID, sample_index, sample_label, discrim_score, group1_prob, group2_prob\n');
for ii=1:length(R.CV.indtest)
   fprintf(fid,'%s, %.4f, %.4f, %.4f, %.4f, %.4f\n', SubjectID{ii}, R.CV.indtest(ii), R.CV.testlabel(ii), R.CV.cls_scores(ii), R.CV.cls_prob(ii,1), R.CV.cls_prob(ii,2));
end
fclose(fid);

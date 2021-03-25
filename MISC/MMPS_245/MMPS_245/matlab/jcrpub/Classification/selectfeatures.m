function [featselect] = selectfeatures(Xmat,class_labels,cls_options)

%fprintf(1,'Got in selectfeatures.\n');

Xcols = size(Xmat,2);

label_list = unique(class_labels);
ng1 = length(find(class_labels == label_list(1)));
ng2 = length(class_labels)-ng1;

max_features = min([Xcols ng1-1 ng2-1]);

feats_best = [];
lambda_best = 1;
while (length(feats_best) < max_features)

    feats_prev = feats_best;
    lambda_prev = lambda_best;

    nfeats = length(feats_best);
    Fv = zeros(nfeats,1);
    lambda = zeros(nfeats,1);

    % first try to remove redundant, previously selected features
    if nfeats > 1,
       df_term = size(Xmat,1)-nfeats-1;  % true for 2-group case only
       for ii = randperm(nfeats),
          feats_curr = setdiff(feats_best,feats_best(ii));
          [D,P,Stats] = manova1(Xmat(:,feats_curr),class_labels);      
          Fv(ii) = df_term*(Stats.lambda - lambda_prev)/lambda_prev;
          lambda(ii) = Stats.lambda;
       end

       [Fmin,min_ind] = min(Fv);
       if Fmin < cls_options.f_thresh_del,
          feats_best = feats_best(setdiff([1:nfeats],min_ind));
          lambda_best = lambda(min_ind);
%          fprintf(1,'Feature removed = %d. F-val = %f\n', feats_prev(min_ind), Fmin);
          continue;     % skip feature-adding section
       end
    end


    if nfeats > 0,
       [D,P,Stats_best] = manova1(Xmat(:,feats_best),class_labels);
       sepcoeff = Stats_best.eigenvec(:,1);
    else
       sepcoeff = []; R = zeros([2 2]);
    end
    R = zeros([2 2]);

    % Try to add best feature
    df_term = size(Xmat,1)-nfeats-2;  % true for 2-group case only
    feats_to_test = setdiff([1:Xcols],feats_best);
    numTestFeats = length(feats_to_test);
    feats_to_test = feats_to_test(randperm(numTestFeats));
    Fv = zeros(numTestFeats,1);
    lambda = zeros(numTestFeats,1);
    for ii = 1:numTestFeats,
       ifeat = feats_to_test(ii);
       feats_curr = [feats_best ifeat];
       [D,P,Stats] = manova1(Xmat(:,feats_curr),class_labels);
       Fv(ii) = df_term*(lambda_prev - Stats.lambda)/Stats.lambda;
       lambda(ii) = Stats.lambda;
       if nfeats > 0,
          R = corrcoef(Xmat(:,ifeat), Xmat(:,feats_best)*sepcoeff);
       end
       tol(ii) = 1 - R(1,2)*R(1,2);
    end

    [Fmax,max_ind] = max(Fv);
    if (Fmax >= cls_options.f_thresh_add && tol(max_ind) >= cls_options.tol_thresh)
       feats_best = [feats_best feats_to_test(max_ind)];
       lambda_best = lambda(max_ind);
%       fprintf(1,'Feature added = %d. F-val = %f\n', feats_to_test(max_ind), Fmax);
    end

    if isequal(feats_prev,feats_best)
        % No features added or removed so we're done. Print informative report now
%        fprintf(1,'We''re done\n');
        break;
    end
end

featselect = feats_best;

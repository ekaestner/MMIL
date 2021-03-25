function [classout, class_prob, coeff] = classifySamples(Xtrain, Xtest, train_labels, cls_options)

classifier = cls_options.classifier;

switch classifier
  case 'DiscrimAnalysis',

    [classout,err,class_prob,logp,coeff] = classify(Xtest, Xtrain, train_labels, cls_options.modeltype);
    coeff = coeff(2,1); 

  case 'Logistic',

    min_label = min(train_labels);
    if min_label < 1,
       train_labels = train_labels + 1 + abs(min_label);  % makes all labels >= 1 for logda()
    end   
    FF = logda(Xtrain, train_labels, 1);
    [classout, class_prob] = classify(FF, Xtest);   % not MATLAB's built-in classify.m!
    coeff.const = FF.coefs(1);
    coeff.linear = FF.coefs(2:end)';

  case 'SVM',

    error(sprintf('classifySamples: SVM classifier coming soon!'));

  otherwise

    error(sprintf('classifySamples: invalid classifier-> %s\n', classifier));

  end

end

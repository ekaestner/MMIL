function [dscores] = discrim_scores(XM,coeff,model_type)

dscores = [];

K = coeff.const;
L = coeff.linear;

if strmatch(model_type,'linear') || strmatch(model_type,'diaglinear')

      dscores = K + XM*L;

elseif strmatch(model_type,'quadratic')

      Q = coeff.quadratic;
      dscores = K + XM*L + XM*Q*XM';

end

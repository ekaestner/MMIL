function z = dotprod(x,y);
%DOTPROD - Dot product of two vectors (deprecated)
% function z = dotprod(x,y);
% return z = x dot y
% Should use Mathworks "dot" function instead



z = ( sum((x.*y)'))';

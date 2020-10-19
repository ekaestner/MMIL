function p=crossprod_fun(x,y)
% function p=crossprod_fun(x,y)
% p(i,:)=x(i,:) x y(i,:)
% x and y have three columns

if size(x,2)~=3 | size(y,2)~=3,
  error(' Must have three columns ')
end
p=[x(:,2).*y(:,3)-x(:,3).*y(:,2),...
   x(:,3).*y(:,1)-x(:,1).*y(:,3),...
   x(:,1).*y(:,2)-x(:,2).*y(:,1)];
   
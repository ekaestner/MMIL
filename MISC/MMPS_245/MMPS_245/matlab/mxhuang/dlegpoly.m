function [P,dP1] = dlegpoly(N,x)
%DLEGPOLY - Evaluate the first N Legendre polynomials
% function [P,dP1] = dlegpoly(N,x)
%*****************************************************************************
%  function [P,dP1,dP2] = dlegpoly(N,x,flag)
%  This function evaluates the first N Legendre polynomials,the first and
%  second derivatives at the vector x.
%
%  Input:
%    N       (scalar)       last polynomial order to evaluate
%    x       (m x 1)        input vector of values to eval. the polys at
%    flag    (scalar)       if flag=1, compute dP2
%
%  Outputs:
%    P       (N x m)      = [P1(x1) P1(x2)    ...    P1(xm);
%                            P2(x1) P2(x2)    ...    P2(xm);
%                             .
%                             .
%                            PN(x1) PN(x2)    ...    PN(xm)]
% 
%    dP1      (N x m)      =  derivative w.r.t. x of P
%*****************************************************************************


if size(x,1)~=1, x=x'; end      % make it a row vector

m=size(x,2);  % # evaluation points

P=zeros(N,m);     % preallocate
dP1=zeros(N,m);

P(1,:)=x;              % initials of P
P(2,:)=1.5*x.^2-0.5;  
dP1(1,:)=ones(1,m);    % initial derivative

% recursively compute P,dP1
for n=3:N,
  P(n,:)=( (2*n-1).*x.*P(n-1,:)-(n-1)*P(n-2,:) )/n;
end

for n=2:N
  dP1(n,:)=x.*dP1(n-1,:)+n.*P(n-1,:);     
end

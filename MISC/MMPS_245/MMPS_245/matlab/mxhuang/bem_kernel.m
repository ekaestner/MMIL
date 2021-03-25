function [Kv,Kb]= bem_kernel(rQ,basis,test,geometry,nodes,cdv,Te,P_wts,Tm,S);
%BEM_KERNEL - computes kernel matrices for EEG and MEG
% function [Kv,Kb]= bem_kernel(rQ,basis,test,geometry,nodes,cdv,Te,P_wts,Tm,S);
%*******************************************************************************
%
% This computes kernel matrices for EEG and MEG
%       Gv=Te*Gq   
%       Gb= Gp+Tm*Gq
%       Te and Tm are computed off-line using bem1.m
%   
%  Input parameters:
%  rQ        : dipole locations,                            no_dipole x 3
%  To get the following, just load the .mat file saved by bem1.m
%  basis     : string --  constant bem or linear bem
%  test      : string --  collocation, galerkin
%  N         : [no of EEG sensors, total no of triangles/ nodes]   2 x1 
%              e.g. [492,1476].
%  geometry,nodes,cdv : see bem1.m for description
%  S         : MEG sensor locations       
%                   no_sensor x 3
%  P_wts     : EEG Triangle Node Weights (Applies to Linear Basis only;
%              (Set to empty if not used)                          (Meeg x 3)    
%  Te,Tm     : transfer matrices
%  Be sure to provide ALL input parameters
%
%  Output Parameters:
%  Kv        : kernel matrix of EEG ,              no_eegsensor x (3*no_dipoles)
%              Kv(i,3*j-2:3*j): EEG due to the x,y,and z component of the j-th
%              dipole at the i-th EEG sensor
%
%  Kb        : kernel matrix of MEG ,              (3*no_megsensor) x (3*no_dipole)
%              Kb(3*i-2:3*i,3*j-2):  MEG due to x component of the j-th dipole
%              Kb(3*i-2:3*i,3*j-1):  MEG due to y component of the j-th dipole
%              Kb(3*i-2:3*i,3*j)  :  MEG due to z component of the j-th dipole
%              
%
% Last Modified by M.X. Huang, PhD, March 2004.
%

if size(rQ,2)~=3,
    error(' rQ must have three columns')
end
%
%%%% This part determines BEM type and EEG and/or MEG Processing %%%%
%
mode = 0;
if ~isempty(Te), mode=1; end
if ~isempty(Tm), mode=2; end
if ~isempty(Te) & ~isempty(Tm), mode=3; end
if mode==0
  error('Both EEG and MEG transfer matrices are empty');
end;

if lower(basis(1)) =='c', basis_opt = 0; end  % constant BEM
if lower(basis(1)) =='l', basis_opt =1; end  % linear BEM

if lower(test(1)) =='c', choice =0; end  % collocation
if lower(test(1)) =='g', choice =1; end  % galerkin
%
%%%% This part preallocates variables and computes triangle vertices %%%%
%
no_dipole =  size(rQ,1);
%
r1=nodes(geometry(:,1),1:3);  % vertices of the ith triangle 
r2=nodes(geometry(:,2),1:3);
r3=nodes(geometry(:,3),1:3);
%
%Gq=zeros(N(2),3*no_dipole); 
if ~isempty(Te)
    Gq=zeros(size(Te,2),3*no_dipole); 
else
    Gq=zeros(size(Tm,2),3*no_dipole); 
end

%
%%% This part performs processing unique to Constant and Linear Basis Functions %%%%
%
if basis_opt == 0,                          % Constant BEM
    whichsurf = nodes(geometry(:,1),4);
    r = (r1+r2+r3)/3;                         % centroid of each triangle
    no_elt = size(geometry,1);
else                                        % Linear BEM
    whichsurf = nodes(:,4); 
    r = nodes(:,1:3);                         % position of each node
    no_node = size(nodes,1);
end
%
cdvsum = cdv(whichsurf,1)+cdv(whichsurf,2);    % conductivity sum (size depends 
%                                                on Basis Function used
%%%% This part computes Gq for Each of the BEM Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if choice ==0,                              % collocation (constant -or- linear)
    c = 2*pi*cdvsum;                         % (vectorized)
    Gq = gterm_constant(r,rQ);
    Gq = Gq./repmat(c,1,3*no_dipole);
    %
else                                        % Galerkin forms 
    %
    if basis_opt ==0,                         % Constant Galerkin (vectorized)
        %
        %      c = 2*pi*cdvsum.*area;   
        c = 4*pi*cdvsum;        % (Residual term after area/2 cancels out in nit_Gq_cg)
        Gq = nit_gq_cg(r1,r2,r3,rQ);
        Gq = Gq./repmat(c,1,3*no_dipole);
        %    
    elseif basis_opt == 1,                   % Linear Galerkin, 4x4 point GQ
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        whichsurf = nodes(geometry(:,1),4);
        cdvsum = cdv(whichsurf,1)+cdv(whichsurf,2);    % conductivity sum
        %
        DET = dotprod(r1,crossprod(r2,r3));  % Nt x 1
        rnxrn1 = [crossprod(r2,r3),crossprod(r3,r1),crossprod(r1,r2)]; % (Nt x 9)
        %
        ptri = zeros(no_node,3);
        tri_data = zeros(6*no_node,15);
        t_cnt = 0;
        %
        for n=1:no_node  % Loop through all nodes and develop "triangle list"     
            ptri = [];    
            for  i=1:3,    % find triangles containing this node
                t =find(geometry(:,i)==n); 
                ptri= [ptri; t,i*ones(size(t))];  % Save record of all triangles
            end
            %
            for i=1:size(ptri,1)  % Loop through all triangles attached to node
                t= ptri(i,1);     % triangle number
                t_cnt = t_cnt + 1; % Running count of all triangles
                tri_data(t_cnt,:) =[n t ptri(i,2) r1(t,:) r2(t,:) r3(t,:) ....
                        rnxrn1(t,3*ptri(i,2)-2:3*ptri(i,2))]; 
            end
        end  % node index
        tri_data = tri_data(1:t_cnt,:);  % Delete excess terms
        %
        %%% This part integrates and performs post-scaling on all triangles %%%%
        %
        Gq_all = nit_gq_lg(tri_data,rQ);  % Integrate all triangles
        term = (2*pi)*cdvsum(tri_data(:,2),1).*DET(tri_data(:,2),1);
        Gq_all = Gq_all./repmat(term,1,3*no_dipole);
        %
        %%%% This part sums together all triangles associated with a given node %%%%
        %
        for n=1:no_node
            indx = find(tri_data(:,1)==n);
            Gq(n,:) = sum(Gq_all(indx,:));
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end      
end    
%
%%%% This part computes the final Kernel for EEG %%%%
%
Kv=[];
if mode==1 | mode==3,
    if basis_opt==1  % Linear Basis (Interpolated Between 3 Triangle Vertices)
        Kv = zeros(size(Te,1)/3,3*no_dipole);
        for m=1:(size(Te,1)/3)
            Kv(m,:) = sum(repmat(P_wts(m,:)',1,3*no_dipole).*(Te(3*m-2:3*m,:)*Gq));         
        end 
    else               
        Kv = Te*Gq;   % Constant Basis  
    end     
end
%
%%%% This part computes the final Kernel for MEG %%%%
%
Kb=[];
if mode==2 | mode ==3
    no_sensor=size(S,1);
    m=3*no_sensor;
    Gp=zeros(m,3*no_dipole);  
    % compute Gp
    for i=1:no_dipole
        S_rQ=[S(:,1)-rQ(i,1),S(:,2)-rQ(i,2),S(:,3)-rQ(i,3)];
        mag3=rownorm(S_rQ).^3;
        t= S_rQ ./[mag3,mag3,mag3];
        Gp(1:3:m,[3*i-1,3*i])=[t(:,3),-t(:,2)];
        Gp(2:3:m,[3*i-2,3*i])=[-t(:,3),t(:,1)];
        Gp(3:3:m,[3*i-2,3*i-1])=[t(:,2),-t(:,1)];
    end
    
    Kb=Gp+Tm*Gq;      % u0/4pi OMITTED HERE
end     


% SUBFUNCTIONS -----------------------------------------------------------------------


function I = nit_gq_lg(tri_data,Rq)
%NIT_GQ_LG : numerical integration over tessel tri's for const galerkin function
% function I = nit_gq_lg(tri_data,Rq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function integrates the vector function gterm_constant 
% (for the constant galerkin option function over triangles 
% with vertices r1,r2,and r3 using 4-point Gauss Quadrature
%
% where  r=(x,y,z) is a point on the triangle, and P1,P2,.... are parameters
% defined in F.
%
% Inputs:
% r1, r2,r3 : specify the vertices of triangles over which function is to be 
%             integrated (Nt x 3) 
%       Rq  : Dipole Positions (P x 3) 
%
% Output:
%         I : integral (Nt x 3P)
%
% %%% John Ermer 04/02/00
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%<autobegin> -------- 20-Nov-2002 14:07:25 ------------------------------
% ---- Automatically Generated Comments Block using auto_comments -----------
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\crossprod.m
%   toolbox\glg1.m
%   toolbox\rownorm.m
%
% At Check-in: $Author: Mosher $  $Revision: 15 $  $Date: 6/14/04 3:37p $
%
% Copyright (c) 2002 BrainStorm MMII
% Source code may not be distributed in original or altered form.
% See bst_splashscreen, http://neuroimage.usc.edu, or email leahy@sipi.usc.edu
%   for license and copyright notices.
%<autoend> ---------- 20-Nov-2002 14:07:25 ------------------------------


%profile on -detail builtin
%
r1 = tri_data(:,4:6);
r2 = tri_data(:,7:9);
r3 = tri_data(:,10:12);
rnxrn1 = tri_data(:,13:15);
%
%-----------------------------------------------------------
% W: weight table, X : abscissa table
X = [0.3399810, 0.8611363];
W = [0.6521452, 0.3478548];
%-----------------------------------------------------------
%
Nt = size(r1,1);
P = size(Rq,1);
r2_r1=r2-r1;   % (Nt x 3)
r3_r1=r3-r1;   % (Nt x 3)
halfarea = rownorm(crossprod(r2_r1,r3_r1))/4;  % (Nt x 1)
%
I = zeros(Nt,3*P);    % dimension of integral is 3
half_r2_r1= r2_r1/2;  % (Nt x 3)
%
%%% This part loops through 16 points on triangle, evaluates them, and sums them
%
for i=1:2
    t1 = 0.5*(1-X(i));    % scalar
    t2 = 0.5*(1+X(i));
    ct2 = r1+r3_r1*t2;   % common terms
    ct1 = r1+r3_r1*t1;
    s1 =zeros(size(I));
    s2 =zeros(size(I));
    for j=1:2
        d = half_r2_r1*t1;
        r = ct2+ d*(1+X(j));
        Ta = glg1(r,Rq,rnxrn1);
        r = ct2+ d*(1-X(j));
        Tb = glg1(r,Rq,rnxrn1);
        %
        d= half_r2_r1*t2;
        r = ct1+ d*(1+X(j));
        Tc = glg1(r,Rq,rnxrn1);
        r = ct1+ d*(1-X(j));
        Td = glg1(r,Rq,rnxrn1);
        %
        s1 = s1+W(j)*(Ta+Tb);
        s2 = s2+W(j)*(Tc+Td);
        %
    end
    I = I + W(i)*(t1*s1 +t2*s2);
end
%
I = I.*repmat(halfarea,1,3*P);


% Subfunctions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = glg1(r,rq,rnxrn1)
%GLG1 (overwrite succinct one line summary here)
% function g = glg1(r,rq,rnxrn1)
% gterm_constant: Integration Function for Linear Galerkin Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%         r : Location Points on Triangles (Nt x 3) 
%        rq : Dipole Positions (P x 3) 
%     rnxrn1: cross-product term (Ntx3
%
% Output:
%         g : Function Result (Nt x 3P)
%
% %%% John Ermer 04/02/00
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%<autobegin> -------- 20-Nov-2002 14:05:57 ------------------------------
% ---- Automatically Generated Comments Block using auto_comments -----------
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\dotprod.m
%   toolbox\rownorm.m
%
% At Check-in: $Author: Mosher $  $Revision: 15 $  $Date: 6/14/04 3:37p $
%
% Copyright (c) 2002 BrainStorm MMII
% Source code may not be distributed in original or altered form.
% See bst_splashscreen, http://neuroimage.usc.edu, or email leahy@sipi.usc.edu
%   for license and copyright notices.
%<autoend> ---------- 20-Nov-2002 14:05:57 ------------------------------


%
Nt = size(r,1);
P = size(rq,1);
%
%%%% Modified version incorporating vectorized processing
%
r_rq = repmat(r,P,1) - reshape(repmat(rq',Nt,1),3,Nt*P)'; % (Nt*P)x(3)
n3 = norlig(r_rq)'; %rownorm(r_rq);                                       % (Nt*P)x(3)
n3 = n3.*n3.*n3;                                          % (Nt*P)x(3)
dprod = repmat(dotprod(rnxrn1,r),1,P); 
%
g = zeros(Nt,3*P);
g(:,[1:3:3*P-2]) = dprod.*reshape(r_rq(:,1)./n3,Nt,P);  % x-comp
g(:,[2:3:3*P-1]) = dprod.*reshape(r_rq(:,2)./n3,Nt,P);  % y-comp
g(:,[3:3:3*P]) = dprod.*reshape(r_rq(:,3)./n3,Nt,P);    % z-comp
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

function I = nit_gq_cg(r1,r2,r3,Rq)
%NIT_GQ_CG : numerical integration over tessel tri's for const galerkin function
% function I = nit_gq_cg(r1,r2,r3,Rq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function integrates the vector function gterm_constant 
% (for the constant galerkin option function over triangles 
% with vertices r1,r2,and r3 using 4-point Gauss Quadrature
%
% where  r=(x,y,z) is a point on the triangle, and P1,P2,.... are parameters
% defined in F.
%
% Inputs:
% r1, r2,r3 : specify the vertices of triangles over which function is to be 
%             integrated (Nt x 3) 
%       Rq  : Dipole Positions (P x 3) 
%
% Output:
%         I : integral (Nt x 3P)
%
% %%% John Ermer 04/02/00
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%<autobegin> -------- 20-Nov-2002 14:07:23 ------------------------------
% ---- Automatically Generated Comments Block using auto_comments -----------
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\gterm_constant.m
%
% At Check-in: $Author: Mosher $  $Revision: 15 $  $Date: 6/14/04 3:37p $
%
% Copyright (c) 2002 BrainStorm MMII
% Source code may not be distributed in original or altered form.
% See bst_splashscreen, http://neuroimage.usc.edu, or email leahy@sipi.usc.edu
%   for license and copyright notices.
%<autoend> ---------- 20-Nov-2002 14:07:23 ------------------------------


%profile on -detail builtin
%
%-----------------------------------------------------------
% W: weight table, X : abscissa table
X = [0.3399810, 0.8611363];
W = [0.6521452, 0.3478548];
%-----------------------------------------------------------
%
Nt = size(r1,1);
P = size(Rq,1);
r2_r1=r2-r1;   % (Nt x 3)
r3_r1=r3-r1;   % (Nt x 3)
%
I = zeros(Nt,3*P);    % dimension of integral is 3
half_r2_r1= r2_r1/2;  % (Nt x 3)
%
%%% This part loops through 16 points on triangle, evaluates them, and sums them
%
for i=1:2
    t1 = 0.5*(1-X(i));    % scalar
    t2 = 0.5*(1+X(i));
    ct2 = r1+r3_r1*t2;   % common terms
    ct1 = r1+r3_r1*t1;
    s1 =zeros(size(I));
    s2 =zeros(size(I));
    for j=1:2
        d = half_r2_r1*t1;
        r = ct2+ d*(1+X(j));
        Ta = gterm_constant(r,Rq);
        r = ct2+ d*(1-X(j));
        Tb = gterm_constant(r,Rq);
        %
        d= half_r2_r1*t2;
        r = ct1+ d*(1+X(j));
        Tc = gterm_constant(r,Rq);
        r = ct1+ d*(1-X(j));
        Td = gterm_constant(r,Rq);
        %
        s1 = s1+W(j)*(Ta+Tb);
        s2 = s2+W(j)*(Tc+Td);
        %
    end
    I = I + W(i)*(t1*s1 +t2*s2);
end
%
%I = I.*repmat(halfarea,1,3*P); % Cancels out in calling program




function g = gterm_constant(r,rq)
%gterm_constant 
% function g = gterm_constant(r,rq)

%<autobegin> -------- 20-Nov-2002 14:06:02 ------------------------------
% ---- Automatically Generated Comments Block using auto_comments -----------
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\rownorm.m
%
% At Check-in: $Author: Mosher $  $Revision: 15 $  $Date: 6/14/04 3:37p $
%
% Copyright (c) 2002 BrainStorm MMII
% Source code may not be distributed in original or altered form.
% See bst_splashscreen, http://neuroimage.usc.edu, or email leahy@sipi.usc.edu
%   for license and copyright notices.
%<autoend> ---------- 20-Nov-2002 14:06:02 ------------------------------


if size(rq,1)  == 1 % Just one dipole
    r_rq= [r(:,1)-rq(1),r(:,2)-rq(2),r(:,3)-rq(3)];
    n = rownorm(r_rq).^3;
    g = r_rq./[n,n,n];
else
    g = zeros(size(r,1),3*size(rq,1));
    isrc = 1;
    for k = 1:size(rq,1)
        r_rq= [r(:,1)-rq(k,1),r(:,2)-rq(k,2),r(:,3)-rq(k,3)];
        n = rownorm(r_rq).^3;
        g(:,3*(isrc-1)+1: 3*isrc) = r_rq./[n,n,n];
        isrc = isrc + 1;
    end
    
end

function return_var = bem_transfer_1shell(bem_input,NVertMax,LU_filename,Verbose);

% function return_var = bem_transfer_1shell(bem_input,NVertMax,LU_filename,Verbose);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EEG/MEG FORWARD MODEL USING BOUNDARY ELEMENT METHODS
%                      - TRANSFER MATRIX CALCULATION (bem_transfer_1shell.m)
%
% This function computes the "transfer matrix" associated with the voltage
% potential  forward gain matrix for an array of EEG electrodes located on the
% outermost layer of a single layer surface,
%                              -AND/OR-
% the radial magnetic field forward gain matrix for an array of MEG sensors
% located on the  outermost layer of a single layer surface.
% 
% Surface data for each layer is defined via a user specified tesselated grid.
% Each region  of the multilayer surface is assumed to be concentric with
% isotropic conductivity.   EEG sensors are assumed to be located on the surface
% of the single layer at one of the surface tesselation points. Should
% specified EEG sensor points not coincide with a point on the tesselated
% outer surface layer grid, the sensor location will be quantized to the
% closest outer surface tesselation point. Dipole generator(s) are assumed to be
% interior  to the single layer surface
%
% By M.X. Huang, PhD March 2004
% 
% INPUTS (Required)
%       bem_input.R_eeg  : EEG sensor locations on scalp (meters)         (Meeg x 3)
%                         (insert dummy value if mode = 2)
%       bem_input.R_meg  : MEG sensor locations on scalp (meters)         (Mmeg x 3)
%                         (insert dummy value if mode = 1)
%       bem_input.O_meg  : MEG sensor orientations (unit-vector)          (Mmeg x 3)
%                         (insert dummy value if mode = 1)
%    bem_input.vertices  : Tessalated Surface layers from INNERMOST skull. 
%                          Packed as Cell Array where index 1 is innermost surface  (1 x 1)-Cell Array
%                          Individual entries are position row vectors in PCS
%       bem_input.faces  : Tesselated Surface Node connections for each triangle    (1 x 1)-Cell Array
%                          Packed as Cell Array
%                          Individual entries are indexes to "vertices"
%       bem_input.sigma  : conductivity from brain tissue        (1 x 1)
%       bem_input.mode   : 1,  compute EEG only
%                          2,  compute MEG only
%                          3,  compute both EEG and MEG
%   bem_input.basis_opt  : 0, constant basis
%                          1, linear
%   bem_input.test_opt   : 0, collocation,
%                          1, Galerkin
%     bem_input.ISA      : 0,  Inhibit Isolated Skull Approach
%                        : 1,  Enable Isolated Skull Approach
%    bem_input.fn_eeg    : EEG transfer matrix filename                   Char String
%                         - EEG Transfer matrices are computed and then saved 
%                          in this "*.mat" file
%
%    bem_input.fn_meg    : MEG transfer matrix filename                   Char String
%                         - MEG Transfer matrices are computed and then saved 
%                         in this "*.mat" file
%
%         WHERE: M = # of sensors; P = # of dipoles;
%               Nn = # of Surface Tesselation Nodes; Nt = # of Surf Tess Triangles 
%
%       LU_filename     : file name containing previously calculated L_bem, U_bem matrices
%                         If file does not exist, L_bem and U_bem will be calculated and 
%                         the newly calculated L_bem and U_bem are save in a file with this
%                         name.
%
% INPUTS (Optional):
%    NVertMax : Maximum Number of vertices for Surface Tesselation  (for a single layer)
%                   - Surface is re-tesselated if # points exceeds max      
%                   - Default Maximum Number of Triangle Faces is 2292    (1 x NL)
%
%    Verbose : Toggles Verbose mode on/off (on = 1);
%
% OUTPUTS: 
%    return_var: structure variable
%           return_var.basis_opt: string --- constant bem or linear bem
%           return_var.test_opt: string -- collocation, glerkin
%           return_var.mode: mode type, 1 for EEG, 2 for MEG, 3 for both
%           return_var.geometry: geometry matrix
%           return_var.nodes: nodes matrix
%           return_var.cdv: conductivity vector
%           return_var.ISA: 1 for ISA, 0 for no ISA
%           return_var.Te: the EEG transfermation matrix
%           return_var.Tm: the MEG transformation matrix
%           return_var.Te_ISA: the EEG transfermation matrix with ISA
%           return_var.Tm_ISA: the MEG transformation matrix with ISA
%           return_var.eegsensor: eegsensor
%           return_var.R_eeg: R_eeg
%           return_var.P_sens: P_sens
%           return_var.R_meg: the MEG sensor locations, no_sensor by 3
%           return_var.O_meg: the MEG sensor orientations, no_sensor by 3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Last modification: M.X. Huang, March 2004
% Adding and modifying return_var to the output, M. Huang
%
% Other people involved in old version: J. Chang, J. Mosher, E. Ermer

% unpacking the input variables

R_eeg=bem_input.R_eeg;
R_meg=bem_input.R_meg;
O_meg=bem_input.O_meg;
vertices=bem_input.vertices;
faces=bem_input.faces;
sigma=bem_input.sigma;
mode=bem_input.mode;
basis_opt=bem_input.basis_opt;
test_opt=bem_input.test_opt;
ISA=bem_input.ISA;
fn_eeg=bem_input.fn_eeg;
fn_meg=bem_input.fn_meg;

clear bem_input

NVertMax_default = 10000;

if Verbose
    fprintf('Start the computation of BEM Transfer Matrices:\n')
end

Meeg = size(R_eeg,1);   % Number of EEG Sensors
Mmeg = size(R_meg,1);   % Number of MEG Sensors
NL = size(vertices,2);  % Number of Surface Layers


%
%%% THIS PART CHECKS INPUT PARAMETERS FOR THEIR DIMENSION AND VALIDITY %%%
%
%
if (size(R_eeg,1)>=1)&(size(R_eeg,2)~=3)
    error('EEG sensor location(s) must have three columns!!!')
end
%
if (size(R_meg,1)>=1)&(size(R_meg,2)~=3)
    error('MEG sensor location(s) must have three columns!!!')
    if max(abs(rownorm(O_meg)-ones(Mmeg,1))>1.0e-3)
        error('MEG Sensor Orientations are not Unity!!!')
    end
    if size(O_meg,1)~=size(R_meg,1)
        error('Number of MEG Sensors must equal number of MEG Sensor Orientations!!!')
    end
end
%
if (Mmeg>=1)&(size(O_meg,2)~=3)
    error('MEG sensor orientation(s) must have three columns!!!')
end
%
%%%% THIS PART CHECKS THE DIMESION OF PARAMETERS ASSOCIATED WITH EACH SURFACE LAYER %%%%
%
for k=1:NL
    if size(vertices{k},2)~=3
        vertices{k} = vertices{k}';
    end
    
    if size(faces{k},2)~=3
        faces{k} = faces{k}';
    end
    
    if size(faces)~=size(vertices)
        error('Variables "vertices" and "faces" must have same number of layers (NL)!!!')
    end
    
    Nv(k) = size(vertices{k},1);      
    % Number of Surface Node (tesselation) Points (for each surface)
    
    Nf(k) = size(faces{k},1);         
    % Number of Surface Triangles (for each surface
    
end
%
%%% THIS PART DOWNSAMPLES THE NUMBER OF TESSELATION NODES IF THE NUMBER OF POINTS EXCEEDS A %%%%%%%%%
%%% USER SPECIFIED -OR- DEFAULT MAXIMUM %%%%%%
%
if ~exist('NVertMax','var') | isempty(NVertMax),
     NVertMax = NVertMax_default;
end
%
for k=1:NL
    fv(k).faces = double(faces{k});
    fv(k).vertices = double(vertices{k});
    if Nf(k) > NVertMax
        if Verbose & k==1
            fprintf('Reducing all %d tessellated layers to %d vertices. . .\n',NL,NVertMax);
        end
        % Reduce tesselation to required maximum number of vertices 
        % (call to reducepatch/qslim necessitates specification of the
        % number of faces though)
        [fv(k)] = reducepatch(fv(k),2*NVertMax-4); 
        %[fv(k)] = qslim(fv(k),'-t',2*NVertMax-4); 

        % Check for triangle orientation (need to point outwards)
        [fv(k), status] = tessellation_outwards(fv(k));
        
        if ~status & Verbose
            fprintf('Surface #%d had wrong triangle orientation; changed all.\n',k);
        end
        
        
        Nf(k) = size(fv(k).faces,1);                  % Reduced number of triangular faces
        Nv(k) = size(fv(k).vertices,1);                   % Reduced number of vertices
        if Verbose & k==NL
            fprintf('Reducing all %d tessellated layers to %d vertices -> DONE\n',NL,NVertMax);
        end    
    end
end

%if Verbose
%        fprintf('Displaying surface nodes. . .\n');
%        bem_mesh_plot(vertices,faces);
%        drawnow
%end

%
%%%% This part formats parameters for pre-existing routines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NOTE: GEOMETRY AND NODES REFORMATTED TO OUTERMOST->INNERMOST FORMAT %%%%%%%%%%%%%%%%%%%%%%%%%%
%
geometry = zeros(sum(Nf),3);   % Triangle Interconnection Matrix for each layer
% specifying nodes for each triangle
nodes = zeros(sum(Nv),4);      % Nodes for each layer; Col 4 = layer # (Outermost=1!!!)
%
try 
    for k=1:NL  % Reverse order from outermost to innermost !!!
        geometry((sum(Nf(NL-k+2:NL))+1):(sum(Nf(NL-k+1:NL))),1:3) = faces{NL-k+1} + repmat(sum(Nv(NL-k+2:NL)),Nf(NL-k+1),3);
        nodes((sum(Nv(NL-k+2:NL))+1):(sum(Nv(NL-k+1:NL))),1:4) = [vertices{NL-k+1} k*ones(Nv(NL-k+1),1)];
    end
catch
    errordlg('Surfaces have probably too few vertices. Try to increase the number of vertices per tessellated envelope.', mfilename)
    return
end

%
r1 = nodes(geometry(:,1),1:3);   % All Layers Vertex 1 
r2 = nodes(geometry(:,2),1:3);   % All Layers Vertex 2
r3 = nodes(geometry(:,3),1:3);   % All Layers Vertex 3
ctrd = (r1+r2+r3)/3;         % Centroid of each triangle (all layers)
%
% compute the normal vectors
tri_xp=(cross((r2-r1)',(r3-r1)'))';  % Triangle cross product
area2=rownorm(tri_xp);            % twice area of triangle
N= tri_xp./ [area2,area2,area2];  % normalization
%
%%%% FOR THE COLLOCATION CASE ONLY, THIS PART FINDS THE TRIANGLE ON THE OUTERMOST LAYER %%%%%%%%%%%
%%%% WHOSE CENTROID IS CLOSEST TO EACH EEG SENSOR. THE SENSOR LOCATION IS THEN QUANTIZED %%%%%%%%%%
%%%% TO THE LOCATION OF THE CLOSEST TRIANGLE CENTROID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
P = zeros(Meeg,1);
P_sens = zeros(Meeg,3);
P_wts = [];
if Meeg >= 1
    for m=1:Meeg
        [temp indx] = min(rownorm(repmat(R_eeg(m,:),Nf(NL),1)- ctrd(1:Nf(NL),:)));
        P(m) = indx;                        % Index of Tesselation Triangle in "geometry"
        P_sens(m,:) = ctrd(indx,:);
        R_eeg_dist(m,:) = temp;
        R_eeg_quant(m,:) = ctrd(indx,:);
    end
else
    P=[];
end
%
%%%% FOR THE LINEAR CASE ONLY, THE THREE VERTICES ASSOCIATED WITH THE TRIANGLE WHOS CENTROID %%%%%%
%%%% (COMPUTED ABOVE)IS CLOSEST TO THE SPECIFIED POSITION ARE FOUND. DURING KERNEL OPERATIONS, %%%%
%%%% THE FINAL SOLUTION IS INTERPOLATED BETWEEN THE 3 VERTICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if basis_opt==1  % Linear Case
    P = zeros(3*Meeg,1);       % Node index list
    P_sens = zeros(3*Meeg,3);    % Node Position List
    P_wts = zeros(Meeg,3);     % Node Weight List (for sensor position lying interior to triangle)
    if Meeg >= 1
        %
        for m=1:Meeg  % Loop thourgh sensors
            %
            %%%% This part finds closest node to sensor position and those triangles containing node %%%%
            %
            [temp indx] = min(rownorm(repmat(R_eeg(m,:),Nv(NL),1)- nodes(1:Nv(NL),1:3)));
            %
            Pt =[];
            for  i=1:3,                            % find triangles adjacent to node p
                t =find(geometry(:,i)==indx);
                Pt= [Pt; t];  % Append to adjacent triangle list
            end
            %
            %%%% This part finds the projection of sensor point onto nearest triangle surface %%%%
            for i=1:size(Pt,1)  % check each triangle to find the one which contains projected sensor 
                Atemp = [r2(Pt(i,1),:)-r1(Pt(i,1),:); r3(Pt(i,1),:)-r1(Pt(i,1),:)]'; % triangle basis
                btemp = R_eeg(m,:)' - r1(Pt(i,1),:)';    % Sensor point
                xtemp = Atemp\btemp;   % Solve for linear combination
                pnt = r1(Pt(i,1),:) + (Atemp*xtemp)';   % Projection of sensor point onto triangle surface
                %%%% This part computes the weights for each vertices; If point outside triangle %%%%%%%%%%%%%%%%%%
                %%%% nearest vertice node is given full weight %%%%  
                %
                P_wts(m,3) = N(Pt(i,1),:)*crossprod(r2(Pt(i,1),:)-r1(Pt(i,1),:),pnt-r1(Pt(i,1),:))'/area2(Pt(i,1),:); %s
                P_wts(m,2) = N(Pt(i,1),:)*crossprod(pnt-r1(Pt(i,1),:),r3(Pt(i,1),:)-r1(Pt(i,1),:))'/area2(Pt(i,1),:); %t
                P_wts(m,1) = 1 - P_wts(m,2) - P_wts(m,3);
                %
                P(3*m-2:3*m,1) = geometry(Pt(i),1:3)';  % Node List (Triangle Vertices)
                P_sens(3*m-2:3*m,:) = nodes(geometry(Pt(i),1:3)',1:3);
                %
                if (min(P_wts(m,:)')>=0)&(max(P_wts(m,:)')<=1) % Stop if pnt falls interior to triangle
                    break
                end     
                %      
            end  % Triangle List Loop
            %
        end  % Sensor Loop
    else
        P=[];
    end
    %
end 
%
%%%% IF NO MEG SENSORS HAVE BEEN SPECIFIED, THIS PART GENERATES A DUMMY SENSOR TO AVOID ERROR %%%%
%%%% STATEMENTS IN EXISTING CODE (TEMPORARY MEASURE) 
%
if Mmeg >= 1
    S = R_meg;
else
    S = [];
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS PART COMPUTES THE BEM TRANSFER MATRICES BASED ON JAMES CHANG'S ORIGINAL ROUTINE %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
t0 = clock;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% THIS PART FORMATS CONDUCTIVITY PARAMETERS %%%%
%
len = length(sigma);
cdvity = [sigma,0];
cdv =zeros(len,2);

for i=len:-1:1,
    cdv(len+1-i,:)=[cdvity(i),cdvity(i+1)];
end

ss=cdv(:,1)+cdv(:,2);   % conductivity sum no_shell x 1
d=cdv(:,1)-cdv(:,2);   % conductivity difference no_shell x 1
%
%%%% THIS PART COMPUTES THE NORMAL VECTORS FOR EACH TRIANGLE %%%%
%
% change directions to be outward
%
% MXH: this seems to be incrrect, so disable them 05/11/2005
%dp=(dot(ctrd',N'))';
%dex=find(dp < 0);
%N(dex,:)= -N(dex,:);
%temp=r2(dex,:);      %  r1,r2,r3 is CCW if looked outside the head
%r2(dex,:)=r3(dex,:);
%r3(dex,:)=temp;
%
%temp=geometry(dex,2);      %  CCW
%geometry(dex,2)=geometry(dex,3);
%geometry(dex,3)=temp;
no_elt = size(geometry,1);
no_node = size(nodes,1);
%
clear temp dp dex 
%
if basis_opt == 0,
    whichsurf = nodes(geometry(:,1),4);
elseif basis_opt == 1 | basis_opt==2,  
    whichsurf = nodes(:,4);
end  
%
cdvsum=ss(whichsurf);
cdvdiff=d(whichsurf);
no_surf =max(whichsurf);
DIM =zeros(no_surf,1);
for i=1: no_surf
    DIM(i) = length(find( whichsurf ==i));
end
%
%%%% THIS PART COMPUTES THE H AND A MATRICES %%%%
%
H = []; A = [];
if mode ==1, S=[]; end
if basis_opt == 0         % CONSTANT COLLOCATION CASE
    if Verbose
        fprintf('BEM with constant basis, construct system matrix\n') % H : no of triangle x no of triangle
    end
    eegsensor = ctrd(P,1:3);     % EEG SENSOR LOCATION SET TO TRIANGLE CENTROID
    %[H,A] = bem_constant(mode,test_opt,cdvdiff,cdvsum,r1,r2,r3,ctrd,area2/2,N,S);
    if exist(LU_filename)==2 % previous calculated LU exist
        eval(['load ',LU_filename,' pre_bem_basis_opt']);
        if pre_bem_basis_opt==0 % if pre-calculated BEM is constant
            if Verbose, fprintf('LU pre-calculated, loading the matrices...'), end
            eval(['load ',LU_filename,' L_bem U_bem']);  
            if Verbose, fprintf('DONE\n'), end
            if Verbose, fprintf('Only construct system matrix A...'), end
            [H,A] = bem_constant(0,test_opt,cdvdiff,cdvsum,r1,r2,r3,ctrd,area2/2,N,S);
            clear H
            if Verbose, fprintf('DONE\n'), end
        end    
    else    % calculate LU
        if Verbose, fprintf('Construct both system matrices H and A...'), end
        [H,A] = bem_constant(mode,test_opt,cdvdiff,cdvsum,r1,r2,r3,ctrd,area2/2,N,S);
        if Verbose, fprintf('DONE\n'), end
        I = speye(no_elt);            
        H = I- H + 1/no_elt;  
        clear I
        if Verbose, fprintf('Performing LU decomposition...\n'), end
        [L_bem,U_bem] = lu(H);
        L_bem=sparse(L_bem);    
        U_bem=sparse(U_bem);   
        clear H 
        if Verbose, fprintf('DONE\n'), end
        pre_bem_basis_opt=basis_opt;
        if Verbose, fprintf('Saving LU in %s...',LU_filename), end
        eval(['save ',LU_filename,' L_bem U_bem pre_bem_basis_opt -V6']);
        if Verbose, fprintf('DONE\n'), end
    end    
elseif basis_opt == 1    % LINEAR COLLOCATION CASE
    if Verbose
        fprintf('BEM with linear basis, construct system matrix\n')  
    end
    eegsensor = nodes(P,1:3);   % EEG SENSOR LOCATION SET TO CLOSEST NODE
    %[H,A,I] = bem_linear(mode,test_opt,cdvdiff,cdvsum,geometry,...
    %    nodes,r1,r2,r3,N,area2,S);
    if exist(LU_filename)==2 % previous calculated LU exist
        eval(['load ',LU_filename,' pre_bem_basis_opt']);
        if pre_bem_basis_opt==1 % if pre-calculated BEM is linear
            if Verbose, fprintf('LU pre-calculated, loading the matrices...'), end
            eval(['load ',LU_filename,' L_bem U_bem']);  
            if Verbose, fprintf('DONE\n'), end
            if Verbose, fprintf('Only construct system matrix A...'), end
            [H,A,I] = bem_linear(0,test_opt,cdvdiff,cdvsum,geometry,...
                nodes,r1,r2,r3,N,area2,S);
            clear H I
            if Verbose, fprintf('DONE\n'), end
        end    
    else    % calculate LU
        if Verbose, fprintf('Construct both system matrices H and A...'), end
        [H,A,I] = bem_linear(mode,test_opt,cdvdiff,cdvsum,geometry,...
            nodes,r1,r2,r3,N,area2,S);
        if Verbose, fprintf('DONE\n'), end
        H = I- H + 1/no_elt;  
        clear I
        if Verbose, fprintf('Performing LU decomposition...'), end
        [L_bem,U_bem] = lu(H);
        L_bem=sparse(L_bem);    
        U_bem=sparse(U_bem);   
        clear H 
        if Verbose, fprintf('DONE\n'), end
        pre_bem_basis_opt=basis_opt;
        if Verbose, fprintf('Saving LU in %s...',LU_filename), end
        eval(['save ',LU_filename,' L_bem U_bem pre_bem_basis_opt -V6']);
        if Verbose, fprintf('DONE\n'), end
    end   
end
%
%%%% This part does other functions 
%
PP =  sparse(length(P),sum(DIM)); % a sparse zero matrix
for i=1:length(P);
    PP(i,P(i))=1;
end
%
if mode ==1, S=[]; end;
if Verbose, fprintf('Calculate transfer matrix...'), end
[Te_ISA,Te,Tm_ISA,Tm]= trans_matrix(mode,test_opt,ISA,DIM,cdv(no_surf,:),L_bem,U_bem,PP,A,size(S,1));

if Verbose, fprintf('DONE\n'), end

%
Te_ISA = full(Te_ISA);
Te = full(Te);
Tm_ISA = full(Tm_ISA);
Tm = full(Tm);
%
basisname = ['constant BEM';...
        'linear   BEM'];
testname = ['collocation  ';...
        'galerkin     '];
basis = basisname(basis_opt+1,:);
test = testname(test_opt+1,:);
%
telap_xfer = etime(clock,t0);
if Verbose
    fprintf('Computed in %3.1f sec.\n',telap_xfer);
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THIS PART STORES FINAL TRANSFER MATRIX RESULTS AND PARAMETERS IN A ".MAT" FILE %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       
eegsensor = [P(:),P_sens];
%if (mode==1)|(mode==3)
%    modality = 1;
%    try
%        eval(['save ', fn_eeg,[' basis_opt test_opt basis test geometry nodes cdv' ...
%                    ' Te Te_ISA eegsensor modality R_eeg P_sens P_wts']]);
%    catch
%        Users = get_user_directory;
%        cd(Users.STUDIES)
%        eval(['save ', fn_eeg,[' basis_opt test_opt basis test geometry nodes cdv' ...
%                    ' Te Te_ISA eegsensor modality R_eeg P_sens P_wts']]);
%    end
%    
%    if Verbose
%        fprintf('EEG Transfer Matrix data stored in file %s.\n',fn_eeg);
%    end
%    
%    %
%end
%
%if  (mode==2)|(mode==3)
%    modality = 2;
%    try 
%        eval(['save ', fn_meg,[' basis_opt test_opt basis test geometry nodes cdv' ...
%                    ' Tm Tm_ISA S modality R_meg O_meg']])
%    catch
%        Users = get_user_directory;
%        cd(Users.STUDIES)
%        eval(['save ', fn_meg,[' basis_opt test_opt basis test geometry nodes cdv' ...
%                    ' Tm Tm_ISA S modality R_meg O_meg']])
%    end
%    
%    if Verbose
%        fprintf('MEG Transfer Matrix data stored in file %s.\n',fn_meg);
%    end
%    
%    %
%end
%
%if Verbose
%    fprintf('Computing BEM Transfer Matrices -> DONE\n');
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct return variables M. Huang
% 
%    return_var: structure variable
%           return_var.basis_opt: string --- constant bem or linear bem
%           return_var.test_opt: string -- collocation, glerkin
%           return_var.mode: mode type
%           return_var.geometry: geometry matrix
%           return_var.nodes: nodes matrix
%           return_var.cdv: conductivity vector
%           return_var.ISA: 1 for ISA, 0 for no ISA
%           return_var.Te: the EEG transfermation matrix
%           return_var.Tm: the MEG transformation matrix
%           return_var.Te_ISA: the EEG transfermation matrix with ISA
%           return_var.Tm_ISA: the MEG transformation matrix with ISA
%           return_var.eegsensor: eegsensor
%           return_var.R_eeg: R_eeg
%           return_var.P_sens: P_sens
%           return_var.P_wts: P_wts weights for triangle points, linear
%           return_var.R_meg: the MEG sensor locations, no_sensor by 3
%           return_var.O_meg: the MEG sensor orientations, no_sensor by 3

return_var.basis_opt=basis;
return_var.test_opt=test;
return_var.mode=mode;
return_var.geometry=geometry;
return_var.nodes=nodes;
return_var.cdv=cdv;
return_var.ISA=ISA;
return_var.Te=Te;
return_var.Tm=Tm;
return_var.Te_ISA=Te_ISA;
return_var.Tm_ISA=Tm_ISA;
return_var.eegsensor=eegsensor;
return_var.R_eeg=R_eeg;
return_var.P_sens=P_sens;
return_var.P_wts=P_wts;
return_var.R_meg=R_meg;
return_var.O_meg=O_meg;


% ----------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% SUBFUNCTIONS


% ----------------------------------------------------------------------------------------------------------------


function [H,A,G] = bem_linear(mode,choice,cdvdiff,cdvsum,geometry,...
    nodes,r1,r2,r3,N,area2,S)
%bem_linear (overwrite succinct one line summary here)
% function [H,A,G] = bem_linear(mode,choice,cdvdiff,cdvsum,geometry,...
%                             nodes,r1,r2,r3,N,area2,S)
%
% mode: 1 for EEG only, 2 for MEG only, 3 for both, 0 for A only
%           (precalculated LU)
% r1(i,:), r2(i,:), r3(i,:) form three vertices of triangle i and ther are
% CCW as seen from outside			
% References:
% (1) de Munck,IEEE Trans. BME, pp 986-990,1992
% (2) Schlitt,et al, IEEE Trans. BME, pp 52-58,1995
%
% %%% (03/12/00)-Added Case where MEG is not computed to prevent 
%                Matlab error from being flagged (John Ermer)
%

%<autobegin> -------- 20-Nov-2002 14:04:25 ------------------------------
% ---- Automatically Generated Comments Block using auto_comments -----------
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\crossprod.m
%   toolbox\dotprod.m
%   toolbox\hrow_linear.m
%   toolbox\nit_gq.m
%   toolbox\rownorm.m
%<autoend> ---------- 20-Nov-2002 14:04:25 ------------------------------

no_node = size(nodes,1);
no_elt=size(geometry,1);
no_eltx2 =2*no_elt;

PP =  sparse(no_node,3*no_elt);
for i=1:no_elt,
    a=  1/area2(i);
    PP(geometry(i,1),i)=a;
    PP(geometry(i,2),i+no_elt)=a;
    PP(geometry(i,3),i+no_eltx2)=a;
end
% work on 'geometry' matrix.
% three vertices of the ith row contribute to the ith triangle i.       


r21 = r2-r1;
r32 = r3-r2;
r13 = r1-r3;
yset = [r21,r32,r13];
yset_mag = [rownorm(r21), rownorm(r32),rownorm(r13)];

dx = [ 2 3 1 2];
dx1= [1 2 3 1];

if mode == 0, % A only no need for H
    H = [];
    G = [];
else    
    H = zeros(no_node);
    %disp('Computing H .... ')
    
    if choice == 0, % collocation
        G = speye(no_node);
        %  t=  clock;
        for  i=1:no_node
            v =1:no_elt;
            x =find(geometry(:,1)==i | geometry(:,2)==i | geometry(:,3)==i);
            % excluding all triangles having node i as a vertex as the integral
            % across these trianglesa are either zeros or singular
            v(x)=[];
            H(i,:) = hrow_linear(nodes(i,1:3),PP,v,yset,yset_mag,geometry,nodes,cdvdiff,N,area2)/(2*pi*cdvsum(i));
            
        end
        
        for i=1:no_node,    % diagonal terms
            H(i,i)= 1-sum(H(i,:));
        end
    else 
        G = sparse(no_node,no_node);
        rnxrn1 = [crossprod(r2,r3),crossprod(r3,r1),crossprod(r1,r2)];
        DET = dotprod(r1,rnxrn1(:,1:3));
        for p=1:no_node, 
            %    tt=  clock;
            %sprintf('node# = %5.0f',p)
            %       profile on -detail builtin
            ptri =[];
            for  i=1:3,   % find triangles adjacent to node p
                t =find(geometry(:,i)==p);
                ptri= [ptri; t,i*ones(size(t))];
            end
            
            for i=1:size(ptri,1)  % compute the Gram matrix  
                t= ptri(i,1);   
                G(p,geometry(t,:)) = G(p,geometry(t,:))+ nit_gq('kl_gram',4,3,r1(t,:),r2(t,:),r3(t,:),rnxrn1(t,:),ptri(i,2))/DET(t)^2;
            end       
            
            for i=1:size(ptri,1),
                t = ptri(i,1);
                v = 1:no_elt;
                v(t)=[];
                H(p,:) = H(p,:)+nit_gq('kl_galerkin',3,no_node,r1(t,:),r2(t,:),...
                    r3(t,:),PP,ptri(i,2),rnxrn1(t,:),v,yset,yset_mag,...
                    geometry,nodes,cdvdiff,N,area2)/DET(t);
            end
            
            H(p,:) = H(p,:)/(2*pi*cdvsum(p));    % cdvdiff is treated in hrow_linea
            %profile report
        end 
    end    % end if  
end
%
% ----    compute the A matrix for MEG exactly   -----------
% ----    see Ferguson's paper for details ------------------
%
if mode~=1,
    %disp('Computing the matrix A for MEG (omitting u0/4pi) .... ')
    % t=clock;
    if length(S)==3, 
        no_megsensor=1;     % only one meg sensor
        S=S(:)';            % into a row vector
    else  no_megsensor=size(S,1);    
    end
    
    A=zeros(3*no_megsensor,no_node);
    cols = [0,no_elt,no_elt*2,no_elt*3];
    for i=1:no_megsensor,
        y1 = [r1(:,1)-S(i,1),r1(:,2)-S(i,2),r1(:,3)-S(i,3)];
        y2 = [r2(:,1)-S(i,1),r2(:,2)-S(i,2),r2(:,3)-S(i,3)];
        y3 = [r3(:,1)-S(i,1),r3(:,2)-S(i,2),r3(:,3)-S(i,3)];
        mag1=rownorm(y1);
        mag2=rownorm(y2);
        mag3=rownorm(y3);
        sangle = 2*atan2(dotprod(y1,crossprod(y2,y3)),mag1.*mag2.*mag3+...
            mag1.*dotprod(y2,y3) + mag2.*dotprod(y1,y3) +...
            mag3.*dotprod(y1,y2) );
        
        c = dotprod(y1,N);
        c = N.*[c,c,c];
        a = [S(i,1)+c(:,1),S(i,2)+c(:,2),S(i,3)+c(:,3)];
        c1 = r1-a;
        c2 = r2-a;  
        c3 = r3-a;
        
        tmp = zeros(no_elt,3);
        tmpsum =zeros(no_elt,1);
        
        for p=1:3,
            i1 =dx1(p+1);
            i2 = p;
            yp1 = eval(['y',num2str(i1)]);
            yp = eval(['y',num2str(i2)]);
            cp1 = eval(['c',num2str(i1)]);
            cp = eval(['c',num2str(i2)]);
            
            yp1_yp = yset(:,3*p-2:3*p);
            mag_yp1_yp = yset_mag(:,p);
            
            w = (mag_yp1_yp.*rownorm(yp1) + dotprod(yp1,yp1_yp))./...
                (mag_yp1_yp.*rownorm(yp) + dotprod(yp,yp1_yp));
            gammap =log(w)./mag_yp1_yp;
            tmp =tmp+crossprod(cp,cp1).*[gammap,gammap,gammap];
        end
        % twice area is absorbed into PP
        tmpsum =( dotprod(N,tmp)-dotprod(N,c).*sangle);
        rows = 3*i-2:3*i;
        for n=1:3,
            col= geometry(:,n);
            yj_yk = -yset(:,3*dx(n)-2:3*dx(n));
            
            a = cdvdiff(col).*tmpsum;
            a = (yj_yk .*[a,a,a]);
            
            A(rows,:)= A(rows,:)+(PP(:,cols(n)+1:cols(n+1))*a)';
        end
        
    end
else
    A = [];
end

% ----------------------------------------------------------------------------------------------------------------

function [H,A] = bem_constant(mode,choice,cdvdiff,cdvsum,r1,r2, ...
    r3,ctrd,area,N,S)
%bem_constant Compute the geometry matrix H for EEG & MEG by centroid approxiamtion
% function [H,A] = bem_constant(mode,choice,cdvdiff,cdvsum,r1,r2, ...
%                  r3,ctrd,area,N,S)
% - compute the geometry matrix H for EEG & MEG by centroid approximation of numerical integration   
% - compute the matrix A analytic
% mode: 1 for EEG only, 2 for MEG only, 3 for both, 0 for A only
%           (precalculated LU)
% See:
% Hamalainen, Sarvas, IEEE Trans. BME, pp.165-171, 1989
% de Munck,IEEE Trans.BME,pp986-990,1992 (for integration of the MEG kernel)
%
% (03/12/00)-Added Case where MEG is not computed to prevent 
%                Matlab error from being flagged (John Ermer)

%<autobegin> -------- 20-Nov-2002 14:04:16 ------------------------------
% ---- Automatically Generated Comments Block using auto_comments -----------
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\dotprod.m
%   toolbox\nit_gq.m
%   toolbox\rownorm.m
%   toolbox\solid_angle2.m
%<autoend> ---------- 20-Nov-2002 14:04:16 ------------------------------


%TH = max(rownorm(r1)); % threshold used in Galerkin method, see below
%tic
no_elt = size(r1,1);

if mode == 0, % A only no need for H
    H = [];
else    
    %%disp(' Computing  H ......')
    H=zeros(no_elt);
    
    for i=1:no_elt
        %sprintf('Node #%4f...',i)
        x=1:no_elt;
        x(i)=[];     % excluding the diagonal terms
        sangle = zeros(no_elt-1,1);
        if choice ==0,    % ctrd approxmation
            sangle = solid_angle2(ctrd(i,:),r1(x,:),r2(x,:),r3(x,:)); 
            
        else  % Galerkin
            sangle =nit_gq('solid_angle2',3,1,r1(i,:),r2(i,:),r3(i,:),r1(x,:),r2(x,:),r3(x,:))/area(i);
            
        end      % end if
        H(i,x) = ( cdvdiff(x).* sangle )'/( 2*pi*cdvsum(i)) ; 
    end;
end    
    %telap_h1 = toc
    %
    %  compute the A matrix for MEG exactly 
    %  see de Munck's paper for details
    %
if mode~=1 , % MEG and Fusion case
    %disp('Computing the matrix A for MEG(omitting u0/4pi) .... ')
% t = clock;
    if length(S)==3, 
        no_megsensor=1;     % only one meg sensor
        S=S(:)';            % into a row vector
    else  no_megsensor=size(S,1);    
    end
    
    A=zeros(3*no_megsensor,no_elt);
    r1_r3= r1-r3;
    r3_r2= r3-r2;
    r2_r1= r2-r1;
    mag_r2_r1=rownorm(r2_r1);
    mag_r3_r2=rownorm(r3_r2);
    mag_r1_r3=rownorm(r1_r3);
    
    for i=1:no_megsensor,
        S_r1=[S(i,1)-r1(:,1),S(i,2)-r1(:,2),S(i,3)-r1(:,3)];
        S_r2=[S(i,1)-r2(:,1),S(i,2)-r2(:,2),S(i,3)-r2(:,3)];
        S_r3=[S(i,1)-r3(:,1),S(i,2)-r3(:,2),S(i,3)-r3(:,3)];
        mag_S_r1=rownorm(S_r1);
        mag_S_r2=rownorm(S_r2);
        mag_S_r3=rownorm(S_r3);
        
        gamma1 = -log(( mag_r2_r1.*mag_S_r2 - dotprod(r2_r1,S_r2)) ...
            ./( mag_r2_r1.*mag_S_r1-dotprod(r2_r1,S_r1) )) ./mag_r2_r1;
        
        gamma2 = -log(( mag_r3_r2.*mag_S_r3 - dotprod(r3_r2,S_r3)) ...
            ./( mag_r3_r2.*mag_S_r2-dotprod(r3_r2,S_r2) )) ./mag_r3_r2;
        
        gamma3 = -log(( mag_r1_r3.*mag_S_r1 - dotprod(r1_r3,S_r1)) ...
            ./( mag_r1_r3.*mag_S_r3-dotprod(r1_r3,S_r3) )) ./mag_r1_r3;
        
        % omit u0/4pi
        A(3*i-2:3*i,:)= (([gamma1,gamma1,gamma1].*r2_r1+...
            [gamma2,gamma2,gamma2].*r3_r2+...
            [gamma3,gamma3,gamma3].*r1_r3) .*[cdvdiff,cdvdiff,cdvdiff])';
    end
else % EEG case
    A = [];
end
%keyboard

% ----------------------------------------------------------------------------------------------------------------

function I = nit_gq(F,N,Idim,r1,r2,r3,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12)
%NIT_GQ (F,N,Idim,r1,r2,r3,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12)
% function I = nit_gq(F,N,Idim,r1,r2,r3,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12)
% This function integrates a scalar/vector(dimension indicated by Idim)
% function over triangles with vertices r1,r2,and r3, using
% the Gauss Quadrature,i.e.,
%    /
%    | F(r,P1,P2,P3,....)da    , 
%    /
% where  r=(x,y,z) is a point on the triangle, and P1,P2,.... are parameters
% defined in F.
%
% Input:
% F     :  A string containing the function to be integrated. It should take r,P1,P2,P3,...
%          as input, where r is no_triangles by 3,and produces no_triangles by
%          Idim funtion valus.
% N     :  use N point Gauss Quadrature          scalar
% r1, r2,r3 :  specify the vertices of triangles over which F is integrated;
%              each is no_triangles x 3 
% P1,P2..  : parameters defined in F. 
%
% Output:
% I : integral                            no_triangles(=size(r1,1)) x Idim
%

%<autobegin> -------- 20-Nov-2002 14:07:22 ------------------------------
% ---- Automatically Generated Comments Block using auto_comments -----------
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\crossprod.m
%   toolbox\rownorm.m
%
%<autoend> ---------- 20-Nov-2002 14:07:22 ------------------------------


% Author: CCH
%


if N < 2 | N > 10,
    error(' The GQ table of your choice is not built. Sorry');
end

%-----------------------------------------------------------
% W: weight table, X : abscissa table
X2 = 0.5773503;
W2 = 1.0;

X3 = [0.0,       0.7745967];
W3 = [0.4444445, 0.5555556];

X4 = [0.3399810, 0.8611363];
W4 = [0.6521452, 0.3478548];

X5 = [0.0,       0.5384693, 0.9061798];
W5 = [0.2844445, 0.4786287, 0.2369269];

X6 = [0.2386192, 0.6612094, 0.9324695];
W6 = [0.4679139, 0.3607616, 0.1713245];

X7 = [0.0,       0.4058452, 0.7415312, 0.9491079];
W7 = [0.2089796, 0.3818301, 0.2797054, 0.1294850];

X8 = [0.1834346, 0.5255324, 0.7966665, 0.9602899];
W8 = [0.3626838, 0.3137066, 0.2223810, 0.1012285];

X9 = [0.0        0.3242534  0.6133714  0.8360311 0.9681602];
W9 = [0.3302394  0.3123471  0.2606107  0.1806482 0.0812744];

X10 = [0.1488743 0.4333954 0.6794096 0.8650634 0.9739065];    
W10 = [0.2955242 0.2692667 0.2190864 0.1494513 0.0666713];
%-----------------------------------------------------------
% retrieve weights and abscissa
W = eval(['W' int2str(N)]);   % weight
X = eval(['X' int2str(N)]);   % abscissa

r2_r1=r2-r1;
r3_r1=r3-r1;
halfarea = rownorm(crossprod(r2_r1,r3_r1))/4;  % size(r1,1) x 1
evalstr = [F,'(r'];
for i=1:nargin - 6,
    evalstr = [evalstr,',P',int2str(i)];
end
evalstr = [evalstr,')'];

LEN = size(W,2);
I = zeros(size(r1,1),Idim);  % Idim is the dimension of integral
half_r2_r1= r2_r1/2;
for i=1:LEN,
    t1 = 0.5*(1-X(i));    % scalar
    t2 = 0.5*(1+X(i));
    ct2 = r1+r3_r1*t2;   % common terms
    ct1 = r1+r3_r1*t1;
    s1 =zeros(size(I));
    s2 =zeros(size(I));
    for j=1:LEN,
        d = half_r2_r1*t1;
        r = ct2+ d*(1+X(j));
        Ta = eval(evalstr);
        r = ct2+ d*(1-X(j));
        Tb = eval(evalstr);
        
        d= half_r2_r1*t2;
        r = ct1+ d*(1+X(j));
        Tc = eval(evalstr);
        r = ct1+ d*(1-X(j));
        Td = eval(evalstr);
        
        s1 = s1+W(j)*(Ta+Tb);
        s2 = s2+W(j)*(Tc+Td);
        
    end
    I = I + W(i)*(t1*s1 +t2*s2);
end

I = I.*(halfarea*ones(1,Idim));
%disp('in nit_gq')
%keyboard

 
function nfv = make_surfaces_coincide(fv,Verbose)

% Align all surfaces in fv to the innermost surface in array patch structure fv.
% [Align in the sense that resulting vertices are aligned (i.e. coincident)]

% Centroid of innermost surface
head_center = mean(fv(1).vertices);
% Initialize new structure;
for k=1:length(fv)
    nfv(k) = fv(1);
end
    

for isurf = 2:length(fv) % for each outermost surface
    if Verbose
        makeuswait('start')
        hw = waitbar(0,sprintf('Aligning BEM envelope #%d on innermost envelope #1',isurf));
    end
    for ivert = 1:size(fv(isurf).vertices,1) % project each vertex onto inner surface
        [intersect,indx,t,u,v] = ray_intersect(fv(isurf),fv(1).vertices(ivert,:),-fv(1).vertices(ivert,:)+head_center,'b');
        if isempty(t)
            errordlg(sprintf('Surfaces #%d and #1 are probably not registered into the same coordinate system. Please check registration using BrainStorm''s alignement tool.',isurf)), 
            return, 
        end
        [tmp, imin] = min(abs(t));
        t = abs(t(imin));
                
        nfv(isurf).vertices(ivert,:) = intersect(:,imin)';
        
        if Verbose & ~rem(ivert,50)
            waitbar(ivert/size(fv(isurf).vertices,1),hw)
        end
    end

    %nfv(isurf).faces = fv(1).faces;
    if Verbose
        makeuswait('stop')
        close(hw)
    end

end


% -----------------------------------------------------------------

function fv = make_surfaces_embedded(fv,tol,Verbose)

% Check that innermost points are within outermost envelopes
for isurf = 1:length(fv)-1 % for each innermost surface
    iter = 0; 
    IS = inpolyhd(fv(isurf).vertices,fv(isurf+1).vertices,fv(isurf+1).faces);
    iWrong = find(IS == 0)'; % index of vertices frominner surface that are either oustide the outersurface (iWrong = 1) 
    % or inside but within tol (meters) of outer surface (iWrong  = .5)
    head_center = mean([fv(isurf).vertices]);
    
    while ~isempty(iWrong) & iter <100 % some supposedly inner points are outside next outermost envelope
        iter = iter+1;   % while on iter avoid infinte loops
        iWrong_out = find(IS == 0)'; % Vertices outside surface 
        for k = iWrong_out
            [intersect,indx,t,u,v] = ray_intersect(fv(isurf+1),fv(isurf).vertices(k,:),head_center-fv(isurf).vertices(k,:),'b');
            [tmp, imin] = min(abs(t));
            t = t(imin);
            intersect = intersect(:,imin)';
            dist = sign(t)*(fv(isurf).vertices(k,:) - intersect);
            dist = dist/abs(t);
            fv(isurf).vertices(k,:) = intersect - tol*dist;
        end
        
        IS = inpolyhd(fv(isurf).vertices,fv(isurf+1).vertices,fv(isurf+1).faces);
        iWrong = find(IS==0)';
    end
end

% Now test that every inner vertices are at least away of outer surface by tol distance 

% Check that innermost points are within outermost envelopes
for isurf = 1:length(fv)-1 % for each innermost surface

    head_center = mean([fv(isurf).vertices]);
    if Verbose
        makeuswait('start')
        hw = waitbar(0,sprintf('Checking for minimum separation between surfaces #%d and #%d',isurf, isurf+1));
    end

    for k = 1:size(fv(isurf).vertices,1) % For each vertex
        [intersect,indx,t,u,v] = ray_intersect(fv(isurf+1),fv(isurf).vertices(k,:),head_center-fv(isurf).vertices(k,:),'b');
        % Find closest intersection
        [tmp, imin] = min(abs(t));
        t = t(imin);
        if abs(t) < tol % Vertex is too close: move away from outer surface
            intersect = intersect(:,imin)';
            dist = sign(t)*(fv(isurf).vertices(k,:) - intersect);
            if t == 0
                t = tol;
                dist = sign(t)*(fv(isurf).vertices(k,:)+eps - intersect);
            end
            dist = dist/abs(t);
            fv(isurf).vertices(k,:) = intersect - tol*dist;
        end
        if Verbose & ~rem(k,50)
            waitbar(k/size(fv(isurf).vertices,1),hw)
        end
    end
end
if Verbose
    makeuswait('stop')
    close(hw)
end
        

% Now that outer surfaces may have moved
% Check that innermost points of surf#1 are still away enough from surf#2
for isurf = 1 % for each innermost surface

    head_center = mean([fv(isurf).vertices]);

    for k = 1:size(fv(isurf).vertices,1) % For each vertex
        [intersect,indx,t,u,v] = ray_intersect(fv(isurf+1),fv(isurf).vertices(k,:),head_center-fv(isurf).vertices(k,:),'b');
        % Find closest intersection
        [tmp, imin] = min(abs(t));
        t = t(imin);
        if abs(t) < tol % Vertex is too close: move away from outer surface
            intersect = intersect(:,imin)';
            dist = sign(t)*(fv(isurf).vertices(k,:) - intersect);
            if t == 0
                t = tol;
                dist = sign(t)*(fv(isurf).vertices(k,:)+eps - intersect);
            end
            dist = dist/abs(t);
            fv(isurf).vertices(k,:) = intersect - tol*dist;
        end
    end
end


% Final check
% Check that innermost points are within outermost envelopes
for isurf = 1:length(fv)-1 % for each innermost surface
    iter = 0; 
    IS = inpolyhd(fv(isurf).vertices,fv(isurf+1).vertices,fv(isurf+1).faces);
    iWrong = find(IS == 0)'; % index of vertices frominner surface that are either oustide the outersurface (iWrong = 1) 
    % or inside but within tol (meters) of outer surface (iWrong  = .5)
    head_center = mean([fv(isurf).vertices]);
    
    while ~isempty(iWrong) & iter <100 % some supposedly inner points are outside next outermost envelope
        iter = iter+1;   % while on iter avoid infinte loops
        iWrong_out = find(IS == 0)'; % Vertices outside surface 
        for k = iWrong_out
            [intersect,indx,t,u,v] = ray_intersect(fv(isurf+1),fv(isurf).vertices(k,:),-head_center+fv(isurf).vertices(k,:),'o');
            intersect = intersect(:,1)';
            dist = (fv(isurf).vertices(k,:) - intersect);
            dist = dist/abs(t(1));
            fv(isurf).vertices(k,:) = intersect - tol*dist;
        end
        
        IS = inpolyhd(fv(isurf).vertices,fv(isurf+1).vertices,fv(isurf+1).faces);
        iWrong = find(IS<1)';
    end
end

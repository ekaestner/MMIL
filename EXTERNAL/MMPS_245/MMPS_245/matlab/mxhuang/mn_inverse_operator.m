function [dSPM_output,J_output,B_predict_output,factor_J]=mn_inverse_operator(data_input,G_orig);
% function [dSPM_output,J_output,B_predict_output,factor_J]=mn_inverse_operator(data_input,G_orig);
% calculate the dSPM using cortical constrained minimum norm solution
% using: Dale et al Neuron, Vol. 26, 55-67, April, 2000,
% Also using MNE manual by M. Hamalainen Chap 4
%
% (c) M. Huang, Ph.D. 05/19/05
% 
% data_input.meg:   1 x n_cond, cell array, each cell is an
%                   m x T MEG data matrix, T/m for gradiometer and T for
%                   magnetometer.
% data_input.t:     1 x T time index
% data_input.ncov:  m x m the raw data covariance matrix
% data_input.navg:  1 x n_cond the number of trials in averaging for each condition 
% data_input.SNR:   1 x 1 signal to noise ratio
% data_input.R:     3p x 3p, source covariance matrix, usually is SPARSE
%                   identity matrix
% G_orig:           m x 3p MEG gain matrix, T/m for gradiometer, T for
%                   magnetometer for A*m dipole moments
%
% dSPM_output:      1 x n_cond, cell array, each cell is a
%                   p x T MEG dSPM F-score matrix
% J_output:         1 x n_cond, cell array, each cell is a
%                   3p x T MEG dipole moment matrix unit nA*m
% B_predict_output: 1 x n_cond, cell array, each cell is a
%                   m x T predicted MEG fields matrix
% factor_J:         3p x m dSPM inverse operator
%
% m is the number of MEG channel
% p is the number of dipole nodes
% T is the number of time points
% n_cond is the number of conditions in the MEG file
%
%------------------------------
% change 01/25/06
% 
% 1) remove baseline, it should be done before calling this program 
% 2) add rank check for Cinv_sqrt_raw
% 3) delect the bad channel remove, should be done before calling


% unpack the data
n_cond=length(data_input.meg);
dSPM_output=cell(1,n_cond);
J_output=cell(1,n_cond);
B_predict_output=cell(1,n_cond);
t=data_input.t;
nt=length(t);
%t_base=data_input.baseline;
%[dum1,id_bmin]=min(abs(t-t_base(1)));
%[dum1,id_bmax]=min(abs(t-t_base(2)));
Craw=data_input.ncov;
SNR=data_input.SNR;
Rraw=data_input.R;
%id_bad=data_input.bad;

%remove the bad channels
G=G_orig;
%G(id_bad,:)=[];              
%Craw(id_bad,:)=[];
%Craw(:,id_bad)=[];

% now calculate the inverse operator
[m,n]=size(G);
[Uc,gammac]=eig(Craw);
rank_Craw=rank(Craw);
sqrt_inv_gammac=zeros(size(gammac));
for i=1:m
    if (m-i+1)<=rank_Craw,
        sqrt_inv_gammac(i,i)=sqrt(1.0/gammac(i,i));
    end
end    
Cinv_sqrt_raw=sqrt_inv_gammac*Uc';

lamda_sq_SNR=1.0/SNR^2;     % SNR for power is needed here
G=Cinv_sqrt_raw*G;              % prewhitening the data with raw covariance matrix 
Rraw=m*Rraw/trace(G*Rraw*G');  % make sure trace(GRG')/trace(I)=1
Rraw_c=chol(Rraw);     
A=G*Rraw_c;                 % A should be independent of the number of trials
[V,lamda,U]=svd(A',0);      % economy SVD proventing memory overflow
lamda=diag(lamda);
gamma=diag(lamda./(lamda.^2+lamda_sq_SNR));
gamma=sparse(gamma);
Vbar_raw=Rraw_c*V;
factor_J=Vbar_raw*gamma*U'*Cinv_sqrt_raw;
Vbar_gamma_raw=Vbar_raw*gamma;
Vbar_gamma_raw_T=Vbar_gamma_raw';
factor_w=zeros(n,1);
for k=1:n
    factor_w(k)=sqrt(Vbar_gamma_raw(k,:)*Vbar_gamma_raw_T(:,k)); % noise normalization factor Dale et al Eq(5)
end    
clear Vbar_raw A G Uc gammac Rraw Rraw_c Craw U V Vbar_gamma_raw Vbar_gamma_raw_T

% processing minimum norm inverse for every conditions
for i=1:n_cond,
    B=data_input.meg{i};
    %B(id_bad,:)=[];             % remove the bas channels
    %B_base=mean(B(:,id_bmin:id_bmax)')';
    %B=B-B_base*ones(1,nt);      % baseline correction
    L=data_input.navg(i);
    J=factor_J*B;       % independent of the L
    w_normal=factor_w/sqrt(L); % noise normalization 
    J_amp=zeros(n/3,nt);
    F_amp=zeros(n/3,nt);
    for j=1:n/3
        id_temp=[-2:0]+j*3;
        J_amp(j,:)=sqrt(J(id_temp(1),:).^2+J(id_temp(2),:).^2+J(id_temp(3),:).^2); % sum over xyz
        F_amp(j,:)=(J_amp(j,:).^2)/sum(w_normal(id_temp).^2);    
    end    
    dSPM_output(i)={F_amp};
    J_output(i)={J};
    B_predict_output(i)={G_orig*J};
end    


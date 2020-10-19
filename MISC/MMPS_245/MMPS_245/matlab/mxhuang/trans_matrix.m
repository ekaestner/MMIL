function [Te_ISA,Te,Tm_ISA,Tm] = trans_matrix(mode,test_opt,ISA,DIM,s,L,U,P,A,no_megsensor)
%TRANS_MATRIX_V2 - calculate BEM transfer matrices for EEG and MEG
% function [Te_ISA,Te,Tm_ISA,Tm] = trans_matrix(mode,test_opt,ISA,DIM,s,L,U,P,A,no_megsensor)
%**************************************************************************
%
% this function calculates transfer matrices for EEG and MEG
% mode        =1 , compute Te
%              2 , compute Tm
%              3 , compute BOTH
% ISA         see bem.m
% DIM         : DIM(i) is the number of triangles(or nodes) on surface i
%             : 3 x 1
% s           : [conductivity of brain, conductivity of skull]
%             : 2 x 1
% P           : EEG scan matrix       no_eegsensor x sum(DIM)
% L,U         : LU decomposition of EEG system matrix H
% A           : MEG system matrix
% no_megsensor: number of MEG sensors
%
% The isolated skull approach is incorporated for THREE shells model.
%*****************************************************************************
% last Modified by M.X. Huang, March 2006
% other people invloved in old version: J. Chang, J. Mosher, E. Ermer.


s3=s(1);
s2=s(2);
no_elt=sum(DIM); 

%disp(' LU decomposition .......')
% Incoporate the isolated skull approach 
% C3: solid angle matrix of brain layer

if ISA, 
   M=DIM(1)+DIM(2);
   b =  M+1:no_elt;  
   I3 = I(b,b);
   C3 =I3 - H(b,b)*(s3+s2)/(s3-s2);

   for i=1:size(C3,1)
     C3(i,i)=0;
     C3(i,i)= -sum(C3(i,:));
   end
   C3 = C3+1/DIM(3);   % deflation 

   [L3,U3]=lu(C3);
   L3=sparse(L3);  U3=sparse(U3);     
   clear C3 b
end   

Te_ISA = [];   %  transfer matrix for EEG, with isolated skull approach(ISA)
Tm_ISA = [];   %  transfer matrix for MEG, with ISA 
Te = [];   % without isa
Tm = [];

if mode==1 | mode==3,
   %disp(' ')
   %disp('>>>> Compute Te ...')
%   t=clock;
   Te = P/U/L; 
   if ISA   % isolated skull approach
      P3 = P(:,M+1:no_elt);    
      Te_ISA = s2/s3*[Te(:,1:M),Te(:,M+1:no_elt)-2*(Te(:,M+1:no_elt)*I3)/U3/L3+(1+s3/s2)*P3/U3/L3];
    end
%   run_Te_in_min=etime(clock,t)/60
%keyboard
   clear P P3
end; 

if mode==2 | mode==3,
%   disp(' ')
%   disp('>>>> Compute Tm .....')
%   t=clock;
   %Tm = A/U/L;
   Tm = A*inv(U)*inv(L); % reduce memary cost
   if ISA
     A3 = A(:,M+1:no_elt);
     clear U L A
     Tm_ISA = s2/s3*[Tm(:,1:M),Tm(:,M+1:no_elt)-2*(Tm(:,M+1:no_elt)*I3)/U3/L3+(1+s3/s2)*A3/U3/L3];
   end
%   run_Tm_in_min=etime(clock,t)/60
end;
%keyboard









clc;
clear all;
format long;
I1=imread('image11.jpg');
I2=double(I1);
A=I2;
clear vision;
k=1;
for i=1:10:(960-96);
for j=1:10:(1280-128);
m1=max(max(A(i:i+95,j:j+127)));
m2=min(min(A(i:i+95,j:j+127)));
vision(k)=(m1-m2)/(m1);
k=k+1
end;
end;

visionk=mean(vision) 
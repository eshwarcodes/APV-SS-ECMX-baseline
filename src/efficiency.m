function [KGE,R2]=efficiency(S,O)
% delete efficiency.m~;
r=min(min(corr(S,O)));
a=std(S)/std(O);
b=mean(S)/mean(O);
R2=r^2;
KGE=1-sqrt((r-1)^2+(a-1)^2+(b-1)^2);

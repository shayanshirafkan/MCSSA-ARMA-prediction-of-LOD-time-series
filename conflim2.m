function [c,c2]=conflim(L,p) 
% CONFLIM - returns empirical confidence limits for the variables in L  >>>2.5 sigma
% Syntax: [c,c2]=conflim(L,p); 
% 
% Given L, an M by N matrix consisting of N observations (each 
% consisting of M variables) of a random (possibly multidimensional) 
% variate, and a significance level p, CONFLIM sorts the matrix L 
% row-by-row and returns the values in each row which correspond to 
% the confidence limits determined by p. A second pass is then made 
% through L to determine how many significant variables exist in  
% any one column (observation).   
% 
% Inputs: L - The M by N matrix of observations. 
%         p - confidence limit, in a two-tailed form, e. g.  
%             for p=.95 the confidence limits will be at  
%             2.5% and 97.5% of the distribution. 
% 
% Output: c - An M by 2 matrix containing upper and lower confidence limits. 
%         c2- A vector where the i-th element is the probability of 
%             realizing i 'significant' values in one observation.  
% 
% Written by Eric Breitenberger.      Version 1/21/96 
% Please send comments and suggestions to eric@gi.alaska.edu        
% 
c2=[];
 if p==0.95
     k=2.5;
 end
[M,N]=size(L);
Lm=mean(L')';
Ls=std(L')';
c=[Lm-k*Ls,Lm+k*Ls];

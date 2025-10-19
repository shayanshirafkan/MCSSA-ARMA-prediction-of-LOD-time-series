function [c,c2]=conflim(L,p) 
% CONFLIM - returns empirical confidence limits for the variables in L. 
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
 
[M,N]=size(L); 
Ls=sort(L'); Ls=Ls'; 
c=zeros(M,2); 
 
% Now pull out the confidence limits: unless the limits 
% fall exactly on an integer, this is done by   
% selecting the values that lie closest to the *outside* 
% of the selected confidence interval. For example, for 
% N=100 and p=.95 (defaults) the 3rd and 98th sorted 
% observation will be selected as the confidence limits.  
 
p=(1-p)/2;  
index=N*p; 
low=floor(index)+1; 
high=N-floor(index); 
c(:,1)=Ls(:,low); 
c(:,2)=Ls(:,high); 
 
% Now do a second Monte Carlo pass to get the distribution  
% of the number of significant variables in each observation. 
% Each observation is compared to the upper confidence limits 
% and the number of significant variables picked off. The 
% distribution of these numbers is then tabulated. Finally, 
% 'c2' is calculated - c2(i) contains the observed probability 
% that a column of L will yield exactly i significant values. 
 
% Set up to find all the values above the upper confidence limit: 
% If you want the ones *below* the *lower* limit, change the next 
% two statements to:  Ls=c(:,1)*ones(1,N); Ls=L<Ls; 
 
Ls=c(:,2)*ones(1,N);  
Ls=L>Ls; 
Ls=sum(Ls);  
 
% Ls is now a 1 by N vector, where Ls(i) is the number of  
% significant EOFs in the i-th realization. 
Ls=sort(Ls); 
Ls=Ls(find(Ls)); 
Ls=[0 Ls]; 
Ld=diff(Ls); 
 
if max(Ld)==1 
  % Parse this next line and I'll buy you a beer! 
  c2=diff([0 find(Ld)-1 length(Ls)-1]); 
  c2=c2(2:length(c2)); 
elseif max(Ld)==0 
  c2=length(Ls)-1; 
else 
  k=sum(Ld); 
  c2=zeros(1,k);   
  for i=1:k 
    c2(i)=length(find(Ls==i)); 
  end 
end 
 
c2=c2/N; 

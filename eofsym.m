
function s=eofsym(E,tol) 
% EOFSYM - check symmetry of EOFs.  
% Syntax: s=eofsym(E); s=eofsym(E,tol); s=eofsym(E,'1 or 0');  
% 
% Given a matrix of EOFs E, eofsym returns a vector 
% s containing the symmetry of the EOF. 
% 
% If the i-th EOF is symmetric,  s(i)=1,  
% if anti-symmetric,             s(i)=0.  
% if neither sym. or anti-sym., s(i)=-1.  
% 
% s(i)=-1 is only possible if a non-Toeplitz (BK type) 
% covariance matrix was used, or if the tolerance 'tol' 
% is not set high enough. 'tol' is set to tol=10^4*eps  
% by default, or it can be specified as the second argument. 
% 
% If 'tol' is set to the string '1 or 0' or ('0 or 1') the  
% output will be forced to give only ones and zeros. EOFSYM  
% will decide whether the EOFs are symmetric or anti-symmetric 
% based on which assumption gives the lowest rms error. 
% 
% Written by Eric Breitenberger.      Version 5/24/96 
% Please send comments and suggestions to eric@gi.alaska.edu        
% 
 
if nargin==1 
  tol=10^6*eps; % Can be adjusted if necessary 
end 
[M,K]=size(E); 
s=zeros(1,K); 
 
if ~isstr(tol) 
  for k=1:K 
    % Pick off left and right halves of EOF: 
    if rem(M,2)==0 
      L=E(1:M/2,k); R=E(M:-1:M/2+1,k); 
    else 
      L=E(1:(M-1)/2,k); R=E(M:-1:(M+1)/2+1,k); 
    end 
   
    if max(abs(L-R))<tol 
      s(k)=1; 
    elseif max(abs(L+R))<tol 
      s(k)=0; 
    else  
      s(k)=-1; 
%      disp(['Warning: EOF ' num2str(k) ' is neither symmetric or antisymmetric.']) 
    end 
  end 
elseif strcmp(tol,'1 or 0') | strcmp(tol,'0 or 1') 
  Esym=sum(sqrt((E(M:-1:1,:)-E).^2)); 
  Easym=sum(sqrt((-E(M:-1:1,:)-E).^2)); 
  s=(sign(Easym-Esym)+1)/2; 
  ties=find(s==1/2); % Resolve any ties by calling these symmetric 
  s(ties)=ones(1,length(ties)); 
else 
  error('EOFSYM: Improper specification of tolerance.') 
end 
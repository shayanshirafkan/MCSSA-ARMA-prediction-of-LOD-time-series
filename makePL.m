function U=makePL(alpha,Year)
Year=round(Year*365.25);
if size(Year,1)==1
    Year=Year';
end

m=Year(end);

Year1(Year,1)=Year;

H(1)=1;
for i=2:m
    H(i,1)=(alpha/2+i-2)*H(i-1)/(i-1);
end

IDX=ismember(1:m,Year1);

U = toeplitz(H);
U= triu(U);
U=U(:,IDX);
U=U'*U;
return

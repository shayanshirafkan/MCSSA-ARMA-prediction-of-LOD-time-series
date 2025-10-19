clc
clear all
close all
an= dlmread('New Text Document.txt');


[m,n]=size(an);
bb=8;

yn=an(:,bb);
%% number of the predictions
EP=0:1:3;
%% train
z1=20455;
z2=21550;

%%
for bu=1:length(EP);   
oi=EP(1,bu);
ys=yn(z1+oi:1:z2+oi);
yg=ys;
[m1,m2]=size(ys);
vv  = LODR(z1+oi,z2+oi );
ys=ys-vv;
  
t=1:1:m1;
time=t./365.25;

AA=make_TTL_A(time);

AT=[AA(:,1:2)];

x=inv(AT'*AT)*AT'*ys;
y=ys-AT*x;
%% prediction length
po=10;
y4=an(z2+oi+1:1:po-1+z2+oi+1,bb);

vv4=LODR(z2+oi+1,po-1+z2+oi+1);

[tw,ter]=size(y4);
t3=m1+1:1:m1+1+tw-1;
t3=t3/365.25;

A3=make_TTL_A(t3);

AT3=[A3(:,1:2) ];
yp=AT3*x;
%% SSA
M=1*365;
[E,V,A,R,C]=ssa(y,M,'unbiased');     % do the SSA of x 
RT=R;
s=eofsym(E); %the difference of different approachs to provide lag covariance matrix is here when in BK approach all the EOFs are not symmetric but in others it is
RT=R;
% po=365;
[m1,m2]=size(y);
s=m1;
n=zeros(po,1);
y1=[y;n];
for i=1:22;
    tt=1;
    pl=1;
    while tt >=10^(-5) && pl<40;
       
       [E,V,A,R,C]=ssa(y1,M,'unbiased');
       RC1=zeros(s+po,1);
       for j=1:i;
           RC1=R(:,j)+RC1;
       end
       n=RC1(s+1:1:s+po,1);
       y2=[y;n];
[E,V,A,R2,C]=ssa(y2,M,'unbiased');
       RC2=zeros(s+po,1);
       for j=1:i;
           RC2=R2(:,j)+RC2;
       end
tt=norm(RC2-RC1)
i
bu
y1=y2;
pl=pl+1;
    end
y1=y2;
end
n=y1(s+1:1:s+po,1);
h=1:po;
y10=n+yp+vv4;

e1(:,bu)=y4-y10;
% max(max(abs(e1)))*1000;
%% ARMA
noise=[];
noise=zeros(m1,1);
for i=23:M;
noise=RT(:,i)+noise;
end
pMax = 1;
qMax = 7;
LogL = zeros(pMax,qMax);
SumPQ = LogL;

for p = 1:pMax
    for q = 1:qMax
        Mdl = arima(1,0,q);
        [~,~,LogL(p,q)] = estimate(Mdl,noise,  'Display','off');
         
        SumPQ(p,q) = p+q;
    end
end
logL = reshape(LogL,pMax*qMax,1);...
    % Elements taken column-wise 
numParams = reshape(SumPQ,pMax*qMax,1) + 4;
aic = aicbic(logL,numParams);
AIC = reshape(aic,pMax,qMax)
minAIC = min(aic)
[bestP,bestQ] = find(AIC == minAIC)
Mdl = arima(bestP,0,bestQ) ;
[EtsMdl,EstParamCov,logL,info] = estimate(Mdl,noise);
[YF,YMSE] = forecast(EtsMdl,po,'Y0',noise);
df=[];
df=y4-(y10+YF);
bu
e11(:,bu)=df;
end
% po
e11=abs(e11);
[l1,l2]=size(e11);
ut=zeros(po,1);
for i=1:l2;
    ut=abs(e11(:,i))+ut;
end

MAE=ut/l2*1000;

plot(MAE)

save('lodpredictions.mat','e11','MAE')


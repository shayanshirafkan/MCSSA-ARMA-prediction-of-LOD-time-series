clc
clear all
close all
an= dlmread('New Text Document.txt');
[m,n]=size(an);
bb=5;
yn=an(:,bb);
z1=19359;
z2=20819;

ys=yn(z1:1:z2);

[m1,m2]=size(ys);

t=1:1:m1;
time=t./365.25;

%% design matrix
freq=[365.25/365.25,365.25/433];

A(:,1)=ones(length(time),1);
A(:,2)=time;

for i=1:size(freq,2)
A(:,2*i+1)=cos((2*pi*freq(i))*time);
A(:,2*i+2)=sin((2*pi*freq(i))*time);
end

%% LS
x=inv(A'*A)*A'*ys;
noise=ys-A*x;
%%
po=365;

y4=an(z2+1:1:z2+po,bb);

[tw,ter]=size(y4);
t3=m1+1:1:m1+1+tw-1;
t3=t3/365.25;


A3(:,1)=ones(length(t3),1);
A3(:,2)=t3;

for i=1:size(freq,2)
A3(:,2*i+1)=cos((2*pi*freq(i))*t3);
A3(:,2*i+2)=sin((2*pi*freq(i))*t3);
end
%% prediction
yp=A3*x;
%% noise
Mdl = arima(1,0,0) ;
[EtsMdl,EstParamCov,logL,info] = estimate(Mdl,noise);
[YF,YMSE] = forecast(EtsMdl,po,'Y0',noise);
df=[];
df=abs(y4-(yp+YF));

figure;
plot(df)
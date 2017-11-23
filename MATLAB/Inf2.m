%% Q1b
X=  [1.1 0.05;
    1.01 0.01;
    0.99 0.01;
    0.98 0.01;
    1 0.02;
    1.3 0.4];

R1=X(:,2).^2;
X2=1./R1;
X3=X(:,1).*X2;
sumX2=sum(X2);
sumX3=sum(X3);
MLE=sumX3/sumX2;

% The varience is given by 1/sumX2
%therefor std =sqrt(1/sumX2)
Stderr=sqrt(1/sumX2);

%% Q2b
N=10;
T=50;
Slim=30;
 X=linspace(0,30,1000000);

likelihood=@(x)poisspdf(T,x.*N);
evidence=(1./Slim).*integral(likelihood,0,Slim);
prior=1/Slim;
posterior= @(x)likelihood(x).*prior./evidence;

hold on

fplot(posterior,[0,Slim]);

%Mode
mode=T/N;
%Mean
expectation=@(x) posterior(x).*x;
mean=integral(expectation,0,Slim);

% Q2c
mu=7.7;
sig=0.3;
normal=@(x)normpdf(x,mu,sig);
plot(X,normal(X),'g')
Ayy=@(x)likelihood(x).*normal(x);
evidence2=(1./Slim).*integral(Ayy,0,Slim);
posterior2=@(x)Ayy(x).*(1/Slim)./evidence2;

fplot(posterior2,[0,Slim],'r');



%???quad=(-(10*sig^2-mu)+sqrt((10*sig^2-mu)^2+4*N))/2%%%%WTF IS GOING ON???
ALL=posterior2(X);
[M,N]=max(ALL)
Xmax=X(N)

%% Q3
% syms r;
rmin=1/sqrt(4*pi);
R1=linspace(rmin,10000000,100);
E=10e44;
n=8;

N=E.*10e9./(4.*pi.*R1.^2);
NUM=R1.^2.*exp(-N).*N.^n;
plot(R1,N,'r');

% DEN=sum(R1);
% POST=NUM./DEN;
%  
% POSTX1=POST(R1);

% N=@(r)E.*10e9./(4.*pi.*r.^2);
% 
% NUM=@(r)r.^2.*exp(-N).*N.^n
% DEN=integral(N,rmin,10);
% POST=@(r)NUM./DEN;
%  POSTX1=POST(R1);

% fplot(POST);



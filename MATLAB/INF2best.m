%% Q1b
X=  [1.1 0.05;
    1.01 0.01;
    0.99 0.01;
    0.98 0.01;
    1 0.02;
    1.3 0.4];

X1=1./(X(:,2).^2);
X2=X(:,1).*X1;
sumX2=sum(X1);
sumX3=sum(X2);
MLE=sumX3/sumX2;

% The varience is given by 1/sumX2
%therefor std =sqrt(1/sumX2)
Stderr=sqrt(1/sumX2);

%% Q2
% Q2b
N=10;
T=50;
Slim=30;

likelihood=@(x)poisspdf(T,x.*N);
evidence=(1./Slim).*integral(likelihood,0,Slim);
prior=1/Slim;
posterior= @(x)likelihood(x).*prior./evidence;

figure
fplot(posterior,[0,Slim]);
xlabel('Density of stars(1/deg^2)','FontSize',16);
ylabel('Posterior distribution','FontSize',16);
ylim([0 0.6]);

%Mode
mode=T/N;
%Mean
expectation=@(x) posterior(x).*x;
mean=integral(expectation,0,Slim);

% Q2c
mu=7.7;
rsig=0.3;
X=linspace(0,30,1000000);
normal=@(x)normpdf(x,mu,rsig);
% plot(X,normal(X),'g')
Ayy=@(x)likelihood(x).*normal(x);
evidence2=(1./Slim).*integral(Ayy,0,Slim);
posterior2=@(x)Ayy(x).*(1/Slim)./evidence2;
figure
hold on
fplot(posterior,[0,Slim]);
fplot(normal,[0 Slim],'black')
fplot(posterior2,[0,Slim],'r');
legend('Poisson','Gaussian','Combined Posterior')
xlabel('Density of stars(1/deg^2)','FontSize',16);
ylabel('Posterior distribution','FontSize',16);
xlim([0 20]);

%MODE ANALYTICAL
%quad=(-(10*sig^2-mu)+sqrt((10*sig^2-mu)^2+4*N))/2%%%%WTF IS GOING ON???

%MODE NUMERICAL
ALL=posterior2(X);
[M,N]=max(ALL)
X=X(N)

%% Q3d
n=8;
E=1e44;
F=sqrt(E*1e9);
spacing=1/50000;
R=0:spacing:0.5;
r=R.*F;
LHOOD=poisspdf(n,1./(4*pi.*R.^2));
POSTUNORM=(R.^2).*LHOOD;

EVID=F.*trapz(R,POSTUNORM);
POSTNORM=POSTUNORM/EVID;

%MODE NUMERICAL
[M,N]=max(POSTUNORM);
Rmode_num=R(N);
rmode_num=Rmode_num.*F;
%MODE ANALYTICAL
rmode_an=sqrt(E*1e9/(2*pi*(2*n-2)));

%MEAN
EXP_R=trapz(r,R.*POSTNORM);
EXP_r=EXP_R.*F;
%% Q3e
Rsig=Rmode_num*(2-2*n+3/(2*pi*Rmode_num^2))^(-0.5);
%% Q3f
figure
% this is all a bit of code to get the vertical lines and relevant labels
Rlower=Rmode_num-Rsig;
Rupper=Rmode_num+Rsig;
rlower=Rlower.*F;
rupper=Rupper.*F;
[a1,b1]=min(abs(R-Rlower));
[a2,b2]=min(abs(R-Rupper));
linemaxpeak=POSTNORM(N);
linemax1=POSTNORM(b1)
linemax2=POSTNORM(b2)
xpeak=[rmode_num,rmode_num];
xlower=[rlower,rlower];
xupper=[rupper,rupper];
linepeak=[0,linemaxpeak];
linelower=[0,linemax1];
lineupper=[0,linemax2];
hold on;
plot(xlower,linelower,'r --');
plot(xupper,lineupper,'r --');
plot(xpeak,linepeak,'r');
plot(r,POSTNORM,'black'); %Here the normalised posterior is plotted
text(rlower,1e-26,'r_{mode}-\sigma','FontSize',10);
text(rlower,0.1e-26,'\leftarrow 2.73','FontSize',10,'Rotation',20);
text(rupper,1e-26,'r_{mode}+\sigma','FontSize',10);
text(rupper,0.1e-26,'\leftarrow 4.01','FontSize',10,'Rotation',20);
text(rmode_num,2e-26,'r_{mode}','FontSize',10);
text(rmode_num,0.1e-26,'\leftarrow 3.37','FontSize',10,'Rotation',20);
xlabel('Radius from a GRB (m)','FontSize',16);
ylabel('Normalized Posterior distribution','FontSize',16);


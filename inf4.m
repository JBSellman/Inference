%% Q1ai
TD=importdata('Tdisp.txt');
A=TD.data;
a1=A(:,1);
barx=mean(a1);
barx2=mean(a1.^2);
varx=barx2-barx^2;

a2=A(:,2);
bary=mean(a2);
bary2=mean(a2.^2);
vary=bary2-bary^2;

barxy=mean(a1.*a2);
covxy=barxy-barx*bary;

ro=covxy/sqrt(vary*varx);
dof=length(A)-2;
tp=ro.*sqrt(dof/(1-ro^2));
%using the student t distribution table
sig1=tpdf(tp,dof)
%therefor 99.9% sure
%% Q1aii
a=size(A,1);
[A1,B1]=sort(A(:,1));
[C1, D1]=sort(B1);

[A2,B2]=sort(A(:,2));
[C2, D2]=sort(B2);

Order=[D1 D2 D1-D2 (D1-D2).^2];

r=1-6*sum(Order(:,4))/(a*(a^2-1));
ts=r.*sqrt(dof/(1-r^2))
sig2=tpdf(ts,dof)
%therefor 99.8% sure
%% Q1b
figure
plot(a1,a2,'o')
hold on
%Comment-Prefer pearson as data looks pretty linear to begin with and gives
%a higher significance
%% Q1c
%for y(x)=ay*x+by
ay=covxy./(barx2-barx^2);
by=bary-ay.*barx;
fity=ay.*a1+by;
plot(a1,fity)

%forx(y)
ax=covxy/(bary2-bary^2);
bx=barx-ax.*bary;
fitx=1./ax.*a1-bx./ax;
plot(a1,fitx,'g')
%COMMENT ON THE DIFFERENCE
%% Q2a
L0=3.846e26;%Watts
Const=1e9*L0;
alpha=-0.7;
expnorm=@(u) (u/14).^alpha.*exp(-u/14)/14;
unorm=integral(expnorm,1,inf);
uarg=@(u) (u/14).^(alpha+1).*exp(-u/14)
Lexp=integral(uarg,1,inf)/unorm%*Const;
figure
fplot(uarg,[1 100])
%% Q2b
Data=[1.39; 1.40; 1.29; 5.95; 2.97; 1.63; 1.17; 2.06; 4.69; 2.48]%*Const;
Cmltv=zeros(length(Data),1);
Cmltv(1)=Data(1);
for i=2:length(Data)
    Cmltv(i)=Data(i)+Cmltv(i-1);
end
Lexpvec=ones(length(Data),1)*Lexp;
Lexpcum=cumsum(Lexpvec);
read=1:10;
figure;
plot(read,Cmltv);
hold on; 
plot(read,Lexpcum);
%%
% barx=mean(Data);
% barx2=mean(Data.^2);
% varx=barx2-barx^2;
% 
% bary=mean(Lexpvec);
% bary2=mean(Lexpvec.^2);
% vary=bary2-bary^2;
% 
% barxy=mean(Lexpvec.*Data);
% covxy=barxy-barx*bary;
% 
% ro=covxy/sqrt(vary*varx);
% dof=length(Data)-2;
% tp=ro.*sqrt(dof/(1-ro^2));
% %using the student t distribution table
% sig1=tpdf(tp,dof)
%% Q3a
V=importdata('inf4q3.txt');
Nv=length(V);
dofv=Nv-1;
vbar=sum(V)./Nv;

intvl=0.1;
intlow=chi2inv(intvl/2,dofv);
inthigh=chi2inv(1-intvl/2,dofv);
vhigh=inthigh/dofv
vlow=intlow/dofv

% chi2=dofv*vbar^2;
% peen=chi2pdf(12.592,6)

%%
X=0:1000:1000;
chi2=dovf.*v






%% Q1a,b
%Here is the transition matrix describing the iteration of probabilities 

Mat=[0.9 0.1; 0.4 0.6];
N=[2 987654321];    %N of days

% using the Markov matrix technique, one sets the matrix power equal to the
% number of days of iteration, then selects the correct element. In this
% case, we started on a sunny day and want to know the probability of
% ending on a sunny day, so we select the component in the first row and
% first column
% prob that its sunny, N days after a sunny day is given by
Pa=Mat^N(1);
Proba=Pa(1,1);

Pb=Mat^N(2);
Probb=Pb(1,1);
Ps=1;
Pr=0;
X=0:20;
for i=1:20
    MatX=Mat^i;
    Ps(i+1)=MatX(1,1);
    Pr(i+1)=MatX(2,1);
end
figure
plot(X,Ps,'b .','MarkerSize',8);
hold on
plot(X,Pr,'r');
legend('First day sunny','First day Rainy');
xlabel('Number of days ellapsed','FontSize',16);
ylabel('Probability of being sunny','FontSize',16);
%% Q1c
%Using the bayesian formula P(a|b)=P(b|a)P(a)/P(b)
%prior P(S1)
prior= Probb; %as this is the convergent probability of any day being sunny

% likelihood P(S2|S1)
likelihood=0.9;

% evidence P(S2)=P(S2,S1)+P(S2,S1-)
evidence=0.9*Probb+0.4*(1-Probb);

% posterior
posterior=likelihood*prior/evidence;%???
%% Q2
N=14;       %number of galaxies
k=3;        %number of spiral galaxies
Psp=0.75;   %expected probabiliity
Pnonsp=(1-Psp);

Binom=binopdf(k,N,Psp);
% C=factorial(N)/(factorial(k)*factorial(N-k));
% Prob=C*Psp^k*Pnonsp^(N-k);%returns the same as MATLAB binomial func
Ptot=0; %Intialising the funstion to sum probabilities

%This for loop calculates the probaility of observing 0,1,2 or 3spirals and
%sums them into Ptot
for i=0:k
    P(i+1)=binopdf(i,N,Psp);
    Ptot=Ptot+binopdf(i,N,Psp);
end
%% Q3a
i=0;        %initialising
P=1;
Mean=6;  %The mean no# counts over 6 hours

%This loop only plots counts which have a 0.0001% probability and fills a
%vector for x(number of counts) and Poiss(prob on x counts)
while P>0.0001
    P=poisspdf(i,Mean);
    Poiss(i+1)=P;
    x(i+1)=i;
    i=i+1;
end
plot(x,Poiss,'o','MarkerSize',8);
xlabel('Probability','FontSize',20);
ylabel('Number of detections','FontSize',20);
% Q3b,c
Prob6=poisspdf(6,6);
%actual poisson calculation 
P6=Mean^6*exp(-Mean)/factorial(6)

%The probability of measuring 10 or more counts is equal to 1 minus the
%prob of 9 or less, P(>=10)=1-P(<10)

Prob9less=0;   %initialisation
for i=0:9
    Prob9less=Prob9less+poisspdf(i,Mean);
end
Prob10=1-Prob9less;
%% Q3d
% Using format long, the probabilities of seeing 5 counts and 6 counts were
% compared
format long
Prob6=poisspdf(6,6)
Prob5=poisspdf(5,6)

%% Q3e
P=1; 
Mean=6;

%when the probaility measured for i counts falls below 0.01 and when i>0
%the minimum count to satisfy the criterea is output
for i=0:1000
    P=poisspdf(i,Mean);
    if P<0.01 && i~=0
       Mincount=i; 
       break
    end
end
Mincount



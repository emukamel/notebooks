clear;

%% step 1 generate data using xiaowei's code and test the DPsearch algorithm.
N = 4;
M = 20;
K = 2;

%P = unifrnd(0, 1, [N, K]);
p = unifrnd(0, 1, [N, 1]);
P = [p, p];
R = mnrnd(1, [0.5 0.5], M)';
X0 = binornd(1, P * R);
%X = [repmat([0; 0; 0; 1], [1, M / 2]), repmat([1; 0; 1; 0], [1, M / 2])];

[g, ~] = find(R == 1);

X = X0';

w = 4;
addpath('C:\Users\hongxiaozhu\Documents\MATLAB\David_Xie\DPsearch_betabernoulli');
beamSize=20; % 20-100
Y= g'; % or Y=[];
heuristic='i';
alpha = 1; 
% G0alphabeta = [2.1,2.2,1.9,1.8;
%                2,  1.8, 1.85, 2.1];
G0alphabeta = 2*ones(2,w);
[r,marginalOrder] = DPsearch_betabernoulli(X,alpha,G0alphabeta,beamSize,Y,heuristic);

[l0, marg_like]= log_likelihood_given_c(X, r.c, G0alphabeta);
[l01, marg_like1]= log_likelihood_given_c(X, g', G0alphabeta);


%% step 2. generate data using DPM method and test DPsearch algorithm.
clear;
alpha = 1;
w = 4;
G0alphabeta = 2*ones(2,w);
N = 20;
addpath('C:\Users\hongxiaozhu\Documents\MATLAB\David_Xie\DPsearch_betabernoulli');
[pX, c, K] = DP_rnd_beta(alpha, G0alphabeta, w, N);
table(c)

mean(pX(c==1,:))
mean(pX(c==2,:))
X = binornd(1,pX);

beamSize=20; % 20-100
Y= c'; % or Y=[];
heuristic='i';
[r,marginalOrder] = DPsearch_betabernoulli(X,alpha,G0alphabeta,beamSize,Y,heuristic);

[l0, marg_like]= log_likelihood_given_c(X, r.c, G0alphabeta);
[l01, marg_like1]= log_likelihood_given_c(X, c, G0alphabeta);




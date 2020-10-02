function [A,b,Q,delta,x0]=GenerateSElasticNet(m,n,s,alpha,sigma1,sigma2)
init=2055615866;
rand('seed',init);
randn('seed',init);
normrnd('seed',init);
[A,Rtmp] = qr(randn(n,m),0);
[B,Rtmp] = qr(randn(n,n),0);
A  = A';
Q  = B'*B;
x0 = zeros(n,1); 
p  = randperm(n);
index=p(1:s);
x0(index)=normrnd(0,sigma1^2,s,1);
epsilon=normrnd(0,sigma2^2,m,1);
b = A*x0+epsilon;
delta = (alpha*norm(x0,1)+(1-alpha)*x0'*Q*x0)/2;
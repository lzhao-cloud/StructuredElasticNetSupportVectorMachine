function [x,h]=VAPP_SElasticNet_Cone(A,b,Q,delta,alpha,iteration,fstar)
[m,n]=size(A);
x=ones(n,1);
p=zeros(n+1,1);
MM=sqrt(n+1)*norm(b,2)^2/(2*delta)+1;
gamma=1;
epsilon=0.2;%2000%0.3%1000;
for i=1:iteration
    xk=x;
    pk=p;
    qk_intern=pk+gamma*[(1-alpha)*xk'*Q*xk-delta;alpha*xk];
    [qkxx,qky]=Proj_infty(qk_intern(2:n+1),qk_intern(1));
    qk=[qky;qkxx];
    qk=min(1,MM/norm(qk,2))*qk;
    x=xk-epsilon*(A'*(A*xk-b)+[((1-alpha)*(Q'+Q)*xk)';alpha*eye(n)]'*qk);
    x=max((-delta/alpha),min((delta/alpha),x));
    p_intern=pk+gamma*[(1-alpha)*x'*Q*x-delta;alpha*x];
    [pxx,py]=Proj_infty(p_intern(2:n+1),p_intern(1));
    p=[py;pxx];
    p=min(1,MM/norm(p,2))*p;
    %save informations;
    h.obj(i)=abs((1/2)*norm(A*xk-b,2)^2-fstar);
    h.constraint(i)=max(0,alpha*norm(xk,1)+(1-alpha)*xk'*Q*xk-delta);
    h.relation(i)=norm(x-xk,2)/max(norm(x,2),1);
    h.plus(i)=h.obj(i)+h.constraint(i);
end
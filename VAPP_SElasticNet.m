function [x,h]=VAPP_SElasticNet(A,b,Q,delta,alpha,iteration,fstar)
[m,n]=size(A);
x=ones(n,1);
p=0;
MM=norm(b,2)^2/(2*delta)+1;
gamma=0.1;%0.1;2000%0.1;%1;
epsilon=0.07;%0.07;2000%0.15;0.13;
for i=1:iteration
    xk=x;
    pk=p;
    qk=min(MM,max(0,pk+gamma*(alpha*norm(xk,1)+(1-alpha)*xk'*Q*xk-delta)));
    zetak=A'*(A*xk-b)+(1-alpha)*qk*(Q+Q')*xk;
    x=sign(xk-epsilon*zetak).*max(0,abs(xk-epsilon*zetak)-epsilon*alpha*qk*ones(n,1));
    %x=max((-delta/alpha),min((delta/alpha),x));
    p=min(MM,max(0,pk+gamma*(alpha*norm(x,1)+(1-alpha)*x'*Q*x-delta)));
    %save informations;
    h.obj(i)=abs((1/2)*norm(A*xk-b,2)^2-fstar);
    h.constraint(i)=max(0,alpha*norm(xk,1)+(1-alpha)*xk'*Q*xk-delta);
    h.relation(i)=norm(x-xk,2)/max(norm(x,2),1);
    h.plus(i)=h.obj(i)+h.constraint(i);
end

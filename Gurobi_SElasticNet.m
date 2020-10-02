function [u_v,Objective_v] = Gurobi_SElasticNet(A,b,Q,delta,alpha)
n=length(Q);
u=sdpvar(n,1,'full');
Objective=(1/2)*norm(A*u-b,2)^2;
CX=[];
CX=[CX;alpha*norm(u,1)+(1-alpha)*u'*Q*u<=delta];
options=sdpsettings('verbose',1,'solver','cplex');
sol=optimize(CX,Objective,options);
Objective_v=value(Objective);
u_v=value(u);
sol.info
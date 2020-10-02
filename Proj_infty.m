function [xx,y]=Proj_infty(aa,b)
d=length(aa);
y=0;
xx=zeros(d,1);
abs_aa=abs(aa);
sort_abs_aa=sort(abs_aa);
%aaa=[0;sort_abs_aa];
if b>=sort_abs_aa(d)
    y=b;
    xx=aa;
    return;
end
%if -b>=norm(aa,1)
%  return;
%end
t=0;
for j=d-1:-1:1
    t=t+sort_abs_aa(j+1);
    yy=(b+t)/(1+d-j);
    if yy>sort_abs_aa(j)&&yy<=sort_abs_aa(j+1)
       y=yy;
       xx=sign(aa).*min(abs_aa,y);
       return;
    end
end
t=t+sort_abs_aa(1);
yy=(b+t)/(1+d);
if yy>0&&yy<=sort_abs_aa(1)
   y=yy;
   xx=sign(aa).*min(abs_aa,y);
   return;
end

%aaaaa=0;
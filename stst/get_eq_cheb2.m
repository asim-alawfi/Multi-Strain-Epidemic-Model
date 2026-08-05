function [yout,ch2]=get_eq_cheb2(r,g,info,onlypos)
for i=2:-1:1
    pol{i}=@(x,y)info.y_eqn_fun{i}(x,y,r,g);
end
ch2=chebfun2v(pol{:},info.dom);
rt=roots(ch2);
y2=info.y2fun(rt(:,1),rt(:,2),r,g);
y3=info.y3fun(rt(:,1),rt(:,2),r,g);
yout=cat(2,rt(:,1),y2,y3,rt(:,2)).';
yout(:,any(abs(imag(yout))>1e-13,1))=NaN;
if nargin<4 ||~onlypos
    return
end
yout(:,any(yout<0,1))=NaN;
end

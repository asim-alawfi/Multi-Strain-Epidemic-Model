function [x,pars,stab,ev]=x_eq_from_y(r_rg,y,A,par0,in,FJ)
%% augment solution y by z and w component to get full equilibrium
allnan=all(isnan(y),[1,2]);
y(:,:,allnan)=[];
for i=numel(r_rg):-1:1
    u=r_rg(i)*double(A)*squeeze(y(:,i,:));
    ws(:,i,:)=u./(u+1);
    zs(:,i,:)=y(:,i,:)*r_rg(i)./(1+r_rg(i)*y(:,i,:));
    pars(:,i)=par0;
    pars(in.r,i)=r_rg(i);
end
x=cat(1,y,zs,ws);
ev=NaN(size(x));
stab=NaN(size(x,[2,3]));
if nargin<=5
    return
end
for i=1:size(x,2)
    for j=1:size(x,3)
        if any(isnan(x(:,i,j)))
            continue
        end
        J=FJ(x(:,i,j),pars(:,i));
        ev(:,i,j)=eig(J);
        stab(i,j)=sum(real(ev(:,i,j))>0);
    end
end
end
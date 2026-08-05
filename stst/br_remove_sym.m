function [x,par,stab]=br_remove_sym(x,par,stab,crit)
sel=squeeze(all(crit(x),1));
[xdim,nr]=size(x,[1,2]);
[x(:,sel),stab(sel)]=deal(NaN);
[x,stab]=deal(reshape(x,xdim,nr,[]),reshape(stab,nr,[]));
allnan=all(isnan(x),[1,2]);
x(:,:,allnan)=[];
stab(:,allnan)=[];
end

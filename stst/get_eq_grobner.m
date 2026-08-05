function yout=get_eq_grobner(r,g,info,onlypos)
for i=length(info.cf):-1:1
    y{info.rootindex,i}=reshape(roots(info.cf{i}(r,g)),1,[]);
    for k=1:length(info.fun)
        y{info.dependent(k),i}=info.fun{k}(y{info.rootindex,i},r,g);
    end
    for k=1:length(info.target)
        y{info.target(k),i}=info.factor(k)*y{info.source(k),i};
    end
    yc{i}=cat(1,y{:,i});
end
yout=cat(2,yc{:});
yout(:,any(abs(imag(yout))>1e-13,1))=NaN;
if nargin<4 ||~onlypos
    return
end
yout(:,any(yout<0,1))=NaN;
end

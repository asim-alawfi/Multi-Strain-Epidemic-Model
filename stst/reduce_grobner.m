function info=reduce_grobner(y_eqn,y,pars)
%% solve nonlinear equations using Grobner basis
for i=length(y_eqn):-1:1 % convert into polynomial equations
    [y_eqn_poly(i),den(i)]=numden(y_eqn(i));
    assert(logical(0==simplify(y_eqn_poly(i)-y_eqn(i)*den(i))));
end
%% Reduce system of polynomial equations to single polynomial using Groebner basis
gb=gbasis(y_eqn_poly,y,'MonomialOrder','lexicographic').';
%% gb(1) to gb(end-1) are expressing y(i) as a function of y(end)
ny=length(y);
for i=ny-1:-1:1
    sol{i}=solve(gb(i),y(i));
    fun{i}=matlabFunction(sol{i},'Vars',[y(ny),pars]);
end
gbf=factor(gb(end));
degs=poly_deg(gbf,y(ny));
gbfpolys=gbf(degs>0);
for i=length(gbfpolys):-1:1
    cfsym{i}=coeffs(gbfpolys(i),y(ny),'All');
    cf{i}=matlabFunction(cfsym{i},'Vars',pars);
end
info=struct('gb',gb,'sol',{sol},...
    'gbf',gbf,'degs',degs,'gbfpolys',gbfpolys,'cfsym',{cfsym},...
    'cf',{cf},'fun',{fun},'rootindex',ny,'dependent',1:ny-1);
end
%%
function dpoly=poly_deg(expr,y)
dpoly = -ones(size(expr));
for i=1:numel(expr)
    try
        dpoly(i)=polynomialDegree(expr(i),y);
    end
end
end



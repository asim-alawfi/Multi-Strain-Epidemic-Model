%%  plot branches of equilibria with stability information
clear
format compact
format short g
load('stst_sym_grobner.mat');
%%
clr=lines();
stst={...
    r00,     x00,    stab00,    1,'trivial';...
    r_rgfine,x13_240,stab13_240,2,'(13)(2=4=0)';...
    r_rgfine,x13_40, stab13_40, 3,'(13)(4=0)';...
    r_rgfine,x13,    stab13,    4,'(13)';...
    r_rgfine,x1234,  stab1234,  5,'(1234)';...
    r_rgfine,x40g,   stab40g,   6,'(4=0)'...
    };
[nbr,spacing]=deal(size(stst,1),40);
lw={'LineWidth',2,'MarkerSize',6};
maxunst=max([stab00(:);stab13_240(:);stab13_40(:);stab13(:);stab40g(:)]);
symtype='o+x*>...';
comp=1;
txtused=false(1,nbr);
figure(1);clf;ax=gca;hold(ax,'on');
for i=1:nbr
    [rg,x,stab,cl,txt]=deal(stst{i,:});
    for k=0:maxunst
        for j=1:size(x,3)
            if all(isnan(x(:,:,j)),'all')
                continue
            end
            eqstab=find(stab(:,j)==k);
            if isempty(eqstab)
                continue
            end
            spaced=eqstab(1:spacing:end);
            plot(ax,rg(eqstab),squeeze(x(comp,eqstab,j)),'.',...
                'Color',clr(k+1,:),lw{:},'HandleVisibility','off');
            plot(ax,rg(spaced),squeeze(x(comp,spaced,j)),...
                symtype(cl),'Color',clr(k+1,:),lw{:},'HandleVisibility','off');
        end
    end
end
for i=1:nbr
    plot(ax,NaN,NaN,['k',symtype(i)],'DisplayName',stst{i,end},lw{:});
end
for k=0:maxunst
    plot(ax,NaN,NaN,'-','Color',clr(k+1,:),'DisplayName',sprintf('#unst=%0d',k),lw{:});
end
legend(ax,'Location','best');
xlabel(ax,'$r$','Interpreter','latex');
ylabel(ax,'$y_1$ (scaled by $e$)','Interpreter','latex');
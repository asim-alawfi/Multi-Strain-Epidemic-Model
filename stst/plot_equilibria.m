%%  plot branches of equilibria with stability information
clear
format compact
format short g
load('stst_sym_grobner.mat'); % needs a run of compute_equilibria_grobner first
%%
clr=lines();
stst={...
    r00,     x00,    stab00,    'o','1=2=3=4=0';...
    r_rgfine,x13_240,stab13_240,'+','(13)(2=4=0)';...
    r_rgfine,x13_40, stab13_40, 'x','(13)(4=0)';...
    r_rgfine,x13,    stab13,    '*','(13)';...
    r_rgfine,x1234,  stab1234,  '>','(1234)'
    };
[nbr,spacing]=deal(size(stst,1),30);
lw={'LineWidth',2};
maxunst=max([stab00(:);stab13_240(:);stab13_40(:);stab13(:);stab40g(:)]);
comp=1;
txtused=false(1,nbr);
figure(1);clf;ax=gca;hold(ax,'on');
for i=1:nbr
    [rg,x,stab,symb,txt]=deal(stst{i,:});
    for j=1:size(x,3)
        plot(ax,rg,squeeze(x(comp,:,j)),'-',...
            'Color','k',lw{:},'HandleVisibility','off');
        for k=0:maxunst
            eqstab=find(stab(:,j)==k);
            spaced=eqstab(1:spacing:end);
            plot(ax,rg(spaced),squeeze(x(comp,spaced,j)),...
                symb,'Color',clr(k+1,:),lw{:},'HandleVisibility','off');
        end
    end
end
%% create labels for legend
for i=1:nbr
    plot(ax,NaN,NaN,['k',stst{i,4}],'DisplayName',stst{i,end},lw{:});
end
for k=0:maxunst
    plot(ax,NaN,NaN,'.','Color',clr(k+1,:),'DisplayName',...
        sprintf('#unst=%0d',k),'MarkerSize',14);
end
legend(ax,'Location','best');
xlabel(ax,'$r$','Interpreter','latex');
ylabel(ax,'$y_1$ (scaled by $e$)','Interpreter','latex');
%%  plot branches of equilibria with stability information
clear
format compact
format short g
load(fullfile('stst_sym_e_0_004.mat')); % needs a run of compute_equilibria_grobner first
%%
clr=lines();
stst={...
    r00,     x00,    stab00,    'o','1=2=3=4=0';...
    r_rgfine,x13_240,stab13_240,'+','(13)(2=4=0)';...
    r_rgfine,x_2340, stab_2340, '.','(1)(2=3=4=0)';...
    r_rgfine,x13_40, stab13_40, 'x','(13)(4=0)';...
    r_rgfine,x13_24,    stab13_24,    '*','(13)(24)';...
    r_rgfine,x1234,  stab1234,  '>','(1234)'
    };
[nbr,spacing]=deal(size(stst,1),30);
lw={'LineWidth',2};
maxunst=max([stab00(:);stab13_240(:);stab13_40(:);stab13_24(:);stab40g(:)]);
comp=[1,2];
txtused=false(1,nbr);
figure(1);clf;tl=tiledlayout(2,1,'TileSpacing','tight');
ax(1)=nexttile(tl,1);ax(2)=nexttile(tl,2);hold(ax(1),'on');hold(ax(2),'on');
for i=1:nbr
    [rg,x,stab,symb,txt]=deal(stst{i,:});
    for j=1:size(x,3)
        for pl=1:length(comp)
            plot(ax(pl),rg,squeeze(x(comp(pl),:,j)),'-',...
                'Color','k',lw{:},'HandleVisibility','off');
        end
        for k=0:maxunst
            eqstab=find(stab(:,j)==k);
            spaced=eqstab(1:spacing:end);
            for pl=1:length(comp)
                plot(ax(pl),rg(spaced),squeeze(x(comp(pl),spaced,j)),...
                    symb,'Color',clr(k+1,:),lw{:},'HandleVisibility','off');
            end
        end
    end
end
%% create labels for legend
for i=1:nbr
    plot(ax(1),NaN,NaN,['k',stst{i,4}],'DisplayName',stst{i,end},lw{:});
end
for k=0:maxunst
    plot(ax(1),NaN,NaN,'.','Color',clr(k+1,:),'DisplayName',...
        sprintf('#unst=%0d',k),'MarkerSize',14);
end
legend(ax(1),'Location','best');
xlabel(ax(2),'$r$','Interpreter','latex');
for pl=1:length(comp)
    ylabel(ax(pl),sprintf('$y_%d$ (scaled by $e$)',pl),'Interpreter','latex');
end
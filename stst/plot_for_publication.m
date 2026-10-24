clear
format compact
format short g
load(fullfile('stst_sym_e_0_004.mat')); 
%% 
stst = { ...
    r00,      x00,      stab00,     'o', '$k=0$'                ; ...  % trivial (1=2=3=4=0)
    r_rgfine, x13_240,  stab13_240, 's', '$k\sim(13)(24)$'      ; ...  % (13)(2=4=0)
    r_rgfine, x_2340,   stab_2340,  'd', '$k\sim(24)$'          ; ...  % (1)(2=3=4=0)   <- adjust if S4
    r_rgfine, x13_40,   stab13_40,  '^', '$k\sim(13),\,y_4{=}0$'; ...  % (13)(4=0)
    r_rgfine, x13_24,      stab13_24,     'p', '$k\sim(13)(24)$'          ; ...  % (13)(24)
    r_rgfine, x1234,    stab1234,   '>', '$k\sim(1234)$'        };     % (1234)

comp     = [1 2];      
spacing  = 40;         
markSize = 6;          
nbr      = size(stst,1);

%
maxunst  = max(cellfun(@(s) max(s(:)), stst(:,3)));

% 
base = [0.12 0.30 0.55; 0.20 0.60 0.72; 0.96 0.76 0.20; ...
        0.90 0.45 0.15; 0.75 0.15 0.15];
if maxunst+1 <= size(base,1)
    stabcol = base(1:maxunst+1,:);
else
    stabcol = turbo(maxunst+1);          
end

%% figure 
figure(2); clf
set(gcf,'Color','w','Units','centimeters','Position',[2 2 15 14])
tl = tiledlayout(2,1,'TileSpacing','tight','Padding','compact');

ax(1) = nexttile(tl); hold(ax(1),'on'); box(ax(1),'on')   % top panel
ax(2) = nexttile(tl); hold(ax(2),'on'); box(ax(2),'on')   % bottom panel

%% plot branches
for i = 1:nbr
    [rg,x,stab,symb] = stst{i,1:4};
    for j = 1:size(x,3)                                       
        
        for pl = 1:numel(comp)
            plot(ax(pl), rg, squeeze(x(comp(pl),:,j)), '-k', ...
                 'LineWidth',0.75, 'HandleVisibility','off')
        end
        % stability markers, coloured by #unstable
        for k = 0:maxunst
            idx = find(stab(:,j)==k);
            idx = idx(1:spacing:end);
            if isempty(idx), continue, end
            for pl = 1:numel(comp)
                plot(ax(pl), rg(idx), squeeze(x(comp(pl),idx,j)), symb, ...
                     'MarkerFaceColor',stabcol(k+1,:), 'MarkerEdgeColor','k', ...
                     'MarkerSize',markSize, 'LineWidth',0.5, ...
                     'HandleVisibility','off')
            end
        end
    end
end

%% permutation
hB = gobjects(1,nbr);
for i = 1:nbr
    hB(i) = plot(ax(1), NaN, NaN, stst{i,4}, 'LineStyle','none', ...
        'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k', ...
        'MarkerSize',markSize, 'DisplayName', stst{i,5});
end
lgB = legend(ax(1), hB, 'Interpreter','latex', 'FontSize',9, ...
             'Location','northwest', 'Box','on', 'Color','w');
title(lgB, 'isotropy $k\sim\sigma$', 'Interpreter','latex')

%% stability colours
hS = gobjects(1,maxunst+1);
for k = 0:maxunst
    hS(k+1) = plot(ax(2), NaN, NaN, 's', 'LineStyle','none', ...
        'MarkerFaceColor',stabcol(k+1,:), 'MarkerEdgeColor','none', ...
        'MarkerSize',markSize+3, 'DisplayName', sprintf('$\\#\\mathrm{unst}=%d$',k));
end
lgS = legend(ax(2), hS, 'Interpreter','latex', 'FontSize',9, ...
             'Location','northwest', 'Box','on', 'Color','w');
title(lgS, 'stability', 'Interpreter','latex')

%% 
set(ax, 'FontSize',11, 'Layer','top', 'LineWidth',0.8, ...
        'TickLabelInterpreter','latex', 'XTick',0.5:0.5:3.5)
linkaxes(ax,'x');  xlim(ax(1),[0.5 3.5])
ax(1).XTickLabel = [];                        
xlabel(ax(2), '$r$', 'Interpreter','latex', 'FontSize',12)
for pl = 1:numel(comp)
    ylabel(ax(pl), sprintf('$y_{%d}\\ (\\mathrm{scaled\\ by}\\ e)$',comp(pl)), ...
           'Interpreter','latex', 'FontSize',12)
end

text(ax(1), 0.985,0.07, '(a)', 'Units','normalized', 'HorizontalAlignment','right', ...
     'Interpreter','latex', 'FontSize',12)
text(ax(2), 0.985,0.93, '(b)', 'Units','normalized', 'HorizontalAlignment','right', ...
     'Interpreter','latex', 'FontSize',12)

%% %export ----------------------------------------------------------
% exportgraphics(tl, 'bifurcation.pdf', 'ContentType','vector')
% exportgraphics(tl, 'bifurcation.png', 'Resolution',600)
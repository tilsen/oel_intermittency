function [] = plot_sim(h)

P = gen_panel('s-osc',cos(h.Sp));

P(end+1) = gen_panel('c-osc',h.Cx.*cos(h.Cp));

P(end+1) = gen_panel('s-exc',[h.Sx h.Ex], ...
    'color',[lines(h.Ns); zeros(h.Ne,3)]', ...
    'linestyle',[repelem({'-' '--'},[h.Ns h.Ne])], ...
    'linewidth',[repelem([2 1],[h.Ns h.Ne])], ...
    'panelheight',2);

P(end+1) = gen_panel('s-fb',[h.Fx]);

P(end+1) = gen_panel({'annealer,','monitor'},[h.An h.Me], ...
    'color',[1 0 0; 0 0 0]', ...
    'linewidth',2);

P = struct2table(P);
pans = repelem([1:height(P)],P.panelheight)';

st = dbstack;
if any(contains({st.name},'LiveEditor'))
    figh = figure('units','pixels','OuterPosition',[1 1 640 480]);
    ax = stf([pans],[0.035 0.05 0.01 0.01],[0 0.01],'parent',figh);
else
    ax = stf([pans],[0.035 0.05 0.01 0.01],[0 0.01]);
end

for i=1:height(P)
    data = P.data{i};
    set(gcf,'currentaxes',ax(i));
    for j=1:size(data,2)
        plot(h.t,data(:,j), ...
            'linew',P.linewidth{i}(j), ...
            'Color',P.color{i}(:,j), ...
            'linestyle',P.linestyle{i}{j}); hold on;
    end
    axis tight;
    text(min(xlim)+0.001*diff(xlim),max(ylim),[P.name{i}], ...
        'verti','top','fontweight','bold','fontsize',12);
end

linkaxes(ax,'x');
axis(ax,'tight');
axrescaley(ax,0.01);

defticks;
set(ax,'XGrid','on','YGrid','on','TickLen',0.001*[1 1]);
set(ax(1:end-1),'XTickLabel',[]);

end

%%
function [P] = gen_panel(name,data,varargin)

p = inputParser;

def_color = [];
def_ls = [];
def_lw = [];
def_panheight = 1;

addRequired(p,'name');
addRequired(p,'data');
addParameter(p,'color',def_color);
addParameter(p,'linestyle',def_ls);
addParameter(p,'linewidth',def_lw);
addParameter(p,'panelheight',def_panheight);

parse(p,name,data,varargin{:});

r = p.Results;

P.name = name;
if ~iscell(name)
    P.name = {P.name};
end
P.data = data;

P.color = r.color;
P.linestyle = r.linestyle;
P.linewidth = r.linewidth;
P.panelheight = r.panelheight;

N = size(P.data,2);

if isempty(r.linewidth)
    P.linewidth = 2*ones(1,N);
end
if numel(P.linewidth)<N, P.linewidth = repmat(P.linewidth,1,N); end

if isempty(r.color)
    P.color = lines(N)';
end

if isempty(r.linestyle)
    P.linestyle = repmat({'-'},1,N);
end
if numel(P.linestyle)<N, P.linestyle = repmat(P.linestyle,1,N); end


end

function [h] = intermittency_model(h)

dbstop if error;

%% default parameters
P = model_default_params;
for i=1:height(P)
    if ismatrix(P.val{i})
        eval([P.name{i} '=' mat2str(P.val{i}) ';']);
    else
        eval([P.name{i} '=[' num2str(P.val{i}) '];']);
    end
end

%% overwrite parameters
if nargin>0 && ~isempty(h)
    ff = fieldnames(h);
    for i=1:length(ff)
        if exist(ff{i},'var')
            eval([ff{i} '=[' num2str(h.(ff{i})) '];']);
        end
    end
end

%% derived parameters:
t = 0:dt:(T-dt);        %time vector
nT = T/dt;              %number of time steps
omega = 2*pi*f;         %oscillator intrinsic frequency (radians/s)

%% coupling parameters with random initialization
%c-s coupling (noun or verb)

C_S = repmat(rand(1,Nc)>C_S_prop,3,1);
C_S(end,:) = ~C_S(end,:);

%c-c coupling
C_C = 2*(rand(Nc)-0.5);
C_C(logical(eye(Nc)))=0;

ix_diff_sign = sign(C_C) ~= sign(C_C');
C_C(tril(ix_diff_sign)) = -C_C(tril(ix_diff_sign));

%%
V_E = @(x,c)-(x-c).*(exp(-((x-c).^2)./(2*sigma_VE^2)));
V_Me = @(x)exp(-(x.^2)./(2*sigme_VMe^2));

%% initialize variables
EV = {}; %event log

Cx = nan(nT,Nc);
Cp = nan(nT,Nc);
Sx = nan(nT,Ns);
Sp = nan(nT,Ns);
Fx = nan(nT,Ns);
Ex = nan(nT,Ns+2);

An = nan(nT,1); %annealer
Me = nan(nT,1); %e-organization monitor

%% initial conditions

%random initial phases
Cp(1,:) = 2*pi*rand(1,Nc);
Sp(1,:) = 2*pi*rand(1,Ns);

%random background noise
Cx(1,:) = 0.01*rand(1,Nc);

%initial S gradient
Sx(1,:) = Sx0;

% initial feedback
Fx(1,:) = 0;

%excitation organization
Ex(1,:) = linspace(0,1,Ne);

%e-organization annealer
An(1) = 0;

%e-organization monitor
Me(1) = 0;

%start in initial e-organization regime
init_eorg = true;

%%
%randomly select two noun concepts and a verb to receive energy from the environment
N_ixs = find(C_S(1,:));
V_ixs = find(C_S(end,:));

rN_ixs = N_ixs(randperm(length(N_ixs),env_num_N));
rV_ixs = V_ixs(randperm(length(V_ixs),env_num_V));

ix_env_C = ismember(1:Nc,[rN_ixs rV_ixs]);

%%
EV(end+1,:) = {t(1), 'sim_init'};

for i=1:nT-1

    %% temporary variables
    PHI_SC = Cp(i,:)-Sp(i,:)';
    PHI_CC = Cp(i,:)-Cp(i,:)';
    PHI_SS = Sp(i,:)-Sp(i,:)';
    PHI_CS = Sp(i,:)-Cp(i,:)';

    wC_S = (C_S.*Cx(i,:))'; %c-s map weighted by concept system activation
    wS_C = (C_S'.*Sx(i,:))'; %s-c map weighted by syntactic system activation

    %% concept system phase
    dCp = eta_Cp*rand(1,Nc) + omega + psi_SC*sum(-sin(PHI_SC).*wS_C) + psi_CC*sum(-sin(PHI_CC).*C_C);

    %% syntactic system phase
    dSp = eta_Sp*rand(1,Ns) + omega + psi_SS*sum(-sin(PHI_SS).*S_S) + psi_CS*sum(-sin(PHI_CS).*wC_S);

    %% concept system activation
    Cx_eq = ix_env_C + eta_Cx*rand(1,Nc); %equilibria for conceptual systems

    dCx = -(Cx(i,:)-Cx_eq) + chi_CC*Cx(i,:)*C_C;

    %% syntactic system activation

    %global activation potential
    gw = gain_gw * (1./(Sx(i,:)-gw_lb)) - ...
        gain_gw * (1./(gw_ub-Sx(i,:)));

    %coupling with concepts
    c_coup =  chi_CS * sum(wC_S) .* (Fx(i,:) < tau_Fx_dem);

    %excitation potential
    VE_S = gain_VE * An(i) * (1/dt) * nansum(V_E(Sx(i,:),Ex(i,:)')); %#ok<NANSUM>

    dSx = eta_Sx*rand(1,Ns) + VE_S + c_coup + gw;

    %% syntactic system feedback
    dFx = Sx(i,:)>tau_Sx_sel;

    %% excitation monitor

    %cost of multiple occupation
    dgsx = (Sx(i,:).*Sx(i,:)') .* V_Me(Sx(i,:)-Sx(i,:)') .* ~eye(Ns);

    %cost of unoccupied level
    dgex = Ex(i,:).*min(sqrt((Ex(i,:)-Sx(i,:)').^2));

    dMe = Me_unocc_cost*nansum(dgex) + Me_multiocc_cost*sum(dgsx(:)) + Me_decay*(-1-Me(i)); %#ok<NANSUM>

    %% excitation organization

    %constant E organization steps
    dEx = 0;

    %% annealer
    dAn = (An(i)<1) * (init_eorg*gain_Ai + (~init_eorg)*gain_Ar);

    %% apply updates
    Cp(i+1,:) = mod(Cp(i,:) + dCp*dt,2*pi);
    Sp(i+1,:) = mod(Sp(i,:) + dSp*dt,2*pi);
    Cx(i+1,:) = Cx(i,:) + dCx*dt;
    Sx(i+1,:) = Sx(i,:) + dSx*dt;
    Ex(i+1,:) = Ex(i,:) + dEx*dt;
    Fx(i+1,:) = Fx(i,:) + dFx*dt;
    An(i+1) = An(i) + dAn*dt;
    Me(i+1) = Me(i) + dMe*dt;

    %%
    ix_sel = Sx(i+1,:)>1;
    ix_dem = ix_sel & Fx(i+1,:)>=tau_Fx_dem;
    ix_sup = ~ix_sel & Fx(i+1,:)>=tau_Fx_dem;

    %%
    % e-organization degeneracy before any reorganization
    if Me(i+1)>reset_Me && ~any(ix_sup)
        init_eorg = true;
        An(i+1) = 0;
        Sx(i+1,:) = Sx0 + eta_Sx0*rand(1,Ns);
        Me(i+1) = 0;
        Cx(i+1,:) = Cx(i+1,:)*eps_Cx;
        EV(end+1,:) = {t(i), 'deg_estate'};  %#ok<*AGROW>

    % e-organization degeneracy after re-organization has occurred
    elseif Me(i+1)>reset_Me && any(ix_sup)
        init_eorg = true;
        An(i+1) = 0;
        Me(i+1) = 0;

        Sx(i+1,ix_sup) = eps_Sx;
        Sx(i+1,~ix_sup) = (sum(~ix_sup):1)/sum(~ix_sup);

        Cx(i+1,:) = Cx(i+1,:)*eps_Cx;
        EV(end+1,:) = {t(i), 'deg_estate'};

    % sequence completion
    elseif Me(i+1)<reorg_Me && all(ix_sup)
        if stop_on_completion
            EV(end+1,:) = {t(i), 'seq_comp'}; %#ok<*NASGU>
            break;
        elseif ~strcmp(EV{end,2},'seq_comp')
            EV(end+1,:) = {t(i), 'seq_comp'};
        end

    % initial coherence
    elseif Me(i+1)<reorg_Me && ~any(ix_sel)
        init_eorg = false;
        An(i+1) = 0;
        Me(i+1) = 0;
        Ex(i+1,:) = e_reorg(Ex(i+1,:));
        EV(end+1,:) = {t(i), 'init_coherence'};
        %plot_sim(getvars(whos));

     % canonical reorganization
    elseif any(ix_dem)
        An(i+1) = 0;
        Me(i+1) = 0;
        Ex(i+1,:) = e_reorg(Ex(i+1,:));
        Sx(i+1,ix_dem) = eps_Sx;
        EV(end+1,:) = {t(i), 'can_reorg'};
        %plot_sim(getvars(whos));

    end

end

h = getvars(whos);

if nargin==0 && nargout==0
    plot_sim(h);
end

end

function [E] = e_reorg(E)

ix_first_nonground = find(E>0,1,'first');
E(ix_first_nonground) = nan;

end

function [h] = getvars(ff)
for i=1:length(ff)
    h.(ff(i).name) = evalin('caller',ff(i).name);
end
end




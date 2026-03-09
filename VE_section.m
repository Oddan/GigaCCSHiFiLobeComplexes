mrstModule add co2lab
gravity on

% load Grid
load('grids')

% convert units
rocks_mid.perm = rocks_mid.perm * milli * darcy; % from mD to Darcy
rocks_mid.permz = rocks_mid.permz * milli * darcy; % from mD to Darcy
rocks_mid.poro = rocks_mid.poro; % porosity is unitless


% Define a simple vertical equilibrium grid
[Gt, G, transMult, discarded_cells] = topSurfaceGrid(Gmid);


% loop over fields in rocks_mid and remove values corresponding to discarded
% cells
rock_fields = fieldnames(rocks_mid);
for k = 1:numel(rock_fields)
    field = rock_fields{k};
    rocks_mid.(field)(discarded_cells) = [];
end

rockVE = averageRock(rocks_mid, Gt)

% handle vanishing permeabilities at either endpoint
rockVE.perm(rockVE.perm(:,1) < 1e-15) = min(rockVE.perm(21:end-1));

%% run 3D simulation



%% run VE simulation

gravity reset on;
g       = gravity;
rhow    = 1000;                                 % water density (kg/m^3)
co2     = CO2props();                           % CO2 property functions
p_ref   = 30 *mega*Pascal;                      % reference pressure
t_ref   = 94+273.15;                            % reference temperature
co2_rho = co2.rho(p_ref, t_ref);                % CO2 density
co2_c   = co2.rhoDP(p_ref, t_ref) / co2_rho;    % CO2 compressibility
wat_c   = 0;                                    % water compressibility
c_rock  = 4.35e-5 / barsa;                      % rock compressibility 
srw     = 0.27;                                 % residual water
src     = 0.20;                                 % residual CO2
pe      = 5 * kilo * Pascal;                    % capillary entry pressure
muw     = 8e-4 * Pascal * second;               % brine viscosity
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity

% fluid
invPc3D = @(pc) (1-srw) .* (pe./max(pc, pe)).^2 + srw;
kr3D    = @(s) max((s-src)./(1-src), 0).^2; % uses CO2 saturation
fluid   = makeVEFluid(Gt, rockVE, ...
               'sharp_interface_simple'           , ...
               'co2_mu_ref'  , muco2                  , ...
               'wat_mu_ref'  , muw                    , ...
               'co2_rho_ref' , co2_rho                , ...
               'wat_rho_ref' , rhow                   , ...
               'co2_rho_pvt' , [co2_c, p_ref]         , ...
               'wat_rho_pvt' , [wat_c, p_ref]         , ...
               'residual'    , [srw, src]             , ...
               'pvMult_p_ref', p_ref                  , ...
               'pvMult_fac'  , c_rock                 , ...
               'invPc3D'     , invPc3D                , ...
               'kr3D'        , kr3D                   , ...
               'dissolution' , false, ... % set to 'true' if you want dissolution
               'dis_rate'    , Inf, ...
               'dis_max'     , 0.05, ...
               'transMult'   , transMult);

% initial condition
initState.pressure = rhow * norm(g) * Gt.cells.z;
initState.s        = repmat([1, 0], Gt.cells.num, 1);
initState.sGmax    = initState.s(:,2);
initState.rs       = zeros(Gt.cells.num, 1);

% well
W = addWellVE([], Gt, rockVE, 152, ...
    'type', 'rate', ...
    'val', 0.01, ...
    'comp_i', [0 1], ...
    'name', 'Injector');

% schedule

% hydrostatic pressure conditions for open boundary faces
bc = pside([], Gt, 'Xmin', initState.pressure(1), 'sat', [1, 0]);
bc = pside(bc, Gt, 'Xmax', initState.pressure(end), 'sat', [1, 0]);

% Setting up two copies of the well and boundary specifications. 
% Modifying the well in the second copy to have a zero flow rate.
schedule.control    = struct('W', W, 'bc', bc);
schedule.control(2) = struct('W', W, 'bc', bc);
schedule.control(2).W.val = 0;

% Specifying length of simulation timesteps
% schedule.step.val = [rampupTimesteps(2*year, year/24, 7); ...
%                      repmat(1 * year, 250, 1)];
schedule.step.val = [rampupTimesteps(1*year, year/24, 7); ...
                     rampupTimesteps(100*year, 5 * year, 7)];
                     
% Specifying which control to use for each timestep.
% The first 100 timesteps will use control 1, the last 100
% timesteps will use control 2.
schedule.step.control = [ones(31, 1); ...
                         ones(27, 1) * 2];
%schedule.step.control = [ones(55, 1)]
%ones(57, 1) * 2];
% schedule.step.control = [ones(79, 1)];
                         

% run simulation
model = CO2VEBlackOilTypeModel(Gt, rockVE, fluid);
model.verbose=true;
nls = NonLinearSolver('maxTimestepCuts', 12);
[wellSol, states] = simulateScheduleAD(initState, model, schedule, ...
                                       'NonLinearSolver', nls);
states = [{initState}; states];

mrstModule add mrst-gui
plotToolbar(Gt, states);

sat_VE3D = {};
% converting saturation result to 3D
for i = 1:numel(states)
    s = states{i}.s(:,2);
    smax = states{i}.sGmax;
    p = states{i}.pressure;
    [h, h_max] = upscaledSat2height(s, smax, Gt, ...
                                    'resSat', [srw, src], ...
                                    'poro', rocks_mid.poro, ...
                                    'rhoW', fluid.rhoW, ...
                                    'rhoG', fluid.rhoG, 'p', p);
    sat_VE3D = [sat_VE3D, {height2finescaleSat(h, h_max, Gt, ...
                                               fluid.res_water, fluid.res_gas, ...
                                               'rhoW', fluid.rhoW(p), ...
                                               'rhoG', fluid.rhoG(p))}];
end
% $$$ figure;
% $$$ plotCellData(G, rocks_mid.poro, 'EdgeColor', 'none');
% $$$ figure;
% $$$ plotCellData(G, rocks_mid.perm, 'EdgeColor', 'none');
figure;
plotToolbar(G, sat_VE3D);

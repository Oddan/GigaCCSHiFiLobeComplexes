mrstModule add ad-core ad-props
mrstModule add co2lab
gravity on

load('middle_grid_and_rock_compressed');

% convert units
rocks.perm = rocks.perm * milli * darcy; % from mD to Darcy
rocks.permz = rocks.permz * milli * darcy; % from mD to Darcy
rocks.poro = rocks.poro; % porosity is unitless

% Define a simple vertical equilibrium grid
[Gt, G, transMult, discarded_cells] = topSurfaceGrid(Gs, 'discard_below_holes', false);

rockVE = averageRock(rocks, Gt)


tsurf = 15 + 273.15; % surface temperature in K
tgrad = 30; % temperature gradient in K/km

rhow = 1050; % density of brine (kg/m^3)

pfun = @(z) 1 * atm + rhow * norm(gravity) * z;

initState.pressure = pfun(Gt.cells.z);
initState.temperature = tsurf + Gt.cells.z * tgrad / 1000; 
initState.s = repmat([1, 0], Gt.cells.num, 1); % all water, no CO2
initState.sGmax    = initState.s(:,2);
initState.rs       = zeros(Gt.cells.num, 1);


%% run VE simulation

gravity reset on;
g       = gravity;

t_ref = mean(initState.temperature);
p_ref = mean(initState.pressure);

co2     = CO2props();                           % CO2 property functions
co2_rho = co2.rho(p_ref, t_ref);                % CO2 density
co2_c   = co2.rhoDP(p_ref, t_ref) / co2_rho;    % CO2 compressibility
wat_c   = 4.3e-5/barsa;                         % water compressibility
c_rock  = 3.0e-5 / barsa;                       % rock compressibility 
srw     = 0.081;                                % residual water
src     = 0.40;                                  % residual CO2
pe      = 5 * kilo * Pascal;                    % capillary entry pressure
muw     = 8e-4 * Pascal * second;               % brine viscosity
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity
n = [2,2] ; % relperm exponents

% fluid
invPc3D = @(pc) (1-srw) .* (pe./max(pc, pe)).^2 + srw;
kr3D    = @(s) max((s-src)./(1-src), 0).^2; % uses CO2 saturation
fluid   = makeVEFluid(Gt, rockVE, ...
               'sharp_interface_integrated'           , ...
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
               'kr3D'        , n                      , ... %kr3D                   , ...
               'dissolution' , false, ... % set to 'true' if you want dissolution
               'dis_rate'    , Inf, ...
               'dis_max'     , 0.05, ...
               'transMult'   , transMult);

% well position, with logical coordiantes i=129, j=10
wcell = find(Gt.cells.ij(:,1) == 129 & Gt.cells.ij(:,2)==10);

inj_rate = 0.75 * mega * 1e3 / year / fluid.rhoGS;
inj_rate = inj_rate * 0.5; % @@ reduce injection rate to avoid convergence issues

W = addWellVE([], Gt, rockVE, wcell, ...
    'type', 'rate', ...
    'val', inj_rate, ...
    'comp_i', [0, 1], ...
              'name', 'Injector');

% schedule

% closed boundaries for now
bc = [];

% Setting up two copies of the well and boundary specifications. 
% Modifying the well in the second copy to have a zero flow rate.
schedule.control    = struct('W', W, 'bc', bc);
schedule.control(2) = struct('W', W, 'bc', bc);
schedule.control(2).W.val = 0;


% Specifying length of simulation timesteps
isteps = rampupTimesteps(25*year, year/2, 7);
msteps = rampupTimesteps(200*year, 5 * year, 7);
schedule.step.val = [isteps; msteps];

schedule.step.control = [ones(numel(isteps), 1); ...
                         2*ones(numel(msteps), 1)];

% run simulation
model = CO2VEBlackOilTypeModel(Gt, rockVE, fluid);
model.verbose=true;
nls = NonLinearSolver('maxTimestepCuts' , 10, ...
                      'maxIterations'     ,16);

[wellsol, states] = simulateScheduleAD(initState, model, schedule, ...
                                       'NonLinearSolver', nls);
states = [{initState}; states];

mrstModule add mrst-gui
plotToolbar(Gt, states); view(0,90);


states3D = VEstates23D(states, Gt, model.fluid, 'poro3D', rocks.poro);

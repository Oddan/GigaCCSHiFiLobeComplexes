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

% @@ Alternative rock  --- @@ Experimental, use only the upper part of the grid to compute rock properties for VE
%upper_ixs = find(Gs.cells.indexMap <= 234*182*130);
upper_ixs = find(Gs.cells.indexMap <= 234*182*(129+3));
GtAlt = topSurfaceGrid(extractSubgrid(Gs, upper_ixs), 'discard_below_holes', false);
rocksAlt.perm = rocks.perm(upper_ixs);
rocksAlt.permz = rocks.permz(upper_ixs);
rocksAlt.poro = rocks.poro(upper_ixs);
rockVEAlt = averageRock(rocksAlt, GtAlt);
rockVE = rockVEAlt; % @@

% Define temperature (from report)
tfun = @(z) 273.15 + 62 + 0.025 * (z-2100); % temperature function in K, with depth z in m

rhow = 1021; % density of brine (kg/m^3), from report

% pressure function in Pa, with depth z in m, reference depth is 2100 m
pfun = @(z) 180 * barsa + rhow * norm(gravity) * (z-2100); % from report

initState.pressure = pfun(Gt.cells.z);
initState.temperature = tfun(Gt.cells.z);
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
                                                %co2_c   = co2.rhoDP(p_ref, t_ref) / co2_rho;    % CO2 compressibility
wat_c   = 5e-5/barsa;                           % water compressibility
c_rock  = 1e-4 / barsa;                         % rock compressibility, confirmed with Trine

srw     = 0.2  ;                                % residual water
src     = 0.15;                                  % residual CO2
muw     = 5.946e-4; %@@8e-4 * Pascal * second;               % brine viscosity
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity (@@ note, value doesn't matter, since sampled table is used!)
prange = [0.1, 400] * mega * Pascal; % pressure range for PVT tables
trange = [  4, 250] + 274; % CO2 default temperature range for PVT tables, in K

% cap pressure and relperm

GSF = [0.000000 0.000000 0.000000
       0.150000 0.000000 0.000000
       0.215000 0.003686 0.018167
       0.280000 0.024204 0.022506
       0.345000 0.067386 0.028690
       0.410000 0.132803 0.037971
       0.475000 0.217969 0.052895
       0.540000 0.320682 0.079363
       0.605000 0.439928 0.133900
       0.670000 0.575816 0.279863
       0.735000 0.728993 0.986900
       0.800000 0.900000 64.931419];

% WSF = [0.200000 0.000000
%        0.330000 0.000800
%        0.395000 0.004050
%        0.460000 0.012800
%        0.525000 0.031250
%        0.590000 0.064800
%        0.655000 0.120050
%        0.720000 0.204800
%        0.785000 0.328050
%        0.850000 0.500000
%        1.000000 0.500000];

kr3D = @(sg) interp1(GSF(:,1), GSF(:,2), sg);
%kr3D = @(sg) 0.9 * interp1(GSF(:,1), GSF(:,2), sg); % @@@  what happens if we reduce relperm

swat = 1-GSF(:,1);
pcval = GSF(:,3) * barsa;
%swat(2) = [];
%pcval(2) = []; % avoid degenerate table (and smooth out the jump in saturation)
swat(1) = [];
pcval(1) = []; % avoid degenerate table (and smooth out the jump in saturation)
invPc3D = @(pc) interp1(pcval, swat, pc, 'linear', 'extrap');

%pe = 2.5e-2 * barsa;
%n = [2,2] ; % relperm exponents
%invPc3D = @(pc) (1-srw) .* (pe./max(pc, pe)).^2 + srw;
%kr3D    = @(s) max((s-src)./(1-src), 0).^2; % uses CO2 saturation

% fluid
fluid   = makeVEFluid(Gt, rockVE, ...
                      'P-scaled table'                       , ...
                      'hyst_model'  , 'fixed residual'       , ...
                      'co2_mu_ref'  , muco2                  , ... % doesn't matter, since sampled table is used
                      'co2_mu_pvt'  , [prange, trange]       , ...
                      'wat_mu_ref'  , muw                    , ...
                      'co2_rho_ref' , co2_rho                , ...
                      'wat_rho_ref' , rhow                   , ...
                      'wat_rho_pvt' , [wat_c, p_ref]         , ...
                      'residual'    , [srw, src]             , ...
                      'pvMult_p_ref', p_ref                  , ...
                      'pvMult_fac'  , c_rock                 , ...
                      'invPc3D'     , invPc3D                , ...
                      'kr3D'        , kr3D                   , ...
                      'reservoirT'  , tfun(Gt.cells.z),  ...
                      'dissolution' , false, ... % set to 'true' if you want dissolution
                      'dis_rate'    , Inf, ...
                      'dis_max'     , 0.05, ...
                      'transMult'   , transMult);

%fluid.muG = @(p) 2 * fluid.muG(p); % @@@
%'kr3D'        , n                      , ... %kr3D                   , ...
%'co2_rho_pvt' , [co2_c, p_ref]         , ... % default is sampled table
%'sharp_interface_integrated'           , ...

% well position, with logical coordiantes i=129, j=10
wcell = find(Gt.cells.ij(:,1) == 129 & Gt.cells.ij(:,2)==10);

inj_rate = 0.75 * mega * 1e3 / year / fluid.rhoGS;
%inj_rate = inj_rate * 0.5; % @@ reduce injection rate to avoid convergence issues

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
%msteps = rampupTimesteps(200*year, 5 * year, 7);
msteps = rampupTimesteps(200*year, 2 * year, 7);

% schedule.step.val = [isteps]; % @@
% schedule.step.control = [ones(numel(isteps), 1)]; % @@ 

%isteps = isteps(1:10); % @@@

schedule.step.val = [isteps; msteps];
schedule.step.control = [ones(numel(isteps), 1); ...
                         2*ones(numel(msteps), 1)];

% schedule.step.val = [isteps; msteps(1:7)];
% schedule.step.control = [ones(numel(isteps), 1); ...
%                          2*ones(7, 1)];


% run simulation
model = CO2VEBlackOilTypeModel(Gt, rockVE, fluid);
model.verbose=true;
nls = NonLinearSolver('maxTimestepCuts' , 10, ...
                      'maxIterations'     ,16);

[wellsol, states] = simulateScheduleAD(initState, model, schedule, ...
                                       'NonLinearSolver', nls);
states = [{initState}; states];

mrstModule add mrst-gui
%plotToolbar(Gt, states); view(0,90);


tic; states3D = VEstates23D(states, Gt, model.fluid, 'poro3D', rocks.poro);toc
colormap(custom_jet_colormap(256))
plotToolbar(Gt.parent, states3D); view(0,90);

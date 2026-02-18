%% 3D simulation of injection at bottom of aquifer, to assess
%% timescale of vertical equilibrium.

mrstModule add co2lab
mrstModule add mrst-gui

gravity on
do_plot=false;

[G, rock, Gt, rockVE, transMult] = setup_testgrid();



%% determine wellcells (mid-bottom)
tmp = zeros(G.cartDims);
tmp(G.cells.indexMap) = 1:G.cells.num;

well_cells = unique(tmp(130,:,206:210));
well_cells = well_cells(well_cells > 0);

if do_plot
    field = zeros(G.cells.num,1);
    field(well_cells) = 1;
    plotCellData(G, field, 'edgecolor', 'k', 'edgealpha', 0.5); colorbar; view(0,0);
    title('well cells');
end

%% Set up initial state

tsurf = 15 + 273.15; % surface temperature in K
tgrad = 30; % temperature gradient in K/km

rhow = 1050; % density of brine (kg/m^3)

pfun = @(z) 1 * atm + rhow * norm(gravity) * z;

initState.pressure = pfun(G.cells.centroids(:,3));
initState.temperature = tsurf + G.cells.centroids(:,3) * tgrad / 1000; 
initState.s = repmat([1, 0], G.cells.num, 1); % all water, no CO2

%% set up fluid
t_ref = mean(initState.temperature);
p_ref = mean(initState.pressure);

co2props = CO2props();
muw = 8e-4 * Pascal * second; % brine viscosity;
muc = co2props.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity

rhoc = co2props.rho(p_ref, t_ref); % density of co2 (kg/m^3)
cf_wat = 0; % water compressibility
cf_co2 = 0; % co2 compressibility
cf_rock = 4.35e-5 / barsa; % rock compressibility
srw = 0.1;
src = 0.1;
%n = [1,1]; % @@ relperm exponents
n = [2,2]; % @@ relperm exponents

fluid = initSimpleADIFluid('phases', 'WG'           , ...
                           'mu'  , [muw, muc]       , ...
                           'rho' , [rhow, rhoc]     , ...
                           'pRef', p_ref            , ...
                           'c'   , [cf_wat, cf_co2] , ...
                           'cR'  , cf_rock          , ...
                           'smin', [srw, src]       , ...
                           'n'   , n); 

%% setup model

% @@ hack: avoid 0 ntg
rock_modif = rock;
rock_modif.ntg = max(rock_modif.ntg, 0.05);
rock_modif.perm = [rock.perm, rock.perm, rock.permz];
%rock_modif.perm = rock_modif.permz; % "Worst case"

if do_plot
    figure;plotCellData(G, rock_modif.perm/darcy, 'edgealpha', 0.1); view(0,0);
    colorbar; title('permeability field');
end

model = TwoPhaseWaterGasModel(G, rock_modif, fluid, tsurf, tgrad);
model.verbose=true;

%% setup schedule
inj_rate = 500 * meter^3/day; %1000 * meter^3/day;
W = addWell([], G, rock_modif, well_cells, 'Type', 'rate', 'Val', inj_rate, 'comp_i', [0, 1], 'name', 'Injector');
open_boundaries = true;
bc = open_boundary_conditions(G, pfun, open_boundaries);
inj_period = 4*year;%60*day;
inj_steps = 400; %40;
migr_period = 10*year;
migr_steps = 100; %100;

schedule = simple_injection_migration_schedule(W, bc, inj_period, inj_steps, ...
                                               migr_period, migr_steps);
%% run simulation

nls = NonLinearSolver('maxTimestepCuts', 14, ...
                      'timeStepSelector', ...
                      IterationCountTimeStepSelector('targetIterationCount', 12));


[wellsol, states] = simulateScheduleAD(initState, model, schedule, ...
                                       'NonLinearSolver', nls);

%% Visualize

% add rock parameters to states for plotting
plotstates = [{initState}; states];
for state_ix = 1:length(plotstates)
    plotstates{state_ix}.poro = rock_modif.poro;
    plotstates{state_ix}.perm = rock_modif.perm;
    plotstates{state_ix}.ntg = rock_modif.ntg;
end

if do_plot
    plotToolbar(G, plotstates); view(0,0);
end

%% running VE simulation

wcell = 130;
field = zeros(Gt.cells.num,1);
field(wcell) = 1;
plotCellData(Gt, field, 'edgecolor', 'k', 'edgealpha', 0.5); colorbar; view(0,90);
title('well cell in VE grid');

%% set up initial state
initStateVE.pressure = pfun(Gt.cells.z);
initStateVE.temperature = tsurf + Gt.cells.z * tgrad / 1000; 
initStateVE.s = repmat([1, 0], Gt.cells.num, 1); % all water, no CO2
initStateVE.sGmax = initStateVE.s(:,2);
initStateVE.rs = zeros(Gt.cells.num, 1);

%% setup fluid
fluidVE = makeVEFluid(Gt, rockVE                     , ...
                      'sharp_interface_integrated'   , ...
                      'co2_mu_ref'  , muc            , ...
                      'wat_mu_ref'  , muw            , ...
                      'co2_rho_ref' , rhoc           , ...
                      'wat_rho_ref' , rhow           , ...
                      'co2_rho_pvt' , [cf_co2, p_ref] , ...
                      'wat_rho_pvt' , [cf_wat, p_ref] , ...
                      'residual'    , [srw, src]     , ...
                      'dissolution' , false          , ... 
                      'dis_rate'    , Inf            , ...
                      'dis_max'     , 0.05           , ...
                      'pvMult_p_ref', p_ref          , ...
                      'pvMult_fac'  , cf_rock        , ...
                      'kr3D'        , n              , ...
                      'transMult'   , transMult);
                      %'invPc3D'     , invPc3D        , ...
                      
%% setup model
modelVE = CO2VEBlackOilTypeModel(Gt, rockVE, fluidVE);
modelVE.verbose=true;

%% setup schedule
WVE = addWellVE([], Gt, rockVE, wcell, 'Type', 'rate', ...
                'Val', inj_rate, 'comp_i', [0, 1], 'name', 'Injector');

% schedule
bcVE = open_boundary_conditions(Gt, pfun, open_boundaries);
scheduleVE = simple_injection_migration_schedule(WVE, bcVE, inj_period, inj_steps, ...
                                               migr_period, migr_steps);

%% run simulation
[wellsolVE, statesVE] = simulateScheduleAD(initStateVE, modelVE, scheduleVE, ...
                                           'NonLinearSolver', nls);

plotstatesVE = [{initStateVE}; statesVE];
sat_VE3D = {};
% converting saturation results to 3D
for i = 1:numel(plotstatesVE)
    s = plotstatesVE{i}.s(:,2);
    smax = plotstatesVE{i}.sGmax;
    p = plotstatesVE{i}.pressure;
    [h, h_max] = upscaledSat2height(s, smax, Gt, ...
                                    'resSat', [srw, src], ...
                                    'poro', rock.poro, ...
                                    'rhoW', fluidVE.rhoW, ...
                                    'rhoG', fluidVE.rhoG, 'p', p);
    sat_VE3D = [sat_VE3D, {height2finescaleSat(h, h_max, Gt, ...
                                               fluidVE.res_water, fluidVE.res_gas, ...
                                               'rhoW', fluidVE.rhoW(p), ...
                                               'rhoG', fluidVE.rhoG(p))}];
end
figure;
plotToolbar(G, sat_VE3D); view(0,0);

%% investigate
% - impact of transmult on correspondence with 3D
% - impact of fluid model (heterogeneity or not)
% - impact of corey exponents in relperm
% - are the pressures comparable

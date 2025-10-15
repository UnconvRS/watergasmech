%% Example: Poroelasticity simulation applied to the Norne case.
% 
% The simulation options are gathered in the opt structure. If opt=[] the
% simulation is run with the default options defined below
%
% **  Summary of the options ** 
%
% option 'norne_case' :
%
%     * 'full'       : 7392 cells
%     * 'mini Norne' :  605 cells
%
% option 'bc_case' :
%
%     * 'no displacement' : All nodes belonging to external faces have displacement
%                           equal to zero
%     * 'bottom fixed'    : The nodes that belong to the bottom have zero
%                           displacement, while a given pressure is imposed on
%                           the external faces that are not bottom faces.
%
% option 'method' :
%
%     * 'fully coupled'          : The mechanical and flow equations are solved fully coupled.
%     * 'fixed stress splitting' : The mechanical and flow equations are solved
%                                  sequentially using a fixed stress splitting
%
% option 'fluid_model' :
%
%     * 'blackoil'  : blackoil model is used for the fluid (gas is injected, see
%                     schedule below)
%     * 'oil water' : Two phase oil-water
%     * 'water'     : water model is used for the fluid

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

clc;clear;close all;
mrstModule add hfm_v2 co2lab-common co2lab-ve co2lab-spillpoint ad-mechanics_v3 ad-core ad-props ad-blackoil vemmech deckformat mrst-gui
mrstModule add coarsegrid    

%% setup default option values


% setup default option values
opt = struct('nonlinearTolerance' , 1e-6            , ...
             'splittingTolerance' , 1e-3            , ...
             'verbose'            , false           , ...
             'splittingVerbose'   , false);  


%% Load Johansen model

grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
fullCartDims = grdecl.cartDims;
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);
grdecl = cutGrdecl(grdecl, [10 25; 35 55; 1 22]); 

grdecl.ACTNUM = ones(size(grdecl.ACTNUM));
G = processGRDECL(grdecl);
G = G(1);
G = computeGeometry(G);

figure,
plotGrid(G, 'facecolor', 'none');  
axis tight

% Perm and poro values

perm = [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ];
rock.perm = perm(G.cells.indexMap, :);
rock.poro = max(grdecl.PORO(G.cells.indexMap), 0.1);


%% Grab all the boundary faces in the reservoir domain

nx = G.cartDims(1); ny = G.cartDims(2); nz = G.cartDims(3);

backdir = searchForBoundaryFaces(G, 'BACK');
leftdir = searchForBoundaryFaces(G, 'LEFT');
rightdir = searchForBoundaryFaces(G, 'RIGHT');
frontdir = searchForBoundaryFaces(G, 'FRONT'); 

% figure;
% plotGrid(G, 'faceColor', 'none')
% plotFaces(G, backdir, 'b');
% plotFaces(G, frontdir, 'b');
% plotFaces(G, leftdir, 'y');
% plotFaces(G, rightdir, 'y');
% plotFaces(G, topfaces, 'r');
% view(63, 47), axis tight off;


% title('3D grid');


%% Initial state
gravity on; % tell MRST to turn on gravity
g = gravity; % get the gravity vector
rhow = 1000; % density of brine corresponding to 94 degrees C and 300 bar
initState.pressure = rhow * g(3) * G.cells.centroids(:,3);
initState.s = repmat([1, 0], G.cells.num, 1);
initState.sGmax = initState.s(:,2);

%% Fluid model
co2     = CO2props(); % load sampled tables of co2 fluid properties
pRef   = 30 * mega * Pascal; % choose reference pressure
t_ref   = 94 + 273.15; % choose reference temperature, in Kelvin
rhoc    = co2.rho(pRef, t_ref); % co2 density at ref. press/temp
cf_co2  = co2.rhoDP(pRef, t_ref) / rhoc; % co2 compressibility
cf_wat  = 0; % brine compressibility (zero)
cf_rock = 4.35e-5 / barsa; % rock compressibility
muw     = 8e-4 * Pascal * second; % brine viscosity
muco2   = co2.mu(pRef, t_ref) * Pascal * second; % co2 viscosity

mrstModule add ad-props; % The module where initSimpleADIFluid is found

% Use function 'initSimpleADIFluid' to make a simple fluid object
fluid = initSimpleADIFluid('phases', 'WG'           , ...
                           'mu'  , [muw, muco2]     , ...
                           'rho' , [rhow, rhoc]     , ...
                           'pRef', pRef            , ...
                           'c'   , [cf_wat, cf_co2] , ...
                           'cR'  , cf_rock          , ...
                           'n'   , [2 2]);                     

p_sc = 1*atm;
T_sc = 288.706;
rhoGSinj= co2.rho(p_sc, T_sc);
fluid.rhoSinj = [rhow,rhoGSinj];

% Change relperm curves
srw = 0.27;
src = 0.20;
fluid.krW = @(s) fluid.krW(max((s-srw)./(1-srw), 0));
fluid.krG = @(s) fluid.krG(max((s-src)./(1-src), 0));

% % Plot relperm curves
% figure; hold on;
% sw = linspace(srw, 1, 200);
% plot(sw, fluid.krW(sw),   'b', 'linewidth', 1.5);
% plot(sw, fluid.krG(1-sw), 'r', 'linewidth', 1.5);
% line([srw, srw], [0 1], 'color', 'black', 'linestyle', ':', 'linewidth', 1);
% xlabel('brine saturation'); ylabel('relative permeability');
% set(gca, 'fontsize', 14, 'xlim', [0 1]);


% Add capillary pressure curve
pe = 5 * kilo * Pascal;
pcWG = @(sw) pe * sw.^(-1/2);
fluid.pcWG = @(sg) pcWG(max((1-sg-srw)./(1-srw), 1e-5));

%% Setup material parameters for Biot and mechanics

E          = 1 * giga * Pascal; % Young's module
nu         = 0.3;               % Poisson's ratio
alpha      = 1;                 % Biot's coefficient

% Transform these global properties (uniform) to cell values.
E          = repmat(E, G.cells.num, 1);
nu         = repmat(nu, G.cells.num, 1);
rock.alpha = repmat(alpha, G.cells.num, 1);

%% Setup boundary conditions for mechanics (no displacement)

nx = G.cartDims(1);
ny = G.cartDims(2);
nz = G.cartDims(3);

% Find the bottom nodes. On these nodes, went impose zero displaceme

c = zeros(nx*ny*nz, 1);
c(G.cells.indexMap) = (1 : numel(G.cells.indexMap))';
bottomcells = c(nx*ny*(nz - 1) +  (1 : (nx*ny))');
bottomcells = bottomcells(bottomcells > 0); % Ensure no invalid indices

facesbottom = G.cells.faces(mcolon(G.cells.facePos(bottomcells), G.cells.facePos(bottomcells ...
                                                  + 1) - 1), :);
bottomfaces = facesbottom(facesbottom(:, 2) == 6  , 1);


indbottom_nodes = mcolon(G.faces.nodePos(bottomfaces), G.faces.nodePos(bottomfaces ...
                                                  + 1) - 1);
bottom_nodes = G.faces.nodes(indbottom_nodes);
isbottom_node = false(G.nodes.num, 1);
isbottom_node(bottom_nodes) = true;
bcnodes = find(isbottom_node);

nn = numel(bcnodes);
u = zeros(nn, 3);
m = ones(nn, 3);

disp_bc = struct('nodes', bcnodes, 'uu', u, 'mask', m);

% Find outer faces that are not at the bottom. On these faces, we impose
% a given pressure.

topcells = c(1 : nx*ny);  % indices of the topmost layer
topcells = topcells(topcells > 0);  % keep only active cells

facestop = G.cells.faces(mcolon(G.cells.facePos(topcells), G.cells.facePos(topcells ...
                                                  + 1) - 1), :);
topfaces = facestop(facestop(:, 2) == 5  , 1);

figure;
plotGrid(G, 'faceColor', 'none')
plotFaces(G, topfaces, 'r');
view(63, 47); axis tight off;

rightfaces = rightdir; leftfaces = leftdir; frontfaces = frontdir; backfaces = backdir;

x_dirfaces = [rightfaces; leftfaces];

y_dirfaces = [frontfaces; backfaces];

top_pressure = 2.2.*pRef; xdir_pressure = 1.9.*pRef; ydir_pressure = 1.6.*pRef; 

signcoefx = (G.faces.neighbors(x_dirfaces, 1) == 0) - (G.faces.neighbors(x_dirfaces, 2) == 0);

signcoefy = (G.faces.neighbors(y_dirfaces, 1) == 0) - (G.faces.neighbors(y_dirfaces, 2) == 0);

signcoeftop = (G.faces.neighbors(topfaces, 1) == 0) - (G.faces.neighbors(topfaces, 2) == 0);

n_xdir = bsxfun(@times, G.faces.normals(x_dirfaces, :), signcoefx./ G.faces.areas(x_dirfaces));
force_xdir = bsxfun(@times, n_xdir, xdir_pressure);

n_ydir = bsxfun(@times, G.faces.normals(y_dirfaces, :), signcoefy./ G.faces.areas(y_dirfaces));
force_ydir = bsxfun(@times, n_ydir, ydir_pressure);

n_top = bsxfun(@times, G.faces.normals(topfaces, :), signcoeftop./G.faces.areas(topfaces));
force_top = bsxfun(@times, n_top, top_pressure);

forceall = [force_top; force_xdir; force_ydir];
outerfaces = [topfaces; x_dirfaces; y_dirfaces];

force_bc = struct('faces', outerfaces, 'force', forceall);

el_bc = struct('disp_bc' , disp_bc, 'force_bc', force_bc);


%% Setup load for mechanics

% In this example we do not impose any volumetric force
loadfun = @(x) (0*x);

%% Gather all the mechanical parameters in a struct

mech = struct('E', E, 'nu', nu, 'el_bc', el_bc, 'load', loadfun);


%% Gravity
% The gravity in this option affects only the fluid behavior
gravity on;


    %% Setup model

 model = MechWaterGasModel(G, rock, fluid, mech, 'verbose', false);

    %% Setup wells

W = [];

wc_global = false(G.cartDims);
wc_global(8, 10, 10:12) = true;
wc = find(wc_global(G.cells.indexMap));


W = addWell(W, G, rock, wc, ...
            'Type'    , 'rate', ...
            'Val'     , 1.4e6/day, ...
            'Sign'    , 1,  ...
            'Comp_i'   , [0, 1], ... % inject gas
            'Name'    , 'inj');

% figure,
% plotGrid(G, 'facecolor', 'none');
% plotWell(G, W, 'color', 'r', 'fontSize', 18);
% colormap(jet);
% colorbar;
% axis tight
% view(63, 47)

facilityModel = FacilityModel(model.fluidModel);
facilityModel = facilityModel.setupWells(W);
model.FacilityModel = facilityModel; 

    %% Setup schedule 

schedule.step.val     = [1*day*ones(1, 1); 5*day*ones(20, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control      = struct('W', W);

initState   = computeInitDisp(model, initState, [], 'pressure', initState.pressure);

    %% Simulate
    
[wellSols, states, schedulereport] = simulateScheduleAD(initState, model, schedule);

%Plot results 
figure,
plotWellSols(wellSols, cumsum(schedule.step.val)); 

for i=1:numel(states)
    states{i}.sThresh = zeros(G.cells.num, 1);
    states{i}.sThresh(states{i}.s(:,2) > 1e-2) = 1.0;
    states{i}.Kx = rock.perm(:, 1)/darcy;
    states{i}.Ky = rock.perm(:, 2)/darcy;
    states{i}.Kz = rock.perm(:, 3)/darcy;
    states{i}.sG = zeros(G.cells.num, 1);
    states{i}.sG(states{i}.s(:,2) > 1e-4) = 1.0;
end

figure;
plotToolbar(G, states, 'outline', true); view(-4, 32); axis tight;
colormap jet;



sat_end = states{end}.s(:,2);  % co2 saturation at end of simulation

plume_cells3 = sat_end > 0.0001;

figure,
plotGrid(G, 'facecolor', 'none','edgeAlpha',0.1);  % plot outline of simulation grid
plotGrid(G, plume_cells3, 'facecolor', 'red'); % plot cells with CO2 in red
plotWellData(G, W,  'TextColor', 'blue', 'color', [0 0 1]);
axis off
view(-75,30)


% %% Calculation of Efficiency
% 
% leakedAmntSC=[];
% 
% for i=1:numel(states)
%     sat_cur = states{i}.s(topLayer,2);
%     p_cur = states{i}.pressure(topLayer);
%     rhoco2RC=co2.rho(p_cur, t_ref);
%     Amntco2Leaked=sum(G.cells.volumes(topLayer).*rock.poro(topLayer).*sat_cur.*rhoco2RC)/fluid.rhoSinj(end);
%     leakedAmntSC=[leakedAmntSC;Amntco2Leaked];
% end
% 
% figure,
% plot(cumsum(schedule.step.val)/year,leakedAmntSC);
% xlabel('Time, years');
% ylabel('Total leaked CO_{2}, s-m^{3}')

% Plot amount of injected gas

injectedCO2SC=[];
for i=1:numel(states)
    curInjection=wellSols{i}.qGs;  
    injectedCO2SC=[injectedCO2SC;curInjection];
end

cumeCO2injectedSC=cumtrapz(cumsum(schedule.step.val),injectedCO2SC);

figure;
plot(cumsum(schedule.step.val)/day,cumeCO2injectedSC);
xlabel('Time, years');
ylabel('Amount Injected, s-m3')
title('Volume at surface condition')

Bg = rhoGSinj/rhoc;
cummCO2injres = cumeCO2injectedSC .* Bg;

figure;
plot(cumsum(schedule.step.val)/day,cummCO2injres);
xlabel('Time, years');
ylabel('Amount Injected, r-m3')
title('Volume at reservoir condition')














function [eqs, names, types, state] = equationsWaterGasMech(p0, sW0, state0, ...
                                                            p, sW, wellVars, ...
                                                            state, model, ...
                                                            dt, mechTerm, ...
                                                            drivingForces, varargin)
%Undocumented Utility Function

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

    opt = struct('Verbose'     , mrstVerbose , ...
                 'reverseMode' , false       , ...
                 'resOnly'     , false       , ...
                 'iteration'   , -1          , ...
                 'stepOptions' , []); % compatibility only
    opt = merge_options(opt, varargin{:});
    
    W  = drivingForces.W;
    bc = drivingForces.bc;
    s  = model.operators;
    f  = model.fluid;
    G  = model.G;
    t  = model.t;
    rock = model.rock;

%     % Extract current and previous values of all variables to solve for
%     [p, sG, wellSol] = model.getProps(state, 'pressure', 'sg', 'wellsol');
%     [p0, sG0, wellSol0] = model.getProps(state0, 'pressure', 'sg', 'wellsol');
%     
%     [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);
% 
%     % ------------------ Initialization of independent variables ------------------
%     
%     if ~opt.resOnly
%         if ~opt.reverseMode
%             [p, sG, wellVars{:}] = model.AutoDiffBackend.initVariablesAD(p, sG, wellVars{:});
%         else
%             wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
%             [p0, sG0, wellVars0{:}] = model.AutoDiffBackend.initVariablesAD(p0, sw0, wellVars0{:}); %#ok
%         end
%     end

    sG = 1-sW;
    sG0 = 1 - sW0;

    satstate = state;
    satstate.s = {sW, sG};
    relperm = model.getProps(satstate, 'RelativePermeability');
    [krW, krG] = deal(relperm{:});

%     [krW, krG] = model.evaluateRelPerm({sW, sG});

    % ----------------------------------------------------------------------------
    
    % Check for p-dependent tran mult:
    trMult = 1;
    if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end;
%     
%     % Check for p-dependent porv mult:
%     pvMult = 1; pvMult0 = 1;
%     if isfield(f, 'pvMultR')
%         pvMult  = f.pvMultR(p);
%         pvMult0 = f.pvMultR(p0);
%     end
%     transMult = 1;
%     if isfield(f, 'transMult')
%         transMult = f.transMult(p);
%     end

%     krW = mobMult.*krW; krG = mobMult.*krG;

    % Multipliers for properties
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

    % Modifiy relperm by mobility multiplier (if any)
    krW = mobMult.*krW; krG = mobMult.*krG;

    trans = s.T .* transMult;
    
    % Check for capillary pressure
    pcWG = 0;
    if isfield(f, 'pcWG')
        pcWG = f.pcWG(sG);
    end
    % ----------------------------------------------------------------------------
%     sW = 1-sG;
%     [krW, krG] = model.evaluateRelPerm({sW, sG});
    
    % computing densities, mobilities and upstream indices
    [bW, mobW, fluxW, vW, upcw] = compMFlux(p       , f.bW, f.muW, f.rhoWS, trMult, krW, s, trans, model);
    [bG, mobG, fluxG, vG, upcg] = compMFlux(p + pcWG, f.bG, f.muG, f.rhoGS, trMult, krG, s, trans, model);
    
    if all(isfinite(t))
        bW0 = f.bW(p0, t);
        bG0 = f.bG(p0 + pcWG, t);
    else
        bW0 = f.bW(p0);
        bG0 = f.bG(p0 + pcWG);
    end
    
    % ----------------------------------------------------------------------------
    % **Mechanical Effects on Porosity**
    % Compute effective pore volume considering mechanical changes
    effPorVol  = rock.poro .* (G.cells.volumes .* pvMult)  + rock.alpha .* mechTerm.new;
    effPorVol0 = rock.poro .* (G.cells.volumes .* pvMult0) + rock.alpha .* mechTerm.old;

    % ----------------------------------------------------------------------------
    % **Modified Mass Conservation Equations**
    
    % Conservation of mass for water (updated with effective porosity)
    water = (1/dt) * (effPorVol .* bW .* sW - effPorVol0 .* bW0 .* sW0) + s.Div(fluxW);
    
    % Conservation of mass for gas (updated with effective porosity)
    gas = (1/dt) * (effPorVol .* bG .* sG - effPorVol0 .* bG0 .* sG0) + s.Div(fluxG);

    eqs = {water, gas};  % List of residual equations
    names = {'water', 'gas'};  % Equation names
    types = {'cell', 'cell'};  % Equation types (cell-based mass balance)
    
    % ---------------------------- Boundary conditions ----------------------------
    if model.outputFluxes
        state = model.storeFluxes(state, vW, [], vG);
    end
    if model.extraStateOutput
        state = model.storebfactors(state, bW, [], bG);
        state = model.storeMobilities(state, mobW, [], mobG);
        state = model.storeUpstreamIndices(state, upcw, [], upcg);
    end
    
    rho = {bW.*f.rhoWS, bG.*f.rhoGS};
    mob = {mobW, mobG};
    sat = {sW, sG};
    if ~isempty(drivingForces.bc) && isempty(drivingForces.bc.sat)
       drivingForces.bc.sat = repmat([1 0], numel(drivingForces.bc.face), 1);
    end

%     primaryVars = {'pressure' , 'sG', wellVarNames{:}};
%     types = {'cell'           , 'cell'};
%     names = {'water'          , 'gas'};

    [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                     {p, p + pcWG}, sat, mob, rho, ...
                                                                     {}, {}, ...
                                                                     drivingForces);
    
    % ------------------------------ Well equations ------------------------------

    wellSol = model.getProp(state, 'wellsol');

    [~, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

    wellSol0 = model.getProp(state0, 'wellsol');

    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, {}, dt, opt);

    % ----------------------------------------------------------------------------
%     problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt); 
end

% ============================= END MAIN FUNCTION =============================

% ----------------------------------------------------------------------------
function [b, mob, fluxS, fluxR, upc] = compMFlux(p, bfun, mufun, rhoS, trMult, kr, s, trans, model)
    if all(isfinite(model.t))
        b   = bfun(p, model.t);
        mob = trMult .* kr ./ mufun(p, model.t);
    else
        b   = bfun(p);
        mob = trMult .* kr ./ mufun(p);
    end

    rhoFace = s.faceAvg(b*rhoS);
    
    dp   = s.Grad(p) - rhoFace .* model.getGravityGradient();
    upc  = (value(dp)<=0);
    fluxR = -s.faceUpstr(upc, mob) .* trans .* dp;
    fluxS = s.faceUpstr(upc, b) .* fluxR;
end

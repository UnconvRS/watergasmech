classdef MechWaterGasModel < MechFluidModel
%
% SYNOPSIS:
%   model = MechWaterGasModel(G, rock, fluid, mech_problem, ...)
%
% DESCRIPTION:
%   Model for fully coupled mechanical fluid simulation. The fluid
%   model is a two-phase water-gas model.
%
% PARAMETERS:
%   G            - Grid structure
%   rock         - Rock properties
%   fluid        - Fluid properties
%   mech_problem - Structure containing mechanical parameters
%
% RETURNS:
%   Class instance
%
% SEE ALSO:
%   MechBlackOilModel, MechWaterModel, MechFluidModel

    methods
        function model = MechWaterGasModel(G, rock, fluid, mech_problem, varargin)
            % Initialize the MechWaterGasModel as a subclass of MechFluidModel
            model = model@MechFluidModel(G, rock, fluid, mech_problem, varargin{:});
        end

        function fluidModel = setupFluidModel(model)
            % Setup the fluid model as a Water-Gas model
            fluidModel = GasWaterFluidModel(model.G, model.rock, model.fluid);
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false, ...
                         'resOnly', false, ...
                         'iteration', -1);  % Compatibility only
            opt = merge_options(opt, varargin{:});

            % Extract variables at the current timestep
            [p, sW, wellSol, xd] = model.getProps(state, 'pressure', ...
                                                         'water', ...
                                                         'wellSol', ...
                                                         'xd');

            % Extract variables at the previous timestep
            [p0, sW0, xd0] = model.getProps(state0, 'pressure', ...
                                                    'water', ...
                                                    'xd');

            % Initialize primary variables
            fluidModel = model.fluidModel; % Shortcut for the fluid model
            mechModel  = model.mechModel;  % Shortcut for the mechanical model

            [wellVars, wellVarNames, wellMap] = fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);

            if ~opt.resOnly
                % Define primary variables and initialize them
                [p, sW, wellVars{:}, xd] = initVariablesADI(p, sW, wellVars{:}, xd);
            end

            % Compute mechanical-fluid coupling terms
            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd, ...
                                                      state0.u, state.u);

            % Call the modified Water-Gas equations with mechanical effects
            [wg_eqs, wg_eqsnames, wg_eqstypes, state] = equationsWaterGasMech( ...
                                                  p0, sW0, state0, ...
                                                  p, sW, wellVars, ...
                                                  state, fluidModel, ...
                                                  dt, mechTerm, ...
                                                  drivingForces, ...
                                                  'iteration', ...
                                                  opt.iteration);

            % Get mechanical equations
            [mech_eqs, mech_eqsnames, mech_eqstypes] = equationsPoroMechanics(xd, ...
                                                              mechModel, ...
                                                              fluidp);

            % Combine fluid and mechanical equations
            eqs = horzcat(wg_eqs, mech_eqs);
            names = {wg_eqsnames{:}, mech_eqsnames{:}};
            types = {wg_eqstypes{:}, mech_eqstypes{:}};

            % 🔹 Correcting Primary Variables (Gas Pressure and Water Saturation)
            primaryVars = {'pressure', 'sW', wellVarNames{:}, 'xd'};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end

        function [problem, state] = getAdjointEquations(model, state0, state, dt, drivingForces, varargin)
            opt = struct('Verbose'    , mrstVerbose, ...
                         'reverseMode', true       , ...
                         'resOnly'    , false      , ...
                         'iteration'  , -1);  % Compatibility only
            opt = merge_options(opt, varargin{:});

            % Extract properties at the current timestep
            [p, sW, wellSol, xd] = model.getProps(state, 'pressure', ...
                                                         'water', ...
                                                         'wellSol', ...
                                                         'xd');

            % Extract properties at the previous timestep
            [p0, sW0, wellSol0, xd0] = model.getProps(state0, 'pressure', ...
                                                              'water', ...
                                                              'wellSol', ...
                                                              'xd');

            fluidModel = model.fluidModel; % Shortcut
            mechModel  = model.mechModel;  % Shortcut

            [wellVars, wellVarNames, wellMap] = fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);
            [wellVars0, ~, ~] = fluidModel.FacilityModel.getAllPrimaryVariables(wellSol0);

            if opt.reverseMode
                [p0, sW0, wellVars0{:}, xd0] = initVariablesADI(p0, sW0, wellVars0{:}, xd0);
            else
                [p, sW, wellVars{:}, xd] = initVariablesADI(p, sW, wellVars{:}, xd);
            end

            % Compute mechanical-fluid coupling terms
            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd, ...
                                                      state0.u, state.u);

            % Compute water-gas equations with mechanical effects
            [wg_eqs, wg_eqsnames, wg_eqstypes, state] = equationsWaterGasMech( ...
                                                  p0, sW0, state0, ...
                                                  p, sW, wellVars, ...
                                                  state, fluidModel, ...
                                                  dt, mechTerm, ...
                                                  drivingForces, ...
                                                  'iteration', ...
                                                  opt.iteration);

            % Compute mechanical equations
            [mech_eqs, mech_eqsnames, mech_eqstypes] = equationsPoroMechanics(xd, ...
                                                              mechModel, ...
                                                              fluidp);

            eqs = horzcat(wg_eqs, mech_eqs);
            names = {wg_eqsnames{:}, mech_eqsnames{:}};
            types = {wg_eqstypes{:}, mech_eqstypes{:}};

            if opt.reverseMode
                primaryVars = {'pressure', 'sW', 'xd'};
            else
                primaryVars = {'pressure', 'sW', 'xd', wellVarNames{:}};
            end

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end
        
        function [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd, u0, u)
            G = model.G;
            op = model.mechModel.operators;
            fluidp = p;  % Using pressure as fluid pressure
            
            if isa(xd, 'ADI')
                u = double2ADI(u, xd);
                u(~op.isdirdofs) = xd;  % Ensure correct derivatives
            end
            
            mechTerm.new = (op.ovol_div * u) ./ G.cells.volumes;
            mechTerm.old = (op.ovol_div * u0) ./ G.cells.volumes;
        end
    end
end

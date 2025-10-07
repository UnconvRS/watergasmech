function schedule = HnPschedule(inj, soak, prod, runtime, W, bc)
%Make a schedule for HnP with user-specified durations for injection,
%soaking, and production
%
% SYNOPSIS:
%   schedule = simpleSchedule(timesteps);
%   schedule = simpleSchedule(timesteps, 'W', W, 'src', src, 'bc', bc);
%
% PARAMETERS:
%   
%   inj     - array of [injection period(day) timestep size(day) rampup counts].
%   soak    - array of [soaking period(day) timestep size(day) rampup counts].
%   prod    - array of [production period(day) timestep size(day) rampup counts].
%   W -  Wells to be used in the schedule. The wells will be active in
%        all timesteps. W is a struct with the following order: injector,
%        soaking, and producer.
%
% RETURNS:
%   schedule - struct suitable for HnP 

    dt_cycle= inj(1) + soak(1) + prod(1);
    cycle_count = floor(runtime/dt_cycle);

    schedule = struct();
    [W_inj, W_soak, W_prod] = deal(W);
    
    W_inj(2).status = false; %producer, which is at row 2, is shut in during injection
    W_soak(1).status = false; W_soak(2).status = false; % injector and producer are shut-in during soaking
    W_prod(1).status = false; %injector, which is at row 1, is shut in during production 
    
    schedule.control = [struct('W', W_inj,'bc', bc);...  % injection
                        struct('W', W_soak,'bc', bc);     % soaking
                        struct('W', W_prod,'bc', bc)];... % production
    dt = [];
    control_id = [];  %injection, soaking, & production control ids are 1, 2, & 3, respectively
    dt_inj = rampupTimesteps(inj(1), inj(2), inj(3));
    dt_inj = [dt_inj(1);dt_inj(1);dt_inj];
    
    dt_soak = rampupTimesteps(soak(1), soak(2), soak(3));
    dt_prod = rampupTimesteps(prod(1), prod(2), prod(3));
    for i = 1:cycle_count
        dt = [dt;dt_inj;dt_soak;dt_prod];
        control_id = [control_id;[2;2;ones(numel(dt_inj(3:end)),1)];2*ones(numel(dt_soak),1);3*ones(numel(dt_prod),1)];
%         control_id = [control_id;ones(numel(dt_inj),1);2*ones(numel(dt_soak),1);3*ones(numel(dt_prod),1)];
    end
    schedule.step.val = dt;
    schedule.step.control = control_id;
end


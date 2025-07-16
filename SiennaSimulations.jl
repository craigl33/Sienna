#!/usr/bin/env julia

"""
SiennaSimulations.jl - Updated for Simplified Config Structure
=============================================================

Updated to work with the cleaned config.toml structure while maintaining all functionality.
"""

using PowerSystems
using PowerSimulations
const PSI = PowerSimulations
using PowerModels
const PM = PowerModels
using HiGHS
using Dates
using Logging
using JSON3
using DataFrames
using CSV
using Statistics: mean, std, median  # Only imports specific functions
using Printf
using Infiltrator

# Import our modules
include("SiennaConfig.jl")
include("SiennaSystem.jl")

"""
    SiennaSimulations

Streamlined simulation manager - updated for simplified config structure.
"""
mutable struct SiennaSimulations
    # Core components
    config::SiennaConfig
    sienna_system::SiennaSystem
    
    # Simulation configuration (parsed from simplified config)
    simulation_type::String
    simulation_name::String
    enable_multistage::Bool
    
    # Simulation parameters
    horizon_hours::Int
    interval_hours::Int
    total_steps::Int
    resolution::Period
    
    # PowerSimulations.jl objects
    templates::Dict{String, ProblemTemplate}
    models::Union{SimulationModels, Nothing}
    sequence::Union{SimulationSequence, Nothing}
    simulation::Union{Simulation, Nothing}
    
    # Results and tracking
    results::Union{SimulationResults, Nothing}
    performance_metrics::Dict{String, Any}
    simulation_status::String
    
    # Timestamps
    created_at::DateTime
    started_at::Union{DateTime, Nothing}
    completed_at::Union{DateTime, Nothing}
end

"""
    SiennaSimulations(config_file::String="config.toml", sienna_system::Union{SiennaSystem, Nothing}=nothing)

Constructor - Updated for simplified config structure.
"""
function SiennaSimulations(config_file::String="config.toml", sienna_system::Union{SiennaSystem, Nothing}=nothing)
    @info "🚀 Initializing SiennaSimulations from: $config_file"
    
    # Load and validate config
    config = SiennaConfig(config_file)
    
    # Handle system - use provided or build from config
    if sienna_system === nothing
        @info "No system provided - building from config"
        sienna_system = SiennaSystem(config)
    else
        @info "Using provided SiennaSystem"
    end
    
    # Ensure system is built
    if !is_system_built(sienna_system)
        @info "System not built - building now..."
        build_system!(sienna_system)
    else
        @info "✓ System already built"
    end
    
    return create_simulation_from_components(config, sienna_system)
end

"""
    SiennaSimulations(config::SiennaConfig, sienna_system::Union{SiennaSystem, Nothing}=nothing)

Alternative constructor from SiennaConfig object.
"""
function SiennaSimulations(config::SiennaConfig, sienna_system::Union{SiennaSystem, Nothing}=nothing)
    if !config.is_validated
        error("❌ Configuration must be validated")
    end
    
    @info "🔧 Building SiennaSimulations from SiennaConfig object"
    
    # Handle system - use provided or build from config
    if sienna_system === nothing
        @info "No system provided - building from config"
        sienna_system = SiennaSystem(config)
    else
        @info "Using provided SiennaSystem"
    end
    
    # Ensure system is built
    if !is_system_built(sienna_system)
        @info "System not built - building now..."
        build_system!(sienna_system)
    else
        @info "✓ System already built"
    end
    
    return create_simulation_from_components(config, sienna_system)
end

"""
    create_simulation_from_components(config::SiennaConfig, sienna_system::SiennaSystem)

Internal helper to create simulation object from validated components.
"""
function create_simulation_from_components(config::SiennaConfig, sienna_system::SiennaSystem)
    # Parse simulation configuration from simplified structure
    sim_config = get(config.config_data, "simulations", Dict())
    
    simulation_type = get(sim_config, "simulation_type", "ST")
    simulation_name = "$(config.project_name)_$(simulation_type)_simulation"
    enable_multistage = get(sim_config, "enable_multistage_simulation", false)
    
    # Parse simulation parameters based on type
    horizon_hours, interval_hours, total_steps, resolution = parse_simulation_parameters(config, simulation_type)
    
    # Create simulation manager
    sim = SiennaSimulations(
        config,
        sienna_system,
        simulation_type,
        simulation_name,
        enable_multistage,
        horizon_hours,
        interval_hours,
        total_steps,
        resolution,
        Dict{String, ProblemTemplate}(),
        nothing, nothing, nothing, nothing,
        Dict{String, Any}(),
        "INITIALIZED",
        now(), nothing, nothing
    )
    
    @info "✅ SiennaSimulations ready:"
    @info "   Type: $simulation_type"
    @info "   System: $(get_name(get_power_system(sienna_system)))"
    @info "   Horizon: $(horizon_hours)h, Interval: $(interval_hours)h, Steps: $total_steps"
    @info "   Multistage: $enable_multistage"
    
    return sim
end

"""
    parse_simulation_parameters(config::SiennaConfig, simulation_type::String)

Parse horizon, interval, steps, and resolution from simplified config structure.
"""
function parse_simulation_parameters(config::SiennaConfig, simulation_type::String)
    sim_config = get(config.config_data, "simulations", Dict())
    
    if simulation_type == "test"
        # Single test problem
        horizon = config.default_horizon_hours  # from simplified config
        interval = horizon  # Single step
        steps = 1
        resolution = Hour(1)
        
    elseif simulation_type == "ST"
        # Short-term: 365-day rolling horizon
        horizon = get(sim_config, "horizon_hours", config.default_horizon_hours)  # 48h
        interval = get(sim_config, "interval_hours", 24)  # 24h
        annual_days = get(sim_config, "annual_simulation_days", 365)
        steps = annual_days  # 365 daily steps
        resolution = Hour(1)
        
    elseif simulation_type == "MT"
        # Medium-term: Annual planning with reduced chronology
        horizon_days = get(sim_config, "mt_horizon_days", 365)
        typical_days = get(sim_config, "mt_typical_days", 12)
        resolution_hours = get(sim_config, "mt_resolution_hours", 4)
        
        horizon = horizon_days * 24  # Convert to hours
        interval = horizon  # Single annual solve
        steps = 1
        resolution = Hour(resolution_hours)
        
    elseif simulation_type == "DA"
        # Day-ahead: 24h horizon, daily steps
        horizon = 24
        interval = 24
        steps = get(sim_config, "annual_simulation_days", 365)
        resolution = Hour(1)
        
    elseif simulation_type == "RT"
        # Real-time: Shorter horizons, hourly steps
        horizon = get(sim_config, "rt_horizon_hours", 4)
        interval = get(sim_config, "rt_interval_hours", 1)
        simulation_hours = get(sim_config, "rt_simulation_hours", 8760)  # Full year
        steps = div(simulation_hours, interval)
        resolution = Hour(1)
        
    else
        error("❌ Unknown simulation type: $simulation_type. Valid: MT, ST, DA, RT, test")
    end
    
    @info "Simulation parameters:"
    @info "  Horizon: $(horizon) hours"
    @info "  Interval: $(interval) hours"
    @info "  Steps: $steps"
    @info "  Resolution: $resolution"
    
    return horizon, interval, steps, resolution
end

"""
    run_simulation!(sim::SiennaSimulations)

Main entry point - run simulation based on simplified configuration.
"""
function run_simulation!(sim::SiennaSimulations)
    @info "🎯 Running $(sim.simulation_type) simulation..."
    
    sim.started_at = now()
    sim.simulation_status = "RUNNING"
    
    try
        if sim.simulation_type == "test"
            results = run_test_simulation!(sim)
        elseif sim.simulation_type == "ST"
            results = run_short_term_simulation!(sim)
        elseif sim.simulation_type == "MT"
            results = run_medium_term_simulation!(sim)
        elseif sim.simulation_type == "DA"
            results = run_day_ahead_simulation!(sim)
        elseif sim.simulation_type == "RT"
            results = run_real_time_simulation!(sim)
        else
            error("❌ Unsupported simulation type: $(sim.simulation_type)")
        end
        
        sim.results = results
        sim.completed_at = now()
        sim.simulation_status = "COMPLETED"
        
        # Calculate performance metrics
        calculate_performance_metrics!(sim)
        
        # Print summary
        print_simulation_summary(sim)
        
        @info "✅ $(sim.simulation_type) simulation completed successfully"
        return results
        
    catch e
        sim.simulation_status = "FAILED"
        sim.completed_at = now()
        @error "❌ $(sim.simulation_type) simulation failed: $e"
        rethrow(e)
    end
end

"""
    create_template_for_formulation(sim::SiennaSimulations, formulation_type::String)

Create PowerSimulations template from simplified config structure.
"""
function create_template_for_formulation(sim::SiennaSimulations, formulation_type::String)
    @info "Creating template for: $formulation_type"
    
    # Check if already created
    if haskey(sim.templates, formulation_type)
        return sim.templates[formulation_type]
    end
    
    template = ProblemTemplate()
    sys = get_power_system(sim.sienna_system)
    config = sim.config
    
    # Network model from simplified config
    network_settings = get_network_settings(config)
    network_model_name = network_settings["model"]
    network_model = get_network_model_from_config(network_model_name)
    set_network_model!(template, NetworkModel(network_model))
    
    # Device formulations from simplified config
    device_formulations = get_device_formulations(config, formulation_type)
    
    # Get system components
    thermal_gens = get_components(ThermalStandard, sys)
    renewable_dispatch = get_components(RenewableDispatch, sys)
    renewable_nondispatch = get_components(RenewableNonDispatch, sys)
    hydro_dispatch = get_components(HydroDispatch, sys)
    loads = get_components(PowerLoad, sys)
    lines = get_components(Line, sys)
    dc_lines = get_components(TwoTerminalHVDCLine, sys)
    
    # Set device models using custom formulations
    if !isempty(thermal_gens)
        thermal_formulation_name = get(device_formulations, "thermal_standard", 
                                     formulation_type == "economic_dispatch" ? "ThermalBasicDispatch" : "ThermalBasicUnitCommitment")
        thermal_formulation = get_device_formulation_object(thermal_formulation_name)
        set_device_model!(template, ThermalStandard, thermal_formulation)
        @info "  ✓ ThermalStandard: $(length(thermal_gens)) units ($thermal_formulation_name)"
    end
    
    if !isempty(renewable_dispatch)
        renewable_formulation_name = get(device_formulations, "renewable_dispatch", "RenewableFullDispatch")
        renewable_formulation = get_device_formulation_object(renewable_formulation_name)
        set_device_model!(template, RenewableDispatch, renewable_formulation)
        @info "  ✓ RenewableDispatch: $(length(renewable_dispatch)) units ($renewable_formulation_name)"
    end
    
    if !isempty(renewable_nondispatch)
        renewable_nd_formulation_name = get(device_formulations, "renewable_nondispatch", "FixedOutput")
        renewable_nd_formulation = get_device_formulation_object(renewable_nd_formulation_name)
        set_device_model!(template, RenewableNonDispatch, renewable_nd_formulation)
        @info "  ✓ RenewableNonDispatch: $(length(renewable_nondispatch)) units ($renewable_nd_formulation_name)"
    end
    
    if !isempty(hydro_dispatch)
        hydro_formulation_name = get(device_formulations, "hydro_dispatch", "HydroDispatchRunOfRiver")
        hydro_formulation = get_device_formulation_object(hydro_formulation_name)
        set_device_model!(template, HydroDispatch, hydro_formulation)
        @info "  ✓ HydroDispatch: $(length(hydro_dispatch)) units ($hydro_formulation_name)"
    end
    
    if !isempty(loads)
        load_formulation_name = get(device_formulations, "power_load", "StaticPowerLoad")
        load_formulation = get_device_formulation_object(load_formulation_name)
        set_device_model!(template, PowerLoad, load_formulation)
        @info "  ✓ PowerLoad: $(length(loads)) units ($load_formulation_name)"
    end
    
    if !isempty(lines)
        line_formulation_name = get(device_formulations, "line", "StaticBranch")
        line_formulation = get_device_formulation_object(line_formulation_name)
        set_device_model!(template, Line, line_formulation)
        @info "  ✓ Line: $(length(lines)) branches ($line_formulation_name)"
    end
    
    if !isempty(dc_lines)
        dc_formulation_name = get(device_formulations, "dc_line", "HVDCTwoTerminalDispatch")
        dc_formulation = get_device_formulation_object(dc_formulation_name)
        set_device_model!(template, TwoTerminalHVDCLine, dc_formulation)
        @info "  ✓ TwoTerminalHVDCLine: $(length(dc_lines)) lines ($dc_formulation_name)"
    end
    
    # Store template
    sim.templates[formulation_type] = template
    
    return template
end

"""
    create_optimizer_from_config(config::SiennaConfig, formulation_type::String)

Create optimizer using simplified solver configuration.
"""
function create_optimizer_from_config(config::SiennaConfig, formulation_type::String)
    solver_settings = get_solver_settings(config, formulation_type)
    
    thread_count = solver_settings["threads"] > 0 ? solver_settings["threads"] : min(Sys.CPU_THREADS, 4)
    
    return optimizer_with_attributes(
        HiGHS.Optimizer,
        "time_limit" => solver_settings["time_limit"],
        "mip_rel_gap" => solver_settings["mip_gap"],
        "threads" => thread_count,
        "output_flag" => solver_settings["output_flag"],
        "presolve" => solver_settings["presolve"],
        "parallel" => solver_settings["parallel"]
    )
end

"""
    get_network_model_from_config(model_name::String)

Get PowerSimulations network model from config string.
"""
function get_network_model_from_config(model_name::String)
    if model_name == "CopperPlatePowerModel"
        return PSI.CopperPlatePowerModel
    elseif model_name == "PTDFPowerModel"
        return PSI.PTDFPowerModel
    elseif model_name == "AreaBalancePowerModel"
        return PSI.AreaBalancePowerModel
    elseif model_name == "AreaPTDFPowerModel"
        return PSI.AreaPTDFPowerModel
    elseif model_name == "DCPPowerModel"
        return PM.DCPPowerModel
    elseif model_name == "DCMPPowerModel"
        return PM.DCMPPowerModel
    elseif model_name == "LPACCPowerModel"
        return PM.LPACCPowerModel
    else
        error("❌ Invalid network model: '$model_name'")
    end
end

"""
    get_device_formulation_object(formulation_name::String)

Get PowerSimulations device formulation object from name.
"""
function get_device_formulation_object(formulation_name::String)
    formulation_symbol = Symbol(formulation_name)
    if isdefined(PSI, formulation_symbol)
        return getfield(PSI, formulation_symbol)
    else
        error("❌ Formulation not found in PowerSimulations: '$formulation_name'")
    end
end

"""
    create_decision_model(sim::SiennaSimulations, template::ProblemTemplate, model_name::String)

Create DecisionModel with simulation parameters.
"""
function create_decision_model(sim::SiennaSimulations, template::ProblemTemplate, model_name::String)
    sys = get_power_system(sim.sienna_system)
    optimizer = create_optimizer_from_config(sim.config, "economic_dispatch")
    
    return DecisionModel(
        template,
        sys;
        optimizer = optimizer,
        horizon = Hour(sim.horizon_hours),
        resolution = sim.resolution,
        name = model_name
    )
end

"""
    create_output_directory(sim::SiennaSimulations)

Create timestamped output directory.
"""
function create_output_directory(sim::SiennaSimulations)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    output_dir = joinpath(
        sim.config.output_directory,
        "$(sim.simulation_type)_simulation_$timestamp"
    )
    mkpath(output_dir)
    return output_dir
end

"""
    get_initial_simulation_time(sim::SiennaSimulations)

Get initial time for simulation based on forecast data using correct PowerSystems.jl API.
"""
function get_initial_simulation_time(sim::SiennaSimulations)
    sys = get_power_system(sim.sienna_system)
    

    @info "🕐 Determining simulation initial time from forecast data..."
    
    # Strategy 1: Use get_time_series_array to access Deterministic forecasts
    try
        loads = get_components(PowerLoad, sys)
        if !isempty(loads)
            @info "  Checking load forecasts for initial time..."
            first_load = first(loads)
            
            # Try to get Deterministic forecast using get_time_series_array
            try
                # Get forecast window - this should work if Deterministic forecasts exist
                forecast_array = get_time_series_array(
                    Deterministic, 
                    first_load, 
                    "max_active_power"
                )
                if forecast_array !== nothing
                    initial_time = first(timestamp(forecast_array))
                    @info "  ✓ Found Deterministic forecast, initial time: $initial_time"
                    return initial_time
                end
            catch e
                @debug "  Deterministic forecast array method failed: $e"
            end
        end
    catch e
        @warn "  Load forecast extraction failed: $e"
    end
    
    # Strategy 2: Try renewable generators
    try
        renewables = get_components(RenewableDispatch, sys)
        if !isempty(renewables)
            @info "  Checking renewable forecasts for initial time..."
            first_renewable = first(renewables)
            
            try
                forecast_array = get_time_series_array(
                    Deterministic, 
                    first_renewable, 
                    "max_active_power"
                )
                if forecast_array !== nothing
                    initial_time = first(timestamp(forecast_array))
                    @info "  ✓ Found renewable forecast, initial time: $initial_time"
                    return initial_time
                end
            catch e
                @debug "  Renewable forecast extraction failed: $e"
            end
        end
    catch e
        @warn "  Renewable forecast check failed: $e"
    end
    
    # Strategy 3: Try SingleTimeSeries if Deterministic failed
    try
        @info "  Checking for SingleTimeSeries data..."
        loads = get_components(PowerLoad, sys)
        if !isempty(loads)
            first_load = first(loads)
            
            try
                # Try SingleTimeSeries format
                ts_array = get_time_series_array(
                    SingleTimeSeries, 
                    first_load, 
                    "max_active_power"
                )
                if ts_array !== nothing
                    initial_time = first(timestamp(ts_array))
                    @info "  ✓ Found SingleTimeSeries, initial time: $initial_time"
                    @warn "  Note: Using SingleTimeSeries data - you may need to transform to forecasts"
                    return initial_time
                end
            catch e
                @debug "  SingleTimeSeries check failed: $e"
            end
        end
    catch e
        @warn "  SingleTimeSeries check failed: $e"
    end
    
    # Strategy 4: Check any component for any time series
    try
        @info "  Searching all components for any time series data..."
        
        all_components = [
            get_components(PowerLoad, sys)...,
            get_components(RenewableDispatch, sys)...,
            get_components(ThermalStandard, sys)...
        ]
        
        for component in all_components
            if has_time_series(component)
                @info "  Component $(get_name(component)) has time series"
                
                # Try both Deterministic and SingleTimeSeries
                for ts_type in [Deterministic, SingleTimeSeries]
                    try
                        ts_array = get_time_series_array(
                            ts_type, 
                            component, 
                            "max_active_power"
                        )
                        if ts_array !== nothing
                            initial_time = first(timestamp(ts_array))
                            @info "  ✓ Found $ts_type on $(get_name(component)), initial time: $initial_time"
                            return initial_time
                        end
                    catch e
                        @debug "  $ts_type failed on $(get_name(component)): $e"
                        continue
                    end
                end
            end
        end
    catch e
        @warn "  System-wide time series search failed: $e"
    end
    
    # If we get here, no forecast data was accessible
    @error "❌ Could not extract initial time from any time series data!"
    @error "   This suggests:"
    @error "   1. No time series data was actually loaded"
    @error "   2. Time series exists but with different parameter names"
    @error "   3. Forecast transformation to Deterministic failed"
    
    @info "💡 DEBUGGING HELP:"
    @info "   1. Run: inspect_time_series_data(sim) to see what time series exist"
    @info "   2. Check your timeseries CSV files - what date range do they cover?"
    @info "   3. Verify load_timeseries = true in config"
    @info "   4. Check that transform_single_time_series! succeeded"
    
    # Use a fallback based on common dataset ranges
    fallback_time = DateTime(2020, 1, 1, 0, 0, 0)  # Common starting year
    @warn "⚠️  Using fallback initial time: $fallback_time"
    @warn "   If this fails, check your actual time series data dates"
    
    return fallback_time
end

# ===== SIMULATION EXECUTION METHODS (unchanged) =====

"""
    run_test_simulation!(sim::SiennaSimulations)

Run single test problem for debugging.
"""
function run_test_simulation!(sim::SiennaSimulations)
    @info "🧪 Running test simulation (single DecisionModel)..."
    
    # Create template
    template = create_template_for_formulation(sim, "economic_dispatch")
    
    # Create optimizer
    optimizer = create_optimizer_from_config(sim.config, "economic_dispatch")
    
    # Create DecisionModel
    sys = get_power_system(sim.sienna_system)
    problem = DecisionModel(
        template,
        sys;
        optimizer = optimizer,
        horizon = Hour(sim.horizon_hours),
        resolution = sim.resolution,
        name = sim.simulation_name
    )
    
    # Create output directory
    output_dir = create_output_directory(sim)
    
    # Build and solve
    @info "Building test problem..."
    build_time = @elapsed build!(problem, output_dir = output_dir)
    
    @info "Solving test problem..."
    solve_time = @elapsed solve!(problem)
    
    # Extract results
    problem_results = OptimizationProblemResults(problem)
    
    # Store performance
    sim.performance_metrics = Dict(
        "build_time_seconds" => build_time,
        "solve_time_seconds" => solve_time,
        "total_time_seconds" => build_time + solve_time,
        "objective_value" => get_objective_value(problem_results)
    )
    
    @info "✓ Test simulation completed in $(round(build_time + solve_time, digits=2))s"
    
    # Wrap in SimulationResults-like structure for consistency
    return TestSimulationResults(problem_results, sim.performance_metrics)
end

"""
    run_short_term_simulation!(sim::SiennaSimulations)

Run 365-day rolling horizon simulation.
"""
function run_short_term_simulation!(sim::SiennaSimulations)
    @info "📅 Running short-term simulation (365-day rolling horizon)..."
    @info "   365 daily steps with $(sim.horizon_hours)h horizon and $(sim.interval_hours)h interval"
    
    # Create template for operational simulation
    template = create_template_for_formulation(sim, "economic_dispatch")
    
    # Create DecisionModel
    decision_model = create_decision_model(sim, template, "short_term_operations")
    
    # Create SimulationModels
    sim.models = SimulationModels(decision_models = [decision_model])
    
    # Create SimulationSequence
    sim.sequence = SimulationSequence(
        models = sim.models,
        ini_cond_chronology = InterProblemChronology()
    )
    
    # Create output directory
    output_dir = create_output_directory(sim)
    

        # 🔧 DEBUG: Check values right before Simulation creation
    @info "DEBUG: Right before Simulation() creation:"
    @info "  sim.simulation_name: $(sim.simulation_name)"
    @info "  sim.total_steps: $(sim.total_steps)"
    @info "  sim.models: $(typeof(sim.models))"
    @info "  sim.sequence: $(typeof(sim.sequence))"
    @info "  output_dir: $output_dir"

    # Check if total_steps got corrupted somehow
    if sim.total_steps <= 0
        @error "❌ sim.total_steps is $(sim.total_steps) right before Simulation creation!"
        
        # Try to fix it
        annual_days = get(sim.config.config_data["simulations"], "annual_simulation_days", 365)
        sim.total_steps = max(1, Int(annual_days))
        @info "🔧 Fixed sim.total_steps to: $(sim.total_steps)"
    end

    # Also check available time series
    sys = get_power_system(sim.sienna_system)
    loads = get_components(PowerLoad, sys)
    if !isempty(loads)
        first_load = first(loads)
        if has_time_series(first_load)
            try
                ts_data = get_time_series_array(Deterministic, first_load, "max_active_power")
                available_windows = length(ts_data) ÷ sim.horizon_hours
                @info "  Available forecast windows: $available_windows"
                @info "  Requested simulation steps: $(sim.total_steps)"
                
                if available_windows < sim.total_steps
                    @error "❌ Available windows ($available_windows) < requested steps ($(sim.total_steps))"
                    @info "💡 This explains the '364 vs 365' error!"
                end
            catch e
                @warn "Could not check time series: $e"
            end
        end
    end
    # Create Simulation
    sim.simulation = Simulation(
        name = sim.simulation_name,
        steps = sim.total_steps,  # 365 days
        models = sim.models,
        sequence = sim.sequence,
        simulation_folder = output_dir,
        initial_time = get_initial_simulation_time(sim)
    )

    @infiltrate
    
    # Build simulation
    @info "Building 365-day simulation..."
    build_time = @elapsed build!(sim.simulation)
    
    # Execute simulation
    @info "Executing 365-day simulation (this may take a while)..."
    solve_time = @elapsed execute!(sim.simulation)
    
    # Extract results
    results = SimulationResults(sim.simulation)
    
    # Store performance
    sim.performance_metrics = Dict(
        "build_time_seconds" => build_time,
        "solve_time_seconds" => solve_time,
        "total_time_seconds" => build_time + solve_time,
        "simulation_type" => "short_term",
        "total_steps" => sim.total_steps,
        "horizon_hours" => sim.horizon_hours,
        "interval_hours" => sim.interval_hours
    )
    
    @info "✓ Short-term simulation completed in $(round(build_time + solve_time, digits=2))s"
    
    return results
end

"""
    run_medium_term_simulation!(sim::SiennaSimulations)

Run medium-term planning simulation.
"""
function run_medium_term_simulation!(sim::SiennaSimulations)
    @info "📊 Running medium-term simulation (annual planning)..."
    
    # Use unit commitment for MT planning
    template = create_template_for_formulation(sim, "unit_commitment")
    
    # Create DecisionModel with longer horizon and coarser resolution
    decision_model = create_decision_model(sim, template, "medium_term_planning")
    
    # Single annual optimization
    sim.models = SimulationModels(decision_models = [decision_model])
    sim.sequence = SimulationSequence(models = sim.models)
    
    output_dir = create_output_directory(sim)
    
    sim.simulation = Simulation(
        name = sim.simulation_name,
        steps = sim.total_steps,  # 1 step for annual
        models = sim.models,
        sequence = sim.sequence,
        simulation_folder = output_dir,
        initial_time = get_initial_simulation_time(sim)
    )
    
    # Build and execute
    @info "Building medium-term planning model..."
    build_time = @elapsed build!(sim.simulation)
    
    @info "Solving annual planning problem..."
    solve_time = @elapsed execute!(sim.simulation)
    
    results = SimulationResults(sim.simulation)
    
    sim.performance_metrics = Dict(
        "build_time_seconds" => build_time,
        "solve_time_seconds" => solve_time,
        "total_time_seconds" => build_time + solve_time,
        "simulation_type" => "medium_term",
        "horizon_hours" => sim.horizon_hours
    )
    
    @info "✓ Medium-term simulation completed in $(round(build_time + solve_time, digits=2))s"
    
    return results
end

"""
    run_day_ahead_simulation!(sim::SiennaSimulations)

Run day-ahead market simulation.
"""
function run_day_ahead_simulation!(sim::SiennaSimulations)
    @info "🌅 Running day-ahead simulation..."
    
    template = create_template_for_formulation(sim, "economic_dispatch")
    decision_model = create_decision_model(sim, template, "day_ahead_market")
    
    sim.models = SimulationModels(decision_models = [decision_model])
    sim.sequence = SimulationSequence(
        models = sim.models,
        ini_cond_chronology = InterProblemChronology()
    )
    
    output_dir = create_output_directory(sim)
    
    sim.simulation = Simulation(
        name = sim.simulation_name,
        steps = sim.total_steps,
        models = sim.models,
        sequence = sim.sequence,
        simulation_folder = output_dir,
        initial_time = get_initial_simulation_time(sim)
    )
    
    build_time = @elapsed build!(sim.simulation)
    solve_time = @elapsed execute!(sim.simulation)
    
    results = SimulationResults(sim.simulation)
    
    sim.performance_metrics = Dict(
        "build_time_seconds" => build_time,
        "solve_time_seconds" => solve_time,
        "total_time_seconds" => build_time + solve_time,
        "simulation_type" => "day_ahead"
    )
    
    return results
end

"""
    run_real_time_simulation!(sim::SiennaSimulations)

Run real-time operations simulation.
"""
function run_real_time_simulation!(sim::SiennaSimulations)
    @info "⚡ Running real-time simulation..."
    
    # Use unit commitment for real-time operations
    template = create_template_for_formulation(sim, "unit_commitment")
    decision_model = create_decision_model(sim, template, "real_time_operations")
    
    sim.models = SimulationModels(decision_models = [decision_model])
    sim.sequence = SimulationSequence(
        models = sim.models,
        ini_cond_chronology = InterProblemChronology()
    )
    
    output_dir = create_output_directory(sim)
    
    sim.simulation = Simulation(
        name = sim.simulation_name,
        steps = sim.total_steps,
        models = sim.models,
        sequence = sim.sequence,
        simulation_folder = output_dir,
        initial_time = get_initial_simulation_time(sim)
    )
    
    build_time = @elapsed build!(sim.simulation)
    solve_time = @elapsed execute!(sim.simulation)
    
    results = SimulationResults(sim.simulation)
    
    sim.performance_metrics = Dict(
        "build_time_seconds" => build_time,
        "solve_time_seconds" => solve_time,
        "total_time_seconds" => build_time + solve_time,
        "simulation_type" => "real_time"
    )
    
    return results
end

# ===== PERFORMANCE AND ANALYSIS METHODS =====

"""
    calculate_performance_metrics!(sim::SiennaSimulations)

Calculate and store performance metrics.
"""
function calculate_performance_metrics!(sim::SiennaSimulations)
    if sim.started_at !== nothing && sim.completed_at !== nothing
        total_runtime = Dates.value(sim.completed_at - sim.started_at) / 1000.0
        sim.performance_metrics["total_runtime_seconds"] = total_runtime
        sim.performance_metrics["simulation_type"] = sim.simulation_type
        sim.performance_metrics["horizon_hours"] = sim.horizon_hours
        sim.performance_metrics["interval_hours"] = sim.interval_hours
        sim.performance_metrics["total_steps"] = sim.total_steps
        sim.performance_metrics["status"] = sim.simulation_status
    end
end

"""
    print_simulation_summary(sim::SiennaSimulations)

Print comprehensive simulation summary.
"""
function print_simulation_summary(sim::SiennaSimulations)
    @info "\n" * "="^70
    @info "SIENNA SIMULATION SUMMARY (SIMPLIFIED CONFIG)"
    @info "="^70
    @info "Type: $(uppercase(sim.simulation_type))"
    @info "Project: $(sim.config.project_name)"
    @info "Status: $(sim.simulation_status)"
    
    if haskey(sim.performance_metrics, "total_runtime_seconds")
        runtime = sim.performance_metrics["total_runtime_seconds"]
        @info "Total Runtime: $(round(runtime, digits=1)) seconds"
    end
    
    @info ""
    @info "Configuration:"
    @info "  Horizon: $(sim.horizon_hours) hours"
    @info "  Interval: $(sim.interval_hours) hours"
    @info "  Steps: $(sim.total_steps)"
    @info "  Network Model: $(sim.config.network_model)"
    @info "  Resolution: $(sim.resolution)"
    
    if haskey(sim.performance_metrics, "build_time_seconds")
        build_time = sim.performance_metrics["build_time_seconds"]
        solve_time = sim.performance_metrics["solve_time_seconds"]
        @info ""
        @info "Performance:"
        @info "  Build Time: $(round(build_time, digits=2)) seconds"
        @info "  Solve Time: $(round(solve_time, digits=2)) seconds"
        @info "  Total Time: $(round(build_time + solve_time, digits=2)) seconds"
    end
    
    if haskey(sim.performance_metrics, "objective_value")
        obj_val = sim.performance_metrics["objective_value"]
        @info "  Objective Value: $(round(obj_val, digits=2))"
    end
    
    @info ""
    @info "System Details:"
    sys_stats = get_component_counts(sim.sienna_system)
    for (component, count) in sys_stats
        @info "  $component: $count"
    end
    
    @info "="^70
end

# ===== HELPER STRUCTURES =====

"""
    TestSimulationResults

Wrapper for single DecisionModel results to match SimulationResults interface.
"""
struct TestSimulationResults
    problem_results::OptimizationProblemResults
    performance_metrics::Dict{String, Any}
end

# ===== CONVENIENCE FUNCTIONS =====

"""
    run_test_simulation(config_file::String="config.toml")

Quick test simulation for debugging.
"""
function run_test_simulation(config_file::String="config.toml")
    @info "🧪 Running quick test simulation..."
    
    # Temporarily modify config for test mode
    config = SiennaConfig(config_file)
    config.config_data["simulations"]["simulation_type"] = "test"
    
    sim = SiennaSimulations(config)
    return run_simulation!(sim)
end

"""
    run_annual_simulation(config_file::String="config.toml")

Run your main 365-day simulation.
"""
function run_annual_simulation(config_file::String="config.toml")
    @info "📅 Running annual 365-day simulation..."
    
    config = SiennaConfig(config_file)
    config.config_data["simulations"]["simulation_type"] = "ST"
    
    sim = SiennaSimulations(config)
    return run_simulation!(sim)
end

"""
    run_medium_term_simulation(config_file::String="config.toml")

Run medium-term planning simulation.
"""
function run_medium_term_simulation(config_file::String="config.toml")
    @info "📊 Running medium-term planning simulation..."
    
    config = SiennaConfig(config_file)
    config.config_data["simulations"]["simulation_type"] = "MT"
    
    sim = SiennaSimulations(config)
    return run_simulation!(sim)
end

# ===== RESULTS ANALYSIS FUNCTIONS =====

"""
    analyze_simulation_results(sim::SiennaSimulations)

Analyze and print detailed results analysis.
"""
function analyze_simulation_results(sim::SiennaSimulations)
    if sim.results === nothing
        @warn "No results available for analysis"
        return
    end
    
    @info "\n" * "="^60
    @info "DETAILED RESULTS ANALYSIS"
    @info "="^60
    
    if sim.simulation_type == "test"
        analyze_test_results(sim)
    else
        analyze_full_simulation_results(sim)
    end
end

"""
    analyze_test_results(sim::SiennaSimulations)

Analyze single DecisionModel test results.
"""
function analyze_test_results(sim::SiennaSimulations)
    test_results = sim.results
    problem_results = test_results.problem_results
    
    @info "Test Simulation Analysis:"
    
    # Objective value
    try
        obj_val = get_objective_value(problem_results)
        @info "  Objective Value: $(round(obj_val, digits=2))"
    catch e
        @warn "  Could not get objective value: $e"
    end
    
    # Variables analysis
    try
        var_names = list_variable_names(problem_results)
        @info "  Available Variables: $(length(var_names))"
        
        # Analyze key variables
        if "ActivePowerVariable__ThermalStandard" in var_names
            thermal_data = read_variable(problem_results, "ActivePowerVariable__ThermalStandard")
            if ncol(thermal_data) > 1
                total_thermal = sum(sum(eachcol(thermal_data[:, 2:end])))
                @info "  Total Thermal Generation: $(round(total_thermal, digits=1)) MWh"
            end
        end
        
        if "ActivePowerVariable__RenewableDispatch" in var_names
            renewable_data = read_variable(problem_results, "ActivePowerVariable__RenewableDispatch")
            if ncol(renewable_data) > 1
                total_renewable = sum(sum(eachcol(renewable_data[:, 2:end])))
                @info "  Total Renewable Generation: $(round(total_renewable, digits=1)) MWh"
            end
        end
        
    catch e
        @warn "  Variable analysis failed: $e"
    end
end

"""
    analyze_full_simulation_results(sim::SiennaSimulations)

Analyze multi-step simulation results.
"""
function analyze_full_simulation_results(sim::SiennaSimulations)
    results = sim.results
    
    @info "$(uppercase(sim.simulation_type)) Simulation Analysis:"
    @info "  Total Steps: $(sim.total_steps)"
    @info "  Horizon per Step: $(sim.horizon_hours) hours"
    @info "  Interval per Step: $(sim.interval_hours) hours"
    
    # Try to get aggregate statistics
    try
        # This would require extracting data from all simulation steps
        # Implementation depends on specific analysis needs
        @info "  Detailed analysis available via PowerAnalytics.jl integration"
        
    catch e
        @warn "  Detailed analysis failed: $e"
    end
end

"""
    save_simulation_results(sim::SiennaSimulations, output_path::String)

Save simulation results with metadata.
"""
function save_simulation_results(sim::SiennaSimulations, output_path::String)
    @info "💾 Saving simulation results to: $output_path"
    
    mkpath(dirname(output_path))
    
    # Create comprehensive summary
    summary = Dict(
        "simulation_metadata" => Dict(
            "type" => sim.simulation_type,
            "name" => sim.simulation_name,
            "status" => sim.simulation_status,
            "created_at" => string(sim.created_at),
            "started_at" => string(sim.started_at),
            "completed_at" => string(sim.completed_at)
        ),
        "configuration" => Dict(
            "horizon_hours" => sim.horizon_hours,
            "interval_hours" => sim.interval_hours,
            "total_steps" => sim.total_steps,
            "resolution" => string(sim.resolution),
            "multistage" => sim.enable_multistage,
            "network_model" => sim.config.network_model,
            "project_name" => sim.config.project_name
        ),
        "performance_metrics" => sim.performance_metrics,
        "system_statistics" => get_component_counts(sim.sienna_system),
        "sienna_version" => "v2.0_simplified",
        "powersimulations_version" => "v0.30.2"
    )
    
    # Save summary
    open(output_path, "w") do f
        JSON3.pretty(f, summary)
    end
    
    @info "✓ Results summary saved"
end

# ===== PUBLIC API =====

"""
    get_simulation_results(sim::SiennaSimulations)

Get simulation results (SimulationResults or TestSimulationResults).
"""
function get_simulation_results(sim::SiennaSimulations)
    if sim.results === nothing
        error("❌ No results available. Run simulation first.")
    end
    return sim.results
end

"""
    get_performance_metrics(sim::SiennaSimulations)

Get performance metrics dictionary.
"""
function get_performance_metrics(sim::SiennaSimulations)
    return copy(sim.performance_metrics)
end

"""
    is_simulation_complete(sim::SiennaSimulations)

Check if simulation completed successfully.
"""
function is_simulation_complete(sim::SiennaSimulations)
    return sim.simulation_status == "COMPLETED"
end

# ===== EXPORTS =====

export SiennaSimulations
export run_simulation!
export run_test_simulation, run_annual_simulation, run_medium_term_simulation
export analyze_simulation_results, save_simulation_results
export get_simulation_results, get_performance_metrics, is_simulation_complete

# ===== TESTING =====

if abspath(PROGRAM_FILE) == @__FILE__
    @info "🧪 Testing updated SiennaSimulations with simplified config..."
    
    try
        if isfile("config.toml")
            @info "Running test simulation with simplified config structure..."
            
            # Test the main constructor pattern
            sim = SiennaSimulations("config.toml")
            @info "✅ Constructor successful"
            @info "   Simulation type: $(sim.simulation_type)"
            @info "   Status: $(sim.simulation_status)"
            
        else
            @info "No config.toml found - create simplified configuration first"
            @info "Expected usage: sim = SiennaSimulations(\"config.toml\")"
        end
    catch e
        @error "❌ Updated SiennaSimulations test failed: $e"
        
        # Print detailed error for debugging
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
end
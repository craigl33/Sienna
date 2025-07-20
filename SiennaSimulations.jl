#!/usr/bin/env julia

"""
SiennaSimulations.jl - COMPLETE Drop-in Replacement for Simulation Execution
===========================================================================

Complete implementation adapted for PowerSimulations.jl 0.30.2 simulation execution.
Includes ALL missing functions from your original code.
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
using TOML
using Statistics: mean

# Import our modules
include("SiennaConfig.jl")
include("SiennaSystem.jl")

"""
    SiennaSimulations

Main simulation manager class adapted for simulation execution.
"""
mutable struct SiennaSimulations
    # Core components
    config::SiennaConfig
    sys::SiennaSystem
    
    # Templates and problems
    templates::Dict{String, ProblemTemplate}
    problems::Dict{String, DecisionModel}
    simulation::Union{Simulation, Nothing}
    models::Union{SimulationModels, Nothing}
    sequence::Union{SimulationSequence, Nothing}
    
    # Results (adapted for simulation)
    results::Union{SimulationResults, Nothing}
    
    # Simulation metadata
    simulation_name::String
    created_timestamp::DateTime

    # Performance tracking
    build_times::Dict{String, Float64}
    solve_times::Dict{String, Float64}
    
    # Status tracking
    templates_created::Bool
    problems_built::Dict{String, Bool}
    simulations_completed::Dict{String, Bool}
    
    # Error handling
    simulation_errors::Dict{String, Vector{String}}
    last_error::Union{String, Nothing}
end

"""
    SiennaSimulations(config_in, system_in)

Constructor - Initialize simulation manager with config and system.
"""
function SiennaSimulations(
    config_in::Union{String, SiennaConfig},
    system_in::Union{Nothing, SiennaSystem} = nothing,
)
    @info "🔍 Creating SiennaSimulations..."
    
    # Handle config
    config = if config_in isa String
        @info "Loading config from: $config_in"
        SiennaConfig(config_in)
    else
        config_in
    end
    
    # Handle system
    sys = if system_in === nothing
        @info "Creating system from config"
        system = SiennaSystem(config)
        if !is_system_built(system)
            build_system!(system)
        end
        system
    else
        if !is_system_built(system_in)
            build_system!(system_in)
        end
        system_in
    end

    if !config.is_validated
        error("Configuration must be validated before creating SiennaSimulations")
    end

    # Initialize with default values
    sienna_sim = SiennaSimulations(
        config,                           # 1. config::SiennaConfig
        sys,                              # 2. sys::SiennaSystem
        Dict{String, ProblemTemplate}(),  # 3. templates
        Dict{String, DecisionModel}(),    # 4. problems
        nothing,                          # 5. simulation
        nothing,                          # 6. models
        nothing,                          # 7. sequence
        nothing,                          # 8. results
        "$(config.project_name)_simulation", # 9. simulation_name
        now(),                            # 10. created_timestamp
        Dict{String, Float64}(),          # 11. build_times
        Dict{String, Float64}(),          # 12. solve_times
        false,                            # 13. templates_created
        Dict{String, Bool}(),             # 14. problems_built
        Dict{String, Bool}(),             # 15. simulations_completed
        Dict{String, Vector{String}}(),   # 16. simulation_errors
        nothing                           # 17. last_error
    )
    
    # Initialize tracking dictionaries
    formulation_types = ["mt", "ed", "uc", "simulation"]
    for formulation in formulation_types
        sienna_sim.problems_built[formulation] = false
        sienna_sim.simulations_completed[formulation] = false
        sienna_sim.simulation_errors[formulation] = String[]
    end
    
    @info "✅ SiennaSimulations initialized successfully"
    @info "   System: $(get_name(get_power_system(sys)))"
    @info "   Network Model: $(config.network_model)"
    
    return sienna_sim
end

# ===== TEMPLATE CREATION (FROM YOUR ORIGINAL) =====

"""
    create_templates!(sienna_sim::SiennaSimulations)

Create PowerSimulations templates for all formulation types.
"""
function create_templates!(sienna_sim::SiennaSimulations)
    @info "📋 Creating PowerSimulations templates..."
    
    sys = get_power_system(sienna_sim.sys)
    config = sienna_sim.config
    
    try
        # Create templates for each formulation type
        if should_run_formulation(config, "ed")
            @info "Creating Economic Dispatch template..."
            sienna_sim.templates["ed"] = create_config_driven_template(
                sys, "ed", config
            )
            @info "✓ Economic Dispatch template created"
        end
        
        if should_run_formulation(config, "uc")
            @info "Creating Unit Commitment template..."
            sienna_sim.templates["uc"] = create_config_driven_template(
                sys, "uc", config
            )
            @info "✓ Unit Commitment template created"
        end
        
        sienna_sim.templates_created = true
        @info "✅ All templates created successfully"
        
    catch e
        sienna_sim.last_error = "Template creation failed: $e"
        @error "❌ Template creation failed: $e"
        rethrow(e)
    end
end

"""
    create_config_driven_template(sys::System, formulation_type::String, config::SiennaConfig)

Create a PowerSimulations template using pure configuration approach.
"""
function create_config_driven_template(sys::System, formulation_type::String, config::SiennaConfig)
    @info "Creating config-driven template for: $formulation_type"
    
    template = ProblemTemplate()
    
    # === NETWORK MODEL FROM CONFIG ===
    network_model_name = config.network_model
    network_model = get_network_model_from_config(network_model_name)
    set_network_model!(template, NetworkModel(network_model))
    @info "  ✓ Network Model: $network_model_name"
    
    # === DEVICE FORMULATIONS FROM CONFIG ===
    device_formulations = get_device_formulations(config, formulation_type)
    @info "  Using device formulations from config: [$formulation_type]"
    
    # Get system components for validation
    thermal_gens = get_components(ThermalStandard, sys)
    renewable_dispatch = get_components(RenewableDispatch, sys)
    renewable_nondispatch = get_components(RenewableNonDispatch, sys)
    hydro_dispatch = get_components(HydroDispatch, sys)
    loads = get_components(PowerLoad, sys)
    lines = get_components(Line, sys)
    dc_lines = get_components(TwoTerminalHVDCLine, sys)
    
    # === THERMAL GENERATORS ===
    if !isempty(thermal_gens)
        thermal_formulation_name = get(device_formulations, "thermal_standard", 
                                     formulation_type == "ed" ? "ThermalBasicDispatch" : "ThermalBasicUnitCommitment")
        thermal_formulation = get_device_formulation_object(thermal_formulation_name)
        set_device_model!(template, ThermalStandard, thermal_formulation)
        @info "  ✓ ThermalStandard: $(length(thermal_gens)) units ($thermal_formulation_name)"
    end
    
    # === RENEWABLE DISPATCH ===
    if !isempty(renewable_dispatch)
        renewable_formulation_name = get(device_formulations, "renewable_dispatch", "RenewableFullDispatch")
        renewable_formulation = get_device_formulation_object(renewable_formulation_name)
        set_device_model!(template, RenewableDispatch, renewable_formulation)
        @info "  ✓ RenewableDispatch: $(length(renewable_dispatch)) units ($renewable_formulation_name)"
    end
    
    # === NON-DISPATCHABLE RENEWABLES ===
    if !isempty(renewable_nondispatch)
        renewable_nd_formulation_name = get(device_formulations, "renewable_nondispatch", "FixedOutput")
        renewable_nd_formulation = get_device_formulation_object(renewable_nd_formulation_name)
        set_device_model!(template, RenewableNonDispatch, renewable_nd_formulation)
        @info "  ✓ RenewableNonDispatch: $(length(renewable_nondispatch)) units ($renewable_nd_formulation_name)"
    end
    
    # === HYDRO GENERATORS ===
    if !isempty(hydro_dispatch)
        hydro_formulation_name = get(device_formulations, "hydro_dispatch", "HydroDispatchRunOfRiver")
        hydro_formulation = get_device_formulation_object(hydro_formulation_name)
        set_device_model!(template, HydroDispatch, hydro_formulation)
        @info "  ✓ HydroDispatch: $(length(hydro_dispatch)) units ($hydro_formulation_name)"
    end
    
    # === STORAGE (if present) ===
    try
        battery_storage = get_components(GenericBattery, sys)
        if !isempty(battery_storage)
            battery_formulation_name = get(device_formulations, "generic_battery", "BookKeeping")
            battery_formulation = get_device_formulation_object(battery_formulation_name)
            set_device_model!(template, GenericBattery, battery_formulation)
            @info "  ✓ GenericBattery: $(length(battery_storage)) units ($battery_formulation_name)"
        end
    catch
        # No storage or different storage type
    end
    
    # === LOADS ===
    if !isempty(loads)
        load_formulation_name = get(device_formulations, "power_load", "StaticPowerLoad")
        load_formulation = get_device_formulation_object(load_formulation_name)
        set_device_model!(template, PowerLoad, load_formulation)
        @info "  ✓ PowerLoad: $(length(loads)) units ($load_formulation_name)"
    end
    
    # === TRANSMISSION COMPONENTS ===
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
    
    return template
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
    elseif model_name == "ACPPowerModel"
        return PM.ACPPowerModel
    elseif model_name == "SOCWRPowerModel"
        return PM.SOCWRPowerModel
    else
        available_models = [
            "CopperPlatePowerModel", "PTDFPowerModel", "AreaBalancePowerModel", "AreaPTDFPowerModel",
            "DCPPowerModel", "DCMPPowerModel", "LPACCPowerModel", "ACPPowerModel", "SOCWRPowerModel"
        ]
        error("❌ Invalid network model: '$model_name'\nAvailable: $(join(available_models, ", "))")
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

# ===== PROBLEM BUILDING (FIXED VERSION) =====

"""
    build_problem!(sienna_sim::SiennaSimulations, formulation_type::String)

Build DecisionModel for the specified formulation type.
"""
function build_problem!(sienna_sim::SiennaSimulations, formulation_type::String)
    @info "🔨 Building $formulation_type problem..."

    # Get simulation parameters
    sim_config = get(sienna_sim.config.config_data, "simulations", Dict())
    interval_hours = get(sim_config, "interval_hours", 24)
    lookahead_hours = get(sim_config, "lookahead_hours", 0)
    horizon_hours = interval_hours + lookahead_hours
    resolution = get(sim_config, "resolution", 1)
    
    # Ensure template exists
    if !sienna_sim.templates_created || !haskey(sienna_sim.templates, formulation_type)
        @info "Creating template for $formulation_type..."
        create_templates!(sienna_sim)
    end
    
    sys = get_power_system(sienna_sim.sys)
    config = sienna_sim.config
    template = sienna_sim.templates[formulation_type]
    
    try
        # Create optimizer from config
        optimizer = create_optimizer_from_config(config, formulation_type)
        
        # Create DecisionModel
        @info "Creating DecisionModel..."
        problem = DecisionModel(
            template,
            sys;
            optimizer = optimizer,
            store_variable_names = true,
            resolution = Hour(resolution),
            horizon = Hour(horizon_hours),
            name = "$(formulation_type)_problem"
        )

        sienna_sim.problems[formulation_type] = problem
        @info "✓ DecisionModel created"
        
        # Build the problem
        @info "Building optimization model..."
        output_dir = create_output_directory(config, formulation_type)
        
        build_time = @elapsed begin
            build!(problem, output_dir = output_dir)
        end
        
        sienna_sim.build_times[formulation_type] = build_time
        sienna_sim.problems_built[formulation_type] = true
        
        @info "✅ Problem built successfully in $(round(build_time, digits=2)) seconds"
        
    catch e
        sienna_sim.last_error = "Problem build failed for $formulation_type: $e"
        push!(sienna_sim.simulation_errors[formulation_type], string(e))
        @error "❌ Problem build failed: $e"
        rethrow(e)
    end
end

"""
    build_sequential_simulation!(sienna_sim::SiennaSimulations)

Build sequential simulation using built problems.
"""
function build_sequential_simulation!(sienna_sim::SiennaSimulations)
    @info "🔨 Building sequential simulation from built problems..."

    # Get simulation parameters
    sim_config = get(sienna_sim.config.config_data, "simulations", Dict())
    total_steps = get(sim_config, "total_steps", 365)
    initial_time_str = get(sim_config, "initial_time", nothing)

    if initial_time_str === nothing
        @error "❌ Initial time not found in config. Please set 'initial_time' in your configuration."
        throw(ArgumentError("Initial time must be specified in the configuration."))
    end

    # Check if simulation problems exist
    if isempty(sienna_sim.problems)
        @error "❌ Simulation problems are empty. Run build_problem! first."
        return
    end

    # Sort formulation types by priority and get available problems
    formulation_order = ["mt", "ed", "uc"]
    problems_in_order = [sienna_sim.problems[f] for f in formulation_order if haskey(sienna_sim.problems, f)]

    @info "Creating SimulationModels with $(length(problems_in_order)) decision models..."
    sienna_sim.models = SimulationModels(decision_models = problems_in_order)
    
    # Create SimulationSequence with chronology
    @info "Creating SimulationSequence with chronology..."
    sienna_sim.sequence = SimulationSequence(
        models = sienna_sim.models,
        ini_cond_chronology = InterProblemChronology()
    )
    
    # Create output directory
    output_dir = create_output_directory(sienna_sim.config)
    
    # Create Simulation object
    @info "Creating Simulation object..."
    sienna_sim.simulation = Simulation(
        name = sienna_sim.simulation_name,
        steps = total_steps,
        models = sienna_sim.models,
        sequence = sienna_sim.sequence,
        simulation_folder = output_dir,
        initial_time = DateTime(initial_time_str, "yyyy-mm-dd")
    )
    
    # Build the simulation
    @info "Building sequential simulation..."
    build_time = @elapsed build!(sienna_sim.simulation)
    
    sienna_sim.build_times["simulation"] = build_time
    sienna_sim.problems_built["simulation"] = true
    
    @info "✅ Sequential simulation built successfully in $(round(build_time, digits=2)) seconds"
    @info "   Models: $(length(problems_in_order))"
    @info "   Output: $output_dir"
end

# ===== SIMULATION EXECUTION =====

"""
    run_simulation!(sienna_sim::SiennaSimulations, formulation_types::Vector{String})

Complete simulation workflow: build problems, create simulation, execute, and extract results.
"""
function run_simulation!(sienna_sim::SiennaSimulations, formulation_types::Vector{String}=["ed"])
    @info "🚀 Running complete simulation workflow..."
    
    start_time = now()
    
    try
        # Step 1: Build individual problems
        for formulation_type in formulation_types
            if !should_run_formulation(sienna_sim.config, formulation_type)
                @warn "$formulation_type is disabled in configuration"
                continue
            end
            
            build_problem!(sienna_sim, formulation_type)
        end

        # Step 2: Build sequential simulation from problems
        build_sequential_simulation!(sienna_sim)
        
        # Step 3: Execute simulation
        @info "⚡ Executing simulation..."
        solve_time = @elapsed execute!(sienna_sim.simulation; enable_progress_bar = true)
        sienna_sim.solve_times["simulation"] = solve_time
        
        # Step 4: Extract results
        @info "📊 Extracting simulation results..."
        sienna_sim.results = SimulationResults(sienna_sim.simulation)
        sienna_sim.simulations_completed["simulation"] = true
        
        # Step 5: Save results
        save_simulation_results(sienna_sim)
        
        # Step 6: Print summary
        print_simulation_summary(sienna_sim)
        
        total_time = Dates.value(now() - start_time) / 1000.0
        @info "✅ Simulation completed in $(round(total_time, digits=1)) seconds"
        
        return sienna_sim.results
        
    catch e
        sienna_sim.last_error = "Simulation execution failed: $e"
        push!(sienna_sim.simulation_errors["simulation"], string(e))
        @error "❌ Simulation execution failed: $e"
        rethrow(e)
    end
end

# ===== OPTIMIZER AND OUTPUT UTILITIES =====

"""
    create_optimizer_from_config(config::SiennaConfig, formulation_type::String)

Create optimizer using configuration solver settings.
"""
function create_optimizer_from_config(config::SiennaConfig, formulation_type::String)
    solver_settings = get_solver_settings(config, formulation_type)
    
    @info "Optimizer settings from config:"
    @info "  Time limit: $(solver_settings["time_limit"])s"
    @info "  MIP gap: $(solver_settings["mip_gap"])"
    @info "  Threads: $(solver_settings["threads"] == 0 ? "auto" : solver_settings["threads"])"
    
    # Determine thread count
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
    create_output_directory(config::SiennaConfig, formulation_type::Union{String, Nothing}=nothing)

Create timestamped output directory for simulation results.
"""
function create_output_directory(config::SiennaConfig, formulation_type::Union{String, Nothing}=nothing)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    if isnothing(formulation_type)
        output_dir = joinpath(config.output_directory, "sim_$(timestamp)")
    else
        output_dir = joinpath(config.output_directory, "$(formulation_type)_$(timestamp)")
    end
    
    if !isdir(output_dir)
        mkpath(output_dir)
        @info "Created output directory: $output_dir"
    end
    
    return output_dir
end

# ===== RESULTS HANDLING =====

"""
    save_simulation_results(sienna_sim::SiennaSimulations)

Save simulation results with configuration information.
"""
function save_simulation_results(sienna_sim::SiennaSimulations)
    @info "💾 Saving simulation results..."
    
    if sienna_sim.results === nothing
        @warn "No results to save"
        return
    end
    
    # Create results directory
    output_dir = create_output_directory(sienna_sim.config, "simulation")
    
    try
        # Save results summary
        summary = create_simulation_results_summary(sienna_sim)
        summary_file = joinpath(output_dir, "simulation_results_summary.json")
        open(summary_file, "w") do f
            JSON3.pretty(f, summary)
        end
        
        # Copy configuration to results
        config_file = joinpath(output_dir, "config_used.toml")
        export_config(sienna_sim.config, config_file)
        
        # Save system statistics
        sys_stats = get_system_statistics(sienna_sim.sys)
        stats_file = joinpath(output_dir, "system_statistics.json")
        open(stats_file, "w") do f
            JSON3.pretty(f, sys_stats)
        end
        
        # Save variables based on config settings
        if sienna_sim.config.save_csv_variables
            save_simulation_variables_csv(sienna_sim, output_dir)
        end
        
        @info "✓ Results saved to: $output_dir"
        
    catch e
        @warn "Failed to save results: $e"
    end
end

"""
    create_simulation_results_summary(sienna_sim::SiennaSimulations)

Create comprehensive simulation results summary.
"""
function create_simulation_results_summary(sienna_sim::SiennaSimulations)
    sys = get_power_system(sienna_sim.sys)
    config = sienna_sim.config
    
    # Get simulation statistics
    simulation_stats = try
        if sienna_sim.results !== nothing
            timestamps = get_realized_timestamps(sienna_sim.results)
            Dict(
                "total_steps" => length(timestamps),
                "decision_models" => length(sienna_sim.models.decision_models),
                "timestamps_start" => string(first(timestamps)),
                "timestamps_end" => string(last(timestamps))
            )
        else
            Dict("error" => "No simulation results available")
        end
    catch e
        Dict("error" => "Failed to get simulation statistics: $e")
    end
    
    return Dict(
        "simulation_metadata" => Dict(
            "simulation_type" => "sequential_simulation",
            "timestamp" => string(now()),
            "simulation_name" => sienna_sim.simulation_name,
            "system_name" => get_name(sys)
        ),
        "performance" => Dict(
            "build_time_seconds" => get(sienna_sim.build_times, "simulation", 0.0),
            "solve_time_seconds" => get(sienna_sim.solve_times, "simulation", 0.0),
            "total_time_seconds" => get(sienna_sim.build_times, "simulation", 0.0) + 
                                  get(sienna_sim.solve_times, "simulation", 0.0)
        ),
        "simulation_results" => simulation_stats,
        "config_summary" => Dict(
            "network_model" => config.network_model,
            "solver" => config.default_solver,
            "transmission_limits" => config.enable_transmission_limits,
            "time_series_loaded" => has_time_series_data(sienna_sim.sys),
            "forecasts_available" => has_forecast_data(sienna_sim.sys)
        ),
        "decision_models" => Dict(
            "available_models" => haskey(sienna_sim, :models) && sienna_sim.models !== nothing ? 
                                [get_name(model) for model in sienna_sim.models.decision_models] : [],
            "total_models" => haskey(sienna_sim, :models) && sienna_sim.models !== nothing ? 
                            length(sienna_sim.models.decision_models) : 0
        ),
        "versions" => Dict(
            "powersimulations" => "v0.30.2",
            "powersystems" => "v4.6.2",
            "sienna_framework" => "v2.0"
        )
    )
end

"""
    save_simulation_variables_csv(sienna_sim::SiennaSimulations, output_dir::String)

Save simulation variables as CSV files.
"""
function save_simulation_variables_csv(sienna_sim::SiennaSimulations, output_dir::String)
    if sienna_sim.results === nothing
        @warn "No simulation results to save"
        return
    end
    
    variables_dir = joinpath(output_dir, "variables")
    mkpath(variables_dir)
    
    try
        # Get available variable names from simulation results
        available_variables = list_variable_names(sienna_sim.results)
        
        if isempty(available_variables)
            @warn "No variables found in simulation results"
            return
        end
        
        @info "Saving $(length(available_variables)) simulation variables..."
        
        for var_name in available_variables
            try
                # Read variable data from simulation results
                var_data = read_variable(sienna_sim.results, var_name)
                var_file = joinpath(variables_dir, "$(var_name).csv")
                CSV.write(var_file, var_data)
                @info "  ✓ Saved: $var_name"
            catch e
                @warn "  ✗ Failed to save variable $var_name: $e"
            end
        end
        
    catch e
        @warn "Failed to save simulation variables: $e"
    end
end

# ===== ANALYSIS AND SUMMARY =====

"""
    print_simulation_summary(sienna_sim::SiennaSimulations)

Print detailed simulation summary.
"""
function print_simulation_summary(sienna_sim::SiennaSimulations)
    @info "\n" * "="^60
    @info "SIMULATION EXECUTION SUMMARY"
    @info "="^60
    
    if sienna_sim.results !== nothing
        try
            timestamps = get_realized_timestamps(sienna_sim.results)
            @info "Simulation Steps: $(length(timestamps))"
            @info "Time Range: $(first(timestamps)) to $(last(timestamps))"
            if sienna_sim.models !== nothing
                @info "Decision Models: $(length(sienna_sim.models.decision_models))"
            end
        catch e
            @warn "Could not retrieve simulation statistics: $e"
        end
    else
        @warn "No simulation results available"
    end
    
    # Performance metrics
    build_time = get(sienna_sim.build_times, "simulation", 0.0)
    solve_time = get(sienna_sim.solve_times, "simulation", 0.0)
    @info "Build Time: $(round(build_time, digits=2)) seconds"
    @info "Solve Time: $(round(solve_time, digits=2)) seconds"
    @info "Total Time: $(round(build_time + solve_time, digits=2)) seconds"
    
    # System info
    @info "System: $(get_name(get_power_system(sienna_sim.sys)))"
    @info "Network Model: $(sienna_sim.config.network_model)"
    
    @info "="^60
end

# ===== PUBLIC INTERFACE =====

"""
    get_results(sienna_sim::SiennaSimulations)

Get simulation results.
"""
function get_results(sienna_sim::SiennaSimulations)
    if sienna_sim.results === nothing
        error("❌ No simulation results available. Run simulation first.")
    end
    return sienna_sim.results
end

"""
    has_results(sienna_sim::SiennaSimulations)

Check if simulation results are available.
"""
function has_results(sienna_sim::SiennaSimulations)
    return sienna_sim.results !== nothing
end

"""
    get_variable_data(sienna_sim::SiennaSimulations, variable_name::String)

Get specific variable data from simulation results.
"""
function get_variable_data(sienna_sim::SiennaSimulations, variable_name::String)
    if sienna_sim.results === nothing
        error("❌ No simulation results available")
    end
    
    try
        return read_variable(sienna_sim.results, variable_name)
    catch e
        error("❌ Failed to read variable '$variable_name': $e")
    end
end

"""
    list_available_variables(sienna_sim::SiennaSimulations)

List all available variables in simulation results.
"""
function list_available_variables(sienna_sim::SiennaSimulations)
    if sienna_sim.results === nothing
        @warn "No simulation results available"
        return String[]
    end
    
    try
        return list_variable_names(sienna_sim.results)
    catch e
        @warn "Failed to get variable names: $e"
        return String[]
    end
end

"""
    get_performance_summary(sienna_sim::SiennaSimulations)

Get performance summary for all completed simulations.
"""
function get_performance_summary(sienna_sim::SiennaSimulations)
    return Dict(
        "build_times" => copy(sienna_sim.build_times),
        "solve_times" => copy(sienna_sim.solve_times),
        "total_build_time" => sum(values(sienna_sim.build_times)),
        "total_solve_time" => sum(values(sienna_sim.solve_times)),
        "simulations_completed" => copy(sienna_sim.simulations_completed),
        "errors" => copy(sienna_sim.simulation_errors)
    )
end

"""
    reset_simulations!(sienna_sim::SiennaSimulations)

Reset all simulation state (useful for re-running with different parameters).
"""
function reset_simulations!(sienna_sim::SiennaSimulations)
    @info "🔄 Resetting simulation state..."
    
    # Clear all simulation state
    sienna_sim.templates = Dict{String, ProblemTemplate}()
    sienna_sim.problems = Dict{String, DecisionModel}()
    sienna_sim.results = nothing
    sienna_sim.simulation = nothing
    sienna_sim.models = nothing
    sienna_sim.sequence = nothing
    sienna_sim.build_times = Dict{String, Float64}()
    sienna_sim.solve_times = Dict{String, Float64}()
    
    # Reset status flags
    sienna_sim.templates_created = false
    for formulation in ["ed", "uc", "simulation"]
        sienna_sim.problems_built[formulation] = false
        sienna_sim.simulations_completed[formulation] = false
        sienna_sim.simulation_errors[formulation] = String[]
    end
    
    sienna_sim.last_error = nothing
    
    @info "✓ Simulation state reset - ready for new simulations"
end

# ===== ANALYSIS HELPERS =====

"""
    analyze_generation_dispatch(sienna_sim::SiennaSimulations)

Analyze generation dispatch from simulation results.
"""
function analyze_generation_dispatch(sienna_sim::SiennaSimulations)
    @info "\n📊 Generation Dispatch Analysis:"
    
    if !has_results(sienna_sim)
        @warn "No simulation results available for analysis"
        return
    end
    
    try
        var_names = list_available_variables(sienna_sim)
        total_generation = Dict{String, Float64}()
        
        # Thermal generation
        if "ActivePowerVariable__ThermalStandard" in var_names
            thermal_data = get_variable_data(sienna_sim, "ActivePowerVariable__ThermalStandard")
            if ncol(thermal_data) > 1
                total_generation["Thermal"] = sum(skipmissing(Matrix(thermal_data[:, 2:end])))
                max_thermal = maximum(skipmissing(Matrix(thermal_data[:, 2:end])))
                @info "  Thermal: $(round(total_generation["Thermal"], digits=1)) MWh total, $(round(max_thermal, digits=1)) MW peak"
            end
        end
        
        # Renewable generation
        if "ActivePowerVariable__RenewableDispatch" in var_names
            renewable_data = get_variable_data(sienna_sim, "ActivePowerVariable__RenewableDispatch")
            if ncol(renewable_data) > 1
                total_generation["Renewable"] = sum(skipmissing(Matrix(renewable_data[:, 2:end])))
                max_renewable = maximum(skipmissing(Matrix(renewable_data[:, 2:end])))
                @info "  Renewable: $(round(total_generation["Renewable"], digits=1)) MWh total, $(round(max_renewable, digits=1)) MW peak"
            end
        end
        
        # Hydro generation
        if "ActivePowerVariable__HydroDispatch" in var_names
            hydro_data = get_variable_data(sienna_sim, "ActivePowerVariable__HydroDispatch")
            if ncol(hydro_data) > 1
                total_generation["Hydro"] = sum(skipmissing(Matrix(hydro_data[:, 2:end])))
                max_hydro = maximum(skipmissing(Matrix(hydro_data[:, 2:end])))
                @info "  Hydro: $(round(total_generation["Hydro"], digits=1)) MWh total, $(round(max_hydro, digits=1)) MW peak"
            end
        end
        
        # Generation mix summary
        if !isempty(total_generation)
            total = sum(values(total_generation))
            @info "  Total Generation: $(round(total, digits=1)) MWh"
            
            for (tech, gen) in total_generation
                if gen > 0
                    pct = (gen / total) * 100
                    @info "  $(tech) share: $(round(pct, digits=1))%"
                end
            end
        else
            @warn "  No generation variables found in results"
        end
        
    catch e
        @warn "Generation analysis failed: $e"
    end
end

"""
    analyze_unit_commitment_decisions(sienna_sim::SiennaSimulations)

Analyze unit commitment decisions from simulation results.
"""
function analyze_unit_commitment_decisions(sienna_sim::SiennaSimulations)
    @info "\n🔧 Unit Commitment Analysis:"
    
    if !has_results(sienna_sim)
        @warn "No simulation results available for analysis"
        return
    end
    
    try
        var_names = list_available_variables(sienna_sim)
        
        if "OnVariable__ThermalStandard" in var_names
            on_status = get_variable_data(sienna_sim, "OnVariable__ThermalStandard")
            thermal_gens = get_components(ThermalStandard, get_power_system(sienna_sim.sys))
            
            if ncol(on_status) > 1
                total_units = length(thermal_gens)
                hourly_committed = [sum(skipmissing(on_status[h, 2:end])) for h in 1:nrow(on_status)]
                
                avg_committed = mean(hourly_committed)
                max_committed = maximum(hourly_committed)
                min_committed = minimum(hourly_committed)
                
                @info "  Total thermal units: $total_units"
                @info "  Average committed: $(round(avg_committed, digits=1)) units ($(round(avg_committed/total_units*100, digits=1))%)"
                @info "  Peak commitment: $max_committed units ($(round(max_committed/total_units*100, digits=1))%)"
                @info "  Minimum commitment: $min_committed units ($(round(min_committed/total_units*100, digits=1))%)"
            end
        end
        
        # Start-up analysis
        if "StartVariable__ThermalStandard" in var_names
            start_data = get_variable_data(sienna_sim, "StartVariable__ThermalStandard")
            if ncol(start_data) > 1
                total_starts = sum(skipmissing(Matrix(start_data[:, 2:end])))
                @info "  Total start-ups: $(Int(total_starts))"
            end
        end
        
    catch e
        @warn "Unit commitment analysis failed: $e"
    end
end

# ===== CONVENIENCE FUNCTIONS =====

"""
    run_economic_dispatch_simulation(config_file::String="config.toml")

Convenience function to run economic dispatch simulation.
"""
function run_economic_dispatch_simulation(config_file::String="config.toml")
    @info "🚀 Running economic dispatch simulation..."
    
    sim = SiennaSimulations(config_file)
    return run_simulation!(sim, ["ed"])
end

"""
    run_unit_commitment_simulation(config_file::String="config.toml")

Convenience function to run unit commitment simulation.
"""
function run_unit_commitment_simulation(config_file::String="config.toml")
    @info "🚀 Running unit commitment simulation..."
    
    sim = SiennaSimulations(config_file)
    return run_simulation!(sim, ["uc"])
end

"""
    run_multi_stage_simulation(config_file::String="config.toml")

Convenience function to run multi-stage simulation.
"""
function run_multi_stage_simulation(config_file::String="config.toml")
    @info "🚀 Running multi-stage simulation..."
    
    sim = SiennaSimulations(config_file)
    return run_simulation!(sim, ["ed", "uc"])
end

"""
    quick_analysis(sienna_sim::SiennaSimulations)

Quick analysis of simulation results.
"""
function quick_analysis(sienna_sim::SiennaSimulations)
    @info "🔍 Quick simulation analysis..."
    
    if !has_results(sienna_sim)
        @warn "No results to analyze. Run simulation first."
        return
    end
    
    # Print basic info
    print_simulation_summary(sienna_sim)
    
    # Analyze generation
    analyze_generation_dispatch(sienna_sim)
    
    # Analyze unit commitment if available
    variables = list_available_variables(sienna_sim)
    if "OnVariable__ThermalStandard" in variables
        analyze_unit_commitment_decisions(sienna_sim)
    end
    
    # List all available variables
    @info "\n📋 Available Variables ($(length(variables))):"
    for (i, var) in enumerate(variables)
        @info "  $i. $var"
        if i >= 10
            @info "  ... and $(length(variables) - 10) more variables"
            break
        end
    end
end

# ===== BACKWARDS COMPATIBILITY FUNCTIONS =====

"""
    timeseries_loaded(sienna_sys::SiennaSystem)

Check if time series are loaded (backwards compatibility).
"""
function timeseries_loaded(sienna_sys::SiennaSystem)
    return has_time_series_data(sienna_sys)
end

"""
    forecasts_available(sienna_sys::SiennaSystem)

Check if forecasts are available (backwards compatibility).
"""
function forecasts_available(sienna_sys::SiennaSystem)
    return has_forecast_data(sienna_sys)
end

# ===== EXPORTS =====

export SiennaSimulations
export create_templates!, build_problem!, build_sequential_simulation!
export run_simulation!
export get_results, has_results, get_variable_data, list_available_variables
export get_performance_summary, reset_simulations!
export analyze_generation_dispatch, analyze_unit_commitment_decisions
export print_simulation_summary, quick_analysis
export run_economic_dispatch_simulation, run_unit_commitment_simulation, run_multi_stage_simulation
export timeseries_loaded, forecasts_available

# ===== TESTING =====
if abspath(PROGRAM_FILE) == @__FILE__
    @info "🧪 Testing complete SiennaSimulations..."
    
    try
        if isfile("config.toml")
            @info "Running simulation execution test..."
            sim = SiennaSimulations("config.toml")
            @info "✅ Complete SiennaSimulations loaded successfully"
            @info "   Available functions: create_templates!, build_problem!, run_simulation!"
            @info "   Analysis functions: quick_analysis, analyze_generation_dispatch"
            @info "   Convenience functions: run_economic_dispatch_simulation"
        else
            @info "No config.toml found for testing"
            @info "Expected usage:"
            @info "  sim = SiennaSimulations(\"config.toml\")"
            @info "  results = run_simulation!(sim, [\"ed\", \"uc\"])"
            @info "  quick_analysis(sim)"
        end
    catch e
        @error "❌ Complete adaptation test failed: $e"
    end
end
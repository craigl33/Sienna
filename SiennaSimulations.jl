#!/usr/bin/env julia

"""
SiennaSimulations.jl - PowerSimulations.jl Framework for Sienna Ecosystem
=========================================================================

Class-based simulation manager that encapsulates PowerSimulations.jl functionality.
Provides clean interface for creating templates, problems, and running simulations.

CRITICAL: Use proper Julia function call syntax, NOT method calls:
✅ CORRECT:                          ❌ WRONG:
get_device_formulations(config, x)   config.get_device_formulations(x)
get_solver_settings(config, x)       config.get_solver_settings(x)
should_run_formulation(config, x)    config.should_run_formulation(x)
get_output_directory(config)         config.get_output_directory()
is_system_built(sienna_system)       sienna_system.is_built()
timeseries_loaded(sienna_system)     sienna_system.has_time_series()
get_system_statistics(sienna_system) sienna_system.get_system_statistics()
export_config(config, file)          config.export_config(file)

Features:
- Config-driven template creation
- Problem management with optimization
- Results handling and analysis
- Multi-formulation simulation suites
- Performance monitoring

Usage:
    config = SiennaConfig("config.toml")
    sienna_sys = SiennaSystem(config)
    build_system!(sienna_sys)
    
    sienna_sim = SiennaSimulations(config, sienna_sys)
    results = run_simulation!(sienna_sim, "economic_dispatch")
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
using Infiltrator

# Import our modules
include("SiennaConfig.jl")
include("SiennaSystem.jl")

"""
    SiennaSimulations

Main simulation manager class that encapsulates PowerSimulations.jl functionality.
"""
mutable struct SiennaSimulations
    # Core components
    config::SiennaConfig
    sienna_system::SiennaSystem
    
    # Templates and problems
    templates::Dict{String, ProblemTemplate}
    problems::Dict{String, DecisionModel}
    results::Dict{String, OptimizationProblemResults}
    
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
    SiennaSimulations(config::SiennaConfig, sienna_system::SiennaSystem)

Constructor - Initialize simulation manager with config and system.
"""
function SiennaSimulations(config::SiennaConfig, sienna_system::SiennaSystem)
    @info "⚡ Initializing SiennaSimulations framework"
    
    # Validate inputs
    if !config.is_validated
        error("❌ Configuration must be validated before creating SiennaSimulations")
    end
    
    if !is_system_built(sienna_system)
        error("❌ SiennaSystem must be built before creating SiennaSimulations")
    end
    
    # Initialize with default values
    sienna_sim = SiennaSimulations(
        config,
        sienna_system,
        Dict{String, ProblemTemplate}(),      # templates
        Dict{String, DecisionModel}(),        # problems
        Dict{String, OptimizationProblemResults}(),  # results
        "$(config.project_name)_simulation",  # simulation name
        now(),                                # timestamp
        Dict{String, Float64}(),              # build times
        Dict{String, Float64}(),              # solve times
        false,                                # templates created
        Dict{String, Bool}(),                 # problems built
        Dict{String, Bool}(),                 # simulations completed
        Dict{String, Vector{String}}(),       # simulation errors
        nothing                               # last error
    )
    
    # Initialize tracking dictionaries
    formulation_types = ["economic_dispatch", "unit_commitment"]
    for formulation in formulation_types
        sienna_sim.problems_built[formulation] = false
        sienna_sim.simulations_completed[formulation] = false
        sienna_sim.simulation_errors[formulation] = String[]
    end
    
    @info "✅ SiennaSimulations initialized successfully"
    @info "   System: $(get_name(get_power_system(sienna_system)))"
    @info "   Network Model: $(config.network_model)"
    
    return sienna_sim
end

"""
    create_templates!(sienna_sim::SiennaSimulations)

Create PowerSimulations templates for all formulation types.
"""
function create_templates!(sienna_sim::SiennaSimulations)
    @info "📋 Creating PowerSimulations templates..."
    
    sys = get_power_system(sienna_sim.sienna_system)
    config = sienna_sim.config
    
    try
        # Create templates for each formulation type
        if should_run_formulation(config, "economic_dispatch")
            @info "Creating Economic Dispatch template..."
            sienna_sim.templates["economic_dispatch"] = create_config_driven_template(
                sys, "economic_dispatch", config
            )
            @info "✓ Economic Dispatch template created"
        end
        
        if should_run_formulation(config, "unit_commitment")
            @info "Creating Unit Commitment template..."
            sienna_sim.templates["unit_commitment"] = create_config_driven_template(
                sys, "unit_commitment", config
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
                                     formulation_type == "economic_dispatch" ? "ThermalBasicDispatch" : "ThermalBasicUnitCommitment")
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

"""
    build_problem!(sienna_sim::SiennaSimulations, formulation_type::String)

Build DecisionModel for the specified formulation type.
"""
function build_problem!(sienna_sim::SiennaSimulations, formulation_type::String)
    @info "🔨 Building $formulation_type problem..."
    
    # Ensure template exists
    if !sienna_sim.templates_created || !haskey(sienna_sim.templates, formulation_type)
        @info "Creating template for $formulation_type..."
        create_templates!(sienna_sim)
    end
    
    sys = get_power_system(sienna_sim.sienna_system)
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
            store_variable_names=true,
            horizon = Hour(config.default_horizon_hours),
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
    solve_problem!(sienna_sim::SiennaSimulations, formulation_type::String)

Solve the optimization problem for the specified formulation type.
"""
function solve_problem!(sienna_sim::SiennaSimulations, formulation_type::String)
    @info "⚡ Solving $formulation_type problem..."
    
    # Ensure problem is built
    if !get(sienna_sim.problems_built, formulation_type, false)
        @info "Building problem for $formulation_type..."
        build_problem!(sienna_sim, formulation_type)
    end
    
    problem = sienna_sim.problems[formulation_type]
    
    try
        # Solve the problem
        solve_time = @elapsed begin
            solve!(problem)
        end
        
        sienna_sim.solve_times[formulation_type] = solve_time
        @info "✓ Problem solved in $(round(solve_time, digits=2)) seconds"
        
        # Extract results
        @info "Extracting results..."
        results = OptimizationProblemResults(problem)
        sienna_sim.results[formulation_type] = results
        
        sienna_sim.simulations_completed[formulation_type] = true
        @info "✅ $formulation_type completed successfully"
        
        return results
        
    catch e
        sienna_sim.last_error = "Solve failed for $formulation_type: $e"
        push!(sienna_sim.simulation_errors[formulation_type], string(e))
        @error "❌ Solve failed: $e"
        rethrow(e)
    end
end

"""
    run_simulation!(sienna_sim::SiennaSimulations, formulation_type::String)

Complete simulation workflow: build problem, solve, and extract results.
"""
function run_simulation!(sienna_sim::SiennaSimulations, formulation_type::String)
    @info "🚀 Running complete $formulation_type simulation..."
    
    start_time = now()
    
    try
        # Validate formulation type
        if !should_run_formulation(sienna_sim.config, formulation_type)
            @warn "$formulation_type is disabled in configuration"
            return nothing
        end
        
        # Build and solve
        build_problem!(sienna_sim, formulation_type)
        results = solve_problem!(sienna_sim, formulation_type)
        
        # Save results
        save_simulation_results(sienna_sim, formulation_type, results)
        
        # Print summary
        print_simulation_summary(sienna_sim, formulation_type, results)
        
        total_time = Dates.value(now() - start_time) / 1000.0
        @info "✅ $formulation_type simulation completed in $(round(total_time, digits=1)) seconds"
        
        return results
        
    catch e
        sienna_sim.last_error = "Simulation failed for $formulation_type: $e"
        push!(sienna_sim.simulation_errors[formulation_type], string(e))
        @error "❌ $formulation_type simulation failed: $e"
        rethrow(e)
    end
end

"""
    run_simulation_suite!(sienna_sim::SiennaSimulations)

Run complete simulation suite for all enabled formulation types.
"""
function run_simulation_suite!(sienna_sim::SiennaSimulations)
    @info "🎯 Running complete simulation suite..."
    
    suite_start = now()
    suite_results = Dict{String, Any}()
    
    # Economic Dispatch
    if should_run_formulation(sienna_sim.config, "economic_dispatch")
        @info "\n>>> Economic Dispatch"
        try
            ed_results = run_simulation!(sienna_sim, "economic_dispatch")
            suite_results["economic_dispatch"] = Dict(
                "status" => "SUCCESS",
                "objective_value" => get_objective_value(ed_results),
                "solve_time" => sienna_sim.solve_times["economic_dispatch"],
                "build_time" => sienna_sim.build_times["economic_dispatch"]
            )
        catch e
            @error "Economic Dispatch failed: $e"
            suite_results["economic_dispatch"] = Dict(
                "status" => "FAILED",
                "error" => string(e)
            )
        end
    end
    
    # Unit Commitment
    if should_run_formulation(sienna_sim.config, "unit_commitment")
        @info "\n>>> Unit Commitment"
        try
            uc_results = run_simulation!(sienna_sim, "unit_commitment")
            suite_results["unit_commitment"] = Dict(
                "status" => "SUCCESS",
                "objective_value" => get_objective_value(uc_results),
                "solve_time" => sienna_sim.solve_times["unit_commitment"],
                "build_time" => sienna_sim.build_times["unit_commitment"]
            )
        catch e
            @error "Unit Commitment failed: $e"
            suite_results["unit_commitment"] = Dict(
                "status" => "FAILED",
                "error" => string(e)
            )
        end
    end
    
    # Save suite summary
    save_suite_summary(sienna_sim, suite_results)
    
    # Print final summary
    print_suite_summary(sienna_sim, suite_results, suite_start)
    
    return suite_results
end

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
    create_output_directory(config::SiennaConfig, formulation_type::String)

Create timestamped output directory for simulation results.
"""
function create_output_directory(config::SiennaConfig, formulation_type::String)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    output_dir = joinpath(config.output_directory, "$(formulation_type)_$(timestamp)")
    mkpath(output_dir)
    return output_dir
end

"""
    save_simulation_results(sienna_sim::SiennaSimulations, formulation_type::String, results::OptimizationProblemResults)

Save simulation results with configuration information.
"""
function save_simulation_results(sienna_sim::SiennaSimulations, formulation_type::String, results::OptimizationProblemResults)
    @info "💾 Saving $formulation_type results..."
    
    # Create results directory
    output_dir = create_output_directory(sienna_sim.config, formulation_type)
    
    try
        # Save results summary
        summary = create_results_summary(sienna_sim, formulation_type, results)
        summary_file = joinpath(output_dir, "results_summary.json")
        open(summary_file, "w") do f
            JSON3.pretty(f, summary)
        end
        
        # Copy configuration to results
        config_file = joinpath(output_dir, "config_used.toml")
        export_config(sienna_sim.config, config_file)
        
        # Save system statistics
        sys_stats = get_system_statistics(sienna_sim.sienna_system)
        stats_file = joinpath(output_dir, "system_statistics.json")
        open(stats_file, "w") do f
            JSON3.pretty(f, sys_stats)
        end
        
        # Save variables based on config settings
        if sienna_sim.config.save_csv_variables
            save_variables_csv(results, output_dir, sienna_sim.config)
        end
        
        @info "✓ Results saved to: $output_dir"
        
    catch e
        @warn "Failed to save results: $e"
    end
end

"""
    create_results_summary(sienna_sim::SiennaSimulations, formulation_type::String, results::OptimizationProblemResults)

Create comprehensive results summary.
"""
function create_results_summary(sienna_sim::SiennaSimulations, formulation_type::String, results::OptimizationProblemResults)
    sys = get_power_system(sienna_sim.sienna_system)
    config = sienna_sim.config
    
    # Get optimizer statistics
    optimizer_stats = try
        stats = get_optimizer_stats(results)
        if nrow(stats) > 0
            Dict(
                "solve_time" => stats.solve_time[1],
                "termination_status" => string(stats.termination_status[1]),
                "primal_status" => string(stats.primal_status[1]),
                "dual_status" => string(stats.dual_status[1])
            )
        else
            Dict("error" => "No optimizer statistics available")
        end
    catch e
        Dict("error" => "Failed to get optimizer statistics: $e")
    end
    
    return Dict(
        "simulation_metadata" => Dict(
            "formulation_type" => formulation_type,
            "timestamp" => string(now()),
            "simulation_name" => sienna_sim.simulation_name,
            "system_name" => get_name(sys),
            "horizon_hours" => config.default_horizon_hours
        ),
        "performance" => Dict(
            "build_time_seconds" => get(sienna_sim.build_times, formulation_type, 0.0),
            "solve_time_seconds" => get(sienna_sim.solve_times, formulation_type, 0.0),
            "total_time_seconds" => get(sienna_sim.build_times, formulation_type, 0.0) + 
                                  get(sienna_sim.solve_times, formulation_type, 0.0)
        ),
        "optimization_results" => Dict(
            "objective_value" => try; get_objective_value(results); catch e; "Error: $e"; end,
            "optimizer_stats" => optimizer_stats
        ),
        "config_summary" => Dict(
            "network_model" => config.network_model,
            "solver" => config.default_solver,
            "transmission_limits" => config.enable_transmission_limits,
            "time_series_loaded" => timeseries_loaded(sienna_sim.sienna_system),
            "forecasts_available" => forecasts_available(sienna_sim.sienna_system)
        ),
        "variable_summary" => Dict(
            "available_variables" => try; list_variable_names(results); catch; []; end,
            "available_parameters" => try; list_parameter_names(results); catch; []; end
        ),
        "versions" => Dict(
            "powersimulations" => "v0.30.2",
            "powersystems" => "v4.6.2",
            "sienna_framework" => "v2.0"
        )
    )
end

"""
    save_variables_csv(results::OptimizationProblemResults, output_dir::String, config::SiennaConfig)

Save optimization variables as CSV files.
"""
function save_variables_csv(results::OptimizationProblemResults, output_dir::String, config::SiennaConfig)
    var_names = try
        list_variable_names(results)
    catch e
        @warn "Failed to get variable names: $e"
        return
    end
    
    if isempty(var_names)
        @warn "No variables found in results"
        return
    end
    
    variables_dir = joinpath(output_dir, "variables")
    mkpath(variables_dir)
    
    if config.save_all_variables
        @info "  Saving all $(length(var_names)) variables..."
        for var_name in var_names
            try
                var_data = read_variable(results, var_name)
                var_file = joinpath(variables_dir, "$(var_name).csv")
                CSV.write(var_file, var_data)
            catch e
                @warn "  Failed to save variable $var_name: $e"
            end
        end
    else
        # Save specific variables (if configured)
        specific_vars = get(config.config_data["output"], "specific_variables", [])
        @info "  Saving $(length(specific_vars)) specific variables..."
        for var_name in specific_vars
            if var_name in var_names
                try
                    var_data = read_variable(results, var_name)
                    var_file = joinpath(variables_dir, "$(var_name).csv")
                    CSV.write(var_file, var_data)
                catch e
                    @warn "  Failed to save variable $var_name: $e"
                end
            else
                @warn "  Requested variable $var_name not found in results"
            end
        end
    end
end

"""
    save_suite_summary(sienna_sim::SiennaSimulations, suite_results::Dict)

Save simulation suite summary.
"""
function save_suite_summary(sienna_sim::SiennaSimulations, suite_results::Dict)
    @info "💾 Saving simulation suite summary..."
    
    output_dir = sienna_sim.config.output_directory
    summary_file = joinpath(output_dir, "simulation_suite_summary.json")
    
    suite_summary = Dict(
        "suite_metadata" => Dict(
            "simulation_name" => sienna_sim.simulation_name,
            "timestamp" => string(now()),
            "system_name" => get_name(get_power_system(sienna_sim.sienna_system)),
            "project_name" => sienna_sim.config.project_name
        ),
        "configuration" => Dict(
            "network_model" => sienna_sim.config.network_model,
            "horizon_hours" => sienna_sim.config.default_horizon_hours,
            "solver" => sienna_sim.config.default_solver,
            "formulations_enabled" => [
                formulation for formulation in ["economic_dispatch", "unit_commitment"]
                if should_run_formulation(sienna_sim.config, formulation)
            ]
        ),
        "results" => suite_results,
        "performance_summary" => Dict(
            "total_build_time" => sum(values(sienna_sim.build_times)),
            "total_solve_time" => sum(values(sienna_sim.solve_times)),
            "build_times" => sienna_sim.build_times,
            "solve_times" => sienna_sim.solve_times
        ),
        "error_summary" => Dict(
            "errors_occurred" => any(!isempty(errors) for errors in values(sienna_sim.simulation_errors)),
            "error_details" => sienna_sim.simulation_errors
        )
    )
    
    try
        open(summary_file, "w") do f
            JSON3.pretty(f, suite_summary)
        end
        @info "✓ Suite summary saved to: $summary_file"
    catch e
        @warn "Failed to save suite summary: $e"
    end
end

# ===== ANALYSIS AND REPORTING METHODS =====

"""
    print_simulation_summary(sienna_sim::SiennaSimulations, formulation_type::String, results::OptimizationProblemResults)

Print detailed simulation summary.
"""
function print_simulation_summary(sienna_sim::SiennaSimulations, formulation_type::String, results::OptimizationProblemResults)
    @info "\n" * "="^60
    @info "$(uppercase(formulation_type)) SIMULATION SUMMARY"
    @info "="^60
    
    # Objective value
    try
        obj_value = get_objective_value(results)
        @info "Objective Value: \$(round(obj_value, digits=2))"
    catch e
        @warn "Could not retrieve objective value: $e"
    end
    
    # Performance metrics
    build_time = get(sienna_sim.build_times, formulation_type, 0.0)
    solve_time = get(sienna_sim.solve_times, formulation_type, 0.0)
    @info "Build Time: $(round(build_time, digits=2)) seconds"
    @info "Solve Time: $(round(solve_time, digits=2)) seconds"
    @info "Total Time: $(round(build_time + solve_time, digits=2)) seconds"
    
    # Optimizer statistics
    try
        stats = get_optimizer_stats(results)
        if nrow(stats) > 0
            @info "Termination Status: $(stats.termination_status[1])"
            @info "Primal Status: $(stats.primal_status[1])"
        end
    catch e
        @warn "Could not retrieve optimizer stats: $e"
    end
    
    # Generation analysis
    analyze_generation_dispatch(results, get_power_system(sienna_sim.sienna_system))
    
    # Unit commitment analysis (if applicable)
    if formulation_type == "unit_commitment"
        analyze_unit_commitment_decisions(results, get_power_system(sienna_sim.sienna_system))
    end
    
    @info "="^60
end

"""
    analyze_generation_dispatch(results::OptimizationProblemResults, sys::System)

Analyze generation dispatch from optimization results.
"""
function analyze_generation_dispatch(results::OptimizationProblemResults, sys::System)
    @info "\nGeneration Dispatch Analysis:"
    
    try
        var_names = list_variable_names(results)
        total_generation = Dict{String, Float64}()
        
        # Thermal generation
        if "ActivePowerVariable__ThermalStandard" in var_names
            thermal_data = read_variable(results, "ActivePowerVariable__ThermalStandard")
            if ncol(thermal_data) > 1
                total_generation["Thermal"] = sum(sum(eachcol(thermal_data[:, 2:end])))
                max_thermal = maximum(sum(eachrow(thermal_data[:, 2:end])))
                @info "  Thermal: $(round(total_generation["Thermal"], digits=1)) MWh total, $(round(max_thermal, digits=1)) MW peak"
            end
        end
        
        # Renewable generation
        if "ActivePowerVariable__RenewableDispatch" in var_names
            renewable_data = read_variable(results, "ActivePowerVariable__RenewableDispatch")
            if ncol(renewable_data) > 1
                total_generation["Renewable"] = sum(sum(eachcol(renewable_data[:, 2:end])))
                max_renewable = maximum(sum(eachrow(renewable_data[:, 2:end])))
                @info "  Renewable: $(round(total_generation["Renewable"], digits=1)) MWh total, $(round(max_renewable, digits=1)) MW peak"
            end
        end
        
        # Hydro generation
        if "ActivePowerVariable__HydroDispatch" in var_names
            hydro_data = read_variable(results, "ActivePowerVariable__HydroDispatch")
            if ncol(hydro_data) > 1
                total_generation["Hydro"] = sum(sum(eachcol(hydro_data[:, 2:end])))
                max_hydro = maximum(sum(eachrow(hydro_data[:, 2:end])))
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
    analyze_unit_commitment_decisions(results::OptimizationProblemResults, sys::System)

Analyze unit commitment decisions.
"""
function analyze_unit_commitment_decisions(results::OptimizationProblemResults, sys::System)
    @info "\nUnit Commitment Analysis:"
    
    try
        var_names = list_variable_names(results)
        
        if "OnVariable__ThermalStandard" in var_names
            on_status = read_variable(results, "OnVariable__ThermalStandard")
            thermal_gens = get_components(ThermalStandard, sys)
            
            if ncol(on_status) > 1
                total_units = length(thermal_gens)
                hourly_committed = [sum(on_status[h, 2:end]) for h in 1:nrow(on_status)]
                
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
            start_data = read_variable(results, "StartVariable__ThermalStandard")
            if ncol(start_data) > 1
                total_starts = sum(sum(eachcol(start_data[:, 2:end])))
                @info "  Total start-ups: $(Int(total_starts))"
            end
        end
        
    catch e
        @warn "Unit commitment analysis failed: $e"
    end
end

"""
    print_suite_summary(sienna_sim::SiennaSimulations, suite_results::Dict, suite_start::DateTime)

Print comprehensive suite summary.
"""
function print_suite_summary(sienna_sim::SiennaSimulations, suite_results::Dict, suite_start::DateTime)
    total_time = Dates.value(now() - suite_start) / 1000.0
    
    @info "\n" * "="^70
    @info "SIENNA SIMULATION SUITE SUMMARY"
    @info "="^70
    @info "Project: $(sienna_sim.config.project_name)"
    @info "System: $(get_name(get_power_system(sienna_sim.sienna_system)))"
    @info "Total Runtime: $(round(total_time, digits=1)) seconds"
    @info "Framework: Config-driven Sienna ecosystem"
    
    @info ""
    @info "Simulation Results:"
    successful = 0
    for (sim_type, result) in suite_results
        if result["status"] == "SUCCESS"
            obj_val = result["objective_value"]
            solve_time = result["solve_time"]
            build_time = result["build_time"]
            @info "  ✅ $(titlecase(replace(sim_type, "_" => " "))): \$(round(obj_val, digits=2))"
            @info "     Build: $(round(build_time, digits=2))s, Solve: $(round(solve_time, digits=2))s"
            successful += 1
        else
            @info "  ❌ $(titlecase(replace(sim_type, "_" => " "))): FAILED"
            @info "     Error: $(result["error"])"
        end
    end
    
    @info ""
    @info "Performance Summary:"
    @info "  Success Rate: $successful/$(length(suite_results)) simulations"
    @info "  Total Build Time: $(round(sum(values(sienna_sim.build_times)), digits=2)) seconds"
    @info "  Total Solve Time: $(round(sum(values(sienna_sim.solve_times)), digits=2)) seconds"
    @info "  Network Model: $(sienna_sim.config.network_model)"
    @info "  Time Series: $(timeseries_loaded(sienna_sim.sienna_system) ? "Loaded" : "Not loaded")"
    @info "="^70
end

# ===== PUBLIC INTERFACE METHODS =====

"""
    get_results(sienna_sim::SiennaSimulations, formulation_type::String)

Get optimization results for specified formulation type.
"""
function get_results(sienna_sim::SiennaSimulations, formulation_type::String)
    if !haskey(sienna_sim.results, formulation_type)
        error("❌ No results available for $formulation_type. Run simulation first.")
    end
    return sienna_sim.results[formulation_type]
end

"""
    has_results(sienna_sim::SiennaSimulations, formulation_type::String)

Check if results are available for specified formulation type.
"""
function has_results(sienna_sim::SiennaSimulations, formulation_type::String)
    return haskey(sienna_sim.results, formulation_type)
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
    sienna_sim.results = Dict{String, OptimizationProblemResults}()
    sienna_sim.build_times = Dict{String, Float64}()
    sienna_sim.solve_times = Dict{String, Float64}()
    
    # Reset status flags
    sienna_sim.templates_created = false
    for formulation in ["economic_dispatch", "unit_commitment"]
        sienna_sim.problems_built[formulation] = false
        sienna_sim.simulations_completed[formulation] = false
        sienna_sim.simulation_errors[formulation] = String[]
    end
    
    sienna_sim.last_error = nothing
    
    @info "✓ Simulation state reset - ready for new simulations"
end

# ===== EXPORTS =====

export SiennaSimulations
export create_templates!, build_problem!, solve_problem!
export run_simulation!, run_simulation_suite!
export get_results, has_results, get_performance_summary
export reset_simulations!

# Test functionality when run directly
if abspath(PROGRAM_FILE) == @__FILE__
    @info "🧪 Testing SiennaSimulations..."
    
    try
        # Test with config and system
        if isfile("config.toml")
            config = SiennaConfig("config.toml")
            sienna_sys = SiennaSystem(config)
            
            if isdir(config.data_directory)
                @info "Building system for simulation test..."
                build_system!(sienna_sys)
                
                @info "Creating SiennaSimulations..."
                sienna_sim = SiennaSimulations(config, sienna_sys)
                
                @info "Testing template creation..."
                create_templates!(sienna_sim)
                
                @info "✅ SiennaSimulations test passed"
                @info "   Templates created: $(sienna_sim.templates_created)"
                @info "   Available templates: $(join(keys(sienna_sim.templates), ", "))"
                
            else
                @info "Data directory not found - skipping full simulation test"
                @info "✅ SiennaSimulations initialization test passed"
            end
        else
            @info "No config.toml found - cannot test SiennaSimulations"
            @info "Run SiennaConfig.jl first to create configuration"
        end
    catch e
        @error "❌ SiennaSimulations test failed: $e"
    end
end
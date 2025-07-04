#!/usr/bin/env julia

"""
Sienna Workflow V1.0 - Modern Config-Driven Power System Simulation
==================================================================

Clean, modern workflow orchestrating proven components:
- Config-driven everything (config.toml as single source of truth)
- PowerSystems.jl 4.6.2 + PowerSimulations.jl 0.30.2
- Default templates first, custom by design choice
- Crystal clear error handling and validation
- Foundation-grade code for production use

Architecture:
    config.toml → validate_configuration() → build_sienna_system() → 
    powersimulations_framework() → save_results() → workflow_summary()

Usage:
    julia sienna_workflow_v1.jl                    # Full config-driven workflow
    
    # Or in Julia REPL:
    julia> include("sienna_workflow_v1.jl")
    julia> run_sienna_workflow()                   # Full workflow
    julia> run_sienna_workflow_test()              # Quick test
"""

# Include proven components (unchanged)
include("config_loader.jl")
include("build_sienna_system.jl") 
include("powersimulations_framework.jl")

using PowerSystems
using PowerSimulations
using Dates
using Logging
using JSON3

# Set up logging
global_logger(ConsoleLogger(stderr, Logging.Info))

"""
    run_sienna_workflow()

Main entry point: Complete config-driven Sienna workflow.
This is the primary function users should call.
"""
function run_sienna_workflow()
    @info "🚀 SIENNA WORKFLOW V1.0 - CONFIG-DRIVEN SIMULATION"
    @info "="^70
    
    workflow_start_time = now()
    
    try
        # Step 1: Validate configuration completely upfront
        @info "\n📋 STEP 1: Configuration Validation"
        @info "─"^50
        config = validate_and_load_configuration()
        
        # Step 2: Build system using proven builder
        @info "\n🏗️  STEP 2: System Building"
        @info "─"^50
        sys = build_system_from_config(config)
        
        # Step 3: Run simulations using proven framework
        @info "\n⚡ STEP 3: Power System Simulations"
        @info "─"^50
        simulation_results = run_simulations_from_config(sys, config)
        
        # Step 4: Save and export results
        @info "\n💾 STEP 4: Results Export"
        @info "─"^50
        export_info = export_workflow_results(sys, simulation_results, config)
        
        # Step 5: Workflow summary
        @info "\n📊 STEP 5: Workflow Summary"
        @info "─"^50
        print_workflow_summary(sys, simulation_results, export_info, workflow_start_time)
        
        @info "✅ SIENNA WORKFLOW COMPLETED SUCCESSFULLY!"
        
        return WorkflowResults(sys, simulation_results, export_info, config)
        
    catch e
        @error "❌ WORKFLOW FAILED: $e"
        @error "Stack trace:" exception=(e, catch_backtrace())
        @error "\n🔧 TROUBLESHOOTING TIPS:"
        @error "  1. Check config.toml exists and has required sections"
        @error "  2. Verify data_directory path is correct and contains CSV files"  
        @error "  3. Ensure system has required components (buses, generators, loads)"
        @error "  4. Check solver availability (HiGHS should be installed)"
        return nothing
    end
end

"""
    run_sienna_workflow_test()

Quick test workflow with reduced scope for validation.
Uses config but with minimal horizon for speed.
"""
function run_sienna_workflow_test()
    @info "🧪 SIENNA WORKFLOW TEST MODE"
    @info "="^50
    
    try
        # Load and validate config
        config = validate_and_load_configuration()
        
        # Modify config for quick test
        test_config = deepcopy(config)
        test_config["simulations"]["default_horizon_hours"] = 3  # Quick test
        test_config["simulations"]["run_unit_commitment"] = false  # ED only for speed
        test_config["output"]["create_timestamped_folders"] = false  # Simpler output
        
        @info "Test configuration:"
        @info "  Project: $(test_config["project"]["name"])"
        @info "  Horizon: 3 hours (ED only)"
        @info "  Data: $(test_config["paths"]["data_directory"])"
        
        # Build system
        sys = build_system_from_config(test_config)
        
        # Run quick simulation
        simulation_results = run_simulations_from_config(sys, test_config)
        
        # Check results
        success = validate_test_results(simulation_results)
        
        if success
            @info "✅ TEST PASSED: Workflow is functional"
            return true
        else
            @warn "⚠️  TEST INCOMPLETE: Some simulations failed"
            return false
        end
        
    catch e
        @error "❌ TEST FAILED: $e"
        return false
    end
end

"""
    validate_and_load_configuration()

Comprehensive upfront configuration validation.
Returns validated config or throws descriptive errors.
"""
function validate_and_load_configuration()
    @info "Validating Sienna configuration..."
    
    # Load configuration
    config = try
        load_config()
    catch e
        error("❌ Configuration loading failed: $e\n" *
              "🔧 Solution: Ensure config.toml exists and is valid TOML format")
    end
    
    # Track all validation issues
    validation_errors = String[]
    validation_warnings = String[]
    
    # Validate required sections
    required_sections = ["project", "paths", "system_building", "simulations", "output"]
    for section in required_sections
        if !haskey(config, section)
            push!(validation_errors, "Missing required config section: [$section]")
        end
    end
    
    # Validate critical paths
    if haskey(config, "paths")
        data_dir = get(config["paths"], "data_directory", "")
        if isempty(data_dir)
            push!(validation_errors, "Missing required config: paths.data_directory")
        elseif !isdir(expanduser(data_dir))
            push!(validation_errors, "Data directory does not exist: $data_dir")
        else
            # Check for required CSV files
            required_files = ["bus.csv", "load.csv"]
            for file in required_files
                file_path = joinpath(expanduser(data_dir), file)
                if !isfile(file_path)
                    push!(validation_errors, "Required data file missing: $file in $data_dir")
                end
            end
            
            # Check for generator files
            gen_file = joinpath(expanduser(data_dir), "gen.csv")
            if !isfile(gen_file)
                push!(validation_warnings, "Generator file missing: gen.csv (system may have no generators)")
            end
        end
        
        # Validate output directory is writable
        output_dir = get(config["paths"], "output_directory", "./sienna_results")
        try
            mkpath(expanduser(output_dir))
        catch e
            push!(validation_errors, "Cannot create output directory: $output_dir ($e)")
        end
    else
        push!(validation_errors, "Missing required config section: [paths]")
    end
    
    # Validate simulation settings
    if haskey(config, "simulations")
        solver = get(config["simulations"], "default_solver", "")
        if isempty(solver)
            push!(validation_warnings, "No default_solver specified, will use HiGHS")
            config["simulations"]["default_solver"] = "HiGHS"
        elseif solver != "HiGHS"
            push!(validation_warnings, "Only HiGHS solver is currently supported, switching to HiGHS")
            config["simulations"]["default_solver"] = "HiGHS"
        end
        
        horizon = get(config["simulations"], "default_horizon_hours", 0)
        if horizon <= 0
            push!(validation_warnings, "Invalid horizon_hours, using default: 24")
            config["simulations"]["default_horizon_hours"] = 24
        end
    end
    
    # Validate system building parameters  
    if haskey(config, "system_building")
        base_power = get(config["system_building"], "base_power", 0.0)
        if base_power <= 0
            push!(validation_warnings, "Invalid base_power, using default: 100.0 MW")
            config["system_building"]["base_power"] = 100.0
        end
    end
    
    # Report validation results
    if !isempty(validation_errors)
        @error "❌ CONFIGURATION VALIDATION FAILED:"
        for error in validation_errors
            @error "  • $error"
        end
        @error "\n🔧 Fix these configuration issues and try again"
        error("Configuration validation failed with $(length(validation_errors)) errors")
    end
    
    if !isempty(validation_warnings)
        @warn "⚠️  Configuration warnings (using defaults):"
        for warning in validation_warnings
            @warn "  • $warning"
        end
    end
    
    @info "✓ Configuration validation passed"
    @info "  Project: $(config["project"]["name"])"
    @info "  Data directory: $(config["paths"]["data_directory"])"
    @info "  Solver: $(config["simulations"]["default_solver"])"
    @info "  Horizon: $(config["simulations"]["default_horizon_hours"]) hours"
    
    return config
end

"""
    build_system_from_config(config::Dict)

Build PowerSystems.jl system using proven build_sienna_system.jl with config parameters.
"""
function build_system_from_config(config::Dict)
    @info "Building PowerSystems.jl system from configuration..."
    
    # Extract parameters from config (no hardcoded fallbacks - config should be validated)
    data_dir = config["paths"]["data_directory"]
    base_power = config["system_building"]["base_power"]
    validate_system = config["system_building"]["validate_system"]
    load_timeseries = config["system_building"]["load_timeseries"]
    
    @info "System building parameters:"
    @info "  Data directory: $data_dir"
    @info "  Base power: $base_power MW"
    @info "  Validation: $validate_system"
    @info "  Time series: $load_timeseries"
    
    # Use proven system builder (unchanged)
    sys = build_sienna_system(
        data_dir,
        validate_system=validate_system,
        load_timeseries=load_timeseries,
        base_power=base_power
    )
    
    if sys === nothing
        error("❌ System building failed\n" *
              "🔧 Check data files in $data_dir and build_sienna_system.jl output above")
    end
    
    @info "✓ System built successfully: $(get_name(sys))"
    return sys
end

"""
    run_simulations_from_config(sys::System, config::Dict)

Run power system simulations using proven powersimulations_framework.jl with config parameters.

"""
function run_simulations_from_config(sys::System, config::Dict)
    @info "Running power system simulations from configuration..."
    
    # Extract simulation parameters from config
    horizon_hours = config["simulations"]["default_horizon_hours"]
    run_ed = config["simulations"]["run_economic_dispatch"]
    run_uc = config["simulations"]["run_unit_commitment"]
    solver_type = config["simulations"]["default_solver"]
    use_defaults = get(config["simulations"], "use_default_templates", true)
    
    # Create output directory
    output_base = config["paths"]["output_directory"]
    if config["output"]["create_timestamped_folders"]
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        output_dir = joinpath(expanduser(output_base), "simulation_$(timestamp)")
    else
        output_dir = expanduser(output_base)
    end
    mkpath(output_dir)
    
    @info "Simulation parameters:"
    @info "  Horizon: $horizon_hours hours"
    @info "  Economic Dispatch: $run_ed"
    @info "  Unit Commitment: $run_uc"
    @info "  Solver: $solver_type"
    @info "  Templates: $(use_defaults ? "Default" : "Custom")"
    @info "  Output: $output_dir"
    
    simulation_results = Dict()
    
    # Economic Dispatch
    if run_ed
        @info "\n>>> Running Economic Dispatch"
        try
            # Get solver settings from config
            solver_settings = get_solver_settings_from_config(config, "economic_dispatch")
            
            # Run using proven framework
            ed_results, ed_problem = run_operations_problem(
                sys, 
                "economic_dispatch",
                horizon_hours=horizon_hours,
                output_dir=joinpath(output_dir, "economic_dispatch"),
                solver_type=solver_type,
                use_default_template=use_defaults,
                time_limit=solver_settings["time_limit"],
                mip_gap=solver_settings["mip_gap"]
            )
            
            simulation_results["economic_dispatch"] = Dict(
                "status" => "SUCCESS",
                "objective_value" => get_objective_value(ed_results),
                "solve_time" => get_optimizer_stats(ed_results).solve_time[1],
                "output_directory" => joinpath(output_dir, "economic_dispatch")
            )
            
            @info "  ✓ Economic Dispatch completed: \$$(round(get_objective_value(ed_results), digits=2))"
            
        catch e
            @error "  ❌ Economic Dispatch failed: $e"
            simulation_results["economic_dispatch"] = Dict(
                "status" => "FAILED", 
                "error" => string(e)
            )
        end
    else
        @info "Economic Dispatch skipped (disabled in config)"
    end
    
    # Unit Commitment
    if run_uc
        @info "\n>>> Running Unit Commitment"
        try
            # Get solver settings from config
            solver_settings = get_solver_settings_from_config(config, "unit_commitment")
            
            # Run using proven framework
            uc_results, uc_problem = run_operations_problem(
                sys, 
                "unit_commitment",
                horizon_hours=horizon_hours,
                output_dir=joinpath(output_dir, "unit_commitment"),
                solver_type=solver_type,
                use_default_template=use_defaults,
                time_limit=solver_settings["time_limit"],
                mip_gap=solver_settings["mip_gap"]
            )
            
            simulation_results["unit_commitment"] = Dict(
                "status" => "SUCCESS",
                "objective_value" => get_objective_value(uc_results),
                "solve_time" => get_optimizer_stats(uc_results).solve_time[1],
                "output_directory" => joinpath(output_dir, "unit_commitment")
            )
            
            @info "  ✓ Unit Commitment completed: \$$(round(get_objective_value(uc_results), digits=2))"
            
        catch e
            @error "  ❌ Unit Commitment failed: $e"
            simulation_results["unit_commitment"] = Dict(
                "status" => "FAILED", 
                "error" => string(e)
            )
        end
    else
        @info "Unit Commitment skipped (disabled in config)"
    end
    
    # Validate we ran at least one simulation
    successful_sims = [k for (k, v) in simulation_results if v["status"] == "SUCCESS"]
    if isempty(successful_sims)
        error("❌ All simulations failed or none were enabled\n" *
              "🔧 Check simulation settings in config and error messages above")
    end
    
    @info "✓ Simulations completed: $(length(successful_sims)) successful"
    return simulation_results
end

"""
    get_solver_settings_from_config(config::Dict, formulation_type::String)

Extract solver settings from config for specific formulation type.
"""
function get_solver_settings_from_config(config::Dict, formulation_type::String)
    solver_name = config["simulations"]["default_solver"]
    
    if haskey(config, "simulations") && 
       haskey(config["simulations"], "solver_settings") &&
       haskey(config["simulations"]["solver_settings"], solver_name)
        
        solver_config = config["simulations"]["solver_settings"][solver_name]
        
        if formulation_type == "economic_dispatch"
            return Dict(
                "time_limit" => get(solver_config, "time_limit_ed", 300.0),
                "mip_gap" => get(solver_config, "mip_gap_ed", 0.01),
                "threads" => get(solver_config, "threads", 0)
            )
        elseif formulation_type == "unit_commitment"
            return Dict(
                "time_limit" => get(solver_config, "time_limit_uc", 600.0),
                "mip_gap" => get(solver_config, "mip_gap_uc", 0.02),
                "threads" => get(solver_config, "threads", 0)
            )
        end
    end
    
    # Fallback defaults with warning
    @warn "No solver settings found in config for $formulation_type, using defaults"
    return Dict(
        "time_limit" => formulation_type == "unit_commitment" ? 600.0 : 300.0,
        "mip_gap" => formulation_type == "unit_commitment" ? 0.02 : 0.01,
        "threads" => 0
    )
end

"""
    export_workflow_results(sys::System, simulation_results::Dict, config::Dict)

Export and organize workflow results based on config settings.
"""
function export_workflow_results(sys::System, simulation_results::Dict, config::Dict)
    @info "Exporting workflow results..."
    
    export_info = Dict()
    
    # Export PowerSystemCaseBuilder case if requested
    if get(get(config, "powersystemcasebuilder", Dict()), "create_psb_case", false)
        @info "Creating PowerSystemCaseBuilder case..."
        
        export_dir = config["paths"]["export_directory"]
        case_name = "$(replace(config["project"]["name"], " " => "_"))_$(Dates.format(now(), "yyyy_mm_dd"))"
        
        # Use proven export function from build_sienna_system.jl
        psb_case_dir, psb_metadata = create_psb_case_from_system(
            sys, 
            case_name, 
            expanduser(export_dir),
            description="Generated by Sienna Workflow V1.0"
        )
        
        export_info["psb_case"] = Dict(
            "directory" => psb_case_dir,
            "metadata" => psb_metadata,
            "case_name" => case_name
        )
        
        @info "  ✓ PowerSystemCaseBuilder case: $psb_case_dir"
    end
    
    # Save workflow summary
    summary_dir = config["paths"]["output_directory"]
    if config["output"]["create_timestamped_folders"]
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        summary_dir = joinpath(expanduser(summary_dir), "workflow_summary_$(timestamp)")
    else
        summary_dir = expanduser(summary_dir)
    end
    mkpath(summary_dir)
    
    # Create comprehensive workflow summary
    workflow_summary = Dict(
        "workflow_version" => "Sienna Workflow V1.0",
        "timestamp" => string(now()),
        "config" => Dict(
            "project_name" => config["project"]["name"],
            "data_directory" => config["paths"]["data_directory"],
            "solver" => config["simulations"]["default_solver"],
            "horizon_hours" => config["simulations"]["default_horizon_hours"]
        ),
        "system" => Dict(
            "name" => get_name(sys),
            "base_power" => get_base_power(sys),
            "component_counts" => Dict(
                "buses" => length(get_components(ACBus, sys)),
                "thermal_generators" => length(get_components(ThermalStandard, sys)),
                "renewable_generators" => length(get_components(RenewableDispatch, sys)),
                "loads" => length(get_components(PowerLoad, sys)),
                "lines" => length(get_components(Line, sys))
            )
        ),
        "simulation_results" => simulation_results,
        "export_info" => export_info
    )
    
    # Save summary as JSON
    summary_file = joinpath(summary_dir, "sienna_workflow_summary.json")
    open(summary_file, "w") do f
        JSON3.pretty(f, workflow_summary)
    end
    
    export_info["workflow_summary"] = summary_file
    @info "  ✓ Workflow summary: $summary_file"
    
    return export_info
end

"""
    print_workflow_summary(sys::System, simulation_results::Dict, export_info::Dict, start_time::DateTime)

Print comprehensive workflow summary.
"""
function print_workflow_summary(sys::System, simulation_results::Dict, export_info::Dict, start_time::DateTime)
    total_time = Dates.value(now() - start_time) / 1000.0  # seconds
    
    @info "\n" * "="^70
    @info "SIENNA WORKFLOW V1.0 - FINAL SUMMARY"
    @info "="^70
    @info "Project: $(get_name(sys))"
    @info "Total Runtime: $(round(total_time, digits=1)) seconds"
    @info ""
    @info "System Composition:"
    @info "  • Buses: $(length(get_components(ACBus, sys)))"
    @info "  • Thermal Generators: $(length(get_components(ThermalStandard, sys)))"
    @info "  • Renewable Generators: $(length(get_components(RenewableDispatch, sys)))"
    @info "  • Loads: $(length(get_components(PowerLoad, sys)))"
    @info "  • AC Lines: $(length(get_components(Line, sys)))"
    
    @info ""
    @info "Simulation Results:"
    successful_count = 0
    total_cost = 0.0
    
    for (sim_type, result) in simulation_results
        if result["status"] == "SUCCESS"
            obj_val = result["objective_value"]
            solve_time = result["solve_time"]
            @info "  ✓ $(titlecase(replace(sim_type, "_" => " "))): \$$(round(obj_val, digits=2)) ($(round(solve_time, digits=1))s)"
            successful_count += 1
            total_cost += obj_val
        else
            @info "  ❌ $(titlecase(replace(sim_type, "_" => " "))): FAILED"
        end
    end
    
    if successful_count > 0
        @info "  Total System Cost: \$$(round(total_cost, digits=2))"
    end
    
    @info ""
    @info "Outputs Generated:"
    for (sim_type, result) in simulation_results
        if result["status"] == "SUCCESS" && haskey(result, "output_directory")
            @info "  • $sim_type: $(result["output_directory"])"
        end
    end
    
    if haskey(export_info, "psb_case")
        @info "  • PowerSystemCaseBuilder Case: $(export_info["psb_case"]["directory"])"
    end
    
    if haskey(export_info, "workflow_summary")
        @info "  • Workflow Summary: $(export_info["workflow_summary"])"
    end
    
    @info ""
    @info "Success Rate: $successful_count/$(length(simulation_results)) simulations completed"
    @info "="^70
end

"""
    validate_test_results(simulation_results::Dict)

Validate that test results meet basic sanity checks.
"""
function validate_test_results(simulation_results::Dict)
    success_count = 0
    total_count = length(simulation_results)
    
    for (sim_type, result) in simulation_results
        if result["status"] == "SUCCESS"
            success_count += 1
            
            # Basic sanity checks
            obj_val = result["objective_value"]
            solve_time = result["solve_time"]
            
            if obj_val <= 0
                @warn "Suspicious objective value for $sim_type: $obj_val"
            end
            
            if solve_time > 60
                @warn "Long solve time for test $sim_type: $(round(solve_time, digits=1))s"
            end
        end
    end
    
    return success_count > 0  # At least one simulation must succeed
end

"""
    WorkflowResults

Container for complete workflow results.
"""
struct WorkflowResults
    system::System
    simulation_results::Dict
    export_info::Dict
    config::Dict
end

"""
    show_workflow_help()

Display help information for the Sienna workflow.
"""
function show_workflow_help()
    @info """
    
🚀 SIENNA WORKFLOW V1.0 - HELP
==============================

Main Functions:
  run_sienna_workflow()        # Complete config-driven workflow
  run_sienna_workflow_test()   # Quick test (3 hours, ED only)
  validate_and_load_configuration()  # Check config without running

Configuration:
  • Edit config.toml for all settings
  • Required: [project], [paths], [system_building], [simulations], [output]
  • Data files: bus.csv, load.csv, gen.csv (in data_directory)

Quick Start:
  1. Ensure config.toml exists and data_directory is correct
  2. Run: julia -e "include(\\"sienna_workflow_v1.jl\\"); run_sienna_workflow()"
  3. Check results in output_directory

Troubleshooting:
  • Config errors: Check TOML syntax and required fields
  • Data errors: Verify CSV files exist and have proper format
  • Simulation errors: Check solver installation and system validity
  
For detailed documentation, see function docstrings.
    """
end

# Export main workflow functions
export run_sienna_workflow, run_sienna_workflow_test, validate_and_load_configuration
export show_workflow_help, WorkflowResults

# Main execution when run directly
function main()
    @info "🚀 SIENNA WORKFLOW V1.0"
    @info "="^40
    
    if length(ARGS) > 0
        if ARGS[1] == "test"
            @info "Running in test mode..."
            run_sienna_workflow_test()
        elseif ARGS[1] == "help"
            show_workflow_help()
        elseif ARGS[1] == "validate"
            @info "Validating configuration only..."
            validate_and_load_configuration()
            @info "✓ Configuration is valid"
        else
            @info "Unknown argument: $(ARGS[1])"
            show_workflow_help()
        end
    else
        @info "Running full workflow..."
        run_sienna_workflow()
    end
end

# Run main if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
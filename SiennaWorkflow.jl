#!/usr/bin/env julia

"""
SiennaWorkflow.jl - Orchestrated Workflow for Sienna Ecosystem
===============================================================

Simplified workflow orchestrator using the new class-based architecture:
- SiennaConfig.jl for configuration management
- SiennaSystem.jl for PowerSystems.jl system building
- SiennaSimulations.jl for PowerSimulations.jl operations

This is now a lightweight orchestrator that coordinates the three main classes.

Features:
- Clean separation of concerns
- Class-based architecture
- Comprehensive error handling
- Performance monitoring
- Flexible workflow options

Usage:
    julia SiennaWorkflow.jl                # Full workflow
    julia SiennaWorkflow.jl test           # Quick test
    julia SiennaWorkflow.jl validate       # Validation only
"""

using Dates
using Logging
using JSON3

# Import our class-based modules
include("SiennaConfig.jl")
include("SiennaSystem.jl")
include("SiennaSimulations.jl")

global_logger(ConsoleLogger(stderr, Logging.Info))

"""
    SiennaWorkflowResults

Container for complete workflow results.
"""
struct SiennaWorkflowResults
    config::SiennaConfig
    sienna_system::SiennaSystem
    sienna_simulations::SiennaSimulations
    workflow_summary::Dict{String, Any}
    success::Bool
    total_runtime::Float64
    timestamp::DateTime
end

"""
    run_sienna_workflow(config_file::String="config.toml"; test_mode::Bool=false)

Main workflow orchestrator using class-based architecture.
"""
function run_sienna_workflow(config_file::String="config.toml"; test_mode::Bool=false)
    @info "🚀 SIENNA WORKFLOW - CLASS-BASED ARCHITECTURE"
    @info "="^70
    
    workflow_start = now()
    workflow_summary = Dict{String, Any}()
    
    try
        # === STEP 1: Configuration Management ===
        @info "\n📋 STEP 1: Configuration Management"
        @info "─"^50
        
        step_start = now()
        config = SiennaConfig(config_file)
        show_summary(config)
        workflow_summary["config_load_time"] = Dates.value(now() - step_start) / 1000.0
        
        # === STEP 2: System Building ===
        @info "\n🏗️  STEP 2: System Building"
        @info "─"^50
        
        step_start = now()
        sienna_system = SiennaSystem(config)
        system = build_system!(sienna_system)
        show_summary(sienna_system)
        workflow_summary["system_build_time"] = Dates.value(now() - step_start) / 1000.0
        
        # === STEP 3: Compatibility Validation ===
        @info "\n🔍 STEP 3: System-Configuration Compatibility"
        @info "─"^50
        
        step_start = now()
        compatibility = analyze_compatibility(sienna_system)
        
        if !compatibility["is_compatible"]
            @warn "⚠️  Compatibility issues found:"
            for issue in compatibility["compatibility_issues"]
                @warn "  • $issue"
            end
        else
            @info "✅ System and configuration are fully compatible"
        end
        
        if !isempty(compatibility["suggestions"])
            @info "💡 Optimization suggestions:"
            for suggestion in compatibility["suggestions"]
                @info "  • $suggestion"
            end
        end
        workflow_summary["compatibility_check_time"] = Dates.value(now() - step_start) / 1000.0
        
        # === STEP 4: Simulation Setup ===
        @info "\n⚡ STEP 4: Simulation Framework Setup"
        @info "─"^50
        
        step_start = now()
        sienna_simulations = SiennaSimulations(config, sienna_system)
        
        # Create templates for all enabled formulations
        create_templates!(sienna_simulations)
        workflow_summary["simulation_setup_time"] = Dates.value(now() - step_start) / 1000.0
        
        # === STEP 5: Run Simulations ===
        @info "\n🎯 STEP 5: Power System Simulations"
        @info "─"^50
        
        step_start = now()
        
        if test_mode
            # Test mode - run only Economic Dispatch with reduced horizon
            @info "🧪 Test Mode: Running Economic Dispatch only..."
            
            # Temporarily modify config for testing
            original_horizon = config.default_horizon_hours
            config.default_horizon_hours = min(3, original_horizon)
            
            if should_run_formulation(config, "economic_dispatch")
                ed_results = run_simulation!(sienna_simulations, "economic_dispatch")
                simulation_results = Dict("economic_dispatch" => Dict(
                    "status" => "SUCCESS",
                    "objective_value" => get_objective_value(ed_results)
                ))
            else
                @warn "Economic Dispatch disabled in config - cannot run test"
                simulation_results = Dict()
            end
            
            # Restore original horizon
            config.default_horizon_hours = original_horizon
        else
            # Full simulation suite
            simulation_results = run_simulation_suite!(sienna_simulations)
        end
        
        workflow_summary["simulation_time"] = Dates.value(now() - step_start) / 1000.0
        
        # === STEP 6: Results Summary ===
        @info "\n📊 STEP 6: Workflow Summary"
        @info "─"^50
        
        # Calculate total runtime
        total_runtime = Dates.value(now() - workflow_start) / 1000.0
        workflow_summary["total_runtime"] = total_runtime
        workflow_summary["simulation_results"] = simulation_results
        workflow_summary["performance"] = sienna_simulations.get_performance_summary()
        
        # Create comprehensive workflow results
        workflow_results = SiennaWorkflowResults(
            config,
            sienna_system,
            sienna_simulations,
            workflow_summary,
            true,  # success
            total_runtime,
            now()
        )
        
        # Print final summary
        print_workflow_summary(workflow_results, test_mode)
        
        # Save workflow summary
        save_workflow_summary(workflow_results)
        
        @info "✅ SIENNA WORKFLOW COMPLETED SUCCESSFULLY!"
        @info "Total Runtime: $(round(total_runtime, digits=1)) seconds"
        
        return workflow_results
        
    catch e
        # Handle workflow failure
        total_runtime = Dates.value(now() - workflow_start) / 1000.0
        workflow_summary["total_runtime"] = total_runtime
        workflow_summary["error"] = string(e)
        
        @error "❌ SIENNA WORKFLOW FAILED: $e"
        @error "Total Runtime: $(round(total_runtime, digits=1)) seconds"
        
        # Create failed workflow results
        failed_results = SiennaWorkflowResults(
            nothing,  # config
            nothing,  # sienna_system  
            nothing,  # sienna_simulations
            workflow_summary,
            false,    # success = false
            total_runtime,
            now()
        )
        
        return failed_results
    end
end

"""
    run_workflow_validation_only(config_file::String="config.toml")

Run only configuration and system validation without simulations.
"""
function run_workflow_validation_only(config_file::String="config.toml")
    @info "🔍 SIENNA WORKFLOW - VALIDATION ONLY"
    @info "="^60
    
    validation_start = now()
    
    try
        # Load and validate configuration
        @info "\n📋 Configuration Validation"
        @info "─"^40
        config = SiennaConfig(config_file)
        @info "✅ Configuration loaded and validated"
        
        # Build and validate system
        @info "\n🏗️  System Validation"
        @info "─"^40
        sienna_system = SiennaSystem(config)
        build_system!(sienna_system)
        @info "✅ System built and validated"
        
        # Compatibility analysis
        @info "\n🔍 Compatibility Analysis"
        @info "─"^40
        compatibility = analyze_compatibility(sienna_system)
        
        if compatibility["is_compatible"]
            @info "✅ System and configuration are fully compatible"
        else
            @warn "⚠️  Compatibility issues found:"
            for issue in compatibility["compatibility_issues"]
                @warn "  • $issue"
            end
        end
        
        if !isempty(compatibility["suggestions"])
            @info "💡 Optimization suggestions:"
            for suggestion in compatibility["suggestions"]
                @info "  • $suggestion"
            end
        end
        
        validation_time = Dates.value(now() - validation_start) / 1000.0
        @info "\n✅ VALIDATION COMPLETED SUCCESSFULLY"
        @info "Validation Time: $(round(validation_time, digits=1)) seconds"
        
        return true
        
    catch e
        validation_time = Dates.value(now() - validation_start) / 1000.0
        @error "❌ VALIDATION FAILED: $e"
        @error "Validation Time: $(round(validation_time, digits=1)) seconds"
        return false
    end
end

"""
    run_workflow_test()

Run quick test workflow with minimal configuration.
"""
function run_workflow_test()
    @info "🧪 SIENNA WORKFLOW - QUICK TEST MODE"
    @info "="^60
    
    return run_sienna_workflow("config.toml", test_mode=true)
end

"""
    print_workflow_summary(results::SiennaWorkflowResults, test_mode::Bool=false)

Print comprehensive workflow summary.
"""
function print_workflow_summary(results::SiennaWorkflowResults, test_mode::Bool=false)
    @info "\n" * "="^80
    @info "SIENNA WORKFLOW SUMMARY - CLASS-BASED ARCHITECTURE"
    @info "="^80
    
    # Project information
    @info "Project: $(results.config.project_name)"
    @info "System: $(get_name(get_power_system(results.sienna_system)))"
    @info "Completed: $(Dates.format(results.timestamp, "yyyy-mm-dd HH:MM:SS"))"
    @info "Mode: $(test_mode ? "Test Mode (Reduced Scope)" : "Full Workflow")"
    
    # Architecture summary
    @info ""
    @info "Architecture Components:"
    @info "  ✓ SiennaConfig: Configuration management"
    @info "  ✓ SiennaSystem: PowerSystems.jl integration"
    @info "  ✓ SiennaSimulations: PowerSimulations.jl integration"
    @info "  ✓ SiennaWorkflow: Orchestration layer"
    
    # Performance breakdown
    summary = results.workflow_summary
    @info ""
    @info "Performance Breakdown:"
    @info "  Config Loading: $(round(get(summary, "config_load_time", 0.0), digits=2))s"
    @info "  System Building: $(round(get(summary, "system_build_time", 0.0), digits=2))s"
    @info "  Compatibility Check: $(round(get(summary, "compatibility_check_time", 0.0), digits=2))s"
    @info "  Simulation Setup: $(round(get(summary, "simulation_setup_time", 0.0), digits=2))s"
    @info "  Simulations: $(round(get(summary, "simulation_time", 0.0), digits=2))s"
    @info "  Total Runtime: $(round(results.total_runtime, digits=2))s"
    
    # System statistics
    if results.sienna_system !== nothing
        sys_stats = get_system_statistics(results.sienna_system)
        @info ""
        @info "System Statistics:"
        @info "  Base Power: $(sys_stats["base_power"]) MW"
        @info "  Total Generation: $(sys_stats["capacity"]["total_generation_mw"]) MW"
        @info "  Total Load: $(sys_stats["capacity"]["total_load_mw"]) MW"
        @info "  Reserve Margin: $(sys_stats["capacity"]["reserve_margin"])%"
        @info "  Network Model: $(results.config.network_model)"
        @info "  AC Lines: $(sys_stats["network"]["ac_lines"])"
        @info "  DC Lines: $(sys_stats["network"]["dc_lines"])"
    end
    
    # Simulation results
    if haskey(summary, "simulation_results") && !isempty(summary["simulation_results"])
        simulation_results = summary["simulation_results"]
        @info ""
        @info "Simulation Results:"
        
        successful_sims = 0
        total_sims = length(simulation_results)
        
        for (sim_type, result) in simulation_results
            if result["status"] == "SUCCESS"
                obj_val = result["objective_value"]
                @info "  ✅ $(titlecase(replace(sim_type, "_" => " "))): \$(round(obj_val, digits=2))"
                successful_sims += 1
            else
                @info "  ❌ $(titlecase(replace(sim_type, "_" => " "))): FAILED"
            end
        end
        
        @info "  Success Rate: $successful_sims/$total_sims"
    end
    
    # Configuration summary
    @info ""
    @info "Configuration Summary:"
    @info "  Data Directory: $(results.config.data_directory)"
    @info "  Output Directory: $(results.config.output_directory)"
    @info "  Horizon: $(results.config.default_horizon_hours) hours"
    @info "  Solver: $(results.config.default_solver)"
    @info "  Time Series: $(results.sienna_system !== nothing && timeseries_loaded(results.sienna_system) ? "Loaded" : "Not loaded")"
    
    @info "="^80
    
    if results.success
        @info "🎉 WORKFLOW ARCHITECTURE: Fully modular and class-based!"
        @info "🎉 PERFORMANCE: $(round(results.total_runtime, digits=1))s total runtime"
        @info "🎉 MAINTAINABILITY: Clean separation of concerns achieved"
    else
        @info "❌ WORKFLOW FAILED - Check error messages above"
    end
end

"""
    save_workflow_summary(results::SiennaWorkflowResults)

Save comprehensive workflow summary to JSON file.
"""
function save_workflow_summary(results::SiennaWorkflowResults)
    @info "💾 Saving workflow summary..."
    
    output_dir = results.config.output_directory
    summary_file = joinpath(output_dir, "sienna_workflow_summary.json")
    
    # Create comprehensive summary
    comprehensive_summary = Dict(
        "workflow_metadata" => Dict(
            "version" => "2.0",
            "architecture" => "class-based",
            "timestamp" => string(results.timestamp),
            "success" => results.success,
            "total_runtime_seconds" => results.total_runtime
        ),
        "components" => Dict(
            "config" => Dict(
                "project_name" => results.config.project_name,
                "network_model" => results.config.network_model,
                "horizon_hours" => results.config.default_horizon_hours,
                "solver" => results.config.default_solver,
                "data_directory" => results.config.data_directory
            ),
            "system" => results.sienna_system !== nothing ? results.sienna_system.get_system_statistics() : nothing,
            "simulations" => results.sienna_simulations !== nothing ? results.sienna_simulations.get_performance_summary() : nothing
        ),
        "workflow_performance" => results.workflow_summary,
        "class_based_benefits" => [
            "Clean separation of concerns",
            "Modular architecture",
            "Reusable components",
            "Better error handling",
            "Easier testing and maintenance",
            "Configuration-driven approach"
        ]
    )
    
    try
        open(summary_file, "w") do f
            JSON3.pretty(f, comprehensive_summary)
        end
        @info "✓ Workflow summary saved to: $summary_file"
    catch e
        @warn "Failed to save workflow summary: $e"
    end
end

"""
    show_help()

Display help information for the workflow.
"""
function show_help()
    @info """
    
🚀 SIENNA WORKFLOW - CLASS-BASED ARCHITECTURE HELP
=================================================

ARCHITECTURE OVERVIEW:
  ✅ SiennaConfig.jl    - Configuration management and validation
  ✅ SiennaSystem.jl    - PowerSystems.jl system building and analysis  
  ✅ SiennaSimulations.jl - PowerSimulations.jl templates and solving
  ✅ SiennaWorkflow.jl  - Lightweight orchestration layer

KEY IMPROVEMENTS:
  • Clean separation of concerns
  • Class-based, object-oriented design
  • Reusable, modular components
  • Configuration-driven throughout
  • Better error handling and validation
  • Easier testing and maintenance

USAGE OPTIONS:
  julia SiennaWorkflow.jl              # Full workflow
  julia SiennaWorkflow.jl test         # Quick test (ED only, 3 hours)
  julia SiennaWorkflow.jl validate     # Validation only (no simulations)
  julia SiennaWorkflow.jl help         # This help message

WORKFLOW STEPS:
  1. Configuration Management (SiennaConfig)
     - Load and validate config.toml
     - Check paths and parameters
     - Set up fuel mappings and defaults
  
  2. System Building (SiennaSystem)
     - Build PowerSystems.jl system from CSV data
     - Load time series with metadata
     - Validate system components
  
  3. Compatibility Analysis
     - Check config vs system compatibility
     - Generate optimization suggestions
     - Validate formulation requirements
  
  4. Simulation Setup (SiennaSimulations)  
     - Create PowerSimulations templates
     - Set up optimization problems
     - Configure solvers from config
  
  5. Execute Simulations
     - Run Economic Dispatch and/or Unit Commitment
     - Save results with full metadata
     - Generate performance reports

CONFIGURATION REQUIREMENTS:
  • config.toml with all required sections
  • CSV data files: bus.csv, load.csv, gen.csv (optional: branch.csv, dc_branch.csv)
  • timeseries_metadata.json (if time series enabled)

BENEFITS OF NEW ARCHITECTURE:
  • Each class handles one responsibility
  • Easy to test individual components  
  • Configuration passed cleanly through all layers
  • Better error messages and debugging
  • Extensible for new features
  • Follows modern Julia best practices
    """
end

"""
    main()

Main execution function with command-line argument handling.
"""
function main()
    @info "🚀 SIENNA WORKFLOW - CLASS-BASED ARCHITECTURE v2.0"
    @info "="^70
    
    if length(ARGS) > 0
        arg = ARGS[1]
        
        if arg == "test"
            @info "Running test mode..."
            result = run_workflow_test()
            if result.success
                @info "✅ TEST PASSED"
                exit(0)
            else
                @error "❌ TEST FAILED"
                exit(1)
            end
            
        elseif arg == "validate"
            @info "Running validation only..."
            success = run_workflow_validation_only()
            if success
                @info "✅ VALIDATION PASSED"
                exit(0)
            else
                @error "❌ VALIDATION FAILED"
                exit(1)
            end
            
        elseif arg == "help"
            show_help()
            exit(0)
            
        else
            @warn "Unknown argument: $arg"
            @info "Valid options: test, validate, help"
            show_help()
            exit(1)
        end
    else
        @info "Running full workflow..."
        result = run_sienna_workflow()
        
        if result.success
            @info "✅ WORKFLOW COMPLETED SUCCESSFULLY"
            exit(0)
        else
            @error "❌ WORKFLOW FAILED"
            exit(1)
        end
    end
end

# Export main functions
export run_sienna_workflow, run_workflow_test, run_workflow_validation_only
export SiennaWorkflowResults, print_workflow_summary, save_workflow_summary

# Run main if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
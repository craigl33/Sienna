#!/usr/bin/env julia

"""
SiennaDebug.jl - Clean Power System Debugging

Simple, class-based debugger with flexible constructors.
"""

using PowerSystems
using PowerSimulations
const PSI = PowerSimulations
using HiGHS
using Dates
using Logging
using JSON3
using DataFrames
using CSV
using Printf
using JuMP
using MathOptInterface
const MOI = MathOptInterface

# Import existing modules
include("SiennaConfig.jl")
include("SiennaSystem.jl")
include("SiennaSimulations.jl")

# ===== MAIN DEBUGGER CLASS =====

"""
    SiennaDebug

Power system debugger with flexible inputs.
"""
mutable struct SiennaDebug
    config::SiennaConfig
    system::SiennaSystem
    simulation::Union{SiennaSimulations, Nothing}
    output_dir::String
    session_name::String
    quick_results::Union{Nothing, NamedTuple}
    lp_exports::Dict{String, String}
    solve_attempts::Dict{String, NamedTuple}
    models::Dict{String, Any}
    created_at::DateTime
    last_operation::String
end

# ===== CONSTRUCTOR =====

"""
    SiennaDebug(config, system=nothing, simulation=nothing)

Create debugger with flexible inputs.
"""
function SiennaDebug(
    config::Union{String, SiennaConfig},
    system::Union{Nothing, SiennaSystem} = nothing,
    simulation::Union{Nothing, SiennaSimulations} = nothing
)
    @info "🔍 Creating SiennaDebug..."
    
    # Handle config
    config_obj = if config isa String
        @info "Loading config from: $config"
        SiennaConfig(config)
    else
        config
    end
    
    # Handle system
    system_obj = if system === nothing
        @info "Creating system from config"
        sys = SiennaSystem(config_obj)
        if !is_system_built(sys)
            build_system!(sys)
        end
        sys
    else
        if !is_system_built(system)
            build_system!(system)
        end
        system
    end
    
    # Create session
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    session_name = "debug_session_$timestamp"
    output_dir = joinpath(config_obj.output_directory, session_name)
    mkpath(output_dir)
    
    debugger = SiennaDebug(
        config_obj,
        system_obj,
        simulation,
        output_dir,
        session_name,
        nothing,
        Dict{String, String}(),
        Dict{String, NamedTuple}(),
        Dict{String, Any}(),
        now(),
        "initialized"
    )
    
    @info "✅ SiennaDebug ready: $session_name"
    return debugger
end

# ===== JUMP MODEL EXTRACTION =====

"""
    get_jump_model(psi_model)

Extract the JuMP model from a PowerSimulations DecisionModel.
"""
function get_jump_model(psi_model)
    try
        # Try accessing through get_optimization_container
        container = PSI.get_optimization_container(psi_model)
        return container.JuMPmodel
    catch e1
        try
            # Try direct access through internal structure
            return psi_model.internal.container.JuMPmodel
        catch e2
            try
                # Alternative field name
                container = PSI.get_optimization_container(psi_model)
                return container.jump_model
            catch e3
                @error "Cannot extract JuMP model from PowerSimulations DecisionModel"
                @error "Error 1: $e1"
                @error "Error 2: $e2" 
                @error "Error 3: $e3"
                
                # Debug info
                @info "Model type: $(typeof(psi_model))"
                if hasfield(typeof(psi_model), :internal)
                    @info "Internal type: $(typeof(psi_model.internal))"
                    if hasfield(typeof(psi_model.internal), :container)
                        @info "Container type: $(typeof(psi_model.internal.container))"
                        @info "Container fields: $(fieldnames(typeof(psi_model.internal.container)))"
                    end
                end
                
                rethrow(e1)
            end
        end
    end
end

# ===== MAIN METHODS =====

"""
    quick_check(debugger::SiennaDebug)

Fast diagnosis of common issues.
"""
function quick_check(debugger::SiennaDebug)
    @info "🔍 Running quick check..."
    
    if debugger.quick_results !== nothing
        @info "Using cached results"
        return debugger.quick_results
    end
    
    issues = String[]
    recommendations = String[]
    
    check_capacity(debugger, issues, recommendations)
    check_solver(debugger, issues, recommendations)
    check_network(debugger, issues, recommendations)
    check_timeseries(debugger, issues, recommendations)
    
    debugger.quick_results = (
        has_issues = !isempty(issues),
        issues = issues,
        recommendations = recommendations,
        timestamp = now()
    )
    
    print_results(debugger.quick_results)
    save_results(debugger)
    
    return debugger.quick_results
end

"""
    export_lp(debugger::SiennaDebug, formulation::String="economic_dispatch")

Export LP/MPS files.
"""
function export_lp(debugger::SiennaDebug, formulation::String="economic_dispatch")
    @info "📁 Exporting LP files for $formulation..."
    
    if haskey(debugger.lp_exports, formulation)
        @info "Using cached files"
        return (success = true, lp_file = debugger.lp_exports[formulation], cached = true)
    end
    
    form_dir = joinpath(debugger.output_dir, "lp_files", formulation)
    mkpath(form_dir)
    
    try
        # Get or create simulation
        if debugger.simulation === nothing
            debugger.simulation = SiennaSimulations(debugger.config, debugger.system)
        end
        
        # Build model
        model_key = "$(formulation)_model"
        if !haskey(debugger.models, model_key)
            @info "Building model..."
            template = create_template_for_formulation(debugger.simulation, formulation)
            model = create_decision_model(debugger.simulation, template, model_key)
            build!(model, output_dir = form_dir)
            debugger.models[model_key] = model
        else
            @info "Using cached model"
            model = debugger.models[model_key]
        end
        
        # Export files - FIXED: Extract JuMP model first
        lp_file = joinpath(form_dir, "$(formulation).lp")
        mps_file = joinpath(form_dir, "$(formulation).mps")
        
        # Extract the JuMP model from PowerSimulations DecisionModel
        jump_model = get_jump_model(model)
        
        # Now use JuMP's write_to_file with the extracted JuMP model
        write_to_file(jump_model, lp_file)
        write_to_file(jump_model, mps_file)
        
        debugger.lp_exports[formulation] = lp_file
        
        @info "✅ Files exported:"
        @info "   LP: $lp_file"
        @info "   MPS: $mps_file"
        
        # Try solve
        solve_result = try_solve(debugger, formulation)
        
        return (
            success = true,
            lp_file = lp_file,
            mps_file = mps_file,
            solve_result = solve_result
        )
        
    catch e
        @error "❌ Export failed: $e"
        return (success = false, error = string(e))
    end
end

"""
    try_solve(debugger::SiennaDebug, formulation::String)

Attempt to solve model.
"""
function try_solve(debugger::SiennaDebug, formulation::String)
    @info "🔍 Attempting solve..."
    
    model_key = "$(formulation)_model"
    if !haskey(debugger.models, model_key)
        return (solved = false, error = "Model not built")
    end
    
    if haskey(debugger.solve_attempts, formulation)
        return debugger.solve_attempts[formulation]
    end
    
    model = debugger.models[model_key]
    
    try
        solve_time = @elapsed solve!(model)
        status = primal_status(model)
        
        if status == MOI.FEASIBLE_POINT
            obj_val = objective_value(model)
            @info "✅ Solved: obj = $(round(obj_val, digits=2))"
            result = (solved = true, objective = obj_val, solve_time = solve_time)
        else
            @warn "❌ Infeasible: $status"
            result = (solved = false, status = string(status), solve_time = solve_time)
        end
        
        debugger.solve_attempts[formulation] = result
        return result
        
    catch e
        @error "❌ Solve failed: $e"
        result = (solved = false, error = string(e))
        debugger.solve_attempts[formulation] = result
        return result
    end
end

"""
    full_debug(debugger::SiennaDebug)

Comprehensive debugging analysis.
"""
function full_debug(debugger::SiennaDebug)
    @info "🔍 Full debug analysis..."
    
    # Run checks
    if debugger.quick_results === nothing
        quick_check(debugger)
    end
    
    # Export both formulations
    export_lp(debugger, "ed")
    export_lp(debugger, "uc")
    
    # Create report
    report_file = create_report(debugger)
    
    @info "✅ Full debug complete: $report_file"
    
    return (
        success = true,
        output_dir = debugger.output_dir,
        report_file = report_file,
        results = debugger.quick_results
    )
end

# ===== HELPER FUNCTIONS =====

function check_capacity(debugger, issues, recommendations)
    @info "🔋 Checking capacity..."
    
    try
        sys = get_power_system(debugger.system)
        
        thermal_cap = sum(get_active_power_limits(gen).max for gen in get_components(ThermalStandard, sys))
        renewable_cap = sum(get_rating(gen) for gen in get_components(RenewableDispatch, sys))
        total_gen = thermal_cap + renewable_cap
        total_load = sum(get_max_active_power(load) for load in get_components(PowerLoad, sys))
        
        @info "  Generation: $(round(total_gen, digits=1)) MW"
        @info "  Load: $(round(total_load, digits=1)) MW"
        
        if total_gen < total_load
            push!(issues, "Insufficient generation: $(round(total_gen, digits=1)) MW < $(round(total_load, digits=1)) MW")
            push!(recommendations, "Add generation capacity")
        else
            @info "  ✅ Capacity adequate"
        end
        
    catch e
        push!(issues, "Capacity check failed: $e")
    end
end

function check_solver(debugger, issues, recommendations)
    @info "⚙️  Checking solver..."
    
    try
        settings = get_solver_settings(debugger.config, "economic_dispatch")
        
        if settings["time_limit"] < 60
            push!(issues, "Short time limit: $(settings["time_limit"])s")
            push!(recommendations, "Increase solver time limit")
        end
        
        @info "  ✅ Solver settings OK"
        
    catch e
        push!(issues, "Solver check failed: $e")
    end
end

function check_network(debugger, issues, recommendations)
    @info "🌐 Checking network..."
    
    try
        sys = get_power_system(debugger.system)
        settings = get_network_settings(debugger.config)
        
        n_lines = length(get_components(Line, sys)) + length(get_components(TwoTerminalHVDCLine, sys))
        model = settings["model"]
        
        if model != "CopperPlatePowerModel" && n_lines == 0
            push!(issues, "No transmission lines with $model")
            push!(recommendations, "Add lines or use CopperPlatePowerModel")
        else
            @info "  ✅ Network consistent"
        end
        
    catch e
        push!(issues, "Network check failed: $e")
    end
end

function check_timeseries(debugger, issues, recommendations)
    @info "📊 Checking timeseries..."
    
    try
        if debugger.config.load_timeseries
            if !has_time_series_data(debugger.system)
                push!(issues, "No timeseries data loaded")
                push!(recommendations, "Check timeseries files")
            else
                @info "  ✅ Timeseries OK"
            end
        else
            @info "  ℹ️  Timeseries disabled"
        end
        
    catch e
        push!(issues, "Timeseries check failed: $e")
    end
end

function print_results(results)
    println("\n" * "="^50)
    println("QUICK CHECK SUMMARY")
    println("="^50)
    
    if !results.has_issues
        println("✅ No critical issues found")
    else
        println("❌ Issues found:")
        for (i, issue) in enumerate(results.issues)
            println("  $i. $issue")
        end
        
        if !isempty(results.recommendations)
            println("\n💡 Recommendations:")
            for (i, rec) in enumerate(results.recommendations)
                println("  $i. $rec")
            end
        end
    end
    
    println("="^50 * "\n")
end

function save_results(debugger)
    if debugger.quick_results !== nothing
        file = joinpath(debugger.output_dir, "quick_check.json")
        open(file, "w") do f
            JSON3.pretty(f, debugger.quick_results)
        end
    end
end

function create_report(debugger)
    report_file = joinpath(debugger.output_dir, "debug_report.txt")
    
    open(report_file, "w") do f
        println(f, "SIENNA DEBUG REPORT")
        println(f, "="^50)
        println(f, "Session: $(debugger.session_name)")
        println(f, "Generated: $(now())")
        println(f, "Project: $(debugger.config.project_name)")
        println(f, "="^50)
        println(f)
        
        if debugger.quick_results !== nothing && debugger.quick_results.has_issues
            println(f, "ISSUES FOUND:")
            for (i, issue) in enumerate(debugger.quick_results.issues)
                println(f, "  $i. $issue")
            end
            println(f)
            
            println(f, "RECOMMENDATIONS:")
            for (i, rec) in enumerate(debugger.quick_results.recommendations)
                println(f, "  $i. $rec")
            end
        else
            println(f, "✅ NO CRITICAL ISSUES FOUND")
        end
        
        println(f)
        println(f, "Output directory: $(debugger.output_dir)")
        println(f, "="^50)
    end
    
    return report_file
end

# ===== USAGE EXAMPLE =====
if abspath(PROGRAM_FILE) == @__FILE__
    println("🧪 Testing SiennaDebug.jl")
    
    if isfile("config.toml")
        try
            debugger = SiennaDebug("config.toml")
            issues = debugger.quick_check()
            debugger.export_lp("economic_dispatch")
            results = debugger.full_debug()
            
            println("✅ Test completed!")
            println("📁 Results: $(results.output_dir)")
            
        catch e
            println("❌ Test failed: $e")
            # Print stack trace for debugging
            for (exc, bt) in Base.catch_stack()
                showerror(stdout, exc, bt)
                println()
            end
        end
    else
        println("Usage:")
        println("  debugger = SiennaDebug(\"config.toml\")")
        println("  issues = debugger.quick_check()")
        println("  debugger.export_lp(\"economic_dispatch\")")
        println("  debugger.full_debug()")
    end
end
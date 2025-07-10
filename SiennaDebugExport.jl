#!/usr/bin/env julia

"""
SiennaDebugExport.jl - Model Export and Debugging for Sienna Ecosystem
======================================================================

Class-based debugging and model export manager that provides comprehensive
model introspection, export capabilities, and constraint analysis for 
PowerSimulations.jl models.

Features:
- LP/MPS model export for external solver analysis
- Constraint breakdown and analysis
- Infeasibility diagnosis tools
- Variable and parameter inspection
- Solver log capture and analysis
- JuMP model structure export
- Conflict analysis for infeasible models

Usage:
    config = SiennaConfig("config.toml")
    sienna_sys = SiennaSystem(config)
    sienna_sim = SiennaSimulations(config, sienna_sys)
    
    # Create debugger
    debugger = SiennaDebugExport(config, sienna_sim)
    
    # Export model for debugging
    export_model_for_debugging!(debugger, "economic_dispatch")
    
    # Analyze infeasibility
    analyze_infeasibility!(debugger, "economic_dispatch")
"""

using PowerSystems
using PowerSimulations
const PSI = PowerSimulations
using JuMP
using HiGHS
using Dates
using Logging
using JSON3
using DataFrames
using CSV
using Printf

# Import our modules
include("SiennaConfig.jl")
include("SiennaSimulations.jl")

"""
    SiennaDebugExport

Main debugging and export manager class for PowerSimulations.jl models.
"""
mutable struct SiennaDebugExport
    # Core components
    config::SiennaConfig
    sienna_sim::SiennaSimulations
    
    # Export settings
    export_directory::String
    export_timestamp::String
    
    # Model introspection data
    model_info::Dict{String, Dict{String, Any}}
    constraint_breakdown::Dict{String, DataFrame}
    variable_breakdown::Dict{String, DataFrame}
    
    # Debugging state
    last_formulation_analyzed::String
    debugging_enabled::Bool
    
    # Export tracking
    exported_files::Dict{String, Vector{String}}
    export_summary::Dict{String, Any}
end

"""
    SiennaDebugExport(config::SiennaConfig, sienna_sim::SiennaSimulations)

Constructor - Initialize debugging manager with config and simulation objects.
"""
function SiennaDebugExport(config::SiennaConfig, sienna_sim::SiennaSimulations)
    @info "🔍 Initializing SiennaDebugExport for model debugging..."
    
    # Create timestamped export directory
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    export_dir = joinpath(config.output_directory, "debug_exports_$timestamp")
    mkpath(export_dir)
    
    # Initialize debugger
    debugger = SiennaDebugExport(
        config,
        sienna_sim,
        export_dir,
        timestamp,
        Dict{String, Dict{String, Any}}(),      # model_info
        Dict{String, DataFrame}(),               # constraint_breakdown
        Dict{String, DataFrame}(),               # variable_breakdown
        "",                                      # last_formulation_analyzed
        get(get(config.config_data, "debug", Dict()), "debug_mode", false),  # debugging_enabled with safe fallback
        Dict{String, Vector{String}}(),          # exported_files
        Dict{String, Any}()                      # export_summary
    )
    
    @info "✅ SiennaDebugExport initialized"
    @info "   Export directory: $(debugger.export_directory)"
    @info "   Debug mode: $(debugger.debugging_enabled)"
    
    return debugger
end

"""
    export_model_for_debugging!(debugger::SiennaDebugExport, formulation_type::String)

Export complete model in multiple formats for debugging analysis.
"""
function export_model_for_debugging!(debugger::SiennaDebugExport, formulation_type::String)
    @info "🔧 Exporting $formulation_type model for debugging analysis..."
    
    # Check if problem exists and is built
    if !haskey(debugger.sienna_sim.problems, formulation_type)
        @error "Problem $formulation_type not found. Build the problem first."
        return false
    end
    
    if !get(debugger.sienna_sim.problems_built, formulation_type, false)
        @error "Problem $formulation_type not built. Build the problem first."
        return false
    end
    
    problem = debugger.sienna_sim.problems[formulation_type]
    debugger.last_formulation_analyzed = formulation_type
    
    # Create formulation-specific directory
    form_dir = joinpath(debugger.export_directory, formulation_type)
    mkpath(form_dir)
    
    debugger.exported_files[formulation_type] = String[]
    
    try
        # 1. Export LP/MPS files using PowerSimulations.jl API
        export_optimization_files!(debugger, problem, form_dir, formulation_type)
        
        # 2. Export JuMP model structure
        if get(debugger.config.config_data["debug"], "export_jump_model", false)
            export_jump_model_structure!(debugger, problem, form_dir, formulation_type)
        end
        
        # 3. Export constraint breakdown
        export_constraint_breakdown!(debugger, problem, form_dir, formulation_type)
        
        # 4. Export variable breakdown  
        export_variable_breakdown!(debugger, problem, form_dir, formulation_type)
        
        # 5. Export model statistics
        export_model_statistics!(debugger, problem, form_dir, formulation_type)
        
        # 6. Export solver settings used
        export_solver_settings!(debugger, form_dir, formulation_type)
        
        # 7. Save system information
        export_system_information!(debugger, form_dir, formulation_type)
        
        @info "✅ Model export completed for $formulation_type"
        @info "   Files exported to: $form_dir"
        @info "   Total files: $(length(debugger.exported_files[formulation_type]))"
        
        return true
        
    catch e
        @error "❌ Model export failed for $formulation_type: $e"
        return false
    end
end

"""
    export_optimization_files!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)

Export LP and MPS files using PowerSimulations.jl built-in functions.
"""
function export_optimization_files!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)
    @info "Exporting optimization model files (LP/MPS)..."
    
    debug_config = get(debugger.config.config_data, "debug", Dict())
    
    try
        # Get the JuMP model from the DecisionModel
        jump_model = PSI.get_jump_model(problem)
        
        if jump_model === nothing
            @warn "No JuMP model found in DecisionModel - model may not be built"
            return
        end
        
        # Export MPS file (universal format for most solvers)
        if get(debug_config, "export_optimization_model", true)
            mps_file = joinpath(form_dir, "$(formulation_type)_model.mps")
            @info "  Exporting MPS file: $mps_file"
            write_to_file(jump_model, mps_file)
            push!(debugger.exported_files[formulation_type], mps_file)
            @info "  ✓ MPS file exported successfully"
        end
        
        # Export LP file (human-readable format)
        if get(debug_config, "export_lp_model", true)
            lp_file = joinpath(form_dir, "$(formulation_type)_model.lp")
            @info "  Exporting LP file: $lp_file"
            write_to_file(jump_model, lp_file)
            push!(debugger.exported_files[formulation_type], lp_file)
            @info "  ✓ LP file exported successfully"
        end
        
        @info "✓ Optimization files exported successfully"
        
    catch e
        @error "Failed to export optimization files: $e"
        # Try alternative approach using JuMP directly
        try_alternative_export!(debugger, problem, form_dir, formulation_type)
    end
end

"""
    try_alternative_export!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)

Alternative export method if main method fails.
"""
function try_alternative_export!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)
    @info "Trying alternative export method..."
    
    try
        # Try to get model through PowerSimulations internal API
        container = PSI.get_optimization_container(problem)
        if container !== nothing
            jump_model = PSI.get_jump_model(container)
            if jump_model !== nothing
                # Export using JuMP directly
                mps_file = joinpath(form_dir, "$(formulation_type)_model_alt.mps")
                write_to_file(jump_model, mps_file)
                push!(debugger.exported_files[formulation_type], mps_file)
                @info "  ✓ Alternative MPS export successful"
            end
        end
    catch e
        @warn "Alternative export also failed: $e"
    end
end

"""
    export_jump_model_structure!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)

Export detailed JuMP model structure for analysis.
"""
function export_jump_model_structure!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)
    @info "Exporting JuMP model structure..."
    
    try
        jump_model = PSI.get_jump_model(problem)
        if jump_model === nothing
            @warn "No JuMP model available for structure export"
            return
        end
        
        # Create comprehensive model structure report
        structure_file = joinpath(form_dir, "$(formulation_type)_model_structure.txt")
        
        open(structure_file, "w") do f
            println(f, "="^80)
            println(f, "JUMP MODEL STRUCTURE REPORT")
            println(f, "="^80)
            println(f, "Generated: $(now())")
            println(f, "Formulation: $formulation_type")
            println(f, "="^80)
            println(f)
            
            # Model summary
            println(f, "MODEL SUMMARY:")
            println(f, "-"^40)
            println(f, "Objective sense: $(objective_sense(jump_model))")
            println(f, "Number of variables: $(num_variables(jump_model))")
            println(f, "Number of constraints: $(num_constraints(jump_model))")
            println(f, "Solver: $(solver_name(jump_model))")
            println(f)
            
            # Variables breakdown
            println(f, "VARIABLES BREAKDOWN:")
            println(f, "-"^40)
            var_types = Dict{String, Int}()
            for (i, var) in enumerate(all_variables(jump_model))
                var_type = string(typeof(var))
                var_types[var_type] = get(var_types, var_type, 0) + 1
                if i <= 20  # Show first 20 variables
                    println(f, "  $(name(var)): $(var_type)")
                elseif i == 21
                    println(f, "  ... (and $(num_variables(jump_model) - 20) more variables)")
                    break
                end
            end
            
            println(f)
            println(f, "Variable types summary:")
            for (vtype, count) in var_types
                println(f, "  $vtype: $count")
            end
            println(f)
            
            # Constraints breakdown
            println(f, "CONSTRAINTS BREAKDOWN:")
            println(f, "-"^40)
            constraint_types = constraint_type_counts(jump_model)
            for (ctype, count) in constraint_types
                println(f, "  $ctype: $count")
            end
            println(f)
            
            # Objective function
            println(f, "OBJECTIVE FUNCTION:")
            println(f, "-"^40)
            obj = objective_function(jump_model)
            obj_str = string(obj)
            if length(obj_str) > 1000
                println(f, obj_str[1:1000] * "... (truncated)")
            else
                println(f, obj_str)
            end
            println(f)
            
            println(f, "="^80)
            println(f, "END OF STRUCTURE REPORT")
            println(f, "="^80)
        end
        
        push!(debugger.exported_files[formulation_type], structure_file)
        @info "  ✓ JuMP model structure exported to: $structure_file"
        
    catch e
        @error "Failed to export JuMP model structure: $e"
    end
end

"""
    constraint_type_counts(model::JuMP.Model)

Get counts of each constraint type in the model.
"""
function constraint_type_counts(model::JuMP.Model)
    counts = Dict{String, Int}()
    
    # Get all constraint types
    for (F, S) in list_of_constraint_types(model)
        constraint_name = "$(F)-in-$(S)"
        constraint_count = num_constraints(model, F, S)
        counts[constraint_name] = constraint_count
    end
    
    return counts
end

"""
    export_constraint_breakdown!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)

Export detailed constraint analysis.
"""
function export_constraint_breakdown!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)
    @info "Exporting constraint breakdown..."
    
    try
        jump_model = PSI.get_jump_model(problem)
        if jump_model === nothing
            @warn "No JuMP model available for constraint analysis"
            return
        end
        
        # Create constraint breakdown DataFrame
        constraint_data = []
        
        for (F, S) in list_of_constraint_types(jump_model)
            constraint_refs = all_constraints(jump_model, F, S)
            constraint_name = "$(F)-in-$(S)"
            constraint_count = length(constraint_refs)
            
            # Sample a few constraints for detailed analysis
            sample_constraints = constraint_refs[1:min(5, length(constraint_refs))]
            sample_names = [name(c) for c in sample_constraints]
            
            push!(constraint_data, Dict(
                "constraint_type" => constraint_name,
                "function_type" => string(F),
                "set_type" => string(S),
                "count" => constraint_count,
                "sample_names" => join(sample_names, "; ")
            ))
        end
        
        # Create DataFrame and save
        constraint_df = DataFrame(constraint_data)
        debugger.constraint_breakdown[formulation_type] = constraint_df
        
        # Save to CSV
        constraint_file = joinpath(form_dir, "$(formulation_type)_constraints.csv")
        CSV.write(constraint_file, constraint_df)
        push!(debugger.exported_files[formulation_type], constraint_file)
        
        # Save detailed constraint report
        detailed_file = joinpath(form_dir, "$(formulation_type)_constraints_detailed.txt")
        export_detailed_constraints!(jump_model, detailed_file)
        push!(debugger.exported_files[formulation_type], detailed_file)
        
        @info "  ✓ Constraint breakdown exported ($(nrow(constraint_df)) constraint types)"
        
    catch e
        @error "Failed to export constraint breakdown: $e"
    end
end

"""
    export_detailed_constraints!(model::JuMP.Model, output_file::String)

Export detailed constraint information to text file.
"""
function export_detailed_constraints!(model::JuMP.Model, output_file::String)
    open(output_file, "w") do f
        println(f, "="^80)
        println(f, "DETAILED CONSTRAINT ANALYSIS")
        println(f, "="^80)
        println(f, "Generated: $(now())")
        println(f, "="^80)
        println(f)
        
        for (F, S) in list_of_constraint_types(model)
            constraint_refs = all_constraints(model, F, S)
            constraint_name = "$(F)-in-$(S)"
            
            println(f, "CONSTRAINT TYPE: $constraint_name")
            println(f, "-"^60)
            println(f, "Function Type: $F")
            println(f, "Set Type: $S")
            println(f, "Count: $(length(constraint_refs))")
            println(f)
            
            # Show first few constraints in detail
            max_show = min(10, length(constraint_refs))
            println(f, "Sample constraints (first $max_show):")
            for (i, cref) in enumerate(constraint_refs[1:max_show])
                println(f, "  [$i] $(name(cref))")
                try
                    # Try to get constraint function and set
                    constraint_obj = constraint_object(cref)
                    func = constraint_obj.func
                    set = constraint_obj.set
                    println(f, "      Function: $(func)")
                    println(f, "      Set: $(set)")
                catch e
                    println(f, "      (Details unavailable: $e)")
                end
                println(f)
            end
            
            if length(constraint_refs) > max_show
                println(f, "  ... and $(length(constraint_refs) - max_show) more constraints")
            end
            
            println(f, "="^60)
            println(f)
        end
    end
end

"""
    export_variable_breakdown!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)

Export detailed variable analysis.
"""
function export_variable_breakdown!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)
    @info "Exporting variable breakdown..."
    
    try
        jump_model = PSI.get_jump_model(problem)
        if jump_model === nothing
            @warn "No JuMP model available for variable analysis"
            return
        end
        
        # Create variable breakdown
        variable_data = []
        variable_types = Dict{String, Int}()
        
        for var in all_variables(jump_model)
            var_name = name(var)
            var_type = string(typeof(var))
            variable_types[var_type] = get(variable_types, var_type, 0) + 1
            
            # Get variable bounds
            lower_bound_val = has_lower_bound(var) ? JuMP.lower_bound(var) : "unbounded"
            upper_bound_val = has_upper_bound(var) ? JuMP.upper_bound(var) : "unbounded"
            is_binary_val = is_binary(var)
            is_integer_val = is_integer(var)
            
            push!(variable_data, Dict(
                "variable_name" => var_name,
                "variable_type" => var_type,
                "lower_bound" => string(lower_bound_val),
                "upper_bound" => string(upper_bound_val),
                "is_binary" => is_binary_val,
                "is_integer" => is_integer_val
            ))
        end
        
        # Create DataFrame and save
        variable_df = DataFrame(variable_data)
        debugger.variable_breakdown[formulation_type] = variable_df
        
        # Save to CSV
        variable_file = joinpath(form_dir, "$(formulation_type)_variables.csv")
        CSV.write(variable_file, variable_df)
        push!(debugger.exported_files[formulation_type], variable_file)
        
        # Save variable types summary
        summary_file = joinpath(form_dir, "$(formulation_type)_variable_summary.txt")
        open(summary_file, "w") do f
            println(f, "VARIABLE TYPES SUMMARY")
            println(f, "="^50)
            println(f, "Total variables: $(num_variables(jump_model))")
            println(f)
            for (vtype, count) in sort(collect(variable_types), by=x->x[2], rev=true)
                println(f, @sprintf("%-40s: %6d", vtype, count))
            end
        end
        push!(debugger.exported_files[formulation_type], summary_file)
        
        @info "  ✓ Variable breakdown exported ($(nrow(variable_df)) variables)"
        
    catch e
        @error "Failed to export variable breakdown: $e"
    end
end

"""
    export_model_statistics!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)

Export comprehensive model statistics.
"""
function export_model_statistics!(debugger::SiennaDebugExport, problem::DecisionModel, form_dir::String, formulation_type::String)
    @info "Exporting model statistics..."
    
    try
        jump_model = PSI.get_jump_model(problem)
        if jump_model === nothing
            @warn "No JuMP model available for statistics"
            return
        end
        
        # Collect comprehensive statistics
        stats = Dict{String, Any}(
            "timestamp" => string(now()),
            "formulation_type" => formulation_type,
            "model_summary" => Dict(
                "num_variables" => num_variables(jump_model),
                "num_constraints" => num_constraints(jump_model),
                "objective_sense" => string(objective_sense(jump_model)),
                "solver_name" => string(solver_name(jump_model))
            ),
            "variable_statistics" => get_variable_statistics(jump_model),
            "constraint_statistics" => get_constraint_statistics(jump_model),
            "model_size_metrics" => get_model_size_metrics(jump_model)
        )
        
        # Save as JSON
        stats_file = joinpath(form_dir, "$(formulation_type)_model_statistics.json")
        open(stats_file, "w") do f
            JSON3.pretty(f, stats)
        end
        push!(debugger.exported_files[formulation_type], stats_file)
        
        # Save human-readable summary
        summary_file = joinpath(form_dir, "$(formulation_type)_model_summary.txt")
        write_model_summary(stats, summary_file)
        push!(debugger.exported_files[formulation_type], summary_file)
        
        # Store in debugger
        debugger.model_info[formulation_type] = stats
        
        @info "  ✓ Model statistics exported"
        
    catch e
        @error "Failed to export model statistics: $e"
    end
end

"""
    get_variable_statistics(model::JuMP.Model)

Get detailed variable statistics.
"""
function get_variable_statistics(model::JuMP.Model)
    binary_count = 0
    integer_count = 0
    continuous_count = 0
    bounded_count = 0
    unbounded_count = 0
    
    for var in all_variables(model)
        if is_binary(var)
            binary_count += 1
        elseif is_integer(var)
            integer_count += 1
        else
            continuous_count += 1
        end
        
        if has_lower_bound(var) || has_upper_bound(var)
            bounded_count += 1
        else
            unbounded_count += 1
        end
    end
    
    return Dict(
        "total" => num_variables(model),
        "binary" => binary_count,
        "integer" => integer_count,
        "continuous" => continuous_count,
        "bounded" => bounded_count,
        "unbounded" => unbounded_count
    )
end

"""
    get_constraint_statistics(model::JuMP.Model)

Get detailed constraint statistics.
"""
function get_constraint_statistics(model::JuMP.Model)
    constraint_types = constraint_type_counts(model)
    
    # Categorize constraints
    equality_count = 0
    inequality_count = 0
    set_count = 0
    
    for (F, S) in list_of_constraint_types(model)
        count = num_constraints(model, F, S)
        set_string = string(S)
        if occursin("EqualTo", set_string)
            equality_count += count
        elseif occursin("LessThan", set_string) || occursin("GreaterThan", set_string)
            inequality_count += count
        else
            set_count += count
        end
    end
    
    return Dict(
        "total" => num_constraints(model),
        "equality" => equality_count,
        "inequality" => inequality_count,
        "set_constraints" => set_count,
        "constraint_types" => constraint_types
    )
end

"""
    get_model_size_metrics(model::JuMP.Model)

Get model size and complexity metrics.
"""
function get_model_size_metrics(model::JuMP.Model)
    # Basic size metrics
    n_vars = num_variables(model)
    n_constraints = num_constraints(model)
    
    # Calculate matrix density (approximate)
    # This is a rough estimate - actual density calculation would require matrix access
    estimated_nonzeros = n_vars * n_constraints * 0.1  # Assume 10% density
    
    return Dict(
        "variables" => n_vars,
        "constraints" => n_constraints,
        "variable_constraint_ratio" => n_vars / max(n_constraints, 1),
        "estimated_matrix_density" => 0.1,  # Placeholder
        "estimated_nonzeros" => Int(round(estimated_nonzeros)),
        "model_complexity" => categorize_model_complexity(n_vars, n_constraints)
    )
end

"""
    categorize_model_complexity(n_vars::Int, n_constraints::Int)

Categorize model complexity based on size.
"""
function categorize_model_complexity(n_vars::Int, n_constraints::Int)
    total_elements = n_vars + n_constraints
    
    if total_elements < 1000
        return "Small"
    elseif total_elements < 10000
        return "Medium"
    elseif total_elements < 100000
        return "Large"
    else
        return "Very Large"
    end
end

"""
    write_model_summary(stats::Dict, output_file::String)

Write human-readable model summary.
"""
function write_model_summary(stats::Dict, output_file::String)
    open(output_file, "w") do f
        println(f, "="^80)
        println(f, "POWERSIMULATIONS.JL MODEL SUMMARY")
        println(f, "="^80)
        println(f, "Generated: $(stats["timestamp"])")
        println(f, "Formulation: $(stats["formulation_type"])")
        println(f, "="^80)
        println(f)
        
        # Model summary
        model_sum = stats["model_summary"]
        println(f, "MODEL OVERVIEW:")
        println(f, "-"^40)
        println(f, @sprintf("%-25s: %s", "Objective Sense", model_sum["objective_sense"]))
        println(f, @sprintf("%-25s: %s", "Solver", model_sum["solver_name"]))
        println(f, @sprintf("%-25s: %d", "Total Variables", model_sum["num_variables"]))
        println(f, @sprintf("%-25s: %d", "Total Constraints", model_sum["num_constraints"]))
        println(f)
        
        # Variable statistics
        var_stats = stats["variable_statistics"]
        println(f, "VARIABLE BREAKDOWN:")
        println(f, "-"^40)
        println(f, @sprintf("%-25s: %d", "Binary Variables", var_stats["binary"]))
        println(f, @sprintf("%-25s: %d", "Integer Variables", var_stats["integer"]))
        println(f, @sprintf("%-25s: %d", "Continuous Variables", var_stats["continuous"]))
        println(f, @sprintf("%-25s: %d", "Bounded Variables", var_stats["bounded"]))
        println(f, @sprintf("%-25s: %d", "Unbounded Variables", var_stats["unbounded"]))
        println(f)
        
        # Constraint statistics
        con_stats = stats["constraint_statistics"]
        println(f, "CONSTRAINT BREAKDOWN:")
        println(f, "-"^40)
        println(f, @sprintf("%-25s: %d", "Equality Constraints", con_stats["equality"]))
        println(f, @sprintf("%-25s: %d", "Inequality Constraints", con_stats["inequality"]))
        println(f, @sprintf("%-25s: %d", "Set Constraints", con_stats["set_constraints"]))
        println(f)
        
        # Model size metrics
        size_metrics = stats["model_size_metrics"]
        println(f, "MODEL SIZE METRICS:")
        println(f, "-"^40)
        println(f, @sprintf("%-25s: %.2f", "Var/Constraint Ratio", size_metrics["variable_constraint_ratio"]))
        println(f, @sprintf("%-25s: %s", "Model Complexity", size_metrics["model_complexity"]))
        println(f, @sprintf("%-25s: %d", "Est. Matrix Nonzeros", size_metrics["estimated_nonzeros"]))
        println(f)
        
        # Top constraint types
        println(f, "TOP CONSTRAINT TYPES:")
        println(f, "-"^40)
        constraint_types = sort(collect(con_stats["constraint_types"]), by=x->x[2], rev=true)
        for (i, (ctype, count)) in enumerate(constraint_types[1:min(10, length(constraint_types))])
            println(f, @sprintf("%-2d. %-35s: %d", i, ctype, count))
        end
        
        println(f)
        println(f, "="^80)
    end
end

function export_solver_settings!(debugger::SiennaDebugExport, form_dir::String, formulation_type::String)
    @info "Exporting solver settings..."
    
    try
        # Get solver settings from config
        solver_settings = get_solver_settings(debugger.config, formulation_type)
        
        # Create comprehensive solver settings report
        settings_file = joinpath(form_dir, "$(formulation_type)_solver_settings.json")
        settings_report = Dict{String, Any}(
            "timestamp" => string(now()),
            "formulation_type" => formulation_type,
            "solver_name" => debugger.config.default_solver,
            "solver_settings" => solver_settings,
            "config_debug_settings" => get(debugger.config.config_data, "debug", Dict()),
            "config_advanced_settings" => get(debugger.config.config_data, "advanced", Dict())
        )
        
        open(settings_file, "w") do f
            JSON3.pretty(f, settings_report)
        end
        push!(debugger.exported_files[formulation_type], settings_file)
        
        @info "  ✓ Solver settings exported"
        
    catch e
        @error "Failed to export solver settings: $e"
    end
end

"""
    export_system_information!(debugger::SiennaDebugExport, form_dir::String, formulation_type::String)

Export PowerSystems.jl system information for debugging context.
"""
function export_system_information!(debugger::SiennaDebugExport, form_dir::String, formulation_type::String)
    @info "Exporting system information..."
    
    try
        # Get system from sienna_sim
        sys = get_power_system(debugger.sienna_sim.sienna_system)
        
        # Collect system information
        system_info = Dict{String, Any}(
            "timestamp" => string(now()),
            "system_name" => get_name(sys),
            "base_power" => get_base_power(sys),
            "units_base" => string(get_units_base(sys)),
            "component_counts" => get_component_counts(debugger.sienna_sim.sienna_system),
            "time_series_loaded" => has_timeseries_data(debugger.sienna_sim.sienna_system),
            "forecast_available" => has_forecast_data(debugger.sienna_sim.sienna_system),
            "network_model" => debugger.config.network_model,
            "horizon_hours" => debugger.config.default_horizon_hours
        )
        
        # Add generation summary
        try
            gen_summary = get_generator_fuel_summary(debugger.sienna_sim.sienna_system)
            system_info["generator_summary"] = gen_summary
        catch e
            @warn "Could not get generator summary: $e"
        end
        
        # Save system info
        info_file = joinpath(form_dir, "$(formulation_type)_system_info.json")
        open(info_file, "w") do f
            JSON3.pretty(f, system_info)
        end
        push!(debugger.exported_files[formulation_type], info_file)
        
        @info "  ✓ System information exported"
        
    catch e
        @error "Failed to export system information: $e"
    end
end

"""
    analyze_infeasibility!(debugger::SiennaDebugExport, formulation_type::String)

Analyze model infeasibility and provide debugging insights.
"""
function analyze_infeasibility!(debugger::SiennaDebugExport, formulation_type::String)
    @info "🔍 Analyzing model infeasibility for $formulation_type..."
    
    # Check if problem exists
    if !haskey(debugger.sienna_sim.problems, formulation_type)
        @error "Problem $formulation_type not found"
        return false
    end
    
    problem = debugger.sienna_sim.problems[formulation_type]
    
    try
        jump_model = PSI.get_jump_model(problem)
        if jump_model === nothing
            @error "No JuMP model available for infeasibility analysis"
            return false
        end
        
        # Create infeasibility analysis directory
        infeas_dir = joinpath(debugger.export_directory, formulation_type, "infeasibility_analysis")
        mkpath(infeas_dir)
        
        # 1. Check current optimization status
        analyze_termination_status!(debugger, jump_model, infeas_dir, formulation_type)
        
        # 2. Run conflict analysis if supported
        run_conflict_analysis!(debugger, jump_model, infeas_dir, formulation_type)
        
        # 3. Check for common infeasibility patterns
        check_common_infeasibility_patterns!(debugger, jump_model, infeas_dir, formulation_type)
        
        # 4. Generate infeasibility report
        generate_infeasibility_report!(debugger, infeas_dir, formulation_type)
        
        @info "✅ Infeasibility analysis completed"
        @info "   Analysis files saved to: $infeas_dir"
        
        return true
        
    catch e
        @error "❌ Infeasibility analysis failed: $e"
        return false
    end
end

"""
    analyze_termination_status!(debugger::SiennaDebugExport, model::JuMP.Model, infeas_dir::String, formulation_type::String)

Analyze and report termination status.
"""
function analyze_termination_status!(debugger::SiennaDebugExport, model::JuMP.Model, infeas_dir::String, formulation_type::String)
    @info "Analyzing termination status..."
    
    try
        status_file = joinpath(infeas_dir, "termination_status.txt")
        
        open(status_file, "w") do f
            println(f, "="^80)
            println(f, "TERMINATION STATUS ANALYSIS")
            println(f, "="^80)
            println(f, "Generated: $(now())")
            println(f, "Formulation: $formulation_type")
            println(f, "="^80)
            println(f)
            
            # Get termination status
            term_status = termination_status(model)
            primal_status = primal_status(model)
            dual_status = dual_status(model)
            
            println(f, "OPTIMIZATION STATUS:")
            println(f, "-"^40)
            println(f, "Termination Status: $term_status")
            println(f, "Primal Status: $primal_status")
            println(f, "Dual Status: $dual_status")
            println(f)
            
            # Interpret status
            println(f, "STATUS INTERPRETATION:")
            println(f, "-"^40)
            
            if term_status == MOI.INFEASIBLE
                println(f, "❌ MODEL IS INFEASIBLE")
                println(f, "   The problem has no feasible solution.")
                println(f, "   Check constraints for conflicts.")
            elseif term_status == MOI.DUAL_INFEASIBLE
                println(f, "❌ MODEL IS DUAL INFEASIBLE (UNBOUNDED)")
                println(f, "   The objective can be improved without bound.")
                println(f, "   Check for missing constraints.")
            elseif term_status == MOI.INFEASIBLE_OR_UNBOUNDED
                println(f, "❌ MODEL IS INFEASIBLE OR UNBOUNDED")
                println(f, "   Need to determine which case applies.")
            elseif term_status == MOI.TIME_LIMIT
                println(f, "⏱️  SOLVER TIME LIMIT REACHED")
                println(f, "   Solution may be suboptimal or infeasible.")
                println(f, "   Consider increasing time limit.")
            elseif term_status == MOI.OPTIMAL
                println(f, "✅ MODEL SOLVED TO OPTIMALITY")
                println(f, "   No infeasibility detected.")
            else
                println(f, "⚠️  UNEXPECTED STATUS: $term_status")
                println(f, "   Manual investigation required.")
            end
            
            println(f)
            
            # Get objective value if available
            try
                obj_val = objective_value(model)
                println(f, "Objective Value: $obj_val")
            catch
                println(f, "Objective Value: Not available")
            end
            
            println(f)
            println(f, "="^80)
        end
        
        push!(debugger.exported_files[formulation_type], status_file)
        @info "  ✓ Termination status analysis saved"
        
    catch e
        @error "Failed to analyze termination status: $e"
    end
end

"""
    run_conflict_analysis!(debugger::SiennaDebugExport, model::JuMP.Model, infeas_dir::String, formulation_type::String)

Run conflict analysis to identify infeasible constraint sets.
"""
function run_conflict_analysis!(debugger::SiennaDebugExport, model::JuMP.Model, infeas_dir::String, formulation_type::String)
    @info "Running conflict analysis..."
    
    try
        # Check if solver supports conflict analysis
        if !MOI.supports(backend(model), MOI.ConflictStatus())
            @warn "Solver does not support conflict analysis"
            return
        end
        
        # Run conflict analysis
        @info "Computing conflict..."
        compute_conflict!(model)
        
        # Get conflict status
        conflict_status = MOI.get(backend(model), MOI.ConflictStatus())
        
        conflict_file = joinpath(infeas_dir, "conflict_analysis.txt")
        open(conflict_file, "w") do f
            println(f, "="^80)
            println(f, "CONFLICT ANALYSIS RESULTS")
            println(f, "="^80)
            println(f, "Generated: $(now())")
            println(f, "Conflict Status: $conflict_status")
            println(f, "="^80)
            println(f)
            
            if conflict_status == MOI.CONFLICT_FOUND
                println(f, "✅ CONFLICT FOUND - IDENTIFYING CONFLICTING CONSTRAINTS")
                println(f, "="^60)
                println(f)
                
                # Get conflicting constraints
                conflicting_constraints = []
                
                for (F, S) in list_of_constraint_types(model)
                    constraints = all_constraints(model, F, S)
                    for constraint in constraints
                        try
                            conflict_status = MOI.get(backend(model), MOI.ConstraintConflictStatus(), constraint)
                            if conflict_status == MOI.IN_CONFLICT
                                push!(conflicting_constraints, constraint)
                                println(f, "CONFLICTING CONSTRAINT: $(name(constraint))")
                                println(f, "  Type: $(F)-in-$(S)")
                                
                                # Try to get constraint details
                                try
                                    constraint_obj = constraint_object(constraint)
                                    println(f, "  Function: $(constraint_obj.func)")
                                    println(f, "  Set: $(constraint_obj.set)")
                                catch
                                    println(f, "  (Details unavailable)")
                                end
                                println(f)
                            end
                        catch e
                            # Some constraint types may not support conflict status
                            continue
                        end
                    end
                end
                
                println(f, "="^60)
                println(f, "SUMMARY:")
                println(f, "Total conflicting constraints: $(length(conflicting_constraints))")
                
                if length(conflicting_constraints) > 0
                    println(f)
                    println(f, "RECOMMENDATIONS:")
                    println(f, "1. Review the conflicting constraints above")
                    println(f, "2. Check if constraint bounds are too restrictive")
                    println(f, "3. Verify data consistency in the input files")
                    println(f, "4. Consider relaxing some constraints")
                end
                
            else
                println(f, "❌ NO CONFLICT FOUND OR CONFLICT ANALYSIS FAILED")
                println(f, "Status: $conflict_status")
                println(f)
                println(f, "POSSIBLE REASONS:")
                println(f, "1. Model may not actually be infeasible")
                println(f, "2. Solver may not support full conflict analysis")
                println(f, "3. Numerical issues preventing conflict detection")
            end
            
            println(f)
            println(f, "="^80)
        end
        
        push!(debugger.exported_files[formulation_type], conflict_file)
        @info "  ✓ Conflict analysis completed"
        
    catch e
        @warn "Conflict analysis failed: $e"
        # Create a file noting the failure
        conflict_file = joinpath(infeas_dir, "conflict_analysis_failed.txt")
        open(conflict_file, "w") do f
            println(f, "CONFLICT ANALYSIS FAILED")
            println(f, "="^40)
            println(f, "Error: $e")
            println(f, "This may be due to:")
            println(f, "1. Solver limitations")
            println(f, "2. Model structure issues")
            println(f, "3. Numerical problems")
        end
        push!(debugger.exported_files[formulation_type], conflict_file)
    end
end

"""
    check_common_infeasibility_patterns!(debugger::SiennaDebugExport, model::JuMP.Model, infeas_dir::String, formulation_type::String)

Check for common patterns that cause infeasibility in power system models.
"""
function check_common_infeasibility_patterns!(debugger::SiennaDebugExport, model::JuMP.Model, infeas_dir::String, formulation_type::String)
    @info "Checking common infeasibility patterns..."
    
    try
        patterns_file = joinpath(infeas_dir, "common_patterns_check.txt")
        
        open(patterns_file, "w") do f
            println(f, "="^80)
            println(f, "COMMON INFEASIBILITY PATTERNS CHECK")
            println(f, "="^80)
            println(f, "Generated: $(now())")
            println(f, "="^80)
            println(f)
            
            # Pattern 1: Check for zero-capacity generation
            println(f, "1. ZERO-CAPACITY GENERATION CHECK:")
            println(f, "-"^50)
            check_zero_capacity_generation!(debugger, f)
            println(f)
            
            # Pattern 2: Check load vs generation capacity
            println(f, "2. LOAD VS GENERATION CAPACITY:")
            println(f, "-"^50)
            check_load_generation_balance!(debugger, f)
            println(f)
            
            # Pattern 3: Check for units base inconsistencies
            println(f, "3. UNITS BASE CONSISTENCY:")
            println(f, "-"^50)
            check_units_base_consistency!(debugger, f)
            println(f)
            
            # Pattern 4: Check for missing time series
            println(f, "4. TIME SERIES DATA:")
            println(f, "-"^50)
            check_time_series_issues!(debugger, f)
            println(f)
            
            # Pattern 5: Check network connectivity
            println(f, "5. NETWORK CONNECTIVITY:")
            println(f, "-"^50)
            check_network_connectivity!(debugger, f)
            println(f)
            
            # Pattern 6: Check constraint bounds
            println(f, "6. CONSTRAINT BOUNDS:")
            println(f, "-"^50)
            check_constraint_bounds!(model, f)
            println(f)
            
            println(f, "="^80)
        end
        
        push!(debugger.exported_files[formulation_type], patterns_file)
        @info "  ✓ Common patterns check completed"
        
    catch e
        @error "Failed to check common patterns: $e"
    end
end

"""
    check_zero_capacity_generation!(debugger::SiennaDebugExport, f::IOStream)

Check for generators with zero capacity.
"""
function check_zero_capacity_generation!(debugger::SiennaDebugExport, f::IOStream)
    try
        sys = get_power_system(debugger.sienna_sim.sienna_system)
        
        # Check thermal generators
        thermal_gens = get_components(ThermalStandard, sys)
        zero_thermal = []
        for gen in thermal_gens
            max_power = get_active_power_limits(gen).max
            if max_power <= 0.0
                push!(zero_thermal, get_name(gen))
            end
        end
        
        # Check renewable generators
        renewable_gens = get_components(RenewableDispatch, sys)
        zero_renewable = []
        for gen in renewable_gens
            rating = get_rating(gen)
            if rating <= 0.0
                push!(zero_renewable, get_name(gen))
            end
        end
        
        if !isempty(zero_thermal) || !isempty(zero_renewable)
            println(f, "⚠️  FOUND ZERO-CAPACITY GENERATORS:")
            if !isempty(zero_thermal)
                println(f, "   Thermal: $(join(zero_thermal, ", "))")
            end
            if !isempty(zero_renewable)
                println(f, "   Renewable: $(join(zero_renewable, ", "))")
            end
            println(f, "   → This can cause infeasibility in unit commitment")
        else
            println(f, "✅ No zero-capacity generators found")
        end
        
    catch e
        println(f, "❌ Could not check generation capacity: $e")
    end
end

"""
    check_load_generation_balance!(debugger::SiennaDebugExport, f::IOStream)

Check if total generation capacity can meet total load.
"""
function check_load_generation_balance!(debugger::SiennaDebugExport, f::IOStream)
    try
        sys = get_power_system(debugger.sienna_sim.sienna_system)
        
        # Calculate total generation capacity
        total_gen_capacity = 0.0
        
        # Thermal generation
        for gen in get_components(ThermalStandard, sys)
            total_gen_capacity += get_active_power_limits(gen).max
        end
        
        # Renewable generation
        for gen in get_components(RenewableDispatch, sys)
            total_gen_capacity += get_rating(gen)
        end
        
        for gen in get_components(RenewableNonDispatch, sys)
            total_gen_capacity += get_rating(gen)
        end
        
        # Calculate total load
        total_load = 0.0
        for load in get_components(PowerLoad, sys)
            total_load += get_max_active_power(load)
        end
        
        reserve_margin = (total_gen_capacity - total_load) / total_load * 100
        
        println(f, "Total Generation Capacity: $(round(total_gen_capacity, digits=1)) MW")
        println(f, "Total Load: $(round(total_load, digits=1)) MW")
        println(f, "Reserve Margin: $(round(reserve_margin, digits=1))%")
        
        if reserve_margin < 0
            println(f, "❌ INSUFFICIENT GENERATION CAPACITY")
            println(f, "   → Model will be infeasible - not enough generation to meet load")
        elseif reserve_margin < 15
            println(f, "⚠️  LOW RESERVE MARGIN")
            println(f, "   → Model may be infeasible with additional constraints")
        else
            println(f, "✅ Adequate generation capacity")
        end
        
    catch e
        println(f, "❌ Could not check load/generation balance: $e")
    end
end

"""
    check_units_base_consistency!(debugger::SiennaDebugExport, f::IOStream)

Check for units base consistency issues.
"""
function check_units_base_consistency!(debugger::SiennaDebugExport, f::IOStream)
    try
        sys = get_power_system(debugger.sienna_sim.sienna_system)
        
        units_base = string(get_units_base(sys))
        base_power = get_base_power(sys)
        
        println(f, "System Units Base: $units_base")
        println(f, "System Base Power: $(base_power) MW")
        
        if units_base != "NATURAL_UNITS"
            println(f, "⚠️  NON-NATURAL UNITS DETECTED")
            println(f, "   → May cause scaling issues in optimization")
            println(f, "   → Consider using NATURAL_UNITS for PowerSimulations")
        else
            println(f, "✅ Using NATURAL_UNITS (recommended)")
        end
        
        if base_power != 100.0
            println(f, "⚠️  NON-STANDARD BASE POWER")
            println(f, "   → Standard is 100 MW for NATURAL_UNITS")
        else
            println(f, "✅ Standard base power (100 MW)")
        end
        
    catch e
        println(f, "❌ Could not check units base: $e")
    end
end

"""
    check_time_series_issues!(debugger::SiennaDebugExport, f::IOStream)

Check for time series related issues.
"""
function check_time_series_issues!(debugger::SiennaDebugExport, f::IOStream)
    try
        has_ts = has_timeseries_data(debugger.sienna_sim.sienna_system)
        has_forecasts = has_forecast_data(debugger.sienna_sim.sienna_system)
        
        println(f, "Time Series Loaded: $has_ts")
        println(f, "Forecasts Available: $has_forecasts")
        
        if !has_ts && debugger.config.load_timeseries
            println(f, "❌ TIME SERIES LOADING REQUESTED BUT NOT LOADED")
            println(f, "   → Check timeseries_metadata.json and data files")
        elseif !has_forecasts && has_ts
            println(f, "⚠️  TIME SERIES LOADED BUT NO FORECASTS")
            println(f, "   → May cause issues in PowerSimulations")
        else
            println(f, "✅ Time series configuration appears correct")
        end
        
    catch e
        println(f, "❌ Could not check time series: $e")
    end
end

"""
    check_network_connectivity!(debugger::SiennaDebugExport, f::IOStream)

Check basic network connectivity.
"""
function check_network_connectivity!(debugger::SiennaDebugExport, f::IOStream)
    try
        sys = get_power_system(debugger.sienna_sim.sienna_system)
        
        n_buses = length(get_components(ACBus, sys))
        n_lines = length(get_components(Line, sys))
        n_dc_lines = length(get_components(TwoTerminalHVDCLine, sys))
        
        println(f, "Buses: $n_buses")
        println(f, "AC Lines: $n_lines") 
        println(f, "DC Lines: $n_dc_lines")
        
        network_model = debugger.config.network_model
        println(f, "Network Model: $network_model")
        
        if network_model == "CopperPlatePowerModel"
            println(f, "✅ Using CopperPlate - network topology irrelevant")
        elseif n_lines == 0 && n_dc_lines == 0
            println(f, "⚠️  NO TRANSMISSION LINES WITH NON-COPPERPLATE MODEL")
            println(f, "   → May cause network flow issues")
        else
            # Basic connectivity check - ensure minimum spanning tree exists
            min_lines_needed = n_buses - 1
            total_lines = n_lines + n_dc_lines
            
            if total_lines < min_lines_needed
                println(f, "❌ INSUFFICIENT LINES FOR CONNECTIVITY")
                println(f, "   → Need at least $min_lines_needed lines for $n_buses buses")
            else
                println(f, "✅ Sufficient lines for basic connectivity")
            end
        end
        
    catch e
        println(f, "❌ Could not check network connectivity: $e")
    end
end

"""
    check_constraint_bounds!(model::JuMP.Model, f::IOStream)

Check for obviously infeasible constraint bounds.
"""
function check_constraint_bounds!(model::JuMP.Model, f::IOStream)
    try
        infeasible_bounds_count = 0
        
        # Check variable bounds
        for var in all_variables(model)
            if has_lower_bound(var) && has_upper_bound(var)
                lb = JuMP.lower_bound(var)
                ub = JuMP.upper_bound(var)
                if lb > ub
                    infeasible_bounds_count += 1
                    if infeasible_bounds_count <= 5  # Show first 5 examples
                        println(f, "❌ Variable $(name(var)): LB=$lb > UB=$ub")
                    end
                end
            end
        end
        
        if infeasible_bounds_count > 0
            println(f, "❌ FOUND $infeasible_bounds_count VARIABLES WITH INFEASIBLE BOUNDS")
            if infeasible_bounds_count > 5
                println(f, "   (showing first 5 examples)")
            end
        else
            println(f, "✅ No obviously infeasible variable bounds found")
        end
        
    catch e
        println(f, "❌ Could not check constraint bounds: $e")
    end
end

"""
    generate_infeasibility_report!(debugger::SiennaDebugExport, infeas_dir::String, formulation_type::String)

Generate comprehensive infeasibility report with recommendations.
"""
function generate_infeasibility_report!(debugger::SiennaDebugExport, infeas_dir::String, formulation_type::String)
    @info "Generating comprehensive infeasibility report..."
    
    try
        report_file = joinpath(infeas_dir, "infeasibility_diagnosis_report.txt")
        
        open(report_file, "w") do f
            println(f, "="^80)
            println(f, "SIENNA INFEASIBILITY DIAGNOSIS REPORT")
            println(f, "="^80)
            println(f, "Generated: $(now())")
            println(f, "Formulation: $formulation_type")
            println(f, "Project: $(debugger.config.project_name)")
            println(f, "="^80)
            println(f)
            
            println(f, "DIAGNOSIS SUMMARY:")
            println(f, "-"^50)
            println(f, "This report analyzes potential causes of model infeasibility")
            println(f, "in your PowerSimulations.jl optimization model.")
            println(f)
            
            println(f, "FILES GENERATED:")
            println(f, "-"^50)
            println(f, "1. termination_status.txt - Solver termination analysis")
            println(f, "2. conflict_analysis.txt - Conflicting constraints (if supported)")
            println(f, "3. common_patterns_check.txt - Common infeasibility patterns")
            println(f, "4. This report - Comprehensive diagnosis")
            println(f)
            
            println(f, "RECOMMENDED DEBUGGING STEPS:")
            println(f, "-"^50)
            println(f, "1. CHECK TERMINATION STATUS:")
            println(f, "   → Review 'termination_status.txt'")
            println(f, "   → Identify if model is infeasible, unbounded, or other")
            println(f)
            
            println(f, "2. EXAMINE CONFLICTING CONSTRAINTS:")
            println(f, "   → Review 'conflict_analysis.txt' if available")
            println(f, "   → Focus on constraints listed as 'IN_CONFLICT'")
            println(f, "   → Check your CSV data for inconsistencies")
            println(f)
            
            println(f, "3. VERIFY COMMON ISSUES:")
            println(f, "   → Review 'common_patterns_check.txt'")
            println(f, "   → Pay attention to any '❌' or '⚠️' warnings")
            println(f, "   → Fix data issues identified")
            println(f)
            
            println(f, "4. ANALYZE EXPORTED MODEL FILES:")
            println(f, "   → Open $(formulation_type)_model.lp in a text editor")
            println(f, "   → Load $(formulation_type)_model.mps in external solver (Gurobi, CPLEX)")
            println(f, "   → Use external solver's infeasibility analysis tools")
            println(f)
            
            println(f, "5. CHECK DATA CONSISTENCY:")
            println(f, "   → Verify bus names match between CSV files")
            println(f, "   → Ensure generator capacities are reasonable")
            println(f, "   → Check that total generation ≥ total load")
            println(f, "   → Validate time series data if used")
            println(f)
            
            println(f, "COMMON FIXES:")
            println(f, "-"^50)
            println(f, "• Increase generator capacities")
            println(f, "• Add slack variables for debugging (set add_slack_variables=true)")
            println(f, "• Use CopperPlatePowerModel for initial testing")
            println(f, "• Disable time series initially to isolate issues")
            println(f, "• Check for typos in CSV files")
            println(f, "• Ensure consistent units (MW, not kW)")
            println(f)
            
            println(f, "EXTERNAL SOLVER ANALYSIS:")
            println(f, "-"^50)
            println(f, "For advanced debugging, load the .mps file in:")
            println(f, "• Gurobi: model.computeIIS() for irreducible infeasible set")
            println(f, "• CPLEX: conflict refiner tool")
            println(f, "• SCIP: infeasibility analysis plugins")
            println(f)
            
            println(f, "SIENNA-SPECIFIC DEBUGGING:")
            println(f, "-"^50)
            println(f, "• Enable debug mode in config.toml")
            println(f, "• Set solver output_flag=true to see solver messages")
            println(f, "• Try simpler formulation first (Economic Dispatch before Unit Commitment)")
            println(f, "• Validate PowerSystems.jl system with validate_system=true")
            println(f)
            
            println(f, "="^80)
            println(f, "END OF DIAGNOSIS REPORT")
            println(f, "="^80)
        end
        
        push!(debugger.exported_files[formulation_type], report_file)
        @info "  ✓ Comprehensive infeasibility report generated"
        
    catch e
        @error "Failed to generate infeasibility report: $e"
    end
end

"""
    capture_solver_logs!(debugger::SiennaDebugExport, formulation_type::String)

Capture and save detailed solver logs for analysis.
"""
function capture_solver_logs!(debugger::SiennaDebugExport, formulation_type::String)
    @info "Setting up solver log capture for $formulation_type..."
    
    try
        # Create logs directory
        logs_dir = joinpath(debugger.export_directory, formulation_type, "solver_logs")
        mkpath(logs_dir)
        
        # Check if solver logging is enabled in config
        debug_config = get(debugger.config.config_data, "debug", Dict())
        save_solver_logs = get(debug_config, "save_solver_logs", true)
        
        if save_solver_logs
            # Configure solver to output to file
            log_file = joinpath(logs_dir, "$(formulation_type)_solver.log")
            
            # Update solver settings to enable logging
            problem = debugger.sienna_sim.problems[formulation_type]
            jump_model = PSI.get_jump_model(problem)
            
            if jump_model !== nothing
                # Enable solver output
                set_attribute(jump_model, "output_flag", true)
                set_attribute(jump_model, "log_to_console", true)
                
                @info "  ✓ Solver logging configured"
                @info "    Log file: $log_file"
            end
        end
        
    catch e
        @warn "Failed to setup solver log capture: $e"
    end
end

"""
    export_debug_summary!(debugger::SiennaDebugExport)

Export a comprehensive summary of all debugging information.
"""
function export_debug_summary!(debugger::SiennaDebugExport)
    @info "💾 Exporting comprehensive debug summary..."
    
    try
        summary_file = joinpath(debugger.export_directory, "debug_export_summary.json")
        
        # Create comprehensive summary
        summary = Dict{String, Any}(
            "export_metadata" => Dict(
                "timestamp" => string(now()),
                "export_directory" => debugger.export_directory,
                "debug_enabled" => debugger.debugging_enabled,
                "sienna_version" => "2.0",
                "config_file" => debugger.config.config_file_path
            ),
            "project_info" => Dict(
                "name" => debugger.config.project_name,
                "description" => debugger.config.project_description,
                "version" => debugger.config.project_version,
                "network_model" => debugger.config.network_model,
                "horizon_hours" => debugger.config.default_horizon_hours
            ),
            "exported_formulations" => collect(keys(debugger.exported_files)),
            "exported_files_by_formulation" => debugger.exported_files,
            "model_statistics" => debugger.model_info,
            "total_files_exported" => sum(length(files) for files in values(debugger.exported_files)),
            "constraint_summaries" => Dict(
                form => Dict(
                    "total_constraint_types" => nrow(df),
                    "total_constraints" => sum(df.count)
                ) for (form, df) in debugger.constraint_breakdown
            ),
            "variable_summaries" => Dict(
                form => Dict(
                    "total_variables" => nrow(df)
                ) for (form, df) in debugger.variable_breakdown
            )
        )
        
        # Add system information
        try
            sys = get_power_system(debugger.sienna_sim.sienna_system)
            summary["system_info"] = Dict(
                "name" => get_name(sys),
                "base_power" => get_base_power(sys),
                "units_base" => string(get_units_base(sys)),
                "component_counts" => get_component_counts(debugger.sienna_sim.sienna_system),
                "time_series_loaded" => has_timeseries_data(debugger.sienna_sim.sienna_system),
                "forecast_available" => has_forecast_data(debugger.sienna_sim.sienna_system)
            )
        catch e
            @warn "Could not add system info to summary: $e"
        end
        
        # Save summary
        open(summary_file, "w") do f
            JSON3.pretty(f, summary)
        end
        
        debugger.export_summary = summary
        
        @info "✅ Debug export summary saved"
        @info "   Summary file: $summary_file"
        @info "   Total formulations exported: $(length(debugger.exported_files))"
        @info "   Total files exported: $(summary["total_files_exported"])"
        
        # Print file listing
        print_export_file_listing(debugger)
        
    catch e
        @error "Failed to export debug summary: $e"
    end
end

"""
    print_export_file_listing(debugger::SiennaDebugExport)

Print a organized listing of all exported files.
"""
function print_export_file_listing(debugger::SiennaDebugExport)
    @info "\n" * "="^80
    @info "SIENNA DEBUG EXPORT FILE LISTING"
    @info "="^80
    @info "Export Directory: $(debugger.export_directory)"
    @info "="^80
    
    for (formulation, files) in debugger.exported_files
        @info "\n$(uppercase(formulation)) FILES ($(length(files)) files):"
        @info "-"^60
        
        # Group files by type
        file_groups = Dict{String, Vector{String}}()
        
        for file in files
            filename = basename(file)
            if endswith(filename, ".mps") || endswith(filename, ".lp")
                group = "Optimization Models"
            elseif endswith(filename, ".csv")
                group = "Data Tables"
            elseif endswith(filename, ".json")
                group = "Structured Data"
            elseif endswith(filename, ".txt")
                group = "Analysis Reports"
            else
                group = "Other Files"
            end
            
            if !haskey(file_groups, group)
                file_groups[group] = String[]
            end
            push!(file_groups[group], filename)
        end
        
        # Print each group
        for (group_name, group_files) in sort(collect(file_groups))
            @info "  $group_name:"
            for file in sort(group_files)
                @info "    • $file"
            end
        end
    end
    
    @info "\n" * "="^80
    @info "Use these files for:"
    @info "• Loading .mps/.lp files in external solvers (Gurobi, CPLEX)"
    @info "• Analyzing constraint/variable breakdowns in Excel/Python"
    @info "• Understanding model structure and infeasibility causes"
    @info "• Debugging PowerSystems.jl data issues"
    @info "="^80
end

# ===== INTEGRATION WITH SIENNA SIMULATIONS =====

"""
    enable_debugging_for_simulation!(sienna_sim::SiennaSimulations, enable_export::Bool=true)

Enable debugging mode for a SiennaSimulations object with automatic export.
"""
function enable_debugging_for_simulation!(sienna_sim::SiennaSimulations, enable_export::Bool=true)
    @info "🔧 Enabling debugging mode for SiennaSimulations..."
    
    # Create debugger
    debugger = SiennaDebugExport(sienna_sim.config, sienna_sim)
    
    if enable_export
        @info "   Auto-export enabled - models will be exported after each build"
    end
    
    return debugger
end

"""
    debug_build_and_export!(debugger::SiennaDebugExport, formulation_type::String)

Build problem with debugging enabled and automatically export for analysis.
"""
function debug_build_and_export!(debugger::SiennaDebugExport, formulation_type::String)
    @info "🔧 Building and exporting $formulation_type with debugging..."
    
    try
        # Build the problem
        build_problem!(debugger.sienna_sim, formulation_type)
        
        # Immediately export for debugging
        export_model_for_debugging!(debugger, formulation_type)
        
        @info "✅ Debug build and export completed for $formulation_type"
        return true
        
    catch e
        @error "❌ Debug build failed for $formulation_type: $e"
        
        # Even if build fails, try to export whatever we can
        try
            @info "Attempting to export partial information for debugging..."
            export_system_information!(debugger, 
                joinpath(debugger.export_directory, formulation_type), 
                formulation_type)
        catch export_e
            @warn "Could not export partial information: $export_e"
        end
        
        return false
    end
end

"""
    debug_solve_with_analysis!(debugger::SiennaDebugExport, formulation_type::String)

Solve problem with debugging and automatically run infeasibility analysis if needed.
"""
function debug_solve_with_analysis!(debugger::SiennaDebugExport, formulation_type::String)
    @info "⚡ Solving $formulation_type with debugging and analysis..."
    
    try
        # Ensure problem is built and exported
        if !get(debugger.sienna_sim.problems_built, formulation_type, false)
            debug_build_and_export!(debugger, formulation_type)
        end
        
        # Setup solver logging
        capture_solver_logs!(debugger, formulation_type)
        
        # Attempt to solve
        results = solve_problem!(debugger.sienna_sim, formulation_type)
        
        @info "✅ Solve completed successfully for $formulation_type"
        return results
        
    catch e
        @error "❌ Solve failed for $formulation_type: $e"
        
        # Run infeasibility analysis
        @info "Running infeasibility analysis due to solve failure..."
        analyze_infeasibility!(debugger, formulation_type)
        
        # Export comprehensive debug information
        export_debug_summary!(debugger)
        
        @info "🔍 Debugging files have been exported for analysis"
        @info "   Check: $(debugger.export_directory)"
        
        rethrow(e)
    end
end

# ===== UTILITY FUNCTIONS =====

"""
    clean_debug_exports!(output_directory::String, keep_recent::Int=5)

Clean old debug export directories, keeping only the most recent ones.
"""
function clean_debug_exports!(output_directory::String, keep_recent::Int=5)
    @info "🧹 Cleaning old debug export directories..."
    
    try
        # Find all debug export directories
        debug_dirs = []
        
        for item in readdir(output_directory, join=true)
            if isdir(item) && occursin("debug_exports_", basename(item))
                push!(debug_dirs, item)
            end
        end
        
        if length(debug_dirs) <= keep_recent
            @info "   Only $(length(debug_dirs)) debug directories found, keeping all"
            return
        end
        
        # Sort by modification time (newest first)
        sort!(debug_dirs, by=dir -> stat(dir).mtime, rev=true)
        
        # Remove old directories
        dirs_to_remove = debug_dirs[(keep_recent+1):end]
        
        for dir in dirs_to_remove
            @info "   Removing old debug directory: $(basename(dir))"
            rm(dir, recursive=true)
        end
        
        @info "✅ Cleaned $(length(dirs_to_remove)) old debug directories"
        @info "   Kept $(keep_recent) most recent directories"
        
    catch e
        @warn "Failed to clean debug exports: $e"
    end
end

"""
    quick_infeasibility_check(config::SiennaConfig, sienna_sys::SiennaSystem)

Quick check for obvious infeasibility issues before building optimization model.
"""
function quick_infeasibility_check(config::SiennaConfig, sienna_sys::SiennaSystem)
    @info "🔍 Running quick infeasibility check..."
    
    issues = String[]
    
    try
        sys = get_power_system(sienna_sys)
        
        # Check 1: Load vs Generation
        total_load = sum(get_max_active_power(load) for load in get_components(PowerLoad, sys))
        total_thermal = sum(get_active_power_limits(gen).max for gen in get_components(ThermalStandard, sys))
        total_renewable = sum(get_rating(gen) for gen in get_components(RenewableDispatch, sys))
        total_generation = total_thermal + total_renewable
        
        if total_generation < total_load
            push!(issues, "Insufficient generation capacity: $(round(total_generation, digits=1)) MW < $(round(total_load, digits=1)) MW")
        end
        
        # Check 2: Zero capacity generators
        zero_cap_gens = []
        for gen in get_components(ThermalStandard, sys)
            if get_active_power_limits(gen).max <= 0
                push!(zero_cap_gens, get_name(gen))
            end
        end
        
        if !isempty(zero_cap_gens)
            push!(issues, "Zero capacity thermal generators: $(join(zero_cap_gens[1:min(3, length(zero_cap_gens))], ", "))$(length(zero_cap_gens) > 3 ? "..." : "")")
        end
        
        # Check 3: Units consistency
        if string(get_units_base(sys)) != "NATURAL_UNITS"
            push!(issues, "Non-natural units may cause scaling issues: $(get_units_base(sys))")
        end
        
        # Check 4: Network model vs transmission
        n_lines = length(get_components(Line, sys))
        if config.network_model != "CopperPlatePowerModel" && n_lines == 0
            push!(issues, "Network model $(config.network_model) specified but no transmission lines found")
        end
        
        # Report results
        if isempty(issues)
            @info "✅ No obvious infeasibility issues detected"
            return true
        else
            @warn "⚠️  Potential infeasibility issues detected:"
            for issue in issues
                @warn "   • $issue"
            end
            return false
        end
        
    catch e
        @warn "Quick infeasibility check failed: $e"
        return false
    end
end

# ===== EXPORTS =====

export SiennaDebugExport
export export_model_for_debugging!, analyze_infeasibility!
export enable_debugging_for_simulation!, debug_build_and_export!, debug_solve_with_analysis!
export export_debug_summary!, clean_debug_exports!, quick_infeasibility_check

# ===== TESTING AND EXAMPLE USAGE =====

if abspath(PROGRAM_FILE) == @__FILE__
    @info "🧪 Testing SiennaDebugExport..."
    
    try
        # Test with existing config and system
        if isfile("config.toml")
            config = SiennaConfig("config.toml")
            
            if isdir(get_data_directory(config))
                @info "Testing with real system..."
                
                # Build system
                sienna_sys = SiennaSystem(config)
                build_system!(sienna_sys)
                
                # Run quick check
                quick_infeasibility_check(config, sienna_sys)
                
                # Create simulations
                sienna_sim = SiennaSimulations(config, sienna_sys)
                
                # Create debugger
                debugger = SiennaDebugExport(config, sienna_sim)
                
                @info "Testing debug export (without solving)..."
                
                # Test build and export without solving
                if debug_build_and_export!(debugger, "economic_dispatch")
                    export_debug_summary!(debugger)
                    @info "✅ SiennaDebugExport test passed"
                    @info "   Check exported files in: $(debugger.export_directory)"
                else
                    @warn "Debug export test had issues but completed"
                end
                
            else
                @info "Data directory not found - creating minimal test"
                
                # Test basic functionality
                debugger = SiennaDebugExport(config, SiennaSimulations(config, SiennaSystem(config)))
                @info "✅ SiennaDebugExport initialization test passed"
            end
        else
            @info "No config.toml found - cannot test SiennaDebugExport"
        end
    catch e
        @error "❌ SiennaDebugExport test failed: $e"
        
        # Print stack trace for debugging
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
end
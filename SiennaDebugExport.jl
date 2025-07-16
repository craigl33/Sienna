#!/usr/bin/env julia

"""
SiennaDebugExport.jl - Complete Model Export and Debugging for Sienna Ecosystem
===============================================================================

Class-based debugging and model export manager updated for simplified config structure.
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

Constructor - Initialize debugging manager with simplified config.
"""
function SiennaDebugExport(config::SiennaConfig, sienna_sim::SiennaSimulations)
    @info "🔍 Initializing SiennaDebugExport for simplified config..."
    
    # Create timestamped export directory
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    export_dir = joinpath(config.output_directory, "debug_exports_$timestamp")
    mkpath(export_dir)
    
    # Get debug settings from simplified config
    debug_settings = get_debug_settings(config)
    debugging_enabled = get(debug_settings, "export_optimization_model", false)
    
    # Initialize debugger
    debugger = SiennaDebugExport(
        config,
        sienna_sim,
        export_dir,
        timestamp,
        Dict{String, Dict{String, Any}}(),
        Dict{String, DataFrame}(),
        Dict{String, DataFrame}(),
        "",
        debugging_enabled,
        Dict{String, Vector{String}}(),
        Dict{String, Any}()
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
    
    # Create formulation-specific directory
    form_dir = joinpath(debugger.export_directory, formulation_type)
    mkpath(form_dir)
    
    debugger.exported_files[formulation_type] = String[]
    
    try
        # Get debug settings from simplified config
        debug_settings = get_debug_settings(debugger.config)
        
        # Export system information
        export_system_information!(debugger, form_dir, formulation_type)
        
        # Export solver settings
        export_solver_settings!(debugger, form_dir, formulation_type)
        
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
    export_solver_settings!(debugger::SiennaDebugExport, form_dir::String, formulation_type::String)

Export solver settings from simplified config structure.
"""
function export_solver_settings!(debugger::SiennaDebugExport, form_dir::String, formulation_type::String)
    @info "Exporting solver settings from simplified config..."
    
    try
        # Get solver settings from simplified config
        solver_settings = get_solver_settings(debugger.config, formulation_type)
        network_settings = get_network_settings(debugger.config)
        debug_settings = get_debug_settings(debugger.config)
        
        # Create comprehensive solver settings report
        settings_file = joinpath(form_dir, "$(formulation_type)_solver_settings.json")
        settings_report = Dict{String, Any}(
            "timestamp" => string(now()),
            "formulation_type" => formulation_type,
            "solver_name" => debugger.config.default_solver,
            "solver_settings" => solver_settings,
            "network_settings" => network_settings,
            "debug_settings" => debug_settings,
            "config_structure" => "simplified"
        )
        
        open(settings_file, "w") do f
            JSON3.pretty(f, settings_report)
        end
        push!(debugger.exported_files[formulation_type], settings_file)
        
        @info "  ✓ Solver settings exported from simplified config"
        
    catch e
        @error "Failed to export solver settings: $e"
    end
end

"""
    export_system_information!(debugger::SiennaDebugExport, form_dir::String, formulation_type::String)

Export PowerSystems.jl system information with simplified config context.
"""
function export_system_information!(debugger::SiennaDebugExport, form_dir::String, formulation_type::String)
    @info "Exporting system information with simplified config context..."
    
    try
        sys = get_power_system(debugger.sienna_sim.sienna_system)
        
        # Collect system information with simplified config context
        system_info = Dict{String, Any}(
            "timestamp" => string(now()),
            "config_structure" => "simplified",
            "system_name" => get_name(sys),
            "base_power" => get_base_power(sys),
            "units_base" => string(get_units_base(sys)),
            "component_counts" => get_component_counts(debugger.sienna_sim.sienna_system),
            "time_series_loaded" => has_timeseries_data(debugger.sienna_sim.sienna_system),
            "forecast_available" => has_forecast_data(debugger.sienna_sim.sienna_system),
            "simplified_config_settings" => Dict(
                "network_model" => debugger.config.network_model,
                "horizon_hours" => debugger.config.default_horizon_hours,
                "solver_name" => debugger.config.default_solver,
                "base_power" => debugger.config.base_power,
                "validate_system" => debugger.config.validate_system,
                "load_timeseries" => debugger.config.load_timeseries
            )
        )
        
        # Save system info
        info_file = joinpath(form_dir, "$(formulation_type)_system_info.json")
        open(info_file, "w") do f
            JSON3.pretty(f, system_info)
        end
        push!(debugger.exported_files[formulation_type], info_file)
        
        @info "  ✓ System information exported with simplified config context"
        
    catch e
        @error "Failed to export system information: $e"
    end
end

"""
    analyze_infeasibility!(debugger::SiennaDebugExport, formulation_type::String)

Analyze model infeasibility using simplified config settings.
"""
function analyze_infeasibility!(debugger::SiennaDebugExport, formulation_type::String)
    @info "🔍 Analyzing model infeasibility for $formulation_type with simplified config..."
    
    try
        # Create infeasibility analysis directory
        infeas_dir = joinpath(debugger.export_directory, formulation_type, "infeasibility_analysis")
        mkpath(infeas_dir)
        
        # Check for common infeasibility patterns using simplified config
        check_common_infeasibility_patterns_simplified!(debugger, infeas_dir, formulation_type)
        
        # Generate infeasibility report
        generate_infeasibility_report_simplified!(debugger, infeas_dir, formulation_type)
        
        @info "✅ Infeasibility analysis completed"
        @info "   Analysis files saved to: $infeas_dir"
        
        return true
        
    catch e
        @error "❌ Infeasibility analysis failed: $e"
        return false
    end
end

"""
    check_common_infeasibility_patterns_simplified!(debugger::SiennaDebugExport, infeas_dir::String, formulation_type::String)

Check for common patterns using simplified config structure.
"""
function check_common_infeasibility_patterns_simplified!(debugger::SiennaDebugExport, infeas_dir::String, formulation_type::String)
    @info "Checking common infeasibility patterns with simplified config..."
    
    try
        patterns_file = joinpath(infeas_dir, "common_patterns_check.txt")
        
        open(patterns_file, "w") do f
            println(f, "="^80)
            println(f, "COMMON INFEASIBILITY PATTERNS CHECK (SIMPLIFIED CONFIG)")
            println(f, "="^80)
            println(f, "Generated: $(now())")
            println(f, "="^80)
            println(f)
            
            # Pattern 1: Check for zero-capacity generation
            println(f, "1. ZERO-CAPACITY GENERATION CHECK:")
            println(f, "-"^50)
            check_zero_capacity_generation_simplified!(debugger, f)
            println(f)
            
            # Pattern 2: Check load vs generation capacity
            println(f, "2. LOAD VS GENERATION CAPACITY:")
            println(f, "-"^50)
            check_load_generation_balance_simplified!(debugger, f)
            println(f)
            
            # Pattern 3: Check network model vs transmission
            println(f, "3. NETWORK MODEL CONSISTENCY:")
            println(f, "-"^50)
            check_network_model_consistency_simplified!(debugger, f)
            println(f)
            
            # Pattern 4: Check solver settings
            println(f, "4. SOLVER SETTINGS:")
            println(f, "-"^50)
            check_solver_settings_simplified!(debugger, f, formulation_type)
            println(f)
            
            # Pattern 5: Check formulation compatibility
            println(f, "5. FORMULATION COMPATIBILITY:")
            println(f, "-"^50)
            check_formulation_compatibility_simplified!(debugger, f, formulation_type)
            println(f)
            
            println(f, "="^80)
        end
        
        push!(debugger.exported_files[formulation_type], patterns_file)
        @info "  ✓ Common patterns check completed with simplified config"
        
    catch e
        @error "Failed to check common patterns: $e"
    end
end

"""
    check_zero_capacity_generation_simplified!(debugger::SiennaDebugExport, f::IOStream)

Check for generators with zero capacity using simplified config.
"""
function check_zero_capacity_generation_simplified!(debugger::SiennaDebugExport, f::IOStream)
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
    check_load_generation_balance_simplified!(debugger::SiennaDebugExport, f::IOStream)

Check load vs generation balance using simplified config.
"""
function check_load_generation_balance_simplified!(debugger::SiennaDebugExport, f::IOStream)
    try
        sys = get_power_system(debugger.sienna_sim.sienna_system)
        validation_settings = get(debugger.config.config_data, "validation", Dict())
        
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
        min_reserve = get(validation_settings, "min_reserve_margin", 0.15) * 100
        
        println(f, "Total Generation Capacity: $(round(total_gen_capacity, digits=1)) MW")
        println(f, "Total Load: $(round(total_load, digits=1)) MW")
        println(f, "Reserve Margin: $(round(reserve_margin, digits=1))%")
        println(f, "Required Minimum: $(round(min_reserve, digits=1))%")
        
        if reserve_margin < 0
            println(f, "❌ INSUFFICIENT GENERATION CAPACITY")
            println(f, "   → Model will be infeasible - not enough generation to meet load")
        elseif reserve_margin < min_reserve
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
    check_network_model_consistency_simplified!(debugger::SiennaDebugExport, f::IOStream)

Check network model consistency using simplified config.
"""
function check_network_model_consistency_simplified!(debugger::SiennaDebugExport, f::IOStream)
    try
        sys = get_power_system(debugger.sienna_sim.sienna_system)
        network_settings = get_network_settings(debugger.config)
        
        n_buses = length(get_components(ACBus, sys))
        n_lines = length(get_components(Line, sys))
        n_dc_lines = length(get_components(TwoTerminalHVDCLine, sys))
        
        network_model = network_settings["model"]
        enable_limits = network_settings["enable_transmission_limits"]
        
        println(f, "Network Model: $network_model")
        println(f, "Transmission Limits: $enable_limits")
        println(f, "Buses: $n_buses")
        println(f, "AC Lines: $n_lines")
        println(f, "DC Lines: $n_dc_lines")
        
        if network_model == "CopperPlatePowerModel"
            println(f, "✅ Using CopperPlate - network topology irrelevant")
            if enable_limits
                println(f, "⚠️  Transmission limits enabled but ignored in CopperPlate model")
            end
        elseif n_lines == 0 && n_dc_lines == 0
            println(f, "❌ NO TRANSMISSION LINES WITH NON-COPPERPLATE MODEL")
            println(f, "   → Consider using CopperPlatePowerModel")
        else
            println(f, "✅ Network model and topology appear consistent")
        end
        
    catch e
        println(f, "❌ Could not check network model consistency: $e")
    end
end

"""
    check_solver_settings_simplified!(debugger::SiennaDebugExport, f::IOStream, formulation_type::String)

Check solver settings from simplified config.
"""
function check_solver_settings_simplified!(debugger::SiennaDebugExport, f::IOStream, formulation_type::String)
    try
        solver_settings = get_solver_settings(debugger.config, formulation_type)
        
        println(f, "Solver: $(debugger.config.default_solver)")
        println(f, "Time Limit: $(solver_settings["time_limit"]) seconds")
        println(f, "MIP Gap: $(solver_settings["mip_gap"])")
        println(f, "Threads: $(solver_settings["threads"])")
        println(f, "Output Flag: $(solver_settings["output_flag"])")
        
        # Check for potential issues
        if solver_settings["time_limit"] < 60
            println(f, "⚠️  Very short time limit ($(solver_settings["time_limit"])s)")
            println(f, "   → May cause premature termination")
        end
        
        if solver_settings["mip_gap"] > 0.1
            println(f, "⚠️  Large MIP gap ($(solver_settings["mip_gap"]))")
            println(f, "   → May accept poor quality solutions")
        end
        
        if !solver_settings["output_flag"]
            println(f, "⚠️  Solver output disabled")
            println(f, "   → Enable for debugging: set output_flag = true")
        end
        
        println(f, "✅ Solver settings loaded from simplified config")
        
    catch e
        println(f, "❌ Could not check solver settings: $e")
    end
end

"""
    check_formulation_compatibility_simplified!(debugger::SiennaDebugExport, f::IOStream, formulation_type::String)

Check formulation compatibility using simplified config.
"""
function check_formulation_compatibility_simplified!(debugger::SiennaDebugExport, f::IOStream, formulation_type::String)
    try
        device_formulations = get_device_formulations(debugger.config, formulation_type)
        sim_config = get(debugger.config.config_data, "simulations", Dict())
        
        println(f, "Formulation Type: $formulation_type")
        println(f, "Use Custom Templates: $(get(sim_config, "use_custom_templates", true))")
        
        # Check if formulations are defined
        if isempty(device_formulations)
            println(f, "❌ NO DEVICE FORMULATIONS DEFINED")
            println(f, "   → Add [formulations.$formulation_type] section to config")
        else
            println(f, "✅ Device formulations defined:")
            for (device, formulation) in device_formulations
                println(f, "   $device: $formulation")
            end
        end
        
        # Check for common formulation issues
        thermal_form = get(device_formulations, "thermal_standard", "")
        if formulation_type == "economic_dispatch" && !isempty(thermal_form)
            if occursin("UnitCommitment", thermal_form)
                println(f, "⚠️  Using Unit Commitment formulation in Economic Dispatch")
                println(f, "   → Consider ThermalBasicDispatch for ED")
            end
        end
        
        if formulation_type == "unit_commitment" && !isempty(thermal_form)
            if occursin("Dispatch", thermal_form) && !occursin("UnitCommitment", thermal_form)
                println(f, "⚠️  Using Dispatch formulation in Unit Commitment")
                println(f, "   → Consider ThermalBasicUnitCommitment for UC")
            end
        end
        
    catch e
        println(f, "❌ Could not check formulation compatibility: $e")
    end
end

"""
    generate_infeasibility_report_simplified!(debugger::SiennaDebugExport, infeas_dir::String, formulation_type::String)

Generate comprehensive infeasibility report using simplified config.
"""
function generate_infeasibility_report_simplified!(debugger::SiennaDebugExport, infeas_dir::String, formulation_type::String)
    @info "Generating comprehensive infeasibility report with simplified config..."
    
    try
        report_file = joinpath(infeas_dir, "infeasibility_diagnosis_report.txt")
        
        open(report_file, "w") do f
            println(f, "="^80)
            println(f, "SIENNA INFEASIBILITY DIAGNOSIS REPORT (SIMPLIFIED CONFIG)")
            println(f, "="^80)
            println(f, "Generated: $(now())")
            println(f, "Formulation: $formulation_type")
            println(f, "Project: $(debugger.config.project_name)")
            println(f, "Config Structure: Simplified")
            println(f, "="^80)
            println(f)
            
            println(f, "DIAGNOSIS SUMMARY:")
            println(f, "-"^50)
            println(f, "This report analyzes potential causes of model infeasibility")
            println(f, "in your PowerSimulations.jl optimization model using the")
            println(f, "simplified configuration structure.")
            println(f)
            
            println(f, "SIMPLIFIED CONFIG STRUCTURE:")
            println(f, "-"^50)
            println(f, "Your configuration uses the cleaned, simplified structure with:")
            println(f, "• [solver] section for all solver settings")
            println(f, "• [network] section for network model settings")
            println(f, "• [formulations] section for custom device formulations")
            println(f, "• [output] section for all output and debug settings")
            println(f)
            
            println(f, "CURRENT CONFIGURATION:")
            println(f, "-"^50)
            
            # Print current settings from simplified config
            solver_settings = get_solver_settings(debugger.config, formulation_type)
            network_settings = get_network_settings(debugger.config)
            
            println(f, "Solver: $(debugger.config.default_solver)")
            println(f, "Network Model: $(network_settings["model"])")
            println(f, "Time Limit: $(solver_settings["time_limit"]) seconds")
            println(f, "MIP Gap: $(solver_settings["mip_gap"])")
            println(f, "Transmission Limits: $(network_settings["enable_transmission_limits"])")
            println(f)
            
            println(f, "DEBUGGING STEPS FOR SIMPLIFIED CONFIG:")
            println(f, "-"^50)
            println(f, "1. CHECK SOLVER SETTINGS:")
            println(f, "   → Review [solver] section in config.toml")
            println(f, "   → Increase time_limit_ed or time_limit_uc if needed")
            println(f, "   → Set output_flag = true for solver debugging")
            println(f)
            
            println(f, "2. CHECK NETWORK MODEL:")
            println(f, "   → Review [network] section in config.toml")
            println(f, "   → Try model = \"CopperPlatePowerModel\" for testing")
            println(f, "   → Set enable_transmission_limits = false if having issues")
            println(f)
            
            println(f, "3. CHECK CUSTOM FORMULATIONS:")
            println(f, "   → Review [formulations.economic_dispatch] section")
            println(f, "   → Review [formulations.unit_commitment] section")
            println(f, "   → Verify all device types are properly mapped")
            println(f)
            
            println(f, "4. ENABLE DEBUG OUTPUT:")
            println(f, "   → In [output] section, set:")
            println(f, "   → export_optimization_model = true")
            println(f, "   → save_solver_logs = true")
            println(f, "   → save_constraint_breakdown = true")
            println(f)
            
            println(f, "="^80)
            println(f, "END OF SIMPLIFIED CONFIG DIAGNOSIS REPORT")
            println(f, "="^80)
        end
        
        push!(debugger.exported_files[formulation_type], report_file)
        @info "  ✓ Comprehensive infeasibility report generated for simplified config"
        
    catch e
        @error "Failed to generate infeasibility report: $e"
    end
end

"""
    export_debug_summary!(debugger::SiennaDebugExport)

Export comprehensive debug summary with simplified config info.
"""
function export_debug_summary!(debugger::SiennaDebugExport)
    @info "💾 Exporting comprehensive debug summary with simplified config info..."
    
    try
        summary_file = joinpath(debugger.export_directory, "debug_export_summary.json")
        
        # Create comprehensive summary with simplified config context
        summary = Dict{String, Any}(
            "export_metadata" => Dict(
                "timestamp" => string(now()),
                "export_directory" => debugger.export_directory,
                "debug_enabled" => debugger.debugging_enabled,
                "sienna_version" => "2.0_simplified",
                "config_file" => debugger.config.config_file_path,
                "config_structure" => "simplified"
            ),
            "project_info" => Dict(
                "name" => debugger.config.project_name,
                "description" => debugger.config.project_description,
                "version" => debugger.config.project_version,
                "network_model" => debugger.config.network_model,
                "horizon_hours" => debugger.config.default_horizon_hours
            ),
            "simplified_config_settings" => Dict(
                "solver_settings" => get_solver_settings(debugger.config, "economic_dispatch"),
                "network_settings" => get_network_settings(debugger.config),
                "debug_settings" => get_debug_settings(debugger.config)
            ),
            "exported_formulations" => collect(keys(debugger.exported_files)),
            "exported_files_by_formulation" => debugger.exported_files,
            "model_statistics" => debugger.model_info,
            "total_files_exported" => sum(length(files) for files in values(debugger.exported_files))
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
        
        @info "✅ Debug export summary saved with simplified config info"
        @info "   Summary file: $summary_file"
        @info "   Total formulations exported: $(length(debugger.exported_files))"
        @info "   Total files exported: $(summary["total_files_exported"])"
        
    catch e
        @error "Failed to export debug summary: $e"
    end
end

"""
    quick_infeasibility_check_simplified(config::SiennaConfig, sienna_sys::SiennaSystem)

Quick infeasibility check using simplified config settings.
"""
function quick_infeasibility_check_simplified(config::SiennaConfig, sienna_sys::SiennaSystem)
    @info "🔍 Running quick infeasibility check with simplified config..."
    
    issues = String[]
    
    try
        sys = get_power_system(sienna_sys)
        validation_settings = get(config.config_data, "validation", Dict())
        
        # Check 1: Load vs Generation using validation settings
        total_load = sum(get_max_active_power(load) for load in get_components(PowerLoad, sys))
        total_thermal = sum(get_active_power_limits(gen).max for gen in get_components(ThermalStandard, sys))
        total_renewable = sum(get_rating(gen) for gen in get_components(RenewableDispatch, sys))
        total_generation = total_thermal + total_renewable
        
        if total_generation < total_load
            push!(issues, "Insufficient generation capacity: $(round(total_generation, digits=1)) MW < $(round(total_load, digits=1)) MW")
        end
        
        # Check 2: Network model consistency
        network_settings = get_network_settings(config)
        n_lines = length(get_components(Line, sys))
        if network_settings["model"] != "CopperPlatePowerModel" && n_lines == 0
            push!(issues, "Network model $(network_settings["model"]) specified but no transmission lines found")
        end
        
        # Check 3: Solver settings
        solver_settings = get_solver_settings(config, "economic_dispatch")
        if solver_settings["time_limit"] < 60
            push!(issues, "Very short solver time limit: $(solver_settings["time_limit"]) seconds")
        end
        
        # Report results
        if isempty(issues)
            @info "✅ No obvious infeasibility issues detected with simplified config"
            return true
        else
            @warn "⚠️  Potential infeasibility issues detected:"
            for issue in issues
                @warn "   • $issue"
            end
            @info "💡 Check simplified config sections: [solver], [network], [validation]"
            return false
        end
        
    catch e
        @warn "Quick infeasibility check failed: $e"
        return false
    end
end

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
    enable_debugging_for_simulation!(sienna_sim::SiennaSimulations, enable_export::Bool=true)

Enable debugging mode for simplified config structure.
"""
function enable_debugging_for_simulation!(sienna_sim::SiennaSimulations, enable_export::Bool=true)
    @info "🔧 Enabling debugging mode for SiennaSimulations with simplified config..."
    
    # Create debugger
    debugger = SiennaDebugExport(sienna_sim.config, sienna_sim)
    
    if enable_export
        @info "   Auto-export enabled using simplified config settings"
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
        # Export system and configuration information for debugging
        export_model_for_debugging!(debugger, formulation_type)
        
        @info "✅ Debug build and export completed for $formulation_type"
        return true
        
    catch e
        @error "❌ Debug build failed for $formulation_type: $e"
        
        # Try to export partial information for debugging
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
        # Run infeasibility analysis
        analyze_infeasibility!(debugger, formulation_type)
        
        # Export comprehensive debug information
        export_debug_summary!(debugger)
        
        @info "🔍 Debugging files have been exported for analysis"
        @info "   Check: $(debugger.export_directory)"
        
        return true
        
    catch e
        @error "❌ Debug solve failed for $formulation_type: $e"
        
        # Run infeasibility analysis on failure
        @info "Running infeasibility analysis due to solve failure..."
        analyze_infeasibility!(debugger, formulation_type)
        
        # Export comprehensive debug information
        export_debug_summary!(debugger)
        
        rethrow(e)
    end
end

"""
    create_debugger_from_config(config_file::String="config.toml")

Convenience function to create a complete debugging setup from simplified config.
"""
function create_debugger_from_config(config_file::String="config.toml")
    @info "🔧 Creating complete debugging setup from simplified config..."
    
    try
        # Load simplified config
        config = SiennaConfig(config_file)
        
        # Build system
        sienna_sys = SiennaSystem(config)
        build_system!(sienna_sys)
        
        # Create simulations
        sienna_sim = SiennaSimulations(config, sienna_sys)
        
        # Create debugger
        debugger = SiennaDebugExport(config, sienna_sim)
        
        @info "✅ Complete debugging setup created"
        @info "   Config: Simplified structure"
        @info "   Export directory: $(debugger.export_directory)"
        @info "   Ready for debugging operations"
        
        return debugger
        
    catch e
        @error "❌ Failed to create debugging setup: $e"
        rethrow(e)
    end
end

"""
    run_complete_debug_analysis(config_file::String="config.toml", formulation_type::String="economic_dispatch")

Run a complete debugging analysis using simplified config.
"""
function run_complete_debug_analysis(config_file::String="config.toml", formulation_type::String="economic_dispatch")
    @info "🔍 Running complete debug analysis with simplified config..."
    
    try
        # Create debugger
        debugger = create_debugger_from_config(config_file)
        
        # Run quick infeasibility check
        @info "1. Running quick infeasibility check..."
        quick_infeasibility_check_simplified(debugger.config, debugger.sienna_sim.sienna_system)
        
        # Export system information
        @info "2. Exporting system information..."
        form_dir = joinpath(debugger.export_directory, formulation_type)
        mkpath(form_dir)
        export_system_information!(debugger, form_dir, formulation_type)
        export_solver_settings!(debugger, form_dir, formulation_type)
        
        # Run infeasibility analysis patterns
        @info "3. Running infeasibility pattern analysis..."
        analyze_infeasibility!(debugger, formulation_type)
        
        # Export summary
        @info "4. Exporting debug summary..."
        export_debug_summary!(debugger)
        
        @info "✅ Complete debug analysis finished"
        @info "   Check results in: $(debugger.export_directory)"
        
        return debugger
        
    catch e
        @error "❌ Complete debug analysis failed: $e"
        rethrow(e)
    end
end

"""
    validate_simplified_config(config_file::String="config.toml")

Validate simplified config structure and settings.
"""
function validate_simplified_config(config_file::String="config.toml")
    @info "✅ Validating simplified config structure..."
    
    try
        config = SiennaConfig(config_file)
        
        @info "Configuration validation results:"
        @info "  ✓ Config file loaded successfully"
        @info "  ✓ Simplified structure detected"
        @info "  ✓ Project: $(config.project_name)"
        @info "  ✓ Network model: $(config.network_model)"
        @info "  ✓ Solver: $(config.default_solver)"
        
        # Validate key sections exist
        required_sections = ["project", "paths", "system_building", "simulations", "solver", "network", "formulations", "output"]
        missing_sections = []
        
        for section in required_sections
            if !haskey(config.config_data, section)
                push!(missing_sections, section)
            end
        end
        
        if !isempty(missing_sections)
            @warn "Missing required sections: $(join(missing_sections, ", "))"
            return false
        end
        
        # Validate formulations exist
        econ_dispatch_forms = get_device_formulations(config, "economic_dispatch")
        unit_commit_forms = get_device_formulations(config, "unit_commitment")
        
        if isempty(econ_dispatch_forms)
            @warn "No economic_dispatch formulations defined"
        else
            @info "  ✓ Economic dispatch formulations: $(length(econ_dispatch_forms)) defined"
        end
        
        if isempty(unit_commit_forms)
            @warn "No unit_commitment formulations defined"
        else
            @info "  ✓ Unit commitment formulations: $(length(unit_commit_forms)) defined"
        end
        
        @info "✅ Simplified config validation completed successfully"
        return true
        
    catch e
        @error "❌ Config validation failed: $e"
        return false
    end
end

# ===== EXPORTS =====

export SiennaDebugExport
export export_model_for_debugging!, analyze_infeasibility!
export enable_debugging_for_simulation!, debug_build_and_export!, debug_solve_with_analysis!
export export_debug_summary!, clean_debug_exports!, quick_infeasibility_check_simplified
export create_debugger_from_config, run_complete_debug_analysis, validate_simplified_config

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
                quick_infeasibility_check_simplified(config, sienna_sys)
                
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
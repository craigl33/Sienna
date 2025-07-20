#!/usr/bin/env julia

"""
SiennaConfig.jl - Updated Configuration Management for Simplified Config Structure
=================================================================================

Updated to handle the cleaned config.toml structure while maintaining all functionality.
"""

using TOML
using Dates
using Logging

"""
    SiennaConfig

Main configuration class - updated for simplified config structure.
"""
mutable struct SiennaConfig
    # Core configuration data
    config_data::Dict{String, Any}
    config_file_path::String
    
    # Parsed and validated paths
    data_directory::String
    output_directory::String
    export_directory::String
    
    # System building parameters
    base_power::Float64
    default_voltage::Float64
    validate_system::Bool
    load_timeseries::Bool
    
    # Project metadata
    project_name::String
    project_description::String
    project_version::String
    project_author::String
    
    # Simulation parameters
    default_horizon_hours::Int
    default_solver::String
    network_model::String
    enable_transmission_limits::Bool
    
    # Fuel mappings
    thermal_fuels::Set{String}
    renewable_fuels::Set{String}
    
    # Default cost parameters
    thermal_defaults::Dict{String, Float64}
    renewable_defaults::Dict{String, Float64}
    
    # Output settings
    save_all_variables::Bool
    save_csv_variables::Bool
    create_timestamped_folders::Bool
    
    # Validation state
    is_validated::Bool
    validation_warnings::Vector{String}
    validation_errors::Vector{String}
end

"""
    SiennaConfig(config_file::String="config.toml")

Constructor - Load and parse simplified configuration.
"""
function SiennaConfig(config_file::String="config.toml")
    @info "🔧 Initializing SiennaConfig from: $config_file"
    
    if !isfile(config_file)
        @error "Configuration file not found: $config_file"
        @info "Creating default configuration file..."
        create_default_config(config_file)
        error("❌ Configuration file created. Please edit $config_file and run again.")
    end
    
    config_data = try
        TOML.parsefile(config_file)
    catch e
        error("❌ Failed to parse configuration file: $e")
    end
    
    config = SiennaConfig(
        config_data,
        abspath(config_file),
        "", "", "",  # paths
        100.0, 400.0, true, false,  # system defaults
        "", "", "", "",  # project metadata
        24, "HiGHS", "DCPPowerModel", false,  # simulation defaults
        Set{String}(), Set{String}(),  # fuel sets
        Dict{String, Float64}(), Dict{String, Float64}(),  # cost defaults
        true, true, true,  # output defaults
        false, String[], String[]  # validation state
    )
    
    parse_and_validate!(config)
    
    @info "✅ SiennaConfig initialized successfully"
    @info "   Project: $(config.project_name)"
    @info "   Network Model: $(config.network_model)"
    @info "   Data Directory: $(config.data_directory)"
    
    return config
end

"""
    parse_and_validate!(config::SiennaConfig)

Parse simplified configuration structure.
"""
function parse_and_validate!(config::SiennaConfig)
    @info "📋 Parsing simplified configuration..."
    
    config.validation_warnings = String[]
    config.validation_errors = String[]
    
    try
        parse_paths!(config)
        parse_project_metadata!(config)
        parse_system_building!(config)
        parse_simulation_parameters!(config)
        parse_fuel_mappings!(config)
        parse_cost_defaults!(config)
        parse_output_settings!(config)
        validate_configuration!(config)
        
        config.is_validated = true
        
        if !isempty(config.validation_warnings)
            @warn "⚠️  Configuration loaded with $(length(config.validation_warnings)) warnings"
            for warning in config.validation_warnings
                @warn "   $warning"
            end
        end
        
        if !isempty(config.validation_errors)
            error("❌ Configuration validation failed with $(length(config.validation_errors)) errors:\n" *
                  join(config.validation_errors, "\n"))
        end
        
        @info "✅ Configuration parsed and validated successfully"
        
    catch e
        config.is_validated = false
        @error "❌ Configuration parsing failed: $e"
        rethrow(e)
    end
end

"""
    parse_paths!(config::SiennaConfig)

Parse paths from simplified structure.
"""
function parse_paths!(config::SiennaConfig)
    paths_section = get(config.config_data, "paths", Dict())
    
    config.data_directory = expanduser(get(paths_section, "data_directory", "./data"))
    config.output_directory = expanduser(get(paths_section, "output_directory", "./sienna_results"))
    config.export_directory = expanduser(get(paths_section, "export_directory", "./sienna_psb_cases"))
    
    if !isdir(config.data_directory)
        push!(config.validation_errors, "Data directory not found: $(config.data_directory)")
    end
    
    for dir in [config.output_directory, config.export_directory]
        if !isdir(dir)
            try
                mkpath(dir)
                @info "📁 Created directory: $dir"
            catch e
                push!(config.validation_warnings, "Could not create directory $dir: $e")
            end
        end
    end
end

"""
    parse_project_metadata!(config::SiennaConfig)

Parse project metadata.
"""
function parse_project_metadata!(config::SiennaConfig)
    project_section = get(config.config_data, "project", Dict())
    
    config.project_name = get(project_section, "name", "Sienna Power System Simulation")
    config.project_description = get(project_section, "description", "Power system simulation using Sienna ecosystem")
    config.project_version = get(project_section, "version", "1.0.0")
    config.project_author = get(project_section, "author", "Sienna User")
end

"""
    parse_system_building!(config::SiennaConfig)

Parse system building parameters from simplified structure.
"""
function parse_system_building!(config::SiennaConfig)
    system_section = get(config.config_data, "system_building", Dict())
    
    config.base_power = Float64(get(system_section, "base_power", 100.0))
    config.default_voltage = Float64(get(system_section, "default_voltage", 400.0))
    config.validate_system = get(system_section, "validate_system", true)
    config.load_timeseries = get(system_section, "load_timeseries", false)
    
    if config.base_power <= 0
        push!(config.validation_errors, "base_power must be positive, got: $(config.base_power)")
    end
    
    if config.default_voltage <= 0
        push!(config.validation_errors, "default_voltage must be positive, got: $(config.default_voltage)")
    end
end

"""
    parse_simulation_parameters!(config::SiennaConfig)

Parse simulation parameters from simplified structure.
"""
function parse_simulation_parameters!(config::SiennaConfig)
    sim_section = get(config.config_data, "simulations", Dict())
    
    config.default_horizon_hours = get(sim_section, "horizon_hours", 24)
    config.default_solver = get(get(config.config_data, "solver", Dict()), "name", "HiGHS")
    
    # Parse network settings from new [network] section
    network_section = get(config.config_data, "network", Dict())
    config.network_model = get(network_section, "model", "DCPPowerModel")
    config.enable_transmission_limits = get(network_section, "enable_transmission_limits", false)
    
    if config.default_horizon_hours <= 0
        push!(config.validation_errors, "horizon_hours must be positive, got: $(config.default_horizon_hours)")
    end
    
    if config.default_solver != "HiGHS"
        push!(config.validation_warnings, "Only HiGHS solver is currently supported, got: $(config.default_solver)")
    end
    
    # Validate network model
    valid_network_models = [
        "CopperPlatePowerModel", "PTDFPowerModel", 
        "AreaBalancePowerModel", "AreaPTDFPowerModel",
        "DCPPowerModel", "DCMPPowerModel", "LPACCPowerModel"
    ]
    
    if !(config.network_model in valid_network_models)
        push!(config.validation_errors, 
              "Invalid network model: $(config.network_model). " *
              "Valid options: $(join(valid_network_models, ", "))")
    end
end

"""
    parse_fuel_mappings!(config::SiennaConfig)

Parse fuel mappings from simplified structure.
"""
function parse_fuel_mappings!(config::SiennaConfig)
    system_section = get(config.config_data, "system_building", Dict())
    fuel_section = get(system_section, "fuel_mapping", Dict())
    
    thermal_fuels = get(fuel_section, "thermal_fuels", ["COAL", "NATURAL_GAS", "DIESEL", "NUCLEAR", "BIOMASS"])
    renewable_fuels = get(fuel_section, "renewable_fuels", ["WIND", "SOLAR", "HYDRO", "GEOTHERMAL"])
    
    config.thermal_fuels = Set(uppercase.(thermal_fuels))
    config.renewable_fuels = Set(uppercase.(renewable_fuels))
    
    overlap = intersect(config.thermal_fuels, config.renewable_fuels)
    if !isempty(overlap)
        push!(config.validation_warnings, 
              "Fuel types appear in both thermal and renewable lists: $(join(overlap, ", "))")
    end
end

"""
    parse_cost_defaults!(config::SiennaConfig)

Parse cost defaults from simplified structure.
"""
function parse_cost_defaults!(config::SiennaConfig)
    system_section = get(config.config_data, "system_building", Dict())
    defaults_section = get(system_section, "defaults", Dict())
    
    config.thermal_defaults = Dict{String, Float64}(
        "variable_cost" => get(defaults_section, "thermal_variable_cost", 50.0),
        "startup_cost" => get(defaults_section, "thermal_startup_cost", 1000.0),
        "shutdown_cost" => get(defaults_section, "thermal_shutdown_cost", 0.0),
        "ramp_rate" => get(defaults_section, "thermal_ramp_rate", 100.0),
        "min_power_fraction" => get(defaults_section, "thermal_min_power_fraction", 0.3)
    )
    
    config.renewable_defaults = Dict{String, Float64}(
        "variable_cost" => get(defaults_section, "renewable_variable_cost", 0.0),
        "min_power_fraction" => get(defaults_section, "renewable_min_power_fraction", 0.0)
    )
    
    # Validate cost parameters
    for (param, value) in config.thermal_defaults
        if value < 0 && param != "shutdown_cost"
            push!(config.validation_warnings, "Negative thermal $param: $value")
        end
    end
    
    for (param, value) in config.renewable_defaults
        if value < 0
            push!(config.validation_warnings, "Negative renewable $param: $value")
        end
    end
end

"""
    parse_output_settings!(config::SiennaConfig)

Parse output settings from simplified structure.
"""
function parse_output_settings!(config::SiennaConfig)
    output_section = get(config.config_data, "output", Dict())
    
    config.save_all_variables = get(output_section, "save_all_variables", true)
    config.save_csv_variables = get(output_section, "save_csv_variables", true)
    config.create_timestamped_folders = get(output_section, "create_timestamped_folders", true)
end

"""
    validate_configuration!(config::SiennaConfig)

Validate the parsed configuration.
"""
function validate_configuration!(config::SiennaConfig)
    required_files = ["bus.csv", "load.csv"]
    optional_files = ["gen.csv", "branch.csv", "dc_branch.csv"]
    
    for file in required_files
        file_path = joinpath(config.data_directory, file)
        if !isfile(file_path)
            push!(config.validation_errors, "Required file missing: $file")
        end
    end
    
    for file in optional_files
        file_path = joinpath(config.data_directory, file)
        if isfile(file_path)
            @info "   ✓ Found optional file: $file"
        end
    end
    
    if config.load_timeseries
        metadata_file = joinpath(config.data_directory, "timeseries_metadata.json")
        if !isfile(metadata_file)
            push!(config.validation_warnings, 
                  "Time series loading enabled but timeseries_metadata.json not found")
        end
    end
end

# ===== UPDATED PUBLIC INTERFACE METHODS =====

"""
    get_solver_settings(config::SiennaConfig, formulation_type::String="ed")

Get solver settings from simplified [solver] section.
"""
function get_solver_settings(config::SiennaConfig, formulation_type::String="ed")
    ensure_validated(config)
    
    solver_section = get(config.config_data, "solver", Dict())
    
    if formulation_type == "ed"
        return Dict{String, Any}(
            "time_limit" => get(solver_section, "time_limit_ed", 300.0),
            "mip_gap" => get(solver_section, "mip_gap_ed", 0.01),
            "threads" => get(solver_section, "threads", 0),
            "output_flag" => get(solver_section, "output_flag", false),
            "presolve" => get(solver_section, "presolve", "on"),
            "parallel" => get(solver_section, "parallel", "on")
        )
    elseif formulation_type == "uc"
        return Dict{String, Any}(
            "time_limit" => get(solver_section, "time_limit_uc", 600.0),
            "mip_gap" => get(solver_section, "mip_gap_uc", 0.02),
            "threads" => get(solver_section, "threads", 0),
            "output_flag" => get(solver_section, "output_flag", false),
            "presolve" => get(solver_section, "presolve", "on"),
            "parallel" => get(solver_section, "parallel", "on")
        )
    else
        error("❌ Unknown formulation type: $formulation_type")
    end
end

"""
    get_device_formulations(config::SiennaConfig, formulation_type::String)

Get device formulations from simplified [formulations] section.
"""
function get_device_formulations(config::SiennaConfig, formulation_type::String)
    ensure_validated(config)
    
    formulations_section = get(config.config_data, "formulations", Dict())
    return get(formulations_section, formulation_type, Dict())
end

"""
    should_run_formulation(config::SiennaConfig, formulation_type::String)

Check if formulation should run based on simplified config.
"""
function should_run_formulation(config::SiennaConfig, formulation_type::String)
    ensure_validated(config)
    
    sim_section = get(config.config_data, "simulations", Dict())
    
    if formulation_type == "ed"
        return get(sim_section, "run_economic_dispatch", true)
    elseif formulation_type == "uc"
        return get(sim_section, "run_unit_commitment", true)
    else
        return false
    end
end

"""
    get_network_settings(config::SiennaConfig)

Get network settings from simplified [network] section.
"""
function get_network_settings(config::SiennaConfig)
    ensure_validated(config)
    
    network_section = get(config.config_data, "network", Dict())
    return Dict{String, Any}(
        "model" => get(network_section, "model", "DCPPowerModel"),
        "enable_transmission_limits" => get(network_section, "enable_transmission_limits", false),
        "use_slacks" => get(network_section, "use_slacks", true),
        "add_slack_variables" => get(network_section, "add_slack_variables", true),
        "slack_penalty" => get(network_section, "slack_penalty", 1000000.0)
    )
end

"""
    get_debug_settings(config::SiennaConfig)

Get debug settings from simplified [output] section.
"""
function get_debug_settings(config::SiennaConfig)
    ensure_validated(config)
    
    output_section = get(config.config_data, "output", Dict())
    return Dict{String, Any}(
        "export_optimization_model" => get(output_section, "export_optimization_model", true),
        "save_solver_logs" => get(output_section, "save_solver_logs", true),
        "save_constraint_breakdown" => get(output_section, "save_constraint_breakdown", true),
        "save_variable_breakdown" => get(output_section, "save_variable_breakdown", true)
    )
end

# ===== UTILITY METHODS (unchanged) =====

function ensure_validated(config::SiennaConfig)
    if !config.is_validated
        error("❌ Configuration has not been validated. Call parse_and_validate! first.")
    end
end

function get_data_directory(config::SiennaConfig)
    ensure_validated(config)
    return config.data_directory
end

function get_output_directory(config::SiennaConfig)
    ensure_validated(config)
    return config.output_directory
end

function get_export_directory(config::SiennaConfig)
    ensure_validated(config)
    return config.export_directory
end

function is_fuel_thermal(config::SiennaConfig, fuel_type::String)
    ensure_validated(config)
    return uppercase(fuel_type) in config.thermal_fuels
end

function is_fuel_renewable(config::SiennaConfig, fuel_type::String)
    ensure_validated(config)
    return uppercase(fuel_type) in config.renewable_fuels
end

function get_thermal_defaults(config::SiennaConfig)
    ensure_validated(config)
    return copy(config.thermal_defaults)
end

function get_renewable_defaults(config::SiennaConfig)
    ensure_validated(config)
    return copy(config.renewable_defaults)
end

function show_summary(config::SiennaConfig)
    ensure_validated(config)
    
    @info "\n" * "="^60
    @info "SIENNA CONFIGURATION SUMMARY (SIMPLIFIED)"
    @info "="^60
    @info "Project: $(config.project_name)"
    @info "Version: $(config.project_version)"
    @info "Author: $(config.project_author)"
    @info ""
    @info "Paths:"
    @info "  Data: $(config.data_directory)"
    @info "  Output: $(config.output_directory)"
    @info "  Export: $(config.export_directory)"
    @info ""
    @info "System Parameters:"
    @info "  Base Power: $(config.base_power) MW"
    @info "  Default Voltage: $(config.default_voltage) kV"
    @info "  Load Time Series: $(config.load_timeseries)"
    @info ""
    @info "Simulation Settings:"
    @info "  Horizon: $(config.default_horizon_hours) hours"
    @info "  Solver: $(config.default_solver)"
    @info "  Network Model: $(config.network_model)"
    @info "  Transmission Limits: $(config.enable_transmission_limits)"
    @info ""
    @info "Fuel Classifications:"
    @info "  Thermal: $(join(sort(collect(config.thermal_fuels)), ", "))"
    @info "  Renewable: $(join(sort(collect(config.renewable_fuels)), ", "))"
    @info ""
    @info "Output Settings:"
    @info "  Save All Variables: $(config.save_all_variables)"
    @info "  Save CSV Variables: $(config.save_csv_variables)"
    @info "  Timestamped Folders: $(config.create_timestamped_folders)"
    @info "="^60
    
    if !isempty(config.validation_warnings)
        @info "\n⚠️  Validation Warnings:"
        for warning in config.validation_warnings
            @info "  • $warning"
        end
    end
end

function reload!(config::SiennaConfig)
    @info "🔄 Reloading configuration from: $(config.config_file_path)"
    config.config_data = TOML.parsefile(config.config_file_path)
    parse_and_validate!(config)
    @info "✅ Configuration reloaded successfully"
end

function export_config(config::SiennaConfig, output_file::String)
    ensure_validated(config)
    
    @info "💾 Exporting configuration to: $output_file"
    
    export_data = deepcopy(config.config_data)
    
    if !haskey(export_data, "export_metadata")
        export_data["export_metadata"] = Dict()
    end
    
    export_data["export_metadata"]["exported_at"] = string(now())
    export_data["export_metadata"]["exported_from"] = config.config_file_path
    export_data["export_metadata"]["sienna_config_version"] = "2.0_simplified"
    
    open(output_file, "w") do f
        TOML.print(f, export_data)
    end
    
    @info "✅ Configuration exported successfully"
end

"""
    create_default_config(config_file::String)

Create simplified default configuration file.
"""
function create_default_config(config_file::String)
    @info "📝 Creating simplified default configuration: $config_file"
    
    default_config = """
# Sienna Power System Configuration v2.0 - Simplified
# Edit this file for your specific project

[project]
name = "Sienna Power System Simulation"
description = "Power system simulation using Sienna ecosystem"
version = "1.0.0"
author = "Sienna User"
created_date = "$(Dates.format(now(), "yyyy-mm-dd"))"

[paths]
data_directory = "./data"
output_directory = "./sienna_results"
export_directory = "./sienna_psb_cases"

[system_building]
base_power = 100.0
default_voltage = 400.0
validate_system = true
load_timeseries = false

[system_building.areas]
area_strategy = "per_bus"
validate_area_assignments = true

[system_building.fuel_mapping]
thermal_fuels = ["COAL", "NATURAL_GAS", "DIESEL", "NUCLEAR", "BIOMASS"]
renewable_fuels = ["WIND", "SOLAR", "HYDRO", "GEOTHERMAL"]

[system_building.defaults]
thermal_variable_cost = 50.0
thermal_startup_cost = 1000.0
renewable_variable_cost = 0.0
thermal_ramp_rate = 100.0

[simulations]
simulation_type = "ST"
horizon_hours = 48
interval_hours = 24
annual_simulation_days = 365
run_economic_dispatch = true
run_unit_commitment = false

[formulations.economic_dispatch]
thermal_standard = "ThermalBasicDispatch"
renewable_dispatch = "RenewableFullDispatch"
renewable_nondispatch = "FixedOutput"
power_load = "StaticPowerLoad"
line = "StaticBranch"
dc_line = "HVDCTwoTerminalDispatch"

[formulations.unit_commitment]
thermal_standard = "ThermalBasicUnitCommitment"
renewable_dispatch = "RenewableFullDispatch"
renewable_nondispatch = "FixedOutput"
power_load = "StaticPowerLoad"
line = "StaticBranch"
dc_line = "HVDCTwoTerminalDispatch"

[network]
model = "DCPPowerModel"
enable_transmission_limits = true

[solver]
name = "HiGHS"
time_limit_ed = 300.0
time_limit_uc = 600.0
mip_gap_ed = 0.01
mip_gap_uc = 0.02
threads = 0
output_flag = false

[output]
save_all_variables = true
save_csv_variables = true
create_timestamped_folders = true

[validation]
check_power_balance = true
max_generation_capacity = 100000.0
max_load = 50000.0
min_reserve_margin = 0.15
"""
    
    open(config_file, "w") do f
        write(f, default_config)
    end
    
    @info "✅ Simplified default configuration created: $config_file"
end

# ===== EXPORTS =====

export SiennaConfig
export get_data_directory, get_output_directory, get_export_directory
export get_solver_settings, get_device_formulations, should_run_formulation
export get_network_settings, get_debug_settings
export is_fuel_thermal, is_fuel_renewable
export get_thermal_defaults, get_renewable_defaults
export show_summary, reload!, export_config
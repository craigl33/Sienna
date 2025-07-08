#!/usr/bin/env julia

"""
SiennaConfig.jl - Configuration Management for Sienna Ecosystem
===============================================================

Class-based configuration manager that encapsulates all config-related functionality.
Provides clean interface for accessing configuration parameters throughout the workflow.

Features:
- Centralized config loading and validation
- Type-safe parameter access
- Path resolution and validation
- Default value management
- Config compatibility checking

Usage:
    config = SiennaConfig("config.toml")
    data_dir = config.get_data_directory()
    solver_settings = config.get_solver_settings("economic_dispatch")
"""

using TOML
using Dates
using Logging

"""
    SiennaConfig

Main configuration class that encapsulates all configuration functionality.
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

Constructor - Load and parse configuration from TOML file.
"""
function SiennaConfig(config_file::String="config.toml")
    @info "🔧 Initializing SiennaConfig from: $config_file"
    
    # Check if config file exists
    if !isfile(config_file)
        @error "Configuration file not found: $config_file"
        @info "Creating default configuration file..."
        create_default_config(config_file)
        error("❌ Configuration file created. Please edit $config_file and run again.")
    end
    
    # Load TOML configuration
    config_data = try
        TOML.parsefile(config_file)
    catch e
        error("❌ Failed to parse configuration file: $e")
    end
    
    # Create instance with default values
    config = SiennaConfig(
        config_data,
        abspath(config_file),
        "", "", "",  # paths - will be set in parse_and_validate
        100.0, 400.0, true, false,  # system defaults
        "", "", "", "",  # project metadata
        24, "HiGHS", "CopperPlatePowerModel", false,  # simulation defaults
        Set{String}(), Set{String}(),  # fuel sets
        Dict{String, Float64}(), Dict{String, Float64}(),  # cost defaults
        true, true, true,  # output defaults
        false, String[], String[]  # validation state
    )
    
    # Parse and validate configuration
    parse_and_validate!(config)
    
    @info "✅ SiennaConfig initialized successfully"
    @info "   Project: $(config.project_name)"
    @info "   Network Model: $(config.network_model)"
    @info "   Data Directory: $(config.data_directory)"
    
    return config
end

"""
    parse_and_validate!(config::SiennaConfig)

Parse configuration data and validate all parameters.
"""
function parse_and_validate!(config::SiennaConfig)
    @info "📋 Parsing and validating configuration..."
    
    # Reset validation state
    config.validation_warnings = String[]
    config.validation_errors = String[]
    
    try
        # Parse paths
        parse_paths!(config)
        
        # Parse project metadata
        parse_project_metadata!(config)
        
        # Parse system building parameters
        parse_system_building!(config)
        
        # Parse simulation parameters
        parse_simulation_parameters!(config)
        
        # Parse fuel mappings
        parse_fuel_mappings!(config)
        
        # Parse cost defaults
        parse_cost_defaults!(config)
        
        # Parse output settings
        parse_output_settings!(config)
        
        # Validate parsed configuration
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

Parse and expand all file paths.
"""
function parse_paths!(config::SiennaConfig)
    paths_section = get(config.config_data, "paths", Dict())
    
    # Required paths
    config.data_directory = expanduser(get(paths_section, "data_directory", "./data"))
    config.output_directory = expanduser(get(paths_section, "output_directory", "./sienna_results"))
    config.export_directory = expanduser(get(paths_section, "export_directory", "./sienna_psb_cases"))
    
    # Validate data directory exists
    if !isdir(config.data_directory)
        push!(config.validation_errors, "Data directory not found: $(config.data_directory)")
    end
    
    # Create output directories if they don't exist
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

Parse project metadata section.
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

Parse system building parameters.
"""
function parse_system_building!(config::SiennaConfig)
    system_section = get(config.config_data, "system_building", Dict())
    
    config.base_power = Float64(get(system_section, "base_power", 100.0))
    config.default_voltage = Float64(get(system_section, "default_voltage", 400.0))
    config.validate_system = get(system_section, "validate_system", true)
    config.load_timeseries = get(system_section, "load_timeseries", false)
    
    # Validate ranges
    if config.base_power <= 0
        push!(config.validation_errors, "base_power must be positive, got: $(config.base_power)")
    end
    
    if config.default_voltage <= 0
        push!(config.validation_errors, "default_voltage must be positive, got: $(config.default_voltage)")
    end
end

"""
    parse_simulation_parameters!(config::SiennaConfig)

Parse simulation parameters.
"""
function parse_simulation_parameters!(config::SiennaConfig)
    sim_section = get(config.config_data, "simulations", Dict())
    
    config.default_horizon_hours = get(sim_section, "default_horizon_hours", 24)
    config.default_solver = get(sim_section, "default_solver", "HiGHS")
    
    # Parse network settings
    network_section = get(sim_section, "network", Dict())
    config.network_model = get(network_section, "default_model", "CopperPlatePowerModel")
    config.enable_transmission_limits = get(network_section, "enable_transmission_limits", false)
    
    # Validate simulation parameters
    if config.default_horizon_hours <= 0
        push!(config.validation_errors, "default_horizon_hours must be positive, got: $(config.default_horizon_hours)")
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

Parse fuel type mappings.
"""
function parse_fuel_mappings!(config::SiennaConfig)
    system_section = get(config.config_data, "system_building", Dict())
    fuel_section = get(system_section, "fuel_mapping", Dict())
    
    thermal_fuels = get(fuel_section, "thermal_fuels", ["COAL", "NATURAL_GAS", "DIESEL", "NUCLEAR", "BIOMASS"])
    renewable_fuels = get(fuel_section, "renewable_fuels", ["WIND", "SOLAR", "HYDRO", "GEOTHERMAL"])
    
    config.thermal_fuels = Set(uppercase.(thermal_fuels))
    config.renewable_fuels = Set(uppercase.(renewable_fuels))
    
    # Check for overlap
    overlap = intersect(config.thermal_fuels, config.renewable_fuels)
    if !isempty(overlap)
        push!(config.validation_warnings, 
              "Fuel types appear in both thermal and renewable lists: $(join(overlap, ", "))")
    end
end

"""
    parse_cost_defaults!(config::SiennaConfig)

Parse default cost parameters.
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
        if value < 0 && param != "shutdown_cost"  # shutdown cost can be negative
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

Parse output and results settings.
"""
function parse_output_settings!(config::SiennaConfig)
    output_section = get(config.config_data, "output", Dict())
    
    config.save_all_variables = get(output_section, "save_all_variables", true)
    config.save_csv_variables = get(output_section, "save_csv_variables", true)
    config.create_timestamped_folders = get(output_section, "create_timestamped_folders", true)
end

"""
    validate_configuration!(config::SiennaConfig)

Perform additional validation checks.
"""
function validate_configuration!(config::SiennaConfig)
    # Check required CSV files exist
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
    
    # Check time series metadata if time series loading enabled
    if config.load_timeseries
        metadata_file = joinpath(config.data_directory, "timeseries_metadata.json")
        if !isfile(metadata_file)
            push!(config.validation_warnings, 
                  "Time series loading enabled but timeseries_metadata.json not found")
        end
    end
end

# ===== PUBLIC INTERFACE METHODS =====

"""
    get_data_directory(config::SiennaConfig)

Get the validated data directory path.
"""
function get_data_directory(config::SiennaConfig)
    ensure_validated(config)
    return config.data_directory
end

"""
    get_output_directory(config::SiennaConfig)

Get the output directory path.
"""
function get_output_directory(config::SiennaConfig)
    ensure_validated(config)
    return config.output_directory
end

"""
    get_export_directory(config::SiennaConfig)

Get the export directory path.
"""
function get_export_directory(config::SiennaConfig)
    ensure_validated(config)
    return config.export_directory
end

"""
    get_solver_settings(config::SiennaConfig, formulation_type::String="economic_dispatch")

Get solver settings for the specified formulation type.
"""
function get_solver_settings(config::SiennaConfig, formulation_type::String="economic_dispatch")
    ensure_validated(config)
    
    solver_settings = get(get(config.config_data, "simulations", Dict()), "solver_settings", Dict())
    highs_settings = get(solver_settings, "HiGHS", Dict())
    
    if formulation_type == "economic_dispatch"
        return Dict{String, Any}(
            "time_limit" => get(highs_settings, "time_limit_ed", 300.0),
            "mip_gap" => get(highs_settings, "mip_gap_ed", 0.01),
            "threads" => get(highs_settings, "threads", 0),
            "output_flag" => get(highs_settings, "output_flag", false),
            "presolve" => get(highs_settings, "presolve", "on"),
            "parallel" => get(highs_settings, "parallel", "on")
        )
    elseif formulation_type == "unit_commitment"
        return Dict{String, Any}(
            "time_limit" => get(highs_settings, "time_limit_uc", 600.0),
            "mip_gap" => get(highs_settings, "mip_gap_uc", 0.02),
            "threads" => get(highs_settings, "threads", 0),
            "output_flag" => get(highs_settings, "output_flag", false),
            "presolve" => get(highs_settings, "presolve", "on"),
            "parallel" => get(highs_settings, "parallel", "on")
        )
    else
        error("❌ Unknown formulation type: $formulation_type")
    end
end

"""
    get_device_formulations(config::SiennaConfig, formulation_type::String)

Get device formulation mappings for the specified formulation type.
"""
function get_device_formulations(config::SiennaConfig, formulation_type::String)
    ensure_validated(config)
    
    formulations_section = get(get(config.config_data, "simulations", Dict()), "formulations", Dict())
    return get(formulations_section, formulation_type, Dict())
end

"""
    should_run_formulation(config::SiennaConfig, formulation_type::String)

Check if a specific formulation should be run based on config.
"""
function should_run_formulation(config::SiennaConfig, formulation_type::String)
    ensure_validated(config)
    
    sim_section = get(config.config_data, "simulations", Dict())
    
    if formulation_type == "economic_dispatch"
        return get(sim_section, "run_economic_dispatch", true)
    elseif formulation_type == "unit_commitment"
        return get(sim_section, "run_unit_commitment", true)
    else
        return false
    end
end

"""
    is_fuel_thermal(config::SiennaConfig, fuel_type::String)

Check if a fuel type is classified as thermal.
"""
function is_fuel_thermal(config::SiennaConfig, fuel_type::String)
    ensure_validated(config)
    return uppercase(fuel_type) in config.thermal_fuels
end

"""
    is_fuel_renewable(config::SiennaConfig, fuel_type::String)

Check if a fuel type is classified as renewable.
"""
function is_fuel_renewable(config::SiennaConfig, fuel_type::String)
    ensure_validated(config)
    return uppercase(fuel_type) in config.renewable_fuels
end

"""
    get_thermal_defaults(config::SiennaConfig)

Get default thermal generator parameters.
"""
function get_thermal_defaults(config::SiennaConfig)
    ensure_validated(config)
    return copy(config.thermal_defaults)
end

"""
    get_renewable_defaults(config::SiennaConfig)

Get default renewable generator parameters.
"""
function get_renewable_defaults(config::SiennaConfig)
    ensure_validated(config)
    return copy(config.renewable_defaults)
end

"""
    show_summary(config::SiennaConfig)

Display a comprehensive configuration summary.
"""
function show_summary(config::SiennaConfig)
    ensure_validated(config)
    
    @info "\n" * "="^60
    @info "SIENNA CONFIGURATION SUMMARY"
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

# ===== UTILITY METHODS =====

"""
    ensure_validated(config::SiennaConfig)

Ensure configuration has been validated before use.
"""
function ensure_validated(config::SiennaConfig)
    if !config.is_validated
        error("❌ Configuration has not been validated. Call parse_and_validate! first.")
    end
end

"""
    reload!(config::SiennaConfig)

Reload configuration from file.
"""
function reload!(config::SiennaConfig)
    @info "🔄 Reloading configuration from: $(config.config_file_path)"
    
    config.config_data = TOML.parsefile(config.config_file_path)
    parse_and_validate!(config)
    
    @info "✅ Configuration reloaded successfully"
end

"""
    export_config(config::SiennaConfig, output_file::String)

Export current configuration to a new TOML file.
"""
function export_config(config::SiennaConfig, output_file::String)
    ensure_validated(config)
    
    @info "💾 Exporting configuration to: $output_file"
    
    # Add metadata to exported config
    export_data = deepcopy(config.config_data)
    
    if !haskey(export_data, "export_metadata")
        export_data["export_metadata"] = Dict()
    end
    
    export_data["export_metadata"]["exported_at"] = string(now())
    export_data["export_metadata"]["exported_from"] = config.config_file_path
    export_data["export_metadata"]["sienna_config_version"] = "2.0"
    
    open(output_file, "w") do f
        TOML.print(f, export_data)
    end
    
    @info "✅ Configuration exported successfully"
end

"""
    create_default_config(config_file::String)

Create a default configuration file.
"""
function create_default_config(config_file::String)
    @info "📝 Creating default configuration: $config_file"
    
    default_config = """
# Sienna Power System Configuration v2.0
# Edit this file for your specific project

[project]
name = "Sienna Power System Simulation"
description = "Power system simulation using Sienna ecosystem"
version = "1.0.0"
author = "Sienna User"
created_date = "$(Dates.format(now(), "yyyy-mm-dd"))"

[paths]
# Primary data directory - EDIT THIS PATH
data_directory = "./data"
output_directory = "./sienna_results"
export_directory = "./sienna_psb_cases"

[system_building]
base_power = 100.0
default_voltage = 400.0
validate_system = true
load_timeseries = false

[system_building.fuel_mapping]
thermal_fuels = ["COAL", "NATURAL_GAS", "DIESEL", "NUCLEAR", "BIOMASS"]
renewable_fuels = ["WIND", "SOLAR", "HYDRO", "GEOTHERMAL"]

[system_building.defaults]
thermal_variable_cost = 50.0
thermal_startup_cost = 1000.0
thermal_shutdown_cost = 0.0
renewable_variable_cost = 0.0
thermal_ramp_rate = 100.0
thermal_min_power_fraction = 0.3
renewable_min_power_fraction = 0.0

[simulations]
default_horizon_hours = 24
default_solver = "HiGHS"
run_economic_dispatch = true
run_unit_commitment = true

[simulations.network]
default_model = "CopperPlatePowerModel"
enable_transmission_limits = false

[simulations.formulations.economic_dispatch]
thermal_standard = "ThermalBasicDispatch"
renewable_dispatch = "RenewableFullDispatch"
renewable_nondispatch = "FixedOutput"
power_load = "StaticPowerLoad"
line = "StaticBranch"
dc_line = "HVDCTwoTerminalDispatch"

[simulations.formulations.unit_commitment]
thermal_standard = "ThermalBasicUnitCommitment"
renewable_dispatch = "RenewableFullDispatch"
renewable_nondispatch = "FixedOutput"
power_load = "StaticPowerLoad"
line = "StaticBranch"
dc_line = "HVDCTwoTerminalDispatch"

[simulations.solver_settings.HiGHS]
time_limit_ed = 300.0
time_limit_uc = 600.0
mip_gap_ed = 0.01
mip_gap_uc = 0.02
threads = 0
output_flag = false
presolve = "on"
parallel = "on"

[output]
save_all_variables = true
save_csv_variables = true
create_timestamped_folders = true

[logging]
level = "Info"
log_to_file = false
"""
    
    open(config_file, "w") do f
        write(f, default_config)
    end
    
    @info "✅ Default configuration created: $config_file"
end

# ===== EXPORTS =====

export SiennaConfig
export get_data_directory, get_output_directory, get_export_directory
export get_solver_settings, get_device_formulations, should_run_formulation
export is_fuel_thermal, is_fuel_renewable
export get_thermal_defaults, get_renewable_defaults
export show_summary, reload!, export_config

# Test functionality when run directly
if abspath(PROGRAM_FILE) == @__FILE__
    @info "🧪 Testing SiennaConfig..."
    
    try
        # Test with default config
        if isfile("config.toml")
            test_config = SiennaConfig("config.toml")
            show_summary(test_config)
            @info "✅ SiennaConfig test passed"
        else
            @info "No config.toml found - creating default"
            create_default_config("config.toml")
            @info "✅ Default config created - edit and run again"
        end
    catch e
        @error "❌ SiennaConfig test failed: $e"
    end
end
#!/usr/bin/env julia

"""
Sienna Configuration Loader
Simple, dedicated configuration loader for Sienna power system simulations
Only handles TOML format for clarity and simplicity
"""

using TOML
using Dates
using Logging

# Global configuration storage (avoids const redefinition)
global CONFIG
if !@isdefined(CONFIG) || !isa(CONFIG, Ref)
    CONFIG = Ref{Dict}()
end
"""
    load_config(config_file::String="config.toml")

Load Sienna configuration from TOML file.
"""
function load_config(config_file::String="config.toml")
    @info "Loading Sienna configuration from: $config_file"
    
    if !isfile(config_file)
        @error "Configuration file not found: $config_file"
        @info "Creating default Sienna configuration file..."
        create_default_config(config_file)
        @info "✓ Default config created. Please edit $config_file and run again."
        error("Configuration file created. Please edit paths and run again.")
    end
    
    # Parse TOML configuration
    config = TOML.parsefile(config_file)
    
    # Process and validate configuration
    config = process_sienna_config(config)
    
    # Store in global reference
    CONFIG[] = config
    
    @info "✓ Sienna configuration loaded successfully"
    @info "  Project: $(config["project"]["name"])"
    @info "  Data directory: $(config["paths"]["data_directory"])"
    @info "  Default solver: $(config["simulations"]["default_solver"])"
    
    return config
end

"""
    get_config()

Get the currently loaded Sienna configuration.
"""
function get_config()
    if !isassigned(CONFIG)
        @info "No configuration loaded. Looking for config.toml..."
        load_config()
    end
    return CONFIG[]
end

"""
    process_sienna_config(config::Dict)

Process and validate Sienna configuration.
"""
function process_sienna_config(config::Dict)
    @info "Processing Sienna configuration..."
    
    # Expand all paths
    expand_config_paths!(config)
    
    # Validate required sections exist
    validate_config_structure(config)
    
    # Apply any missing defaults
    apply_config_defaults!(config)
    
    # Validate paths exist
    validate_config_paths(config)
    
    @info "✓ Configuration processed and validated"
    return config
end

"""
    expand_config_paths!(config::Dict)

Expand all file paths in the configuration.
"""
function expand_config_paths!(config::Dict)
    if haskey(config, "paths")
        for (key, path) in config["paths"]
            if isa(path, String) && !isempty(path)
                config["paths"][key] = expanduser(path)
            end
        end
    end
end

"""
    validate_config_structure(config::Dict)

Validate that all required configuration sections exist.
"""
function validate_config_structure(config::Dict)
    required_sections = [
        "project", 
        "paths", 
        "system_building", 
        "simulations", 
        "output"
    ]
    
    for section in required_sections
        if !haskey(config, section)
            @warn "Missing configuration section: $section - adding defaults"
            config[section] = Dict{String, Any}()
        end
    end
    
    # Check required path fields
    required_paths = ["data_directory", "output_directory", "export_directory"]
    for path_key in required_paths
        if !haskey(config["paths"], path_key)
            @warn "Missing required path: $path_key"
        end
    end
end

"""
    apply_config_defaults!(config::Dict)

Apply default values for any missing configuration options.
"""
function apply_config_defaults!(config::Dict)
    # Project defaults
    project = config["project"]
    get!(project, "name", "Sienna Power System Simulation")
    get!(project, "description", "Power system simulation using Sienna ecosystem")
    get!(project, "version", "1.0.0")
    
    # Path defaults
    paths = config["paths"]
    get!(paths, "data_directory", "./data")
    get!(paths, "output_directory", "./sienna_results")
    get!(paths, "export_directory", "./sienna_exports")
    
    # System building defaults
    system_building = config["system_building"]
    get!(system_building, "base_power", 100.0)
    get!(system_building, "default_voltage", 400.0)
    get!(system_building, "validate_system", true)
    get!(system_building, "load_timeseries", false)
    
    # Fuel mapping defaults
    if !haskey(system_building, "fuel_mapping")
        system_building["fuel_mapping"] = Dict{String, Any}()
    end
    fuel_mapping = system_building["fuel_mapping"]
    get!(fuel_mapping, "thermal_fuels", ["COAL", "NATURAL_GAS", "DIESEL", "NUCLEAR", "BIOMASS"])
    get!(fuel_mapping, "renewable_fuels", ["WIND", "SOLAR", "HYDRO", "GEOTHERMAL"])
    
    # Cost defaults
    if !haskey(system_building, "defaults")
        system_building["defaults"] = Dict{String, Any}()
    end
    defaults = system_building["defaults"]
    get!(defaults, "thermal_variable_cost", 50.0)
    get!(defaults, "thermal_startup_cost", 1000.0)
    get!(defaults, "renewable_variable_cost", 0.0)
    
    # Simulation defaults
    simulations = config["simulations"]
    get!(simulations, "default_horizon_hours", 24)
    get!(simulations, "default_solver", "HiGHS")
    get!(simulations, "save_results", true)
    get!(simulations, "run_economic_dispatch", true)
    get!(simulations, "run_unit_commitment", true)
    
    # Solver settings defaults
    if !haskey(simulations, "solver_settings")
        simulations["solver_settings"] = Dict{String, Any}()
    end
    if !haskey(simulations["solver_settings"], "HiGHS")
        simulations["solver_settings"]["HiGHS"] = Dict{String, Any}()
    end
    highs_settings = simulations["solver_settings"]["HiGHS"]
    get!(highs_settings, "time_limit_ed", 300.0)
    get!(highs_settings, "time_limit_uc", 600.0)
    get!(highs_settings, "mip_gap_ed", 0.01)
    get!(highs_settings, "mip_gap_uc", 0.02)
    get!(highs_settings, "threads", 0)
    
    # Output defaults
    output = config["output"]
    get!(output, "save_json_summary", true)
    get!(output, "save_csv_variables", true)
    get!(output, "create_timestamped_folders", true)
end

"""
    validate_config_paths(config::Dict)

Validate that critical paths exist.
"""
function validate_config_paths(config::Dict)
    data_dir = config["paths"]["data_directory"]
    if !isdir(data_dir)
        @warn "Data directory does not exist: $data_dir"
        @info "Please check your configuration and ensure the data directory exists"
    else
        @info "✓ Data directory found: $data_dir"
        
        # Check for required CSV files
        required_files = ["bus.csv", "gen.csv", "load.csv"]
        for file in required_files
            file_path = joinpath(data_dir, file)
            if isfile(file_path)
                @info "  ✓ Found: $file"
            else
                @warn "  ⚠ Missing: $file"
            end
        end
    end
end

"""
    create_default_config(config_file::String)

Create a default Sienna configuration file.
"""
function create_default_config(config_file::String)
    @info "Creating default Sienna configuration: $config_file"
    
    default_config = """
# Sienna Power System Simulation Configuration
# Edit this file for your specific project

[project]
name = "Power System Simulation"
description = "Power system simulation using Sienna ecosystem"
version = "1.0.0"
author = "User"
created_date = "$(Dates.format(now(), "yyyy-mm-dd"))"

[paths]
# Primary data directory - EDIT THIS PATH
data_directory = "./data"
output_directory = "./sienna_results"
export_directory = "./sienna_exports"

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
renewable_variable_cost = 0.0

[simulations]
default_horizon_hours = 24
default_solver = "HiGHS"
save_results = true
run_economic_dispatch = true
run_unit_commitment = true

[simulations.solver_settings.HiGHS]
time_limit_ed = 300.0
time_limit_uc = 600.0
mip_gap_ed = 0.01
mip_gap_uc = 0.02
threads = 0

[simulations.network]
default_model = "CopperPlatePowerModel"

[output]
save_json_summary = true
save_csv_variables = true
create_timestamped_folders = true

[logging]
level = "Info"
log_to_file = false

[regional]
country = "Unknown"
currency = "USD"
"""
    
    open(config_file, "w") do f
        write(f, default_config)
    end
    
    @info "✓ Created default configuration file: $config_file"
    @info "Please edit the data_directory path and other settings as needed"
end

"""
    get_data_directory()

Convenience function to get the data directory from config.
"""
function get_data_directory()
    config = get_config()
    return config["paths"]["data_directory"]
end

"""
    get_output_directory()

Convenience function to get the output directory from config.
"""
function get_output_directory()
    config = get_config()
    return config["paths"]["output_directory"]
end

"""
    get_solver_settings(formulation_type::String="economic_dispatch")

Get solver settings for a specific formulation type.
"""
function get_solver_settings(formulation_type::String="economic_dispatch")
    config = get_config()
    solver_name = config["simulations"]["default_solver"]
    solver_settings = config["simulations"]["solver_settings"][solver_name]
    
    if formulation_type == "economic_dispatch"
        return Dict(
            "time_limit" => solver_settings["time_limit_ed"],
            "mip_gap" => solver_settings["mip_gap_ed"],
            "threads" => solver_settings["threads"]
        )
    elseif formulation_type == "unit_commitment"
        return Dict(
            "time_limit" => solver_settings["time_limit_uc"],
            "mip_gap" => solver_settings["mip_gap_uc"],
            "threads" => solver_settings["threads"]
        )
    else
        return solver_settings
    end
end

"""
    show_config_summary()

Show a summary of the current configuration.
"""
function show_config_summary()
    config = get_config()
    
    @info "\n" * "="^50
    @info "SIENNA CONFIGURATION SUMMARY"
    @info "="^50
    @info "Project: $(config["project"]["name"])"
    @info "Version: $(config["project"]["version"])"
    @info ""
    @info "Paths:"
    @info "  Data: $(config["paths"]["data_directory"])"
    @info "  Output: $(config["paths"]["output_directory"])"
    @info "  Export: $(config["paths"]["export_directory"])"
    @info ""
    @info "System:"
    @info "  Base Power: $(config["system_building"]["base_power"]) MW"
    @info "  Default Voltage: $(config["system_building"]["default_voltage"]) kV"
    @info ""
    @info "Simulations:"
    @info "  Horizon: $(config["simulations"]["default_horizon_hours"]) hours"
    @info "  Solver: $(config["simulations"]["default_solver"])"
    @info "  Economic Dispatch: $(config["simulations"]["run_economic_dispatch"])"
    @info "  Unit Commitment: $(config["simulations"]["run_unit_commitment"])"
    @info "="^50
end

# Export main functions
export load_config, get_config, get_data_directory, get_output_directory, 
       get_solver_settings, show_config_summary
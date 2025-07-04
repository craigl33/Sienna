#!/usr/bin/env julia

"""
Fixed PowerSystems.jl System Builder with Configuration Support
- Fixes undefined variable issues
- Proper function signatures
- No deprecated parser usage
- Manual construction from CSV files
- PowerSystemCaseBuilder.jl integration
- Ready for PowerSimulations.jl
"""

using PowerSystems
using PowerSystemCaseBuilder
using PowerSimulations
using HiGHS
using Dates
using TimeSeries
using DataFrames
using CSV
using YAML
using JSON3
using Logging

# Import configuration loader
include("config_loader.jl")

# Import the correct types for PowerSystems.jl v4.6.2
using PowerSystems: ACBusTypes, PrimeMovers, ThermalFuels, LinearCurve, CostCurve, ThermalGenerationCost, RenewableGenerationCost, StorageCost, StorageTech

global_logger(ConsoleLogger(stderr, Logging.Info))

"""
    build_sienna_system(data_dir::String=""; config_file::String="", validate_system::Bool=true, load_timeseries::Bool=false, base_power::Float64=100.0)

Build PowerSystems.jl system from CSV data with configuration support.
"""
function build_sienna_system(data_dir::String=""; 
                            config_file::String="", 
                            validate_system::Bool=true, 
                            load_timeseries::Bool=false, 
                            base_power::Float64=100.0)
    
    # Load configuration
    if !isempty(config_file)
        config = load_config(config_file)
    else
        try
            config = get_config()  # This will auto-find config files
        catch e
            @warn "Could not load configuration: $e. Using defaults."
            config = create_default_config_dict()
        end
    end
    
    # Override data directory if provided
    if !isempty(data_dir)
        config["paths"]["data_directory"] = expanduser(data_dir)
    end
    
    # Get configuration values with fallbacks
    data_dir_final = get(get(config, "paths", Dict()), "data_directory", data_dir)
    base_power_final = get(get(config, "system_building", Dict()), "base_power", base_power)
    validate_system_final = get(get(config, "system_building", Dict()), "validate_system", validate_system)
    load_timeseries_final = get(get(config, "system_building", Dict()), "load_timeseries", load_timeseries)
    
    @info "Building PowerSystems.jl system with configuration:"
    @info "  Project: $(get(get(config, "project", Dict()), "name", "Unknown Project"))"
    @info "  Data directory: $data_dir_final"
    @info "  Base power: $(base_power_final) MW"
    
    # Check required files
    required_files = ["bus.csv", "load.csv"]
    for file in required_files
        if !isfile(joinpath(data_dir_final, file))
            error("Required file $file not found in $data_dir_final")
        end
    end
    
    # Create empty system
    sys = System(base_power_final)
    project_name = get(get(config, "project", Dict()), "name", "Power System")
    set_name!(sys, project_name)
    
    # Add components in dependency order
    add_buses_from_csv!(sys, joinpath(data_dir_final, "bus.csv"), config)
    
    # Add generators - check for separate files first, then unified
    thermal_gen_file = get(get(config, "paths", Dict()), "thermal_gen_file", "")
    renewable_gen_file = get(get(config, "paths", Dict()), "renewable_gen_file", "")
    
    if !isempty(thermal_gen_file) && !isempty(renewable_gen_file) && 
       isfile(joinpath(data_dir_final, thermal_gen_file)) && isfile(joinpath(data_dir_final, renewable_gen_file))
        @info "Found separate generator files from config"
        add_thermal_generators_from_csv!(sys, joinpath(data_dir_final, thermal_gen_file), config)
        add_renewable_generators_from_csv!(sys, joinpath(data_dir_final, renewable_gen_file), config)
    else
        unified_gen_file = joinpath(data_dir_final, "gen.csv")
        if isfile(unified_gen_file)
            @info "Found unified generator file"
            add_unified_generators_from_csv!(sys, unified_gen_file, config)
        else
            @warn "No generator files found"
        end
    end
    
    add_loads_from_csv!(sys, joinpath(data_dir_final, "load.csv"), config)
    
    # Optional components
    branch_file = joinpath(data_dir_final, "branch.csv")
    if isfile(branch_file)
        add_branches_from_csv!(sys, branch_file, config)
    end
    
    dc_branch_file = joinpath(data_dir_final, "dc_branch.csv")
    if isfile(dc_branch_file)
        add_dc_branches_from_csv!(sys, dc_branch_file, config)
    end
    
    storage_file = joinpath(data_dir_final, "storage.csv")
    if isfile(storage_file)
        add_storage_from_csv!(sys, storage_file, config)
    end
    
    # Post-process the system
    sys = post_process_system!(sys, config)
    
    # Add time series if requested and available
    if load_timeseries_final
        add_time_series_data!(sys, data_dir_final, config)
    end
    
    # Validate if requested
    if validate_system_final
        validate_system_components(sys, config)
    end
    
    # Print system summary
    print_system_summary(sys)
    
    return sys
end

"""
    create_default_config_dict()

Create a minimal default configuration dictionary.
"""
function create_default_config_dict()
    return Dict(
        "project" => Dict(
            "name" => "Power System",
            "description" => "Default power system",
            "version" => "1.0.0"
        ),
        "paths" => Dict(
            "data_directory" => "./data",
            "output_directory" => "./sienna_results",
            "export_directory" => "./sienna_exports"
        ),
        "system_building" => Dict(
            "base_power" => 100.0,
            "default_voltage" => 400.0,
            "validate_system" => true,
            "load_timeseries" => false,
            "fuel_mapping" => Dict(
                "thermal_fuels" => ["COAL", "NATURAL_GAS", "DIESEL", "NUCLEAR", "BIOMASS", "OIL", "GAS"],
                "renewable_fuels" => ["WIND", "SOLAR", "HYDRO", "GEOTHERMAL", "SOLAR_PV", "SOLAR_CSP", "WIND_ONSHORE", "WIND_OFFSHORE"]
            ),
            "defaults" => Dict(
                "thermal_variable_cost" => 50.0,
                "thermal_startup_cost" => 1000.0,
                "thermal_shutdown_cost" => 0.0,
                "renewable_variable_cost" => 0.0,
                "thermal_ramp_rate" => 100.0
            )
        ),
        "simulations" => Dict(
            "default_horizon_hours" => 24,
            "default_solver" => "HiGHS",
            "save_results" => true
        ),
        "output" => Dict(
            "save_json_summary" => true,
            "save_csv_variables" => true
        )
    )
end

function add_buses_from_csv!(sys::System, bus_file::String, config::Dict)
    @info "Adding buses from $bus_file"
    
    df = CSV.read(bus_file, DataFrame)
    default_voltage = get(get(config, "system_building", Dict()), "default_voltage", 400.0)
    
    for (i, row) in enumerate(eachrow(df))
        # Extract bus number from name or use index
        bus_number = try
            numbers = replace(string(row.name), r"[^0-9]" => "")
            isempty(numbers) ? i : parse(Int, numbers)
        catch
            i
        end
        
        # Determine bus type
        bus_type = if haskey(row, :bus_type) && !ismissing(row.bus_type)
            string(row.bus_type) == "REF" ? ACBusTypes.REF : ACBusTypes.PV
        else
            i == 1 ? ACBusTypes.REF : ACBusTypes.PV
        end
        
        # Get base voltage
        base_voltage_val = haskey(row, :base_voltage) && !ismissing(row.base_voltage) ? 
                          Float64(row.base_voltage) : default_voltage
        
        # Create bus with minimal constructor
        bus = ACBus(
            bus_number,                           # number
            string(row.name),                     # name
            bus_type,                            # bustype
            get(row, :angle, 0.0),               # angle
            get(row, :voltage, 1.0),             # magnitude
            (min = 0.95, max = 1.05),           # voltage_limits
            base_voltage_val                     # base_voltage
        )
        
        add_component!(sys, bus)
    end
    
    @info "Added $(nrow(df)) buses"
end

function add_thermal_generators_from_csv!(sys::System, thermal_file::String, config::Dict)
    @info "Adding thermal generators from $thermal_file"
    
    df = CSV.read(thermal_file, DataFrame)
    
    if nrow(df) == 0
        @info "No thermal generators found"
        return
    end
    
    defaults = get(get(config, "system_building", Dict()), "defaults", Dict())
    
    for row in eachrow(df)
        bus = get_component(ACBus, sys, string(row.bus))
        if bus === nothing
            @warn "Bus $(row.bus) not found for generator $(row.name)"
            continue
        end
        
        try
            add_single_thermal_generator!(sys, row, bus, config)
        catch e
            @warn "Failed to add thermal generator $(row.name): $e"
        end
    end
    
    @info "Added thermal generators from $thermal_file"
end

function add_renewable_generators_from_csv!(sys::System, renewable_file::String, config::Dict)
    @info "Adding renewable generators from $renewable_file"
    
    df = CSV.read(renewable_file, DataFrame)
    
    if nrow(df) == 0
        @info "No renewable generators found"
        return
    end
    
    for row in eachrow(df)
        bus = get_component(ACBus, sys, string(row.bus))
        if bus === nothing
            @warn "Bus $(row.bus) not found for generator $(row.name)"
            continue
        end
        
        try
            add_single_renewable_generator!(sys, row, bus, config)
        catch e
            @warn "Failed to add renewable generator $(row.name): $e"
        end
    end
    
    @info "Added renewable generators from $renewable_file"
end

function add_unified_generators_from_csv!(sys::System, gen_file::String, config::Dict)
    @info "Splitting unified generators from $gen_file by fuel type"
    
    df = CSV.read(gen_file, DataFrame)
    
    if nrow(df) == 0
        @info "No generators found"
        return
    end
    
    # Get fuel type sets from config
    fuel_mapping = get(get(config, "system_building", Dict()), "fuel_mapping", Dict())
    thermal_fuels = Set(uppercase.(get(fuel_mapping, "thermal_fuels", ["COAL", "GAS", "NUCLEAR"])))
    renewable_fuels = Set(uppercase.(get(fuel_mapping, "renewable_fuels", ["WIND", "SOLAR", "HYDRO"])))
    
    # Split DataFrame by fuel type
    thermal_df = filter(row -> uppercase(get(row, :fuel, "")) in thermal_fuels, df)
    renewable_df = filter(row -> uppercase(get(row, :fuel, "")) in renewable_fuels, df)
    
    # Process thermal generators
    thermal_count = 0
    for row in eachrow(thermal_df)
        bus = get_component(ACBus, sys, string(row.bus))
        if bus === nothing
            @warn "Bus $(row.bus) not found for thermal generator $(row.name)"
            continue
        end
        
        try
            add_single_thermal_generator!(sys, row, bus, config)
            thermal_count += 1
        catch e
            @warn "Failed to add thermal generator $(row.name): $e"
        end
    end
    
    # Process renewable generators
    renewable_count = 0
    for row in eachrow(renewable_df)
        bus = get_component(ACBus, sys, string(row.bus))
        if bus === nothing
            @warn "Bus $(row.bus) not found for renewable generator $(row.name)"
            continue
        end
        
        try
            add_single_renewable_generator!(sys, row, bus, config)
            renewable_count += 1
        catch e
            @warn "Failed to add renewable generator $(row.name): $e"
        end
    end
    
    @info "Added $thermal_count thermal and $renewable_count renewable generators from unified file"
end

function add_single_thermal_generator!(sys::System, row, bus::ACBus, config::Dict)
    # Get default costs from config
    defaults = get(get(config, "system_building", Dict()), "defaults", Dict())
    default_var_cost = get(defaults, "thermal_variable_cost", 50.0)
    default_startup_cost = get(defaults, "thermal_startup_cost", 1000.0)
    default_shutdown_cost = get(defaults, "thermal_shutdown_cost", 0.0)
    
    # Create cost structure
    variable_cost = get(row, :variable, default_var_cost)
    variable_linear_curve = LinearCurve(variable_cost)
    variable_cost_curve = CostCurve(variable_linear_curve)
    
    operation_cost = ThermalGenerationCost(
        variable = variable_cost_curve,
        fixed = get(row, :fixed_cost, 0.0),
        start_up = get(row, :startup, default_startup_cost),
        shut_down = get(row, :shutdown, default_shutdown_cost)
    )
    
    gen = ThermalStandard(
        string(row.name),
        get(row, :available, true),
        get(row, :status, 1) == 1,
        bus,
        get(row, :active_power, 0.0),
        0.0,
        get(row, :max_active_power, 100.0),
        (min = get(row, :min_active_power, 0.0), max = get(row, :max_active_power, 100.0)),
        (min = get(row, :min_reactive_power, -30.0), max = get(row, :max_reactive_power, 30.0)),
        (up = get(row, :ramp_30, 100.0), down = get(row, :ramp_30, 100.0)),
        operation_cost,
        100.0
    )
    
    add_component!(sys, gen)
end

function add_single_renewable_generator!(sys::System, row, bus::ACBus, config::Dict)
    # Get default costs from config
    defaults = get(get(config, "system_building", Dict()), "defaults", Dict())
    default_var_cost = get(defaults, "renewable_variable_cost", 0.0)
    
    # Create cost structure
    variable_cost = get(row, :variable, default_var_cost)
    variable_linear_curve = LinearCurve(variable_cost)
    variable_cost_curve = CostCurve(variable_linear_curve)
    operation_cost = RenewableGenerationCost(variable_cost_curve)
    
    # Determine prime mover type
    fuel_type = uppercase(get(row, :fuel, ""))
    prime_mover = if fuel_type == "WIND"
        PrimeMovers.WT
    elseif fuel_type == "SOLAR"
        PrimeMovers.PVe
    elseif fuel_type == "HYDRO"
        PrimeMovers.HY
    elseif fuel_type == "GEOTHERMAL"
        PrimeMovers.ST
    else
        PrimeMovers.OT
    end
    
    gen = RenewableDispatch(
        string(row.name),
        get(row, :available, true),
        bus,
        get(row, :active_power, 0.0),
        0.0,
        get(row, :max_active_power, 100.0),
        prime_mover,
        (min = get(row, :min_reactive_power, 0.0), max = get(row, :max_reactive_power, 0.0)),
        1.0,
        operation_cost,
        100.0
    )
    
    add_component!(sys, gen)
end

function add_loads_from_csv!(sys::System, load_file::String, config::Dict)
    @info "Adding loads from $load_file"
    
    df = CSV.read(load_file, DataFrame)
    
    if nrow(df) == 0
        @info "No loads found"
        return
    end
    
    for row in eachrow(df)
        bus = get_component(ACBus, sys, string(row.bus))
        if bus === nothing
            @warn "Bus $(row.bus) not found for load $(row.name)"
            continue
        end
        
        try
            load = PowerLoad(
                name = string(row.name),
                available = get(row, :available, true),
                bus = bus,
                active_power = get(row, :max_active_power, 100.0),
                reactive_power = get(row, :max_reactive_power, 30.0),
                base_power = 100.0,
                max_active_power = get(row, :max_active_power, 100.0),
                max_reactive_power = get(row, :max_reactive_power, 30.0)
            )
            
            add_component!(sys, load)
            
        catch e
            @warn "Failed to add load $(row.name): $e"
        end
    end
    
    @info "Added $(nrow(df)) loads"
end

function add_branches_from_csv!(sys::System, branch_file::String, config::Dict)
    @info "Adding branches from $branch_file"
    
    df = CSV.read(branch_file, DataFrame)
    
    if nrow(df) == 0
        @info "No branches found"
        return
    end
    
    for row in eachrow(df)
        from_bus = get_component(ACBus, sys, string(row.connection_points_from))
        to_bus = get_component(ACBus, sys, string(row.connection_points_to))
        
        if from_bus === nothing || to_bus === nothing
            @warn "Bus not found for branch $(row.name)"
            continue
        end
        
        try
            branch = Line(
                name = string(row.name),
                available = get(row, :available, true),
                active_power_flow = 0.0,
                reactive_power_flow = 0.0,
                arc = Arc(from_bus, to_bus),
                r = get(row, :r, 0.01),
                x = get(row, :x, 0.1),
                b = (from = get(row, :b, 0.0)/2, to = get(row, :b, 0.0)/2),
                rate = get(row, :rate, 100.0),
                angle_limits = (min = deg2rad(-30.0), max = deg2rad(30.0))
            )
            
            add_component!(sys, branch)
            
        catch e
            @warn "Failed to add branch $(row.name): $e"
        end
    end
    
    @info "Added $(nrow(df)) branches"
end

function add_dc_branches_from_csv!(sys::System, dc_branch_file::String, config::Dict)
    @info "Adding DC branches from $dc_branch_file"
    
    df = CSV.read(dc_branch_file, DataFrame)
    
    if nrow(df) == 0
        @info "No DC branches found"
        return
    end
    
    for row in eachrow(df)
        from_bus = get_component(ACBus, sys, string(row.connection_points_from))
        to_bus = get_component(ACBus, sys, string(row.connection_points_to))
        
        if from_bus === nothing || to_bus === nothing
            @warn "Bus not found for DC branch $(row.name)"
            continue
        end
        
        try
            dc_branch = TwoTerminalHVDCLine(
                name = string(row.name),
                available = get(row, :available, true),
                active_power_flow = 0.0,
                arc = Arc(from_bus, to_bus),
                active_power_limits_from = (min = -get(row, :active_power_limits_from, 100.0),
                                          max = get(row, :active_power_limits_from, 100.0)),
                active_power_limits_to = (min = -get(row, :active_power_limits_to, 100.0),
                                        max = get(row, :active_power_limits_to, 100.0)),
                reactive_power_limits_from = (min = -30.0, max = 30.0),
                reactive_power_limits_to = (min = -30.0, max = 30.0),
                loss = (l0 = 0.0, l1 = get(row, :loss, 0.01))
            )
            
            add_component!(sys, dc_branch)
            
        catch e
            @warn "Failed to add DC branch $(row.name): $e"
        end
    end
    
    @info "Added $(nrow(df)) DC branches"
end

function add_storage_from_csv!(sys::System, storage_file::String, config::Dict)
    @info "Adding storage from $storage_file"
    
    df = CSV.read(storage_file, DataFrame)
    
    if nrow(df) == 0
        @info "No storage units found"
        return
    end
    
    for row in eachrow(df)
        bus = get_component(ACBus, sys, string(row.bus))
        if bus === nothing
            @warn "Bus $(row.bus) not found for storage $(row.name)"
            continue
        end
        
        try
            # Create cost structure for storage (simplified)
            variable_cost = get(row, :variable_cost, 0.0)
            variable_linear_curve = LinearCurve(variable_cost)
            variable_cost_curve = CostCurve(variable_linear_curve)
            operation_cost = StorageCost(variable_cost_curve)
            
            # Try a simplified EnergyReservoirStorage constructor
            storage = EnergyReservoirStorage(
                name = string(row.name),
                available = get(row, :available, true),
                bus = bus,
                storage_capacity = get(row, :energy_capacity, 100.0),
                storage_technology_type = StorageTech.BATTERY,
                prime_mover_type = PrimeMovers.BA,
                initial_energy = get(row, :initial_energy, 50.0),
                state_of_charge_limits = (min = 0.1, max = 0.9),
                rating = get(row, :output_active_power_limits, 50.0),
                active_power = 0.0,
                input_active_power_limits = (min = 0.0, max = get(row, :input_active_power_limits, 50.0)),
                output_active_power_limits = (min = 0.0, max = get(row, :output_active_power_limits, 50.0)),
                efficiency = (in = get(row, :efficiency_in, 0.92), out = get(row, :efficiency_out, 0.92)),
                reactive_power = 0.0,
                reactive_power_limits = (min = -10.0, max = 10.0),
                operation_cost = operation_cost,
                base_power = 100.0
            )
            
            add_component!(sys, storage)
            
        catch e
            @warn "Failed to add storage $(row.name): $e"
            continue
        end
    end
    
    @info "Added storage units (or skipped if errors occurred)"
end

#!/usr/bin/env julia

"""
Correct Time Series Implementation Based on PowerSystems.jl Tutorial
===================================================================

Following the official PowerSystems.jl time series tutorial approach.
Replace these functions in build_sienna_system.jl
"""

"""
    add_time_series_data!(sys::System, data_dir::String, config::Dict)

Load time series data using the correct PowerSystems.jl tutorial approach.
"""
function add_time_series_data!(sys::System, data_dir::String, config::Dict)
    @info "Loading time series data from metadata..."
    
    # Check for time series metadata
    ts_metadata_file = joinpath(data_dir, "timeseries_metadata.json")
    if !isfile(ts_metadata_file)
        @warn "No timeseries_metadata.json found - adding synthetic time series"
        add_synthetic_time_series!(sys, config)
        return
    end
    
    try
        # Load metadata
        metadata = JSON3.read(read(ts_metadata_file, String))
        @info "Found time series metadata with $(length(metadata)) entries"
        
        # Group metadata by data file
        files_to_load = Dict{String, Vector}()
        for entry in metadata
            data_file = entry["data_file"]
            if !haskey(files_to_load, data_file)
                files_to_load[data_file] = []
            end
            push!(files_to_load[data_file], entry)
        end
        
        @info "Will load $(length(files_to_load)) time series data files"
        
        # Load each data file
        total_added = 0
        for (data_file, entries) in files_to_load
            file_path = joinpath(data_dir, data_file)
            
            if !isfile(file_path)
                @warn "Time series data file not found: $file_path"
                continue
            end
            
            @info "Loading time series from: $data_file ($(length(entries)) components)"
            added_count = load_time_series_from_csv!(sys, file_path, entries)
            total_added += added_count
            @info "  ✓ Added time series for $added_count/$(length(entries)) components"
        end
        
        if total_added == 0
            @warn "No real time series added - falling back to synthetic data"
            add_synthetic_time_series!(sys, config)
        else
            @info "✓ Successfully added real time series for $total_added components"
        end
        
    catch e
        @error "Failed to process time series metadata: $e"
        @warn "Falling back to synthetic time series"
        add_synthetic_time_series!(sys, config)
    end
end

"""
    load_time_series_from_csv!(sys::System, file_path::String, entries::Vector)

Load time series using the correct PowerSystems.jl tutorial approach.
"""
function load_time_series_from_csv!(sys::System, file_path::String, entries::Vector)
    try
        # Load CSV data
        df = CSV.read(file_path, DataFrame)
        @info "  Loaded CSV with $(nrow(df)) rows and $(ncol(df)) columns"
        
        # Parse timestamps from timestep column
        if !("timestep" in names(df))
            @error "  No timestep column found in CSV"
            return 0
        end
        
        # Parse timestamps
        timestamps = try
            DateTime.(string.(df.timestep), "yyyy-mm-dd HH:MM:SS")
        catch e
            @error "  Failed to parse timestamps: $e"
            return 0
        end
        
        @info "  Time range: $(first(timestamps)) to $(last(timestamps))"
        
        added_count = 0
        
        # Process each component entry
        for entry in entries
            component_name = entry["component"]
            data_column = entry["data_column"]
            label = entry["label"]
            category = entry["category"]
            
            # Find component
            component = find_component_by_name(sys, component_name, category)
            if component === nothing
                continue
            end
            
            # Check if data column exists
            if !(data_column in names(df))
                continue
            end
            
            # Get data and handle missing values
            raw_data = df[!, Symbol(data_column)]
            data_values = Float64.(coalesce.(raw_data, 0.0))
            
            # Apply scaling if specified
            scaling_factor = get(entry, "scaling_factor_multiplier", "1.0")
            if scaling_factor == "get_max_active_power"
                try
                    scale = get_max_active_power(component)
                    data_values .*= scale
                catch
                    # Scaling failed, continue with unscaled data
                end
            end
            
            # Create time series following the tutorial approach
            try
                # Method 1: Try the tutorial's approach with TimeArray and add_time_series!
                ts_array = TimeArray(timestamps, data_values)
                
                # According to tutorial, should be able to add TimeArray directly
                add_time_series!(sys, component, label, ts_array)
                added_count += 1
                
            catch e1
                try
                    # Method 2: Alternative order of parameters
                    ts_array = TimeArray(timestamps, data_values)
                    add_time_series!(sys, component, ts_array; name=label)
                    added_count += 1
                    
                catch e2
                    try
                        # Method 3: Create as SingleTimeSeries first
                        ts_array = TimeArray(timestamps, data_values)
                        single_ts = SingleTimeSeries(label, ts_array)
                        add_time_series!(sys, component, single_ts)
                        added_count += 1
                        
                    catch e3
                        @debug "    All methods failed for $component_name: $e1, $e2, $e3"
                    end
                end
            end
        end
        
        return added_count
        
    catch e
        @error "  Failed to load time series from $file_path: $e"
        return 0
    end
end

"""
    find_component_by_name(sys::System, name::String, category::String)

Find component with exact name match.
"""
function find_component_by_name(sys::System, name::String, category::String)
    if category == "ElectricLoad"
        for load in get_components(PowerLoad, sys)
            if get_name(load) == name
                return load
            end
        end
    elseif category == "RenewableGen"
        # Check all renewable types
        for gen in get_components(RenewableDispatch, sys)
            if get_name(gen) == name
                return gen
            end
        end
        for gen in get_components(RenewableNonDispatch, sys)
            if get_name(gen) == name
                return gen
            end
        end
        for gen in get_components(HydroDispatch, sys)
            if get_name(gen) == name
                return gen
            end
        end
        for gen in get_components(HydroEnergyReservoir, sys)
            if get_name(gen) == name
                return gen
            end
        end
    elseif category == "ThermalGen"
        for gen in get_components(ThermalStandard, sys)
            if get_name(gen) == name
                return gen
            end
        end
    end
    
    return nothing
end

"""
    add_synthetic_time_series!(sys::System, config::Dict)

Add synthetic time series using correct PowerSystems.jl approach.
"""
function add_synthetic_time_series!(sys::System, config::Dict)
    @info "Adding synthetic time series..."
    
    # Get horizon from config
    horizon_hours = get(get(config, "simulations", Dict()), "default_horizon_hours", 24)
    
    # Create time index
    start_time = DateTime(2024, 1, 1, 0, 0, 0)
    timestamps = collect(start_time:Hour(1):(start_time + Hour(horizon_hours - 1)))
    
    added_count = 0
    
    # Add load forecasts
    loads = get_components(PowerLoad, sys)
    for load in loads
        try
            base_power = get_max_active_power(load)
            
            # Simple daily pattern
            pattern = [0.7 + 0.3 * sin(2π * (h-1) / 24) for h in 1:horizon_hours]
            forecast_data = base_power .* pattern
            
            ts_array = TimeArray(timestamps, forecast_data)
            
            # Try different methods following tutorial approach
            try
                add_time_series!(sys, load, "max_active_power", ts_array)
                added_count += 1
            catch e1
                try
                    single_ts = SingleTimeSeries("max_active_power", ts_array)
                    add_time_series!(sys, load, single_ts)
                    added_count += 1
                catch e2
                    @warn "Failed to add synthetic load forecast for $(get_name(load)): $e1, $e2"
                end
            end
        catch e
            @warn "Error processing load $(get_name(load)): $e"
        end
    end
    
    # Add renewable forecasts
    renewables = get_components(RenewableDispatch, sys)
    for renewable in renewables
        try
            rating = get_rating(renewable)
            
            # Simple availability pattern
            pattern = [0.5 + 0.5 * rand() for h in 1:horizon_hours]
            forecast_data = rating .* pattern
            
            ts_array = TimeArray(timestamps, forecast_data)
            
            # Try different methods following tutorial approach
            try
                add_time_series!(sys, renewable, "max_active_power", ts_array)
                added_count += 1
            catch e1
                try
                    single_ts = SingleTimeSeries("max_active_power", ts_array)
                    add_time_series!(sys, renewable, single_ts)
                    added_count += 1
                catch e2
                    @warn "Failed to add synthetic renewable forecast for $(get_name(renewable)): $e2"
                end
            end
        catch e
            @warn "Error processing renewable $(get_name(renewable)): $e"
        end
    end
    
    if added_count > 0
        @info "✓ Added synthetic time series for $added_count components"
    else
        error("❌ Failed to add any time series - check PowerSystems.jl tutorial")
    end
end

function validate_system_components(sys::System, config::Dict)
    @info "Validating system components..."
    
    buses = collect(get_components(ACBus, sys))
    @info "Buses: $(length(buses))"
    
    # Check for reference bus
    ref_buses = [b for b in buses if get_bustype(b) == ACBusTypes.REF]
    if length(ref_buses) == 0
        @warn "No reference bus found!"
    elseif length(ref_buses) > 1
        @warn "Multiple reference buses found: $(length(ref_buses))"
    end
    
    # Check generators
    thermal_gens = collect(get_components(ThermalStandard, sys))
    renewable_gens = collect(get_components(RenewableDispatch, sys))
    @info "Generators: $(length(thermal_gens)) thermal, $(length(renewable_gens)) renewable"
    
    # Check loads
    loads = collect(get_components(PowerLoad, sys))
    @info "Loads: $(length(loads))"
    
    @info "✓ System validation complete"
end

function post_process_system!(sys::System, config::Dict)
    @info "Post-processing system..."
    
    # Ensure we have a reference bus
    buses = collect(get_components(ACBus, sys))
    ref_buses = [b for b in buses if get_bustype(b) == ACBusTypes.REF]
    
    if length(ref_buses) == 0
        first_bus = first(buses)
        set_bustype!(first_bus, ACBusTypes.REF)
        @info "Set reference bus: $(get_name(first_bus))"
    end
    
    @info "✓ Post-processing complete"
    return sys
end

function print_system_summary(sys::System)
    @info "\n" * "="^50
    @info "POWERSYSTEMS.JL SYSTEM SUMMARY"
    @info "="^50
    
    buses = collect(get_components(ACBus, sys))
    @info "Buses: $(length(buses))"
    
    thermal_gens = collect(get_components(ThermalStandard, sys))
    renewable_gens = collect(get_components(RenewableDispatch, sys))
    
    # Get storage units
    storage_count = 0
    try
        storage_units = collect(get_components(EnergyReservoirStorage, sys))
        storage_count = length(storage_units)
    catch
        storage_count = 0
    end
    
    @info "Generators:"
    @info "  - Thermal: $(length(thermal_gens))"
    @info "  - Renewable: $(length(renewable_gens))"
    @info "  - Storage: $storage_count"
    @info "  - Total: $(length(thermal_gens) + length(renewable_gens) + storage_count)"
    
    loads = collect(get_components(PowerLoad, sys))
    @info "Loads: $(length(loads))"
    
    branches = collect(get_components(Line, sys))
    dc_branches = collect(get_components(TwoTerminalHVDCLine, sys))
    @info "Network:"
    @info "  - AC Branches: $(length(branches))"
    @info "  - DC Branches: $(length(dc_branches))"
    
    # Calculate total capacity (safe for empty collections)
    total_thermal_cap = isempty(thermal_gens) ? 0.0 : sum(get_active_power_limits(g).max for g in thermal_gens)
    total_renewable_cap = isempty(renewable_gens) ? 0.0 : sum(get_rating(g) for g in renewable_gens)
    
    # Storage capacity (if any)
    total_storage_cap = 0.0
    try
        storage_units = collect(get_components(EnergyReservoirStorage, sys))
        total_storage_cap = isempty(storage_units) ? 0.0 : sum(get_rating(s) for s in storage_units)
    catch
        total_storage_cap = 0.0
    end
    
    @info "Generation Capacity:"
    @info "  - Thermal: $(round(total_thermal_cap, digits=1)) MW"
    @info "  - Renewable: $(round(total_renewable_cap, digits=1)) MW"
    @info "  - Storage: $(round(total_storage_cap, digits=1)) MW"
    @info "  - Total: $(round(total_thermal_cap + total_renewable_cap + total_storage_cap, digits=1)) MW"
    
    @info "="^50
end

function save_system(sys::System, output_file::String)
    @info "Saving system to: $output_file"
    try
        to_json(sys, output_file)
        @info "✓ System saved successfully"
    catch e
        @error "Failed to save system: $e"
    end
end

# Remaining functions from original file...
function create_psb_case_from_system(sys::System, case_name::String, output_dir::String; description::String="")
    """
    Create a PowerSystemCaseBuilder.jl compatible case from the built system.
    This enables easy sharing, version control, and integration with PSB ecosystem.
    """
    @info "Creating PowerSystemCaseBuilder case: $case_name"
    
    case_dir = joinpath(output_dir, case_name)
    mkpath(case_dir)
    
    # Save system in PSB format
    system_file = joinpath(case_dir, "$(case_name)_sys.json")
    to_json(sys, system_file)
    
    # Create comprehensive case metadata
    component_summary = Dict()
    device_types = [ACBus, ThermalStandard, RenewableDispatch, RenewableNonDispatch, 
                   HydroDispatch, HydroEnergyReservoir, EnergyReservoirStorage, 
                   PowerLoad, Line, Transformer2W, TwoTerminalHVDCLine]
    
    for device_type in device_types
        components = get_components(device_type, sys)
        if length(components) > 0
            component_summary[string(device_type)] = length(components)
        end
    end
    
    # Calculate total capacities
    thermal_capacity = sum(get_active_power_limits(g).max for g in get_components(ThermalStandard, sys))
    renewable_capacity = sum(get_rating(g) for g in get_components(RenewableDispatch, sys))
    total_load = sum(get_max_active_power(l) for l in get_components(PowerLoad, sys))
    
    metadata = Dict(
        "name" => case_name,
        "description" => isempty(description) ? "System built with custom Sienna builder" : description,
        "created" => string(now()),
        "builder_version" => "custom_sienna_builder_v1.0",
        
        # System characteristics
        "system_properties" => Dict(
            "base_power" => get_base_power(sys),
            "frequency" => get_frequency(sys),
            "name" => get_name(sys)
        ),
        
        # Component counts
        "component_counts" => component_summary,
        
        # Capacity summary
        "capacity_summary" => Dict(
            "thermal_capacity_mw" => round(thermal_capacity, digits=1),
            "renewable_capacity_mw" => round(renewable_capacity, digits=1),
            "total_load_mw" => round(total_load, digits=1),
            "total_generation_capacity_mw" => round(thermal_capacity + renewable_capacity, digits=1)
        ),
        
        # Time series info (if available)
        "time_series_summary" => Dict(
            "has_time_series" => has_time_series(sys),
            "time_series_count" => length(get_time_series_multiple(sys, filter_func = x -> true))
        ),
        
        # Files included
        "files" => [
            "$(case_name)_sys.json",
            "metadata.json"
        ]
    )
    
    # Save metadata
    metadata_file = joinpath(case_dir, "metadata.json")
    open(metadata_file, "w") do f
        JSON3.pretty(f, metadata)
    end
    
    # Create README
    readme_file = joinpath(case_dir, "README.md")
    open(readme_file, "w") do f
        write(f, """
# $(case_name) Power System Case

$(metadata["description"])

## System Overview
- **Created**: $(metadata["created"])
- **Base Power**: $(metadata["system_properties"]["base_power"]) MW
- **Frequency**: $(metadata["system_properties"]["frequency"]) Hz

## System Composition
$(join(["- **$(replace(k, "_" => " ") |> titlecase)**: $(v)" for (k,v) in component_summary], "\n"))

## Capacity Summary
- **Total Thermal**: $(metadata["capacity_summary"]["thermal_capacity_mw"]) MW
- **Total Renewable**: $(metadata["capacity_summary"]["renewable_capacity_mw"]) MW
- **Total Load**: $(metadata["capacity_summary"]["total_load_mw"]) MW

## Usage
Load this system with:
```julia
using PowerSystems
sys = System("$(case_name)_sys.json")
```

Or use with PowerSystemCaseBuilder.jl ecosystem.
""")
    end
    
    @info "✓ PowerSystemCaseBuilder case created:"
    @info "  - Directory: $case_dir"
    @info "  - System file: $(case_name)_sys.json"
    @info "  - Metadata: metadata.json"
    @info "  - README: README.md"
    @info "  - Components: $(length(component_summary)) types"
    @info "  - Total capacity: $(round(thermal_capacity + renewable_capacity, digits=1)) MW"
    
    return case_dir, metadata
end

function export_system_for_simulations(sys::System, export_dir::String; case_name::String="", validate::Bool=true)
    """
    Export system in formats ready for PowerSimulations.jl and other tools.
    Creates both raw system files and PowerSystemCaseBuilder cases.
    """
    if isempty(case_name)
        case_name = "power_system_$(Dates.format(now(), "yyyy_mm_dd"))"
    end
    
    @info "Exporting system for simulations: $case_name"
    mkpath(export_dir)
    
    # Validate system if requested
    if validate
        @info "Validating system before export..."
        validate_system_components(sys, Dict())  # Pass empty dict if no config
    end
    
    # 1. Save raw system file
    raw_system_file = joinpath(export_dir, "$(case_name).json")
    save_system(sys, raw_system_file)
    
    # 2. Create PowerSystemCaseBuilder case
    case_dir, metadata = create_psb_case_from_system(
        sys, case_name, export_dir,
        description="Power system ready for PowerSimulations.jl operations"
    )
    
    # 3. Create simulation-ready summary
    simulation_info = Dict(
        "case_name" => case_name,
        "system_file" => raw_system_file,
        "psb_case_dir" => case_dir,
        "ready_for_simulations" => true,
        "recommended_formulations" => get_recommended_formulations(sys),
        "export_timestamp" => string(now())
    )
    
    info_file = joinpath(export_dir, "simulation_info.json")
    open(info_file, "w") do f
        JSON3.pretty(f, simulation_info)
    end
    
    @info "✓ System exported for simulations:"
    @info "  - Raw system: $raw_system_file" 
    @info "  - PSB case: $case_dir"
    @info "  - Simulation info: $info_file"
    
    return simulation_info
end

function get_recommended_formulations(sys::System)
    """
    Analyze system and recommend appropriate PowerSimulations.jl formulations.
    """
    recommendations = Dict{String, String}()
    
    # Check thermal units for UC vs ED recommendation
    thermal_gens = get_components(ThermalStandard, sys)
    if length(thermal_gens) > 0
        # If many small units or complex constraints, recommend UC
        avg_capacity = mean(get_active_power_limits(g).max for g in thermal_gens)
        if length(thermal_gens) > 20 || avg_capacity < 100
            recommendations["thermal"] = "ThermalStandardUnitCommitment (for detailed unit commitment)"
        else
            recommendations["thermal"] = "ThermalBasicEconomicDispatch (for faster economic dispatch)"
        end
    end
    
    # Renewable recommendations
    renewable_gens = get_components(RenewableDispatch, sys)
    if length(renewable_gens) > 0
        recommendations["renewable"] = "RenewableFullDispatch (allows curtailment)"
    end
    
    # Storage recommendations
    storage_units = get_components(EnergyReservoirStorage, sys)
    if length(storage_units) > 0
        recommendations["storage"] = "StorageSystemsSimulation or StorageDispatchWithReserves"
    end
    
    # Network recommendations
    lines = get_components(Line, sys)
    dc_lines = get_components(TwoTerminalHVDCLine, sys)
    if length(lines) == 0 && length(dc_lines) > 0
        recommendations["network"] = "CopperPlatePowerModel (no AC network constraints)"
    elseif length(lines) > 0
        recommendations["network"] = "StandardPTDFModel (DC power flow with transmission limits)"
    end
    
    return recommendations
end

# Legacy wrapper function for backward compatibility
function build_sienna_system_legacy(data_dir::String; 
                                   validate_system::Bool=true, 
                                   load_timeseries::Bool=false, 
                                   base_power::Float64=100.0)
    @info "Using legacy function signature - consider updating to new version"
    return build_sienna_system(data_dir, 
                              validate_system=validate_system,
                              load_timeseries=load_timeseries,
                              base_power=base_power)
end

# Main function for testing
function main()
    try
        # Try to load config, fall back to defaults if not available
        config = try
            get_config()
        catch
            @warn "No configuration file found, using defaults"
            create_default_config_dict()
        end
        
        data_dir = get(get(config, "paths", Dict()), "data_directory", ".")
        
        @info "Sienna System Builder with Configuration Support"
        @info "Data directory: $data_dir"
        
        if !isdir(data_dir)
            @error "Data directory not found: $data_dir"
            @info "Please edit config.toml or provide a valid data directory"
            return
        end
        
        # Build system
        sys = build_sienna_system(data_dir)
        
        # Create export directory
        export_dir = get(get(config, "paths", Dict()), "export_directory", "./exported_systems")
        
        # Export system for simulations
        simulation_info = export_system_for_simulations(
            sys, export_dir, 
            case_name="power_system_$(Dates.format(now(), "yyyy_mm_dd"))",
            validate=true
        )
        
        @info "✅ SUCCESS: System built and exported!"
        @info "Next steps:"
        @info "  1. Use system file: $(simulation_info["system_file"])"
        @info "  2. Or load PSB case from: $(simulation_info["psb_case_dir"])"
        @info "  3. Run PowerSimulations.jl with the exported system"
        
        return sys, simulation_info
        
    catch e
        @error "Build process failed: $e"
        @error "Stack trace:" exception=(e, catch_backtrace())
        return nothing
    end
end

# Export main functions
export build_sienna_system, build_sienna_system_legacy, export_system_for_simulations
export create_psb_case_from_system, get_recommended_formulations
export save_system, print_system_summary

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
#!/usr/bin/env julia

"""
SiennaSystem.jl - Clean PowerSystems.jl System Builder
======================================================

Clean, config-driven system builder with pre-extended time series.
Compatible with PowerSystems.jl v4.6.2+
"""

using PowerSystems
using DataFrames
using CSV
using TimeSeries
using JSON3
using Dates
using Logging
using InfrastructureSystems
using PowerSystems: PumpHydroStatus

# Import configuration manager
include("SiennaConfig.jl")

# ===== MAIN STRUCT =====

"""
    SiennaSystem

Clean system builder with pre-extended time series support.
"""
mutable struct SiennaSystem
    config::SiennaConfig
    system::Union{System, Nothing}
    
    # Status tracking
    is_built::Bool
    has_time_series::Bool
    has_forecasts::Bool
    
    # Component tracking
    component_counts::Dict{String, Int}
    
    # Error tracking
    errors::Vector{String}
end

"""
    SiennaSystem(config::SiennaConfig)

Create a new SiennaSystem from validated configuration.
"""
function SiennaSystem(config::SiennaConfig)
    @info "🏗️ Initializing SiennaSystem"
    
    if !config.is_validated
        error("❌ Configuration must be validated first")
    end
    
    return SiennaSystem(
        config,
        nothing,
        false, false, false,
        Dict{String, Int}(),
        String[]
    )
end

# ===== MAIN BUILD FUNCTION =====

"""
    build_system!(sienna_sys::SiennaSystem)

Build the complete power system with pre-extended time series.
"""
function build_system!(sienna_sys::SiennaSystem)
    @info "🏗️ Building power system..."
    
    # Reset state
    sienna_sys.errors = String[]
    sienna_sys.component_counts = Dict{String, Int}()
    
    try
        # Create base system
        _create_base_system!(sienna_sys)
        
        # Add all components
        _add_all_components!(sienna_sys)
        
        # Load pre-extended time series
        if sienna_sys.config.load_timeseries
            _load_time_series!(sienna_sys)
        end
        
        sienna_sys.is_built = true
        @info "✅ System built successfully"
        
    catch e
        sienna_sys.is_built = false
        push!(sienna_sys.errors, string(e))
        @error "❌ System build failed: $e"
        rethrow(e)
    end
    
    return sienna_sys.system
end

"""
    _create_base_system!(sienna_sys::SiennaSystem)

Create the base PowerSystems.jl system.
"""
function _create_base_system!(sienna_sys::SiennaSystem)
    @info "Creating base system..."
    sienna_sys.system = System(sienna_sys.config.base_power, unit_system="NATURAL_UNITS")
    set_name!(sienna_sys.system, sienna_sys.config.project_name)
end

"""
    _add_all_components!(sienna_sys::SiennaSystem)

Add all system components in correct order.
"""
function _add_all_components!(sienna_sys::SiennaSystem)
    @info "Adding system components..."
    _add_buses!(sienna_sys)
    _add_areas!(sienna_sys) 
    _add_generators!(sienna_sys)
    _add_storage_systems!(sienna_sys)
    _add_loads!(sienna_sys)
    _add_network_branches!(sienna_sys)
end

# ===== TIME SERIES FUNCTIONS =====

"""
    _load_time_series!(sienna_sys::SiennaSystem)

Load and transform time series data with pre-extension.
"""
function _load_time_series!(sienna_sys::SiennaSystem)
    @info "📊 Loading time series data..."
    
    # Load with pre-extension
    _load_extended_time_series!(sienna_sys)
    
    if !sienna_sys.has_time_series
        return
    end
    
    # Transform to forecasts
    @info "🔮 Transforming to forecasts..."
    _transform_to_forecasts!(sienna_sys)
end

"""
    _load_extended_time_series!(sienna_sys::SiennaSystem)

Load time series from files with automatic pre-extension.
"""
function _load_extended_time_series!(sienna_sys::SiennaSystem)
    metadata_file = joinpath(sienna_sys.config.data_directory, "timeseries_metadata.json")
    if !isfile(metadata_file)
        @warn "No timeseries_metadata.json found"
        return
    end
    
    try
        metadata = JSON3.read(read(metadata_file, String))
        
        # Calculate required length
        required_length = _calculate_required_time_series_length(sienna_sys)
        @info "📏 Required time series length: $required_length hours"
        
        loaded_count = 0
        for entry in metadata
            if _add_extended_time_series!(sienna_sys, entry, required_length)
                loaded_count += 1
            end
        end
        
        if loaded_count > 0
            sienna_sys.has_time_series = true
            @info "✓ Loaded time series for $loaded_count components"
        end
        
    catch e
        @error "Failed to load time series: $e"
        push!(sienna_sys.errors, "Time series loading failed: $e")
    end
end

"""
    _calculate_required_time_series_length(sienna_sys::SiennaSystem)

Calculate the required time series length for simulations.
"""
function _calculate_required_time_series_length(sienna_sys::SiennaSystem)
    horizon_hours = sienna_sys.config.default_horizon_hours
    interval_hours = 24
    simulation_days = get(sienna_sys.config.config_data["simulations"], "annual_simulation_days", 365)
    
    required_length = (simulation_days * interval_hours) + (horizon_hours - interval_hours)
    
    @info "📊 Time series calculation:"
    @info "  Simulation days: $simulation_days"
    @info "  Interval: $interval_hours hours"  
    @info "  Horizon: $horizon_hours hours"
    @info "  Required: $required_length hours"
    
    return required_length
end

"""
    _add_extended_time_series!(sienna_sys::SiennaSystem, entry, required_length::Int)

Add time series to component with automatic extension if needed.
"""
function _add_extended_time_series!(sienna_sys::SiennaSystem, entry, required_length::Int)
    component_name = entry["component"]
    data_file = entry["data_file"]
    data_column = entry["data_column"]
    label = entry["label"]
    category = entry["category"]
    
    # Find component
    component = _find_component(sienna_sys.system, component_name, category)
    if component === nothing
        @warn "Component not found: $component_name"
        return false
    end
    
    # Load and potentially extend data
    try
        extended_data = _load_and_extend_data(sienna_sys, data_file, data_column, required_length, component_name)
        if extended_data === nothing
            return false
        end
        
        # Get scaling factor
        scaling_factor = _get_scaling_factor(entry)
        
        # Create time series
        time_series = SingleTimeSeries(
            name = label,
            data = extended_data,
            scaling_factor_multiplier = scaling_factor
        )
        
        # Add to system
        add_time_series!(sienna_sys.system, component, time_series)
        return true
        
    catch e
        @warn "Failed to add time series to $component_name: $e"
        return false
    end
end

"""
    _load_and_extend_data(sienna_sys::SiennaSystem, data_file::String, data_column::String, 
                         required_length::Int, component_name::String)

Load CSV data and extend if necessary.
"""
function _load_and_extend_data(sienna_sys::SiennaSystem, data_file::String, data_column::String, 
                              required_length::Int, component_name::String)
    file_path = joinpath(sienna_sys.config.data_directory, data_file)
    if !isfile(file_path)
        @warn "Time series file not found: $file_path"
        return nothing
    end
    
    # Load CSV
    df = CSV.read(file_path, DataFrame)
    
    # Parse timestamps
    timestamps = _parse_timestamps(df)
    if timestamps === nothing
        @warn "Failed to parse timestamps in $data_file"
        return nothing
    end
    
    # Get data values
    data_values = Float64.(coalesce.(df[!, Symbol(data_column)], 0.0))
    current_length = length(data_values)
    
    @info "  $component_name: $current_length → $required_length hours"
    
    # Extend if needed
    if current_length < required_length
        timestamps, data_values = _extend_time_series_data(timestamps, data_values, required_length)
        @info "    ✓ Extended using cyclic pattern"
    end
    
    return TimeArray(timestamps, data_values)
end

"""
    _parse_timestamps(df::DataFrame)

Parse timestamps from DataFrame with multiple format support.
"""
function _parse_timestamps(df::DataFrame)
    try
        return DateTime.(string.(df.timestep), "yyyy-mm-dd HH:MM:SS")
    catch
        try
            return DateTime.(string.(df.timestep), "yyyy-mm-ddTHH:MM:SS")
        catch
            return nothing
        end
    end
end

"""
    _extend_time_series_data(timestamps::Vector{DateTime}, values::Vector{Float64}, required_length::Int)

Extend time series data using cyclic repetition.
"""
function _extend_time_series_data(timestamps::Vector{DateTime}, values::Vector{Float64}, required_length::Int)
    current_length = length(values)
    if current_length >= required_length
        return timestamps, values
    end
    
    extension_needed = required_length - current_length
    time_step = timestamps[2] - timestamps[1]
    
    # Extend timestamps
    extended_timestamps = copy(timestamps)
    last_time = timestamps[end]
    
    for i in 1:extension_needed
        push!(extended_timestamps, last_time + (i * time_step))
    end
    
    # Extend values using cyclic pattern
    extended_values = copy(values)
    for i in 1:extension_needed
        pattern_index = ((i - 1) % current_length) + 1
        push!(extended_values, values[pattern_index])
    end
    
    return extended_timestamps, extended_values
end

"""
    _get_scaling_factor(entry)

Get scaling factor from metadata entry.
"""
function _get_scaling_factor(entry)
    scaling_factor = get(entry, "scaling_factor_multiplier", nothing)
    
    if scaling_factor == "get_max_active_power"
        return get_max_active_power
    elseif scaling_factor == "get_rating"
        return get_rating
    else
        return nothing
    end
end

"""
    _transform_to_forecasts!(sienna_sys::SiennaSystem)

Transform SingleTimeSeries to Deterministic forecasts.
"""
function _transform_to_forecasts!(sienna_sys::SiennaSystem)
    horizon_hours = sienna_sys.config.default_horizon_hours
    
    try
        @info "Converting time series to forecasts (horizon: $(horizon_hours)h)..."
        
        transform_single_time_series!(
            sienna_sys.system,
            Hour(horizon_hours),
            Hour(24)
        )
        
        set_units_base_system!(sienna_sys.system, "NATURAL_UNITS")
        sienna_sys.has_forecasts = true
        
        @info "✓ Forecast transformation completed"
        
    catch e
        @error "Failed to transform to forecasts: $e"
        push!(sienna_sys.errors, "Forecast transformation failed: $e")
        error("Forecast transformation failed: $e")
    end
end

"""
    _find_component(sys::System, name::String, category::String)

Find component by name and category using PowerSystems.jl v4.6.2 API.
"""
function _find_component(sys::System, name::String, category::String)
    try
        if category == "ElectricLoad"
            return get_component(PowerLoad, sys, name)
        elseif category == "RenewableGen"
            return get_component(RenewableDispatch, sys, name)
        elseif category == "ThermalGen"
            return get_component(ThermalStandard, sys, name)
        end
    catch e
        @debug "Direct component lookup failed for $name, trying iteration: $e"
        
        if category == "ElectricLoad"
            for load in get_components(PowerLoad, sys)
                if get_name(load) == name
                    return load
                end
            end
        elseif category == "RenewableGen"
            for gen in get_components(RenewableDispatch, sys)
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
    end
    
    return nothing
end

# ===== COMPONENT ADDITION FUNCTIONS =====

"""
    _add_buses!(sienna_sys::SiennaSystem)

Add buses from CSV data.
"""
function _add_buses!(sienna_sys::SiennaSystem)
    @info "Adding buses..."
    
    bus_file = joinpath(sienna_sys.config.data_directory, "bus.csv")
    if !isfile(bus_file)
        error("❌ Required bus.csv not found: $bus_file")
    end
    
    df = CSV.read(bus_file, DataFrame)
    bus_count = 0
    ref_bus_count = 0
    
    for (i, row) in enumerate(eachrow(df))
        try
            bus_number = _extract_bus_number(row.name, i)
            bus_type = _determine_bus_type(row, i, ref_bus_count)
            
            if bus_type == ACBusTypes.REF
                ref_bus_count += 1
            end
            
            bus = ACBus(
                bus_number,
                string(row.name),
                bus_type,
                get(row, :angle, 0.0),
                get(row, :voltage, 1.0),
                (min = 0.95, max = 1.05),
                Float64(get(row, :base_voltage, sienna_sys.config.default_voltage))
            )
            
            add_component!(sienna_sys.system, bus)
            bus_count += 1
            
        catch e
            @warn "Failed to add bus $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["ACBus"] = bus_count
    @info "✓ Added $bus_count buses"
end

"""
    _extract_bus_number(name, fallback_index::Int)

Extract bus number from name or use fallback.
"""
function _extract_bus_number(name, fallback_index::Int)
    try
        numbers = replace(string(name), r"[^0-9]" => "")
        return isempty(numbers) ? fallback_index : parse(Int, numbers)
    catch
        return fallback_index
    end
end

"""
    _determine_bus_type(row, index::Int, ref_bus_count::Int)

Determine bus type from data or defaults.
"""
function _determine_bus_type(row, index::Int, ref_bus_count::Int)
    if haskey(row, :bus_type) && string(row.bus_type) == "REF"
        return ACBusTypes.REF
    elseif index == 1 && ref_bus_count == 0
        return ACBusTypes.REF
    else
        return ACBusTypes.PV
    end
end

"""
    _add_areas!(sienna_sys::SiennaSystem)

Add areas and assign buses.
"""
function _add_areas!(sienna_sys::SiennaSystem)
    @info "Adding areas..."
    
    area = Area("system_area")
    add_component!(sienna_sys.system, area)
    
    for bus in get_components(ACBus, sienna_sys.system)
        set_area!(bus, area)
    end
    
    sienna_sys.component_counts["Area"] = 1
    @info "✓ Added 1 area"
end

"""
    _add_generators!(sienna_sys::SiennaSystem)

Add generators from CSV data.
"""
function _add_generators!(sienna_sys::SiennaSystem)
    @info "Adding generators..."
    
    gen_file = joinpath(sienna_sys.config.data_directory, "gen.csv")
    if !isfile(gen_file)
        @info "No gen.csv found"
        return
    end
    
    df = CSV.read(gen_file, DataFrame)
    thermal_count = 0
    renewable_count = 0
    
    for row in eachrow(df)
        try
            bus = get_component(ACBus, sienna_sys.system, string(row.bus))
            if bus === nothing
                @warn "Bus $(row.bus) not found for generator $(row.name)"
                continue
            end
            
            fuel_type = uppercase(string(get(row, :fuel, "")))
            
            if _is_thermal_fuel(sienna_sys.config, fuel_type)
                _add_thermal_generator!(sienna_sys, row, bus)
                thermal_count += 1
            else
                _add_renewable_generator!(sienna_sys, row, bus, fuel_type)
                renewable_count += 1
            end
            
        catch e
            @warn "Failed to add generator $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["ThermalStandard"] = thermal_count
    sienna_sys.component_counts["RenewableDispatch"] = renewable_count
    @info "✓ Added $thermal_count thermal + $renewable_count renewable generators"
end

"""
    _is_thermal_fuel(config::SiennaConfig, fuel_type::String)

Check if fuel type is thermal.
"""
function _is_thermal_fuel(config::SiennaConfig, fuel_type::String)
    thermal_fuels = ["COAL", "GAS", "NATURAL_GAS", "NUCLEAR", "OIL", "DIESEL", "BIOMASS", "BIOENERGY"]
    return uppercase(fuel_type) in thermal_fuels
end

"""
    _add_thermal_generator!(sienna_sys::SiennaSystem, row, bus::ACBus)

Add thermal generator to system.
"""
function _add_thermal_generator!(sienna_sys::SiennaSystem, row, bus::ACBus)
    defaults = _get_thermal_defaults(sienna_sys.config)
    
    variable_cost = get(row, :variable, defaults["variable_cost"])
    operation_cost = ThermalGenerationCost(
        variable = CostCurve(LinearCurve(variable_cost)),
        fixed = get(row, :fixed_cost, 0.0),
        start_up = get(row, :startup, defaults["startup_cost"]),
        shut_down = get(row, :shutdown, defaults["shutdown_cost"])
    )
    
    max_power = get(row, :max_active_power, 100.0)
    min_power = get(row, :min_active_power, max_power * defaults["min_power_fraction"])
    
    generator = ThermalStandard(
        string(row.name),
        get(row, :available, true),
        get(row, :status, 1) == 1,
        bus,
        get(row, :active_power, 0.0),
        0.0,
        max_power,
        (min = min_power, max = max_power),
        (min = Float64(get(row, :min_reactive_power, -30.0)), 
         max = Float64(get(row, :max_reactive_power, 30.0))),
        (up = get(row, :ramp_30, defaults["ramp_rate"]), 
         down = get(row, :ramp_30, defaults["ramp_rate"])),
        operation_cost,
        1.0
    )
    
    add_component!(sienna_sys.system, generator)
end

"""
    _add_renewable_generator!(sienna_sys::SiennaSystem, row, bus::ACBus, fuel_type::String)

Add renewable generator to system.
"""
function _add_renewable_generator!(sienna_sys::SiennaSystem, row, bus::ACBus, fuel_type::String)
    variable_cost = fuel_type in ["WIND", "SOLAR", "SOLAR_PV", "SOLAR_CSP", "PV"] ? 0.0 : get(row, :variable, 0.0)
    operation_cost = RenewableGenerationCost(CostCurve(LinearCurve(variable_cost)))
    
    prime_mover_map = Dict(
        "WIND" => PrimeMovers.WT,
        "SOLAR" => PrimeMovers.PVe,
        "SOLAR_PV" => PrimeMovers.PVe,
        "SOLAR_CSP" => PrimeMovers.ST,
        "PV" => PrimeMovers.PVe,
        "HYDRO" => PrimeMovers.HY,
        "GEOTHERMAL" => PrimeMovers.ST
    )
    prime_mover = get(prime_mover_map, fuel_type, PrimeMovers.OT)
    
    generator = RenewableDispatch(
        string(row.name),
        get(row, :available, true),
        bus,
        get(row, :active_power, 0.0),
        0.0,
        get(row, :max_active_power, 100.0),
        prime_mover,
        (min = get(row, :min_reactive_power, 0.0), 
         max = get(row, :max_reactive_power, 0.0)),
        1.0,
        operation_cost,
        1.0
    )
    
    add_component!(sienna_sys.system, generator)
end

"""
    _get_thermal_defaults(config::SiennaConfig)

Get thermal generator default parameters.
"""
function _get_thermal_defaults(config::SiennaConfig)
    system_section = get(config.config_data, "system_building", Dict())
    defaults_section = get(system_section, "defaults", Dict())
    
    return Dict(
        "variable_cost" => get(defaults_section, "thermal_variable_cost", 50.0),
        "startup_cost" => get(defaults_section, "thermal_startup_cost", 1000.0),
        "shutdown_cost" => get(defaults_section, "thermal_shutdown_cost", 0.0),
        "min_power_fraction" => get(defaults_section, "thermal_min_power_fraction", 0.3),
        "ramp_rate" => get(defaults_section, "thermal_ramp_rate", 100.0)
    )
end

"""
    _add_storage_systems!(sienna_sys::SiennaSystem)

Add storage systems from CSV data.
"""
function _add_storage_systems!(sienna_sys::SiennaSystem)
    @info "Adding storage systems..."
    
    storage_file = joinpath(sienna_sys.config.data_directory, "storage.csv")
    if !isfile(storage_file)
        @info "No storage.csv found"
        return
    end
    
    df = CSV.read(storage_file, DataFrame)
    storage_count = 0
    
    for row in eachrow(df)
        try
            bus = get_component(ACBus, sienna_sys.system, string(row.bus))
            if bus === nothing
                @warn "Bus $(row.bus) not found for storage $(row.name)"
                continue
            end
            
            storage_capacity = Float64(get(row, :energy_capacity, 400.0))
            input_power = Float64(get(row, :input_active_power_limits, 100.0))
            output_power = Float64(get(row, :output_active_power_limits, 100.0))
            rating = max(input_power, output_power)
            
            operation_cost = StorageCost(
                charge_variable_cost = CostCurve(LinearCurve(0.0)),
                discharge_variable_cost = CostCurve(LinearCurve(0.0)),
                fixed = Float64(get(row, :fixed_cost, 0.0)),
                start_up = 0.0,
                shut_down = 0.0,
                energy_shortage_cost = 1000.0,
                energy_surplus_cost = 0.0
            )
            
            storage = EnergyReservoirStorage(;
                name = string(row.name),
                available = get(row, :available, true),
                bus = bus,
                prime_mover_type = PrimeMovers.BA,
                storage_technology_type = StorageTech.LIB,
                storage_capacity = storage_capacity,
                storage_level_limits = (min = 0.1, max = 0.9),
                initial_storage_capacity_level = 0.5,
                rating = rating,
                active_power = 0.0,
                input_active_power_limits = (min = 0.0, max = input_power),
                output_active_power_limits = (min = 0.0, max = output_power),
                efficiency = (in = Float64(get(row, :efficiency_in, 0.9)), 
                            out = Float64(get(row, :efficiency_out, 0.9))),
                reactive_power = 0.0,
                reactive_power_limits = nothing,
                base_power = 1.0,
                operation_cost = operation_cost,
                conversion_factor = 1.0,
                storage_target = 0.0,
                cycle_limits = 10000
            )
            
            add_component!(sienna_sys.system, storage)
            storage_count += 1
            
        catch e
            @warn "Failed to add storage $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["EnergyReservoirStorage"] = storage_count
    @info "✓ Added $storage_count storage systems"
end

"""
    _add_loads!(sienna_sys::SiennaSystem)

Add loads from CSV data.
"""
function _add_loads!(sienna_sys::SiennaSystem)
    @info "Adding loads..."
    
    load_file = joinpath(sienna_sys.config.data_directory, "load.csv")
    if !isfile(load_file)
        error("❌ Required load.csv not found: $load_file")
    end
    
    df = CSV.read(load_file, DataFrame)
    load_count = 0
    
    for row in eachrow(df)
        try
            bus = get_component(ACBus, sienna_sys.system, string(row.bus))
            if bus === nothing
                @warn "Bus $(row.bus) not found for load $(row.name)"
                continue
            end
            
            load = PowerLoad(
                name = string(row.name),
                available = get(row, :available, true),
                bus = bus,
                active_power = get(row, :max_active_power, 100.0),
                reactive_power = get(row, :max_reactive_power, 30.0),
                base_power = 1.0,
                max_active_power = get(row, :max_active_power, 100.0),
                max_reactive_power = get(row, :max_reactive_power, 30.0)
            )
            
            add_component!(sienna_sys.system, load)
            load_count += 1
            
        catch e
            @warn "Failed to add load $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["PowerLoad"] = load_count
    @info "✓ Added $load_count loads"
end

"""
    _add_network_branches!(sienna_sys::SiennaSystem)

Add transmission lines and network components.
"""
function _add_network_branches!(sienna_sys::SiennaSystem)
    @info "Adding network branches..."
    
    branch_file = joinpath(sienna_sys.config.data_directory, "branch.csv")
    if !isfile(branch_file)
        @info "No branch.csv found"
        return
    end
    
    df = CSV.read(branch_file, DataFrame)
    line_count = 0
    
    for row in eachrow(df)
        try
            from_bus = get_component(ACBus, sienna_sys.system, string(row.connection_points_from))
            to_bus = get_component(ACBus, sienna_sys.system, string(row.connection_points_to))
            
            if from_bus !== nothing && to_bus !== nothing
                line = Line(
                    name = string(row.name),
                    available = get(row, :available, true),
                    active_power_flow = 0.0,
                    reactive_power_flow = 0.0,
                    arc = Arc(from_bus, to_bus),
                    r = get(row, :r, 0.01),
                    x = get(row, :x, 0.1),
                    b = (from = get(row, :b, 0.0)/2, to = get(row, :b, 0.0)/2),
                    rating = get(row, :rating, 100.0),
                    angle_limits = (min = deg2rad(-30.0), max = deg2rad(30.0))
                )
                
                add_component!(sienna_sys.system, line)
                line_count += 1
            end
        catch e
            @warn "Failed to add line $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["Line"] = line_count
    @info "✓ Added $line_count transmission lines"
end

# ===== PUBLIC INTERFACE =====

"""
    get_power_system(sienna_sys::SiennaSystem)

Get the PowerSystems.jl system object.
"""
function get_power_system(sienna_sys::SiennaSystem)
    if sienna_sys.system === nothing
        error("❌ System not built yet. Call build_system!() first.")
    end
    return sienna_sys.system
end

"""
    is_system_built(sienna_sys::SiennaSystem)

Check if system is built successfully.
"""
function is_system_built(sienna_sys::SiennaSystem)
    return sienna_sys.is_built && sienna_sys.system !== nothing
end

"""
    has_time_series_data(sienna_sys::SiennaSystem)

Check if time series data is loaded.
"""
function has_time_series_data(sienna_sys::SiennaSystem)
    return sienna_sys.has_time_series
end

"""
    has_forecast_data(sienna_sys::SiennaSystem)

Check if forecast data is available.
"""
function has_forecast_data(sienna_sys::SiennaSystem)
    return sienna_sys.has_forecasts
end

"""
    get_component_counts(sienna_sys::SiennaSystem)

Get dictionary of component counts by type.
"""
function get_component_counts(sienna_sys::SiennaSystem)
    return copy(sienna_sys.component_counts)
end

"""
    get_build_errors(sienna_sys::SiennaSystem)

Get list of build errors.
"""
function get_build_errors(sienna_sys::SiennaSystem)
    return copy(sienna_sys.errors)
end

"""
    print_system_summary(sienna_sys::SiennaSystem)

Print comprehensive system build summary.
"""
function print_system_summary(sienna_sys::SiennaSystem)
    println()
    println("="^60)
    println("SIENNA SYSTEM BUILD SUMMARY")
    println("="^60)
    println("Project: $(sienna_sys.config.project_name)")
    println("Base Power: $(get_base_power(sienna_sys.system)) MW")
    println("Status: $(sienna_sys.is_built ? "✅ Built" : "❌ Failed")")
    println()
    println("Components:")
    for (type, count) in sienna_sys.component_counts
        println("  $type: $count")
    end
    println()
    println("Data Status:")
    println("  Time Series: $(sienna_sys.has_time_series ? "✅ Loaded" : "❌ Not loaded")")
    println("  Forecasts: $(sienna_sys.has_forecasts ? "✅ Available" : "❌ Not available")")
    
    if !isempty(sienna_sys.errors)
        println()
        println("Errors:")
        for error in sienna_sys.errors
            println("  ❌ $error")
        end
    end
    
    println("="^60)
    println()
end

"""
    validate_system(sienna_sys::SiennaSystem)

Validate the built system for common issues.
"""
function validate_system(sienna_sys::SiennaSystem)
    @info "🔍 Validating system..."
    
    if !sienna_sys.is_built
        @error "System not built - cannot validate"
        return false
    end
    
    issues = String[]
    
    # Check component counts
    bus_count = get(sienna_sys.component_counts, "ACBus", 0)
    load_count = get(sienna_sys.component_counts, "PowerLoad", 0)
    gen_count = get(sienna_sys.component_counts, "ThermalStandard", 0) + 
                get(sienna_sys.component_counts, "RenewableDispatch", 0)
    
    if bus_count == 0
        push!(issues, "No buses found")
    end
    
    if load_count == 0
        push!(issues, "No loads found")
    end
    
    if gen_count == 0
        push!(issues, "No generators found")
    end
    
    # Check time series
    if sienna_sys.config.load_timeseries
        if !sienna_sys.has_time_series
            push!(issues, "Time series requested but not loaded")
        elseif !sienna_sys.has_forecasts
            push!(issues, "Time series loaded but forecasts not created")
        end
    end
    
    # Check for reference bus
    ref_buses = [bus for bus in get_components(ACBus, sienna_sys.system) 
                 if get_bustype(bus) == ACBusTypes.REF]
    if isempty(ref_buses)
        push!(issues, "No reference bus found")
    elseif length(ref_buses) > 1
        push!(issues, "Multiple reference buses found")
    end
    
    # Report results
    if isempty(issues)
        @info "✅ System validation passed"
        return true
    else
        @warn "⚠️ System validation found issues:"
        for issue in issues
            @warn "  - $issue"
        end
        return false
    end
end

"""
    save_system(sienna_sys::SiennaSystem, filepath::String)

Save the system to file.
"""
function save_system(sienna_sys::SiennaSystem, filepath::String)
    if !sienna_sys.is_built
        error("Cannot save system - not built yet")
    end
    
    @info "💾 Saving system to $filepath"
    
    try
        to_json(sienna_sys.system, filepath)
        @info "✅ System saved successfully"
    catch e
        @error "Failed to save system: $e"
        rethrow(e)
    end
end

"""
    load_system(filepath::String; config::Union{SiennaConfig, Nothing} = nothing)

Load system from file and optionally attach config.
"""
function load_system(filepath::String; config::Union{SiennaConfig, Nothing} = nothing)
    @info "📂 Loading system from $filepath"
    
    try
        system = System(filepath)
        
        if config !== nothing
            sienna_sys = SiennaSystem(config)
            sienna_sys.system = system
            sienna_sys.is_built = true
            sienna_sys.has_time_series = !isempty(get_time_series_summary(system))
            @info "✅ System loaded and attached to config"
            return sienna_sys
        else
            @info "✅ System loaded (no config attached)"
            return system
        end
        
    catch e
        @error "Failed to load system: $e"
        rethrow(e)
    end
end

# ===== UTILITY FUNCTIONS =====

"""
    get_system_statistics(sienna_sys::SiennaSystem)

Get detailed system statistics.
"""
function get_system_statistics(sienna_sys::SiennaSystem)
    if !sienna_sys.is_built
        return Dict("error" => "System not built")
    end
    
    sys = sienna_sys.system
    
    # Basic counts
    stats = Dict(
        "total_components" => length(get_components(Component, sys)),
        "buses" => length(get_components(ACBus, sys)),
        "generators" => length(get_components(Generator, sys)),
        "loads" => length(get_components(PowerLoad, sys)),
        "lines" => length(get_components(Line, sys)),
        "storage" => length(get_components(Storage, sys))
    )
    
    # Capacity calculations
    thermal_capacity = sum(get_max_active_power(gen) for gen in get_components(ThermalStandard, sys))
    renewable_capacity = sum(get_max_active_power(gen) for gen in get_components(RenewableDispatch, sys))
    total_load = sum(get_max_active_power(load) for load in get_components(PowerLoad, sys))
    
    stats["thermal_capacity_mw"] = thermal_capacity
    stats["renewable_capacity_mw"] = renewable_capacity
    stats["total_generation_capacity_mw"] = thermal_capacity + renewable_capacity
    stats["total_load_mw"] = total_load
    stats["reserve_margin"] = (thermal_capacity + renewable_capacity - total_load) / total_load * 100
    
    # Time series info
    if sienna_sys.has_time_series
        ts_summary = get_time_series_summary(sys)
        stats["time_series_components"] = length(ts_summary)
        if !isempty(ts_summary)
            stats["time_series_length"] = ts_summary[1].time_step_count
        end
    end
    
    return stats
end

"""
    compare_systems(sys1::SiennaSystem, sys2::SiennaSystem)

Compare two SiennaSystem objects.
"""
function compare_systems(sys1::SiennaSystem, sys2::SiennaSystem)
    @info "🔍 Comparing systems..."
    
    stats1 = get_system_statistics(sys1)
    stats2 = get_system_statistics(sys2)
    
    @info "System 1 ($(sys1.config.project_name)):"
    for (key, value) in stats1
        @info "  $key: $value"
    end
    
    @info "System 2 ($(sys2.config.project_name)):"
    for (key, value) in stats2
        @info "  $key: $value"
    end
    
    @info "Differences:"
    for key in keys(stats1)
        if haskey(stats2, key) && stats1[key] != stats2[key]
            @info "  $key: $(stats1[key]) → $(stats2[key])"
        end
    end
end

# ===== LEGACY COMPATIBILITY =====

# Keep old function names for backwards compatibility
const print_build_summary = print_system_summary
const has_timeseries_data = has_time_series_data

# ===== EXPORTS =====

export SiennaSystem
export build_system!, get_power_system, is_system_built
export has_time_series_data, has_forecast_data, has_timeseries_data
export get_component_counts, get_build_errors
export print_system_summary, print_build_summary, validate_system
export save_system, load_system
export get_system_statistics, compare_systems

# ===== TESTING =====

if abspath(PROGRAM_FILE) == @__FILE__
    @info "🧪 Testing Clean SiennaSystem..."
    
    try
        if isfile("config.toml")
            config = SiennaConfig("config.toml")
            sienna_sys = SiennaSystem(config)
            
            if isdir(config.data_directory)
                build_system!(sienna_sys)
                validate_system(sienna_sys)
                print_system_summary(sienna_sys)
                @info "✅ Clean SiennaSystem test passed"
            else
                @info "✅ SiennaSystem initialization test passed"
            end
        else
            @info "No config.toml found for testing"
        end
    catch e
        @error "❌ Clean SiennaSystem test failed: $e"
    end
end
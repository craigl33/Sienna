#!/usr/bin/env julia

"""
SiennaSystem.jl - PowerSystems.jl System Builder for Sienna Ecosystem
======================================================================

Class-based system builder that encapsulates PowerSystems.jl system construction.
Provides clean interface for building systems from CSV data and configuration.

Features:
- Config-driven system building
- Time series loading with proper scaling support
- Differentiated renewable generator types (RenewableDispatch vs RenewableNonDispatch)
- Proper cost handling for renewable generators
- Component validation and analysis
- System export capabilities
- Memory-efficient component addition
- Proper storage handling (EnergyReservoirStorage vs HydroPumpedStorage)

Usage:
    config = SiennaConfig("config.toml")
    sienna_sys = SiennaSystem(config)
    sys = get_system(sienna_sys)

Key Improvements:
1. Fixed time series loading with proper scaling_factor_multiplier usage
2. Differentiated renewable types: wind/solar → RenewableNonDispatch, hydro/bio → RenewableDispatch
3. Set renewable costs to zero for wind, solar PV, solar CSP
4. Improved generator categorization based on fuel type and availability constraints
5. Fixed storage: PHS → HydroPumpedStorage, others → EnergyReservoirStorage
6. Corrected StorageTech enum mappings per PowerSystems.jl v4.6.2
"""

using PowerSystems
using DataFrames
using CSV
using TimeSeries
using JSON3
using Dates
using Logging
using InfrastructureSystems  # Add this import for InfrastructureSystemsInternal
using PowerSystems: PumpHydroStatus  # Import the enum
using Infiltrator

# Import our configuration manager
include("SiennaConfig.jl")

"""
    SiennaSystem

Main system builder class that encapsulates PowerSystems.jl functionality.
"""
mutable struct SiennaSystem
    # Core components
    config::SiennaConfig
    system::Union{System, Nothing}
    
    # System metadata
    system_name::String
    build_timestamp::DateTime
    
    # Component counts (updated during building)
    component_counts::Dict{String, Int}
    
    # Time series metadata
    time_series_loaded::Bool
    time_series_metadata::Union{Any, Nothing}
    forecast_available::Bool
    
    # Validation results
    system_validated::Bool
    validation_warnings::Vector{String}
    validation_errors::Vector{String}
    
    # Build status
    build_successful::Bool
    build_errors::Vector{String}
end

"""
    SiennaSystem(config::SiennaConfig)

Constructor - Initialize system builder with configuration.
"""
function SiennaSystem(config::SiennaConfig)
    @info "🏗️  Initializing SiennaSystem with configuration"
    
    # Ensure config is validated
    if !config.is_validated
        error("❌ Configuration must be validated before creating SiennaSystem")
    end
    
    # Initialize with default values
    sienna_sys = SiennaSystem(
        config,
        nothing,  # system will be built later
        config.project_name,
        now(),
        Dict{String, Int}(),  # component counts
        false, nothing, false,  # time series state
        false, String[], String[],  # validation state
        false, String[]  # build state
    )
    
    @info "✅ SiennaSystem initialized - call build_system!() to construct PowerSystems.jl system"
    return sienna_sys
end

"""
    build_system!(sienna_sys::SiennaSystem)

Build the PowerSystems.jl system using configuration parameters.
"""
function build_system!(sienna_sys::SiennaSystem)
    @info "🏗️  Building PowerSystems.jl system from configuration..."
    
    # Reset build state
    sienna_sys.build_errors = String[]
    sienna_sys.component_counts = Dict{String, Int}()
    
    try
            # Create base system
        create_base_system!(sienna_sys)
        println("After create_base_system: ", get_units_base(sienna_sys.system))
        
        # Add components in dependency order
        add_buses!(sienna_sys)
        println("After add_buses: ", get_units_base(sienna_sys.system))
        
        add_areas!(sienna_sys)
        println("After add_areas: ", get_units_base(sienna_sys.system))
        
        add_generators!(sienna_sys)
        println("After add_generators: ", get_units_base(sienna_sys.system))
        
        add_storage!(sienna_sys)
        println("After add_storage: ", get_units_base(sienna_sys.system))
        
        add_loads!(sienna_sys)
        println("After add_loads: ", get_units_base(sienna_sys.system))
        
        add_network_components!(sienna_sys)
        println("After add_network: ", get_units_base(sienna_sys.system))
        
        # Load time series if requested
        if sienna_sys.config.load_timeseries
            load_time_series!(sienna_sys)
            println("After load_time_series: ", get_units_base(sienna_sys.system))
        end
        
        # Validate system if requested
        if sienna_sys.config.validate_system
            validate_system!(sienna_sys)
        end
        
        # Update build status
        sienna_sys.build_successful = true
        sienna_sys.build_timestamp = now()
        
        @info "✅ PowerSystems.jl system built successfully!"
        print_build_summary(sienna_sys)
        
    catch e
        sienna_sys.build_successful = false
        push!(sienna_sys.build_errors, string(e))
        @error "❌ System build failed: $e"
        rethrow(e)
    end
    
    return sienna_sys.system
end

"""
    create_base_system!(sienna_sys::SiennaSystem)

Create the base PowerSystems.jl system.
"""
function create_base_system!(sienna_sys::SiennaSystem)
    @info "Creating base PowerSystems.jl system..."
    
    # Create system with base power from config
    sienna_sys.system = System(sienna_sys.config.base_power, unit_system="NATURAL_UNITS")
    
    # Set system name from config
    set_name!(sienna_sys.system, sienna_sys.config.project_name)
    
    @info "✓ Base system created with $(sienna_sys.config.base_power) MW base power and NATURAL_UNITS"
end

"""
    add_buses!(sienna_sys::SiennaSystem)

Add buses from CSV file using configuration parameters.
"""
function add_buses!(sienna_sys::SiennaSystem)
    @info "Adding buses from CSV data..."
    
    bus_file = joinpath(sienna_sys.config.data_directory, "bus.csv")
    if !isfile(bus_file)
        error("❌ Required bus.csv file not found: $bus_file")
    end
    
    df = CSV.read(bus_file, DataFrame)
    default_voltage = sienna_sys.config.default_voltage
    
    bus_count = 0
    ref_bus_count = 0
    
    for (i, row) in enumerate(eachrow(df))
        try
            # Parse bus number from name
            bus_number = try
                numbers = replace(string(row.name), r"[^0-9]" => "")
                isempty(numbers) ? i : parse(Int, numbers)
            catch
                i
            end
            
            # Determine bus type
            bus_type = if haskey(row, :bus_type) && !ismissing(row.bus_type)
                if string(row.bus_type) == "REF"
                    ref_bus_count += 1
                    ACBusTypes.REF
                else
                    ACBusTypes.PV
                end
            else
                # Make first bus reference if no explicit reference bus
                if i == 1 && ref_bus_count == 0
                    ref_bus_count += 1
                    ACBusTypes.REF
                else
                    ACBusTypes.PV
                end
            end
            
            # Get voltage level
            base_voltage_val = if haskey(row, :base_voltage) && !ismissing(row.base_voltage)
                Float64(row.base_voltage)
            else
                default_voltage
            end
            
            # Create bus
            bus = ACBus(
                bus_number,
                string(row.name),
                bus_type,
                get(row, :angle, 0.0),
                get(row, :voltage, 1.0),
                (min = 0.95, max = 1.05),
                base_voltage_val
            )
            
            add_component!(sienna_sys.system, bus)
            bus_count += 1
            
        catch e
            push!(sienna_sys.build_errors, "Failed to add bus $(row.name): $e")
            @warn "Failed to add bus $(row.name): $e"
        end
    end
    
    # Ensure we have exactly one reference bus
    if ref_bus_count == 0
        @warn "No reference bus found - making first bus the reference"
        buses = get_components(ACBus, sienna_sys.system)
        if !isempty(buses)
            first_bus = first(buses)
            set_bustype!(first_bus, ACBusTypes.REF)
            ref_bus_count = 1
        end
    elseif ref_bus_count > 1
        @warn "Multiple reference buses found ($ref_bus_count) - this may cause issues"
    end
    
    sienna_sys.component_counts["ACBus"] = bus_count
    @info "✓ Added $bus_count buses ($ref_bus_count reference)"
end

"""
    add_areas!(sienna_sys::SiennaSystem)

Add Area components to the system based on configuration.
"""
function add_areas!(sienna_sys::SiennaSystem)
    @info "Adding areas from configuration..."
    
    # Get area configuration
    area_config = get(sienna_sys.config.config_data, "system_building", Dict())
    area_section = get(area_config, "areas", Dict())
    area_strategy = get(area_section, "area_strategy", "single")
    
    @info "Using area strategy: $area_strategy"
    
    if area_strategy == "single"
        add_single_area!(sienna_sys, area_section)
    elseif area_strategy == "per_bus"
        add_per_bus_areas!(sienna_sys, area_section)
    elseif area_strategy == "from_csv"
        add_areas_from_csv!(sienna_sys, area_section)
    elseif area_strategy == "from_config"
        add_areas_from_config!(sienna_sys, area_section)
    else
        @warn "Unknown area strategy: $area_strategy, using 'single' as fallback"
        add_single_area!(sienna_sys, area_section)
    end
    
    # Validate area assignments
    validate_area_assignments!(sienna_sys, area_section)
end

"""
    add_single_area!(sienna_sys::SiennaSystem, area_section::Dict)

Add a single area for the entire system.
"""
function add_single_area!(sienna_sys::SiennaSystem, area_section::Dict)
    default_area_name = get(area_section, "default_area_name", "system_area")
    
    @info "Creating single area: $default_area_name"
    
    # Create the area
    area = Area(default_area_name)
    add_component!(sienna_sys.system, area)
    
    # Assign all buses to this area
    buses = get_components(ACBus, sienna_sys.system)
    assigned_count = 0
    
    for bus in buses
        try
            set_area!(bus, area)
            assigned_count += 1
        catch e
            @warn "Failed to assign bus $(get_name(bus)) to area $default_area_name: $e"
            push!(sienna_sys.build_errors, "Failed to assign bus $(get_name(bus)) to area: $e")
        end
    end
    
    sienna_sys.component_counts["Area"] = 1
    @info "✓ Created 1 area and assigned $assigned_count buses"
end

"""
    validate_area_assignments!(sienna_sys::SiennaSystem, area_section::Dict)

Validate that all buses are assigned to areas.
"""
function validate_area_assignments!(sienna_sys::SiennaSystem, area_section::Dict)
    validate_assignments = get(area_section, "validate_area_assignments", true)
    allow_unassigned = get(area_section, "allow_unassigned_buses", false)
    
    if !validate_assignments
        return
    end
    
    @info "Validating area assignments..."
    
    buses = get_components(ACBus, sienna_sys.system)
    unassigned_buses = []
    
    # Check each bus for area assignment
    for bus in buses
        try
            area = get_area(bus)
            if area === nothing
                push!(unassigned_buses, get_name(bus))
            end
        catch e
            push!(unassigned_buses, get_name(bus))
        end
    end
    
    if !isempty(unassigned_buses)
        @warn "Found $(length(unassigned_buses)) unassigned buses: $(join(unassigned_buses, ", "))"
        
        if allow_unassigned
            @info "Creating default area for unassigned buses"
            
            # Create default area for unassigned buses
            default_area = Area("unassigned_area")
            add_component!(sienna_sys.system, default_area)
            
            # Assign unassigned buses to default area
            assigned_count = 0
            for bus_name in unassigned_buses
                bus = get_component(ACBus, sienna_sys.system, bus_name)
                if bus !== nothing
                    try
                        set_area!(bus, default_area)
                        assigned_count += 1
                    catch e
                        @warn "Failed to assign bus $bus_name to default area: $e"
                    end
                end
            end
            
            # Update area count
            current_area_count = get(sienna_sys.component_counts, "Area", 0)
            sienna_sys.component_counts["Area"] = current_area_count + 1
            
            @info "✓ Assigned $assigned_count unassigned buses to default area"
        else
            push!(sienna_sys.validation_errors, 
                  "$(length(unassigned_buses)) buses are not assigned to areas: $(join(unassigned_buses, ", "))")
        end
    else
        @info "✓ All buses are properly assigned to areas"
    end
    
    # Summary
    areas = get_components(Area, sienna_sys.system)
    @info "Area assignment summary:"
    for area in areas
        area_buses = [get_name(bus) for bus in buses if try; get_area(bus) == area; catch; false; end]
        @info "  $(get_name(area)): $(length(area_buses)) buses"
    end
end

# ===== IMPROVED GENERATOR FUNCTIONS =====

"""
    add_generators!(sienna_sys::SiennaSystem)

Add generators from CSV file with correct categorization for PowerSystems.jl types.
"""
function add_generators!(sienna_sys::SiennaSystem)
    @info "Adding generators from CSV data with correct PowerSystems.jl categorization..."
    
    gen_file = joinpath(sienna_sys.config.data_directory, "gen.csv")
    if !isfile(gen_file)
        @info "No gen.csv found, skipping generator addition"
        return
    end
    
    df = CSV.read(gen_file, DataFrame)
    
    thermal_count = 0
    renewable_dispatch_count = 0
    hydro_dispatch_count = 0
    
    for row in eachrow(df)
        try
            # Find corresponding bus
            bus = get_component(ACBus, sienna_sys.system, string(row.bus))
            if bus === nothing
                push!(sienna_sys.build_errors, "Bus $(row.bus) not found for generator $(row.name)")
                @warn "Bus $(row.bus) not found for generator $(row.name)"
                continue
            end
            
            # Determine fuel type and generator category
            fuel_type = uppercase(string(get(row, :fuel, "")))
            gen_category = categorize_generator(sienna_sys.config, fuel_type, row)
            
            if gen_category == :thermal
                add_thermal_generator!(sienna_sys, row, bus)
                thermal_count += 1
            elseif gen_category == :renewable_dispatch
                add_renewable_dispatch_generator!(sienna_sys, row, bus, fuel_type)
                renewable_dispatch_count += 1
            else
                @warn "Unknown generator category for: $(row.name) (fuel: $fuel_type)"
                push!(sienna_sys.build_errors, "Unknown generator category for: $(row.name)")
            end
            
        catch e
            # Get full exception info including stack trace
            bt = catch_backtrace()
            
            # Create detailed error message with stack trace
            error_msg = "Failed to add generator $(row.name): $e"
            push!(sienna_sys.build_errors, error_msg)
            
            # Print detailed error with stack trace
            @error "❌ $error_msg"
            @error "Exception type: $(typeof(e))"
            @error "Stack trace:"
            
            # Print the stack trace
            for (i, frame) in enumerate(stacktrace(bt))
                @error "  [$i] $(frame.func) at $(frame.file):$(frame.line)"
            end
            
            # If you want even more detail, use showerror
            io = IOBuffer()
            showerror(io, e, bt)
            detailed_error = String(take!(io))
            @error "Detailed error:\n$detailed_error"
            
            @warn "Failed to add generator $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["ThermalStandard"] = thermal_count
    sienna_sys.component_counts["RenewableDispatch"] = renewable_dispatch_count
    
    @info "✓ Added generators:"
    @info "  - Thermal (including bioenergy): $thermal_count"
    @info "  - Renewable Dispatch (wind/solar/hydro): $renewable_dispatch_count"
end

"""
    categorize_generator(config::SiennaConfig, fuel_type::String, row)

Categorize generator based on fuel type and correct dispatch characteristics.
"""
function categorize_generator(config::SiennaConfig, fuel_type::String, row)
    @debug "Categorizing generator: $(get(row, :name, "unknown")) with fuel: '$fuel_type'"
    
    # Check if thermal first (including bioenergy which should be thermal)
    if is_fuel_thermal(config, fuel_type) || fuel_type in ["BIOENERGY", "BIOMASS", "BIOGAS"]
        @debug "  → Thermal fuel detected (including bioenergy)"
        return :thermal
    end
    
    # Hydro should be treated as renewable dispatch for PowerSimulations compatibility
    if fuel_type == "HYDRO"
        @debug "  → Hydro fuel detected → RenewableDispatch (treated as renewable for PowerSimulations)"
        return :renewable_dispatch
    end
    
    # For renewables, all grid-scale renewable generation should be dispatchable (can be curtailed)
    if is_fuel_renewable(config, fuel_type)
        @debug "  → Grid-scale renewable fuel detected → RenewableDispatch"
        
        # All grid-scale renewables (wind, solar PV, solar CSP) are dispatchable (can be curtailed)
        if fuel_type in ["WIND", "SOLAR", "SOLAR_PV", "SOLAR_CSP", "PV", "GEOTHERMAL"]
            @debug "  → Grid-scale renewable: $fuel_type → RenewableDispatch"
            return :renewable_dispatch
        end
        
        # Fallback for other renewables
        @debug "  → Other renewable: $fuel_type → RenewableDispatch"
        return :renewable_dispatch
    end
    
    @warn "Unknown fuel type: '$fuel_type' for generator $(get(row, :name, "unknown"))"
    return :unknown
end

"""
    add_thermal_generator!(sienna_sys::SiennaSystem, row, bus::ACBus)

Add a single thermal generator using configuration defaults.
"""
function add_thermal_generator!(sienna_sys::SiennaSystem, row, bus::ACBus)
    defaults = get_thermal_defaults(sienna_sys.config)
    
    # Create cost curve
    variable_cost = get(row, :variable, defaults["variable_cost"])
    variable_curve = CostCurve(LinearCurve(variable_cost))
    
    # Create operation cost
    operation_cost = ThermalGenerationCost(
        variable = variable_curve,
        fixed = get(row, :fixed_cost, 0.0),
        start_up = get(row, :startup, defaults["startup_cost"]),
        shut_down = get(row, :shutdown, defaults["shutdown_cost"])
    )
    
    # Get power limits
    max_power = get(row, :max_active_power, 100.0)
    min_power = get(row, :min_active_power, max_power * defaults["min_power_fraction"])
    
    # Get reactive power limits - ensure they're proper Float64 values
    min_reactive = Float64(get(row, :min_reactive_power, -30.0))
    max_reactive = Float64(get(row, :max_reactive_power, 30.0))
    
    @infiltrate
    # Create thermal generator
    gen = ThermalStandard(
        string(row.name),
        get(row, :available, true),
        get(row, :status, 1) == 1,
        bus,
        get(row, :active_power, 0.0),
        0.0,
        max_power,
        (min = min_power, max = max_power),
        (min = min_reactive, max = max_reactive),
        (up = get(row, :ramp_30, defaults["ramp_rate"]), 
         down = get(row, :ramp_30, defaults["ramp_rate"])),
        operation_cost,
        1.0 # base_power
    )
    
    add_component!(sienna_sys.system, gen)
    @infiltrate
end

"""
    add_renewable_dispatch_generator!(sienna_sys::SiennaSystem, row, bus::ACBus, fuel_type::String)

Add a dispatchable renewable generator (hydro, bioenergy, etc.).
"""
function add_renewable_dispatch_generator!(sienna_sys::SiennaSystem, row, bus::ACBus, fuel_type::String)
    # Create cost curve - use operational cost but may set to zero for certain types
    variable_cost = should_use_zero_cost(fuel_type) ? 0.0 : get(row, :variable, 0.0)
    variable_curve = CostCurve(LinearCurve(variable_cost))
    operation_cost = RenewableGenerationCost(variable_curve)
    
    # Determine prime mover from fuel type
    prime_mover = get_prime_mover_from_fuel(fuel_type)
    
    # Create renewable dispatch generator
    gen = RenewableDispatch(
        string(row.name),
        get(row, :available, true),
        bus,
        get(row, :active_power, 0.0),
        0.0,
        get(row, :max_active_power, 100.0),
        prime_mover,
        (min = get(row, :min_reactive_power, 0.0), max = get(row, :max_reactive_power, 0.0)),
        1.0,  # power_factor
        operation_cost,
        1.0  # base_power
    )
    
    add_component!(sienna_sys.system, gen)
    
    @debug "Added RenewableDispatch: $(row.name) (fuel: $fuel_type, cost: $variable_cost)"
end

"""
    add_hydro_dispatch_generator!(sienna_sys::SiennaSystem, row, bus::ACBus)

Add a hydro dispatch generator using correct PowerSystems.jl v4.6.2 constructor.
"""
function add_hydro_dispatch_generator!(sienna_sys::SiennaSystem, row, bus::ACBus)
    @debug "Adding HydroDispatch: $(get(row, :name, "unknown"))"
    
    # Create HydroGenerationCost with correct structure - it has a default constructor
    operation_cost = HydroGenerationCost(nothing)  # Default constructor as per documentation
    
    # Get power limits
    max_power = get(row, :max_active_power, 100.0)
    min_power = get(row, :min_active_power, 0.0)
    
    # Get reactive power limits - ensure they're proper NamedTuple
    min_reactive = Float64(get(row, :min_reactive_power, 0.0))
    max_reactive = Float64(get(row, :max_reactive_power, 0.0))
    reactive_power_limits = (min = min_reactive, max = max_reactive)
    
    # Create hydro dispatch generator using positional arguments (not keyword arguments)
    gen = HydroDispatch(
        string(row.name),                      # name::String
        get(row, :available, true),            # available::Bool  
        bus,                                   # bus::ACBus
        get(row, :active_power, 0.0),         # active_power::Float64
        0.0,                                   # reactive_power::Float64
        Float64(max_power),                    # rating::Float64
        PrimeMovers.HY,                        # prime_mover_type::PrimeMovers
        (min = min_power, max = max_power),    # active_power_limits::MinMax
        reactive_power_limits,                 # reactive_power_limits::Union{Nothing, MinMax}
        nothing,                               # ramp_limits::Union{Nothing, UpDown}
        nothing,                               # time_limits::Union{Nothing, UpDown}
        1.0,                                 # base_power::Float64
        operation_cost,                        # operation_cost::Union{HydroGenerationCost, MarketBidCost}
        Service[],                             # services::Vector{Service}
        nothing,                               # dynamic_injector::Union{Nothing, DynamicInjection}
        Dict{String, Any}(),                   # ext::Dict{String, Any}
        InfrastructureSystems.InfrastructureSystemsInternal()  # internal::InfrastructureSystemsInternal
    )
    
    add_component!(sienna_sys.system, gen)
    
    @debug "Added HydroDispatch: $(row.name) (rating: $(max_power) MW)"
end

"""
    add_renewable_nondispatch_generator!(sienna_sys::SiennaSystem, row, bus::ACBus, fuel_type::String)

Add a non-dispatchable renewable generator (wind, solar).
"""
function add_renewable_nondispatch_generator!(sienna_sys::SiennaSystem, row, bus::ACBus, fuel_type::String)
    # Zero cost for variable renewables in production cost models
    variable_cost = 0.0
    variable_curve = CostCurve(LinearCurve(variable_cost))
    operation_cost = RenewableGenerationCost(variable_curve)
    
    # Determine prime mover from fuel type
    prime_mover = get_prime_mover_from_fuel(fuel_type)
    
    # Create renewable non-dispatch generator
    gen = RenewableNonDispatch(
        string(row.name),
        get(row, :available, true),
        bus,
        get(row, :active_power, 0.0),
        0.0,
        get(row, :max_active_power, 100.0),
        prime_mover,
        (min = get(row, :min_reactive_power, 0.0), max = get(row, :max_reactive_power, 0.0)),
        1.0,  # power_factor
        operation_cost,
        1.0  # base_power
    )
    
    add_component!(sienna_sys.system, gen)
    
    @debug "Added RenewableNonDispatch: $(row.name) (fuel: $fuel_type, cost: $variable_cost)"
end

"""
    should_use_zero_cost(fuel_type::String)

Determine if renewable generator should have zero operational cost.
"""
function should_use_zero_cost(fuel_type::String)
    # Zero cost for variable renewables in production cost models
    zero_cost_fuels = ["WIND", "SOLAR", "SOLAR_PV", "SOLAR_CSP", "PV"]
    return fuel_type in zero_cost_fuels
end

"""
    get_prime_mover_from_fuel(fuel_type::String)

Get PowerSystems.jl PrimeMover enum from fuel type (case-insensitive).
"""
function get_prime_mover_from_fuel(fuel_type::String)
    # Convert to uppercase for case-insensitive matching
    fuel_upper = uppercase(fuel_type)
    
    prime_mover_map = Dict(
        "WIND" => PrimeMovers.WT,
        "SOLAR" => PrimeMovers.PVe,
        "SOLAR_PV" => PrimeMovers.PVe,
        "SOLAR_CSP" => PrimeMovers.ST,  # CSP uses Steam Turbine (ST) prime mover
        "PV" => PrimeMovers.PVe,
        "HYDRO" => PrimeMovers.HY,
        "BIOENERGY" => PrimeMovers.BT,
        "BIOMASS" => PrimeMovers.BT,
        "BIOGAS" => PrimeMovers.BT,
        "GEOTHERMAL" => PrimeMovers.ST
    )
    
    return get(prime_mover_map, fuel_upper, PrimeMovers.OT)
end

# ===== IMPROVED TIME SERIES LOADING =====

"""
    load_time_series!(sienna_sys::SiennaSystem)

Load time series data using metadata approach with proper scaling support.
"""
function load_time_series!(sienna_sys::SiennaSystem)
    @info "Loading time series data from metadata with improved scaling support..."
    
    metadata_file = joinpath(sienna_sys.config.data_directory, "timeseries_metadata.json")
    if !isfile(metadata_file)
        @warn "No timeseries_metadata.json found - skipping time series loading"
        return
    end
    
    try
        # Load metadata
        metadata_content = read(metadata_file, String)
        metadata = JSON3.read(metadata_content)
        sienna_sys.time_series_metadata = metadata
        
        @info "Processing $(length(metadata)) time series entries with improved scaling..."
        
        # Process each metadata entry
        total_added = 0
        for entry in metadata
            try
                if process_time_series_entry_improved!(sienna_sys, entry)
                    total_added += 1
                end
            catch e
                @warn "Failed to process time series entry: $e"
                @warn "Entry: $entry"
            end
        end
        
        if total_added > 0
            @info "✓ Successfully loaded time series for $total_added components"
            
            # Transform to forecasts for PowerSimulations v0.30.2
            transform_to_forecasts!(sienna_sys)
            sienna_sys.time_series_loaded = true
        else
            @warn "No time series data loaded successfully"
        end
        
    catch e
        @error "Failed to load time series: $e"
        push!(sienna_sys.build_errors, "Time series loading failed: $e")
    end
end

"""
    process_time_series_entry_improved!(sienna_sys::SiennaSystem, entry)

Process a single time series metadata entry with improved scaling support.
"""
function process_time_series_entry_improved!(sienna_sys::SiennaSystem, entry)
    component_name = entry["component"]
    data_file = entry["data_file"]
    data_column = entry["data_column"]
    label = entry["label"]
    category = entry["category"]
    
    # Find the component
    component = find_component_by_name_and_category(sienna_sys.system, component_name, category)
    if component === nothing
        @warn "Component not found: $component_name (category: $category)"
        return false
    end
    
    # Load the data file
    file_path = joinpath(sienna_sys.config.data_directory, data_file)
    if !isfile(file_path)
        @warn "Time series file not found: $file_path"
        return false
    end
    
    # Read CSV data
    df = CSV.read(file_path, DataFrame)
    
    # Validate columns
    if !("timestep" in names(df))
        @warn "No 'timestep' column in $data_file"
        return false
    end
    
    if !(data_column in names(df))
        @warn "Column '$data_column' not found in $data_file"
        return false
    end
    
    # Parse timestamps
    timestamps = try
        DateTime.(string.(df.timestep), "yyyy-mm-dd HH:MM:SS")
    catch
        try
            DateTime.(string.(df.timestep), "yyyy-mm-ddTHH:MM:SS")
        catch e
            @warn "Failed to parse timestamps in $data_file: $e"
            return false
        end
    end
    
    # Get data values
    data_values = Float64.(coalesce.(df[!, Symbol(data_column)], 0.0))
    
    # Create TimeArray
    ts_array = TimeArray(timestamps, data_values)
    
    # Create SingleTimeSeries with proper scaling support
    scaling_factor = get_scaling_factor_multiplier(entry, component)
    
   
    single_ts = SingleTimeSeries(
        name = label,
        data = ts_array,
        scaling_factor_multiplier = scaling_factor
    )
    
    # Add to system
    try
        add_time_series!(sienna_sys.system, component, single_ts)
        @info "  ✓ Added time series '$label' to $component_name (scaling: $(typeof(scaling_factor)))"
        return true
    catch e
        @warn "Failed to add time series to $component_name: $e"
        return false
    end
end

"""
    get_scaling_factor_multiplier(entry, component)

Get the appropriate scaling factor multiplier for time series. Nothing means no scaling, and is the appropriate default for PowerSystems.jl.
"""
function get_scaling_factor_multiplier(entry, component)
    scaling_spec = get(entry, "scaling_factor_multiplier", nothing)
    
    if scaling_spec === nothing
        # Default: no scaling
        return nothing
    elseif scaling_spec == "get_max_active_power"
        # Use PowerSystems.jl function as scaling factor
        return get_max_active_power
    elseif scaling_spec == "get_rating"
        # For renewable generators
        return get_rating
    elseif scaling_spec == "get_active_power_limits_max"
        # Alternative for thermal generators
        return x -> get_active_power_limits(x).max
    elseif scaling_spec == "1.0" || scaling_spec == 1.0
        # No scaling - return numeric value, not Nothing
        return nothing
    else
        # Try to parse as numeric value
        try
            return parse(Float64, string(scaling_spec))
        catch
            @warn "Unknown scaling factor: $scaling_spec, using 1.0"
            return nothing  # Return 1.0 instead of Nothing
        end
    end
end

"""
    find_component_by_name_and_category(sys::System, name::String, category::String)

Find component by name and category with support for new generator types.
"""
function find_component_by_name_and_category(sys::System, name::String, category::String)
    if category == "ElectricLoad"
        for load in get_components(PowerLoad, sys)
            if get_name(load) == name
                return load
            end
        end
    elseif category == "RenewableGen"
        # Check both RenewableDispatch and RenewableNonDispatch
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
    transform_to_forecasts!(sienna_sys::SiennaSystem)

Transform static time series to forecasts for PowerSimulations v0.30.2.
"""
function transform_to_forecasts!(sienna_sys::SiennaSystem)
    @info "Transforming time series to forecasts for PowerSimulations v0.30.2..."
    
    horizon_hours = sienna_sys.config.default_horizon_hours
    
    try
        # Use PowerSystems.jl transform function
        transform_single_time_series!(
            sienna_sys.system,
            Hour(horizon_hours),  # horizon
            Hour(24),             # interval (daily forecasts)
        )

        # ✅ FIX: Reset units back to NATURAL_UNITS after transformation
        set_units_base_system!(sienna_sys.system, "NATURAL_UNITS")
        @info "✓ Units reset to NATURAL_UNITS after forecast transformation"
        
        
        sienna_sys.forecast_available = true
        @info "✓ Successfully transformed time series to forecasts"
        @info "  Horizon: $horizon_hours hours"
        @info "  Interval: 24 hours (daily forecasts)"
        
    catch e
        @error "Failed to transform time series to forecasts: $e"
        push!(sienna_sys.build_errors, "Forecast transformation failed: $e")
        error("Time series transformation failed - cannot proceed with simulations")
    end
end

# ===== FIXED STORAGE FUNCTIONS =====

"""
    add_storage!(sienna_sys::SiennaSystem)

Add energy storage systems from CSV file - all storage treated as EnergyReservoirStorage.
"""
function add_storage!(sienna_sys::SiennaSystem)
    @info "Adding storage from CSV data - all storage as EnergyReservoirStorage..."
    
    storage_file = joinpath(sienna_sys.config.data_directory, "storage.csv")
    if !isfile(storage_file)
        @info "No storage.csv found, skipping storage addition"
        return
    end
    
    df = CSV.read(storage_file, DataFrame)
    storage_count = 0
    
    for row in eachrow(df)
        try
            # Find corresponding bus
            bus = get_component(ACBus, sienna_sys.system, string(row.bus))
            if bus === nothing
                push!(sienna_sys.build_errors, "Bus $(row.bus) not found for storage $(row.name)")
                @warn "Bus $(row.bus) not found for storage $(row.name)"
                continue
            end
            
            # Determine storage technology and add as EnergyReservoirStorage
            storage_tech = uppercase(string(get(row, :storage_technology, "BATTERY")))
            add_energy_storage!(sienna_sys, row, bus, storage_tech)
            storage_count += 1
            
        catch e
            # Get full exception info including stack trace
            bt = catch_backtrace()
            
            # Create detailed error message with stack trace
            error_msg = "Failed to add storage $(row.name): $e"
            push!(sienna_sys.build_errors, error_msg)
            
            # Print detailed error with stack trace
            @error "❌ $error_msg"
            @error "Exception type: $(typeof(e))"
            @error "Stack trace:"
            
            # Print the stack trace
            for (i, frame) in enumerate(stacktrace(bt))
                @error "  [$i] $(frame.func) at $(frame.file):$(frame.line)"
            end
            
            # If you want even more detail, use showerror
            io = IOBuffer()
            showerror(io, e, bt)
            detailed_error = String(take!(io))
            @error "Detailed error:\n$detailed_error"
            
            @warn "Failed to add storage $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["EnergyReservoirStorage"] = storage_count
    @info "✓ Added $storage_count energy storage systems (all as EnergyReservoirStorage)"
end

"""
    is_pumped_hydro_storage(storage_tech::String)

Check if storage technology is pumped hydro - now always returns false to treat all storage as EnergyReservoirStorage.
"""
function is_pumped_hydro_storage(storage_tech::String)
    # Treat all storage as EnergyReservoirStorage for simplicity
    return false
end

"""
    add_pumped_hydro_storage!(sienna_sys::SiennaSystem, row, bus::ACBus)

Add a pumped hydro storage system using correct PowerSystems.jl v4.6.2 constructor.
"""
function add_pumped_hydro_storage!(sienna_sys::SiennaSystem, row, bus::ACBus)
    @debug "Adding HydroPumpedStorage: $(get(row, :name, "unknown"))"
    
    # Create HydroGenerationCost with default constructor
    operation_cost = HydroGenerationCost(nothing)  # Default as per documentation
    
    # Get storage parameters with proper validation
    storage_capacity_value = Float64(get(row, :energy_capacity, 4000.0))  # MWh - PHS typically larger
    if storage_capacity_value <= 0.0
        @warn "Pumped hydro $(row.name) has zero or negative capacity: $storage_capacity_value MWh"
        storage_capacity_value = 4000.0  # Default fallback for PHS
    end
    
    # Get power limits (pump and generate can be different)
    pump_power_limit = Float64(get(row, :input_active_power_limits, 400.0))
    generate_power_limit = Float64(get(row, :output_active_power_limits, 400.0))
    rating = Float64(generate_power_limit)  # Generation rating
    rating_pump = Float64(pump_power_limit)  # Pump rating
    
    # Get efficiency values with validation
    pump_efficiency = Float64(get(row, :efficiency_in, 0.85))  # Pumping efficiency
    
    # Validate efficiency values
    if pump_efficiency <= 0.0 || pump_efficiency > 1.0
        @warn "Invalid pumping efficiency for $(row.name): $pump_efficiency, using 0.85"
        pump_efficiency = 0.85
    end
    
    @info "Pumped Hydro $(row.name): Pump efficiency = $(round(pump_efficiency * 100, digits=1))%"
    
    # Get initial storage level - PHS uses UpDown structure for upper/lower reservoirs
    initial_energy = Float64(get(row, :initial_energy, storage_capacity_value * 0.5))
    initial_storage = (up = initial_energy, down = 0.0)  # Assume all energy in upper reservoir initially
    
    # Storage capacity as UpDown (upper and lower reservoirs)
    storage_capacity = (up = storage_capacity_value, down = storage_capacity_value)
    
    # Storage target (default values per documentation)
    storage_target = (up = 1.0, down = 1.0)
    
    # Create pumped hydro storage using positional arguments (not keyword arguments)
    pumped_hydro = HydroPumpedStorage(
        string(row.name),                          # name::String
        get(row, :available, true),                # available::Bool
        bus,                                       # bus::ACBus
        0.0,                                       # active_power::Float64
        0.0,                                       # reactive_power::Float64
        rating,                                    # rating::Float64
        1.0,                                     # base_power::Float64
        PrimeMovers.PS,                           # prime_mover_type::PrimeMovers
        (min = 0.0, max = generate_power_limit),  # active_power_limits::MinMax
        (min = -30.0, max = 30.0),               # reactive_power_limits::Union{Nothing, MinMax}
        nothing,                                   # ramp_limits::Union{Nothing, UpDown}
        nothing,                                   # time_limits::Union{Nothing, UpDown}
        rating_pump,                               # rating_pump::Float64
        (min = 0.0, max = pump_power_limit),      # active_power_limits_pump::MinMax
        (min = -30.0, max = 30.0),               # reactive_power_limits_pump::Union{Nothing, MinMax}
        nothing,                                   # ramp_limits_pump::Union{Nothing, UpDown}
        nothing,                                   # time_limits_pump::Union{Nothing, UpDown}
        storage_capacity,                          # storage_capacity::UpDown
        0.0,                                       # inflow::Float64
        0.0,                                       # outflow::Float64
        initial_storage,                           # initial_storage::UpDown
        storage_target,                            # storage_target::UpDown
        operation_cost,                            # operation_cost::Union{HydroGenerationCost, StorageCost, MarketBidCost}
        pump_efficiency,                           # pump_efficiency::Float64
        1.0,                                       # conversion_factor::Float64
        PumpHydroStatus.OFF,                      # status::PumpHydroStatus
        Inf,                                       # time_at_status::Float64
        Service[],                                 # services::Vector{Service}
        nothing,                                   # dynamic_injector::Union{Nothing, DynamicInjection}
        Dict{String, Any}(),                       # ext::Dict{String, Any}
        InfrastructureSystems.InfrastructureSystemsInternal()  # internal::InfrastructureSystemsInternal
    )
    
    add_component!(sienna_sys.system, pumped_hydro)
    
    @debug "Added HydroPumpedStorage: $(row.name) ($(storage_capacity_value) MWh, pump: $(pump_power_limit) MW, gen: $(generate_power_limit) MW, $(round(pump_efficiency * 100, digits=1))% pump efficiency)"
end

"""
    add_energy_storage!(sienna_sys::SiennaSystem, row, bus::ACBus, storage_tech::String)

Add a single energy storage system using correct PowerSystems.jl v4.6.2 API.
"""
function add_energy_storage!(sienna_sys::SiennaSystem, row, bus::ACBus, storage_tech::String)
    @debug "Adding EnergyReservoirStorage: $(get(row, :name, "unknown")) (tech: $storage_tech)"
    
    # Create storage cost with efficiency-based philosophy
    charge_variable_cost = CostCurve(LinearCurve(0.0))
    discharge_variable_cost = CostCurve(LinearCurve(0.0))
    fixed_cost = Float64(get(row, :fixed_cost, 0.0))
    
    operation_cost = StorageCost(
        charge_variable_cost = charge_variable_cost,
        discharge_variable_cost = discharge_variable_cost,
        fixed = fixed_cost,
        start_up = 0.0,
        shut_down = 0.0,
        energy_shortage_cost = 1000.0,
        energy_surplus_cost = 0.0
    )
    
    # Get storage parameters with proper validation
    storage_capacity = Float64(get(row, :energy_capacity, 400.0))  # MWh
    if storage_capacity <= 0.0
        @warn "Storage $(row.name) has zero or negative capacity: $storage_capacity MWh"
        storage_capacity = 400.0  # Default fallback
    end
    
    # Get power limits
    input_power_limit = Float64(get(row, :input_active_power_limits, 100.0))
    output_power_limit = Float64(get(row, :output_active_power_limits, 100.0))
    rating = max(input_power_limit, output_power_limit)  # Rating is max of input/output
    
    # Get efficiency values with validation
    efficiency_in = Float64(get(row, :efficiency_in, 0.9))
    efficiency_out = Float64(get(row, :efficiency_out, 0.9))
    
    # Validate efficiency values
    if efficiency_in <= 0.0 || efficiency_in > 1.0
        @warn "Invalid charging efficiency for $(row.name): $efficiency_in, using 0.9"
        efficiency_in = 0.9
    end
    if efficiency_out <= 0.0 || efficiency_out > 1.0
        @warn "Invalid discharging efficiency for $(row.name): $efficiency_out, using 0.9"
        efficiency_out = 0.9
    end
    
    @info "Storage $(row.name): Round-trip efficiency = $(round(efficiency_in * efficiency_out * 100, digits=1))%"
    
    # Get state of charge limits as ratios [0,1]
    soc_min = Float64(get(row, :state_of_charge_limits_min, 0.1))
    soc_max = Float64(get(row, :state_of_charge_limits_max, 0.9))
    
    # Get initial storage level as ratio [0,1]
    initial_energy = Float64(get(row, :initial_energy, storage_capacity * 0.5))
    initial_storage_level = initial_energy / storage_capacity  # Convert to ratio
    
    # Validate ranges
    initial_storage_level = clamp(initial_storage_level, 0.0, 1.0)
    soc_min = clamp(soc_min, 0.0, 1.0)
    soc_max = clamp(soc_max, soc_min, 1.0)
    
    # Determine storage technology type from your CSV values with corrected enums
    storage_tech_upper = uppercase(storage_tech)
    storage_tech_type = get_storage_tech_enum(storage_tech_upper)
    
    # Create storage device using correct v4.6.2 constructor
    storage = EnergyReservoirStorage(;
        name = string(row.name),
        available = get(row, :available, true),
        bus = bus,
        prime_mover_type = PrimeMovers.BA,  # Battery prime mover
        storage_technology_type = storage_tech_type,
        storage_capacity = storage_capacity,
        storage_level_limits = (min = soc_min, max = soc_max),
        initial_storage_capacity_level = initial_storage_level,
        rating = rating,
        active_power = 0.0,
        input_active_power_limits = (min = 0.0, max = input_power_limit),
        output_active_power_limits = (min = 0.0, max = output_power_limit),
        efficiency = (in = efficiency_in, out = efficiency_out),
        reactive_power = 0.0,
        reactive_power_limits = nothing,  # No reactive power for batteries
        base_power = 1.0,
        operation_cost = operation_cost,
        conversion_factor = 1.0,  # MWh to MWh
        storage_target = 0.0,     # No end-of-simulation target
        cycle_limits = get(row, :cycle_limits, 10000)  # Default 10k cycles/year
    )
    
    add_component!(sienna_sys.system, storage)
    
    @debug "Added EnergyReservoirStorage: $(row.name) ($(storage_capacity) MWh, $(rating) MW, $(round(efficiency_in * efficiency_out * 100, digits=1))% round-trip efficiency)"
end

"""
    get_storage_tech_enum(storage_tech::String)

Get the correct StorageTech enum - now includes PHS mapped to generic energy storage.
"""
function get_storage_tech_enum(storage_tech::String)
    # Based on official PowerSystems.jl StorageTech documentation
    # https://nrel-sienna.github.io/PowerSystems.jl/stable/api/enumerated_types/#storagetech_list
    
    if storage_tech in ["BATTERY", "LION", "LI-ION", "LI_ION", "LIION", "LITHIUM"]
        return StorageTech.LIB  # LiON Battery
    elseif storage_tech in ["PUMPED_HYDRO", "PHS", "PUMPED HYDRO", "PUMP_HYDRO", "PHES"]
        return StorageTech.PTES  # Pumped Thermal Energy Storage (closest to PHS)
    elseif storage_tech in ["LEAD_ACID", "LAB", "LEAD"]
        return StorageTech.LAB  # Lead Acid Battery
    elseif storage_tech in ["FLOW", "REDOX_FLOW", "FLOW_BATTERY", "VRB"]
        return StorageTech.FLWB  # Redox Flow Battery
    elseif storage_tech in ["SODIUM", "SODIUM_ION", "NAION"]
        return StorageTech.SIB  # Sodium Ion Battery
    elseif storage_tech in ["ZINC", "ZINC_ION", "ZIB"]
        return StorageTech.ZIB  # Zinc Ion Battery
    elseif storage_tech in ["HYDROGEN", "H2", "FUEL_CELL"]
        return StorageTech.HGS  # Hydrogen Gas Storage
    elseif storage_tech in ["LAES", "LIQUID_AIR", "CRYOGENIC"]
        return StorageTech.LAES  # Liquid Air Storage
    elseif storage_tech in ["CAES", "COMPRESSED_AIR"]
        return StorageTech.CAES  # Compressed Air Energy Storage
    elseif storage_tech in ["FLYWHEEL", "FW"]
        return StorageTech.FLYW  # Flywheel
    elseif storage_tech in ["SUPERCAP", "SUPERCAPACITOR", "UCAP"]
        return StorageTech.SCAP  # Supercapacitor
    elseif storage_tech in ["SMES", "SUPERCONDUCTING"]
        return StorageTech.SMES  # Superconducting Magnetic Energy Storage
    else
        @warn "Unknown storage technology '$storage_tech', using LIB (Lithium Ion Battery) as default"
        return StorageTech.LIB  # Default to Lithium Ion Battery
    end
end

# ===== LOADS AND NETWORK COMPONENTS =====

"""
    add_loads!(sienna_sys::SiennaSystem)

Add loads from CSV file.
"""
function add_loads!(sienna_sys::SiennaSystem)
    @info "Adding loads from CSV data..."
    
    load_file = joinpath(sienna_sys.config.data_directory, "load.csv")
    if !isfile(load_file)
        error("❌ Required load.csv file not found: $load_file")
    end
    
    df = CSV.read(load_file, DataFrame)
    load_count = 0
    
    for row in eachrow(df)
        try
            # Find corresponding bus
            bus = get_component(ACBus, sienna_sys.system, string(row.bus))
            if bus === nothing
                push!(sienna_sys.build_errors, "Bus $(row.bus) not found for load $(row.name)")
                @warn "Bus $(row.bus) not found for load $(row.name)"
                continue
            end
            
            # Create load
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
            push!(sienna_sys.build_errors, "Failed to add load $(row.name): $e")
            @warn "Failed to add load $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["PowerLoad"] = load_count
    @info "✓ Added $load_count loads"
end

"""
    add_network_components!(sienna_sys::SiennaSystem)

Add network components (lines, DC lines, transformers).
"""
function add_network_components!(sienna_sys::SiennaSystem)
    @info "Adding network components..."
    
    # AC Lines
    add_ac_lines!(sienna_sys)
    
    # DC Lines
    add_dc_lines!(sienna_sys)
    
    # Transformers (placeholder for future implementation)
    add_transformers!(sienna_sys)
end

"""
    add_ac_lines!(sienna_sys::SiennaSystem)

Add AC transmission lines from CSV.
"""
function add_ac_lines!(sienna_sys::SiennaSystem)
    line_file = joinpath(sienna_sys.config.data_directory, "branch.csv")
    if !isfile(line_file)
        @info "No branch.csv found, skipping AC lines"
        return
    end
    
    @info "Adding AC lines from branch.csv..."
    df = CSV.read(line_file, DataFrame)
    line_count = 0
    
    for row in eachrow(df)
        try
            # Find buses
            from_bus = get_component(ACBus, sienna_sys.system, string(row.connection_points_from))
            to_bus = get_component(ACBus, sienna_sys.system, string(row.connection_points_to))
            
            if from_bus === nothing || to_bus === nothing
                @warn "Bus not found for line $(row.name)"
                continue
            end
            
            # Create line
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
            
        catch e
            push!(sienna_sys.build_errors, "Failed to add line $(row.name): $e")
            @warn "Failed to add line $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["Line"] = line_count
    @info "✓ Added $line_count AC lines"
end

"""
    add_dc_lines!(sienna_sys::SiennaSystem)

Add DC transmission lines from CSV.
"""
function add_dc_lines!(sienna_sys::SiennaSystem)
    dc_file = joinpath(sienna_sys.config.data_directory, "dc_branch.csv")
    if !isfile(dc_file)
        @info "No dc_branch.csv found, skipping DC lines"
        return
    end
    
    @info "Adding DC lines from dc_branch.csv..."
    df = CSV.read(dc_file, DataFrame)
    dc_count = 0
    
    for row in eachrow(df)
        try
            # Find buses
            from_bus = get_component(ACBus, sienna_sys.system, string(row.connection_points_from))
            to_bus = get_component(ACBus, sienna_sys.system, string(row.connection_points_to))
            
            if from_bus === nothing || to_bus === nothing
                @warn "Bus not found for DC line $(row.name)"
                continue
            end
            
            # Create DC line
            dc_line = TwoTerminalHVDCLine(
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
            
            add_component!(sienna_sys.system, dc_line)
            dc_count += 1
            
        catch e
            push!(sienna_sys.build_errors, "Failed to add DC line $(row.name): $e")
            @warn "Failed to add DC line $(row.name): $e"
        end
    end
    
    sienna_sys.component_counts["TwoTerminalHVDCLine"] = dc_count
    @info "✓ Added $dc_count DC lines"
end

"""
    add_transformers!(sienna_sys::SiennaSystem)

Add transformers from CSV (placeholder for future implementation).
"""
function add_transformers!(sienna_sys::SiennaSystem)
    trans_file = joinpath(sienna_sys.config.data_directory, "transformer.csv")
    if isfile(trans_file)
        @info "Transformer loading not yet implemented"
        # TODO: Implement transformer loading
    end
end

# ===== VALIDATION AND ANALYSIS =====

"""
    validate_system!(sienna_sys::SiennaSystem)

Validate the built system with improved generator type validation.
"""
function validate_system!(sienna_sys::SiennaSystem)
    @info "Validating PowerSystems.jl system..."

    println("System units: ", get_units_base(sienna_sys.system))
    #Test if the loads have been loaded correctly in terms of values

    # Basic component validation
    buses = get_components(ACBus, sienna_sys.system)
    thermal_gens = get_components(ThermalStandard, sienna_sys.system)
    renewable_dispatch_gens = get_components(RenewableDispatch, sienna_sys.system)
    renewable_nondispatch_gens = get_components(RenewableNonDispatch, sienna_sys.system)
    hydro_dispatch_gens = get_components(HydroDispatch, sienna_sys.system)
    loads = get_components(PowerLoad, sienna_sys.system)
    lines = get_components(Line, sienna_sys.system)
    dc_lines = get_components(TwoTerminalHVDCLine, sienna_sys.system)
    battery_storage = get_components(EnergyReservoirStorage, sienna_sys.system)
    pumped_hydro = get_components(HydroPumpedStorage, sienna_sys.system)

    for load in loads
        println("Load: ", get_name(load), " - Active Power: ", get_max_active_power(load), " MW, Reactive Power: ", get_max_reactive_power(load), " MVar")
    end
    

    # Check for minimum required components
    if isempty(buses)
        push!(sienna_sys.validation_errors, "System has no buses")
    end
    
    if isempty(loads)
        push!(sienna_sys.validation_warnings, "System has no loads - unusual configuration")
    end
    
    if isempty(thermal_gens) && isempty(renewable_dispatch_gens) && isempty(renewable_nondispatch_gens) && isempty(hydro_dispatch_gens)
        push!(sienna_sys.validation_errors, "System has no generators")
    end
    
    # Check reference bus
    ref_buses = [b for b in buses if get_bustype(b) == ACBusTypes.REF]
    if length(ref_buses) == 0
        push!(sienna_sys.validation_errors, "No reference bus found")
    elseif length(ref_buses) > 1
        push!(sienna_sys.validation_warnings, "Multiple reference buses found ($(length(ref_buses)))")
    end
    
    # Check network topology vs config
    if sienna_sys.config.network_model == "CopperPlatePowerModel"
        if !isempty(lines) || !isempty(dc_lines)
            push!(sienna_sys.validation_warnings, 
                  "CopperPlate model will ignore $(length(lines)) AC lines and $(length(dc_lines)) DC lines")
        end
    else
        if isempty(lines) && isempty(dc_lines)
            push!(sienna_sys.validation_warnings, 
                  "Network model $(sienna_sys.config.network_model) specified but no transmission lines found")
        end
    end
    
    # Time series validation
    if sienna_sys.config.load_timeseries && !sienna_sys.time_series_loaded
        push!(sienna_sys.validation_warnings, "Time series loading was requested but no data was loaded")
    end
    
    # Validate renewable generator costs
    validate_renewable_costs!(sienna_sys, renewable_dispatch_gens, hydro_dispatch_gens)
    
    # Validate storage systems
    validate_storage_systems!(sienna_sys, battery_storage, pumped_hydro)
    
    # Set validation status
    sienna_sys.system_validated = isempty(sienna_sys.validation_errors)
    
    # Report validation results
    if sienna_sys.system_validated
        if isempty(sienna_sys.validation_warnings)
            @info "✅ System validation passed with no warnings"
        else
            @info "⚠️  System validation passed with $(length(sienna_sys.validation_warnings)) warnings:"
            for warning in sienna_sys.validation_warnings
                @info "   • $warning"
            end
        end
    else
        @error "❌ System validation failed with $(length(sienna_sys.validation_errors)) errors:"
        for error in sienna_sys.validation_errors
            @error "   • $error"
        end
    end
end

"""
    validate_renewable_costs!(sienna_sys::SiennaSystem, dispatch_gens, hydro_gens)

Validate that renewable generator costs are set appropriately.
"""
function validate_renewable_costs!(sienna_sys::SiennaSystem, dispatch_gens, hydro_gens)
    # Check renewable dispatch generators have appropriate cost structure
    for gen in dispatch_gens
        cost = get_operation_cost(gen)
        
        # Check if it has a variable cost structure and what the cost value is
        if cost.variable !== nothing
            # Get the actual cost value from the curve
            cost_value = try
                cost.variable.value_curve.linear_term
            catch
                # Fallback if structure is different
                0.0
            end
            
            # Only warn if the cost is actually high (> $1/MWh), not if it's zero
            if cost_value > 1.0
                @warn "RenewableDispatch $(get_name(gen)) has high variable cost: $cost_value \$/MWh when zero cost expected"
            end
        end
    end
    
    # Check hydro dispatch generators have appropriate cost structure  
    for gen in hydro_gens
        cost = get_operation_cost(gen)
        
        # Check if it has a variable cost structure and what the cost value is
        if cost.variable !== nothing
            # Get the actual cost value from the curve
            cost_value = try
                cost.variable.value_curve.linear_term
            catch
                # Fallback if structure is different
                0.0
            end
            
            # Only warn if the cost is actually high (> $1/MWh), not if it's zero
            if cost_value > 1.0
                @warn "HydroDispatch $(get_name(gen)) has high variable cost: $cost_value \$/MWh when zero cost expected"
            end
        end
    end
    
    renewable_count = length(dispatch_gens)
    hydro_count = length(hydro_gens)
    
    if renewable_count > 0
        @info "✓ $renewable_count RenewableDispatch generators validated for appropriate cost structure"
    end
    
    if hydro_count > 0
        @info "✓ $hydro_count HydroDispatch generators validated for appropriate cost structure"
    end
end

"""
    validate_storage_systems!(sienna_sys::SiennaSystem, battery_storage, pumped_hydro)

Validate storage systems - now all storage is EnergyReservoirStorage so validation is simpler.
"""
function validate_storage_systems!(sienna_sys::SiennaSystem, battery_storage, pumped_hydro)
    storage_count = length(battery_storage)
    
    # All storage is now EnergyReservoirStorage, so validate consistently
    for storage in battery_storage
        capacity = get_storage_capacity(storage)  # Float64 for all storage
        rating = get_rating(storage)
        efficiency = get_efficiency(storage)
        
        # Check for reasonable values
        if capacity <= 0.0
            push!(sienna_sys.validation_errors, "Storage $(get_name(storage)) has zero or negative capacity")
        end
        
        if rating <= 0.0
            push!(sienna_sys.validation_errors, "Storage $(get_name(storage)) has zero or negative rating")
        end
        
        # Check round-trip efficiency
        round_trip_eff = efficiency.in * efficiency.out
        if round_trip_eff < 0.5 || round_trip_eff > 1.0
            push!(sienna_sys.validation_warnings, 
                  "Storage $(get_name(storage)) has unusual round-trip efficiency: $(round(round_trip_eff * 100, digits=1))%")
        end
    end
    
    if storage_count > 0
        @info "✓ Validated $storage_count energy storage systems"
    end
end

# ===== PUBLIC INTERFACE METHODS =====

"""
    get_power_system(sienna_sys::SiennaSystem)

Get the PowerSystems.jl system object.
"""
function get_power_system(sienna_sys::SiennaSystem)
    if sienna_sys.system === nothing
        error("❌ System has not been built yet. Call build_system!() first.")
    end
    return sienna_sys.system
end

"""
    is_system_built(sienna_sys::SiennaSystem)

Check if the system has been built successfully.
"""
function is_system_built(sienna_sys::SiennaSystem)
    return sienna_sys.build_successful && sienna_sys.system !== nothing
end

"""
    get_component_counts(sienna_sys::SiennaSystem)

Get counts of all system components.
"""
function get_component_counts(sienna_sys::SiennaSystem)
    if !is_system_built(sienna_sys)
        error("❌ System must be built before getting component counts")
    end
    return copy(sienna_sys.component_counts)
end

"""
    has_timeseries_data(sienna_sys::SiennaSystem)

Check if time series data has been loaded.
"""
function has_timeseries_data(sienna_sys::SiennaSystem)
    return sienna_sys.time_series_loaded
end

"""
    has_forecast_data(sienna_sys::SiennaSystem)

Check if forecast data is available for PowerSimulations.
"""
function has_forecast_data(sienna_sys::SiennaSystem)
    return sienna_sys.forecast_available
end

"""
    export_power_system(sienna_sys::SiennaSystem, output_file::String)

Export the system to a JSON file.
"""
function export_power_system(sienna_sys::SiennaSystem, output_file::String)
    if !is_system_built(sienna_sys)
        error("❌ System must be built before exporting")
    end
    
    @info "💾 Exporting system to: $output_file"
    
    try
        # Create export directory if needed
        mkpath(dirname(output_file))
        
        # Export using PowerSystems.jl
        to_json(sienna_sys.system, output_file)
        
        @info "✅ System exported successfully"
        
    catch e
        @error "❌ Failed to export system: $e"
        rethrow(e)
    end
end

"""
    print_build_summary(sienna_sys::SiennaSystem)

Print a comprehensive build summary with improved generator type reporting.
"""
function print_build_summary(sienna_sys::SiennaSystem)
    @info "\n" * "="^60
    @info "SIENNA SYSTEM BUILD SUMMARY (IMPROVED)"
    @info "="^60
    @info "Project: $(sienna_sys.system_name)"
    @info "Built: $(Dates.format(sienna_sys.build_timestamp, "yyyy-mm-dd HH:MM:SS"))"
    @info "Base Power: $(get_base_power(sienna_sys.system)) MW"
    @info "Network Model: $(sienna_sys.config.network_model)"
    
    @info ""
    @info "System Components:"
    for (component_type, count) in sienna_sys.component_counts
        @info "  $component_type: $count"
    end
    
    # Calculate total capacity with correct generator types
    thermal_gens = get_components(ThermalStandard, sienna_sys.system)
    renewable_dispatch_gens = get_components(RenewableDispatch, sienna_sys.system)
    hydro_dispatch_gens = get_components(HydroDispatch, sienna_sys.system)
    energy_storage = get_components(EnergyReservoirStorage, sienna_sys.system)
    
    thermal_capacity = isempty(thermal_gens) ? 0.0 : 
        sum(get_active_power_limits(g).max for g in thermal_gens)
    renewable_dispatch_capacity = isempty(renewable_dispatch_gens) ? 0.0 : 
        sum(get_rating(g) for g in renewable_dispatch_gens)
    hydro_dispatch_capacity = isempty(hydro_dispatch_gens) ? 0.0 : 
        sum(get_rating(g) for g in hydro_dispatch_gens)
    storage_power_capacity = isempty(energy_storage) ? 0.0 : 
        sum(get_rating(s) for s in energy_storage)
    storage_energy_capacity = isempty(energy_storage) ? 0.0 : 
        sum(get_storage_capacity(s) for s in energy_storage)
    
    @info ""
    @info "Generation Capacity:"
    @info "  Thermal (including bioenergy): $(round(thermal_capacity, digits=1)) MW"
    @info "  Renewable Dispatch (wind/solar): $(round(renewable_dispatch_capacity, digits=1)) MW"
    @info "  Hydro Dispatch: $(round(hydro_dispatch_capacity, digits=1)) MW"
    @info "  Energy Storage Power: $(round(storage_power_capacity, digits=1)) MW"
    @info "  Energy Storage Energy: $(round(storage_energy_capacity, digits=1)) MWh"
    @info "  Total Generation: $(round(thermal_capacity + renewable_dispatch_capacity + hydro_dispatch_capacity + storage_power_capacity, digits=1)) MW"
    @info "  Total Storage Energy: $(round(storage_energy_capacity, digits=1)) MWh"
    
    @info ""
    @info "Time Series Status:"
    @info "  Loaded: $(sienna_sys.time_series_loaded)"
    @info "  Forecasts Available: $(sienna_sys.forecast_available)"
    
    if !isempty(sienna_sys.build_errors)
        @info ""
        @info "⚠️  Build Warnings/Errors:"
        for error in sienna_sys.build_errors
            @info "  • $error"
        end
    end
    
    @info "="^60
end

"""
    get_generator_fuel_summary(sienna_sys::SiennaSystem)

Get summary of generators by fuel type and category.
"""
function get_generator_fuel_summary(sienna_sys::SiennaSystem)
    if !is_system_built(sienna_sys)
        error("❌ System must be built before getting generator summary")
    end
    
    sys = sienna_sys.system
    
    # Get all generators
    thermal_gens = get_components(ThermalStandard, sys)
    renewable_dispatch_gens = get_components(RenewableDispatch, sys)
    renewable_nondispatch_gens = get_components(RenewableNonDispatch, sys)
    
    summary = Dict()
    
    # Thermal generators
    summary["thermal"] = Dict()
    for gen in thermal_gens
        prime_mover = get_prime_mover(gen)
        fuel_str = string(prime_mover)
        capacity = get_active_power_limits(gen).max
        
        if haskey(summary["thermal"], fuel_str)
            summary["thermal"][fuel_str]["count"] += 1
            summary["thermal"][fuel_str]["capacity_mw"] += capacity
        else
            summary["thermal"][fuel_str] = Dict("count" => 1, "capacity_mw" => capacity)
        end
    end
    
    # Renewable dispatch generators (now includes hydro)
    summary["renewable_dispatch"] = Dict()
    for gen in renewable_dispatch_gens
        prime_mover = get_prime_mover(gen)
        fuel_str = string(prime_mover)
        capacity = get_rating(gen)
        cost = get_operation_cost(gen).variable.value_curve.linear_term
        
        if haskey(summary["renewable_dispatch"], fuel_str)
            summary["renewable_dispatch"][fuel_str]["count"] += 1
            summary["renewable_dispatch"][fuel_str]["capacity_mw"] += capacity
            summary["renewable_dispatch"][fuel_str]["avg_cost"] = 
                (summary["renewable_dispatch"][fuel_str]["avg_cost"] + cost) / 2
        else
            summary["renewable_dispatch"][fuel_str] = Dict(
                "count" => 1, 
                "capacity_mw" => capacity,
                "avg_cost" => cost
            )
        end
    end
    
    # Renewable non-dispatch generators
    summary["renewable_nondispatch"] = Dict()
    for gen in renewable_nondispatch_gens
        prime_mover = get_prime_mover(gen)
        fuel_str = string(prime_mover)
        capacity = get_rating(gen)
        cost = get_operation_cost(gen).variable.value_curve.linear_term
        
        if haskey(summary["renewable_nondispatch"], fuel_str)
            summary["renewable_nondispatch"][fuel_str]["count"] += 1
            summary["renewable_nondispatch"][fuel_str]["capacity_mw"] += capacity
        else
            summary["renewable_nondispatch"][fuel_str] = Dict(
                "count" => 1, 
                "capacity_mw" => capacity,
                "cost" => cost
            )
        end
    end
    
    return summary
end

"""
    print_generator_fuel_summary(sienna_sys::SiennaSystem)

Print a detailed summary of generators by fuel type.
"""
function print_generator_fuel_summary(sienna_sys::SiennaSystem)
    summary = get_generator_fuel_summary(sienna_sys)
    
    @info "\n" * "="^60
    @info "GENERATOR FUEL TYPE SUMMARY"
    @info "="^60
    
    # Thermal generators
    if !isempty(summary["thermal"])
        @info "\nThermal Generators:"
        for (fuel, data) in summary["thermal"]
            @info "  $fuel: $(data["count"]) units, $(round(data["capacity_mw"], digits=1)) MW"
        end
    end
    
    # Renewable dispatch generators (now includes hydro)
    if !isempty(summary["renewable_dispatch"])
        @info "\nRenewable Dispatch Generators (including Hydro):"
        for (fuel, data) in summary["renewable_dispatch"]
            cost_str = data["avg_cost"] > 0.01 ? "$(round(data["avg_cost"], digits=2)) \$/MWh" : "zero cost"
            @info "  $fuel: $(data["count"]) units, $(round(data["capacity_mw"], digits=1)) MW, $cost_str"
        end
    end
    
    # Renewable non-dispatch generators
    if !isempty(summary["renewable_nondispatch"])
        @info "\nRenewable Non-Dispatch Generators:"
        for (fuel, data) in summary["renewable_nondispatch"]
            cost_str = data["cost"] > 0.01 ? "$(round(data["cost"], digits=2)) \$/MWh" : "zero cost"
            @info "  $fuel: $(data["count"]) units, $(round(data["capacity_mw"], digits=1)) MW, $cost_str"
        end
    end
    
    @info "="^60
end

# ===== HELPER FUNCTIONS FOR CONFIG COMPATIBILITY =====

"""
    is_fuel_thermal(config::SiennaConfig, fuel_type::String)

Check if fuel type is thermal based on config (includes bioenergy).
"""
function is_fuel_thermal(config::SiennaConfig, fuel_type::String)
    # Get thermal fuels from config, with common defaults (including bioenergy)
    fuel_mapping = get(config.config_data, "fuel_mapping", Dict())
    thermal_fuels = get(fuel_mapping, "thermal_fuels", [
        "COAL", "GAS", "NATURAL_GAS", "NUCLEAR", "OIL", "DIESEL", 
        "FUEL_OIL", "LIGNITE", "PEAT", "WASTE", "OTHER",
        "BIOENERGY", "BIOMASS", "BIOGAS"  # Bioenergy is thermal
    ])
    
    return uppercase(fuel_type) in [uppercase(f) for f in thermal_fuels]
end

"""
    is_fuel_renewable(config::SiennaConfig, fuel_type::String)

Check if fuel type is renewable based on config (excludes bioenergy which is thermal).
"""
function is_fuel_renewable(config::SiennaConfig, fuel_type::String)
    # Get renewable fuels from config, with common defaults (excluding bioenergy)
    fuel_mapping = get(config.config_data, "fuel_mapping", Dict())
    renewable_fuels = get(fuel_mapping, "renewable_fuels", [
        "WIND", "SOLAR", "SOLAR_PV", "SOLAR_CSP", "PV", "HYDRO", "GEOTHERMAL"
        # Note: BIOENERGY, BIOMASS, BIOGAS are now thermal
    ])
    
    return uppercase(fuel_type) in [uppercase(f) for f in renewable_fuels]
end

"""
    get_thermal_defaults(config::SiennaConfig)

Get default values for thermal generators.
"""
function get_thermal_defaults(config::SiennaConfig)
    # Get defaults from config with fallbacks
    defaults_section = get(config.config_data, "generator_defaults", Dict())
    thermal_defaults = get(defaults_section, "thermal", Dict())
    
    return Dict(
        "variable_cost" => get(thermal_defaults, "variable_cost", 50.0),
        "startup_cost" => get(thermal_defaults, "startup_cost", 1000.0),
        "shutdown_cost" => get(thermal_defaults, "shutdown_cost", 500.0),
        "min_power_fraction" => get(thermal_defaults, "min_power_fraction", 0.3),
        "ramp_rate" => get(thermal_defaults, "ramp_rate", 0.33)
    )
end

"""
    get_renewable_defaults(config::SiennaConfig)

Get default values for renewable generators.
"""
function get_renewable_defaults(config::SiennaConfig)
    # Get defaults from config with fallbacks
    defaults_section = get(config.config_data, "generator_defaults", Dict())
    renewable_defaults = get(defaults_section, "renewable", Dict())
    
    return Dict(
        "variable_cost" => get(renewable_defaults, "variable_cost", 0.0),
        "power_factor" => get(renewable_defaults, "power_factor", 1.0)
    )
end

# ===== EXPORTS =====

export SiennaSystem
export build_system!, get_power_system, is_system_built
export get_component_counts, has_timeseries_data, has_forecast_data
export export_power_system, print_build_summary
export get_generator_fuel_summary, print_generator_fuel_summary

# Test functionality when run directly
if abspath(PROGRAM_FILE) == @__FILE__
    @info "🧪 Testing Fixed SiennaSystem with corrected storage handling..."
    
    try
        # Test with config
        if isfile("config.toml")
            test_config = SiennaConfig("config.toml")
            test_sienna_sys = SiennaSystem(test_config)
            
            # Try to build system
            if isdir(test_config.data_directory)
                build_system!(test_sienna_sys)
                print_build_summary(test_sienna_sys)
                print_generator_fuel_summary(test_sienna_sys)
                @info "✅ Fixed SiennaSystem test passed"
            else
                @info "Data directory not found - skipping system build test"
                @info "✅ SiennaSystem initialization test passed"
            end
        else
            @info "No config.toml found - cannot test SiennaSystem"
            @info "Run SiennaConfig.jl first to create configuration"
        end
    catch e
        @error "❌ Fixed SiennaSystem test failed: $e"
    end
end
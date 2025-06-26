#!/usr/bin/env julia

"""
Custom PowerSystems.jl System Builder - No Deprecated Parser
Uses only manual construction from CSV files
"""

using PowerSystems
using PowerSimulations
using HiGHS
using Dates
using TimeSeries
using DataFrames
using CSV
using YAML
using JSON3
using Logging

# Import the correct types for PowerSystems.jl v4.6.2
using PowerSystems: ACBusTypes, PrimeMovers, ThermalFuels, LinearCurve, CostCurve, ThermalGenerationCost, RenewableGenerationCost

global_logger(ConsoleLogger(stderr, Logging.Info))

function build_sienna_system(data_dir::String; validate_system::Bool=true, load_timeseries::Bool=true, base_power::Float64=100.0)
    @info "Building PowerSystems.jl system from: $data_dir using custom parser"
    
    # Check required files
    required_files = ["bus.csv", "load.csv"]
    for file in required_files
        if !isfile(joinpath(data_dir, file))
            error("Required file $file not found in $data_dir")
        end
    end
    
    # Create empty system
    sys = System(base_power)
    
    # Add components in dependency order
    add_buses_from_csv!(sys, joinpath(data_dir, "bus.csv"))
    
    # Add generators - check for separate files first, then unified
    thermal_gen_file = joinpath(data_dir, "thermal_gen.csv")
    renewable_gen_file = joinpath(data_dir, "renewable_gen.csv")
    unified_gen_file = joinpath(data_dir, "gen.csv")
    
    if isfile(thermal_gen_file) && isfile(renewable_gen_file)
        @info "Found separate generator files - using dedicated thermal and renewable CSVs"
        add_thermal_generators_from_csv!(sys, thermal_gen_file)
        add_renewable_generators_from_csv!(sys, renewable_gen_file)
    elseif isfile(unified_gen_file)
        @info "Found unified generator file - splitting by fuel type"
        add_unified_generators_from_csv!(sys, unified_gen_file)
    else
        @warn "No generator files found"
    end
    
    add_loads_from_csv!(sys, joinpath(data_dir, "load.csv"))
    
    # Optional components
    branch_file = joinpath(data_dir, "branch.csv")
    if isfile(branch_file)
        add_branches_from_csv!(sys, branch_file)
    end
    
    dc_branch_file = joinpath(data_dir, "dc_branch.csv")
    if isfile(dc_branch_file)
        add_dc_branches_from_csv!(sys, dc_branch_file)
    end
    
    storage_file = joinpath(data_dir, "storage.csv")
    if isfile(storage_file)
        add_storage_from_csv!(sys, storage_file)
    end
    
    # Post-process the system
    sys = post_process_system!(sys)
    
    # Add time series if requested and available
    if load_timeseries
        add_time_series_data!(sys, data_dir)
    end
    
    # Validate if requested
    if validate_system
        validate_system_components(sys)
    end
    
    # Print system summary
    print_system_summary(sys)
    
    return sys
end

function add_buses_from_csv!(sys::System, bus_file::String)
    @info "Adding buses from $bus_file"
    
    df = CSV.read(bus_file, DataFrame)
    
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
                          Float64(row.base_voltage) : 138.0
        
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

function add_thermal_generators_from_csv!(sys::System, thermal_file::String)
    @info "Adding thermal generators from $thermal_file"
    
    df = CSV.read(thermal_file, DataFrame)
    
    if nrow(df) == 0
        @info "No thermal generators found"
        return
    end
    
    for row in eachrow(df)
        bus = get_component(ACBus, sys, string(row.bus))
        if bus === nothing
            @warn "Bus $(row.bus) not found for generator $(row.name)"
            continue
        end
        
        try
            # Create cost structure
            variable_cost = get(row, :variable, 50.0)
            variable_linear_curve = LinearCurve(variable_cost)
            variable_cost_curve = CostCurve(variable_linear_curve)
            
            operation_cost = ThermalGenerationCost(
                variable = variable_cost_curve,
                fixed = get(row, :fixed_cost, 0.0),
                start_up = get(row, :startup, 1000.0),
                shut_down = get(row, :shutdown, 0.0)
            )
            
            # Create thermal generator
            gen = ThermalStandard(
                string(row.name),                                    # name
                get(row, :available, true),                         # available
                get(row, :status, 1) == 1,                         # status
                bus,                                                 # bus
                get(row, :active_power, 0.0),                       # active_power
                0.0,                                                 # reactive_power
                get(row, :max_active_power, 100.0),                 # rating
                (min = get(row, :min_active_power, 0.0),             # active_power_limits
                 max = get(row, :max_active_power, 100.0)),
                (min = get(row, :min_reactive_power, -30.0),         # reactive_power_limits
                 max = get(row, :max_reactive_power, 30.0)),
                (up = get(row, :ramp_30, 100.0),                    # ramp_limits
                 down = get(row, :ramp_30, 100.0)),
                operation_cost,                                      # operation_cost
                100.0                                                # base_power
            )
            
            add_component!(sys, gen)
            
        catch e
            @warn "Failed to add thermal generator $(row.name): $e"
        end
    end
    
    @info "Added $(nrow(df)) thermal generators"
end

function add_renewable_generators_from_csv!(sys::System, renewable_file::String)
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
            # Create cost structure
            variable_cost = get(row, :variable, 0.0)
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
            
            # Create renewable generator
            gen = RenewableDispatch(
                string(row.name),                                    # name
                get(row, :available, true),                         # available
                bus,                                                 # bus
                get(row, :active_power, 0.0),                       # active_power
                0.0,                                                 # reactive_power
                get(row, :max_active_power, 100.0),                 # rating
                prime_mover,                                         # prime_mover_type
                (min = get(row, :min_reactive_power, 0.0),           # reactive_power_limits
                 max = get(row, :max_reactive_power, 0.0)),
                1.0,                                                 # power_factor
                operation_cost,                                      # operation_cost
                100.0                                                # base_power
            )
            
            add_component!(sys, gen)
            
        catch e
            @warn "Failed to add renewable generator $(row.name): $e"
        end
    end
    
    @info "Added $(nrow(df)) renewable generators"
end

function add_unified_generators_from_csv!(sys::System, gen_file::String)
    @info "Splitting unified generators from $gen_file by fuel type"
    
    df = CSV.read(gen_file, DataFrame)
    
    if nrow(df) == 0
        @info "No generators found"
        return
    end
    
    # Define fuel type sets
    thermal_fuels = Set(["COAL", "NATURAL_GAS", "DIESEL", "NUCLEAR", "BIOMASS"])
    renewable_fuels = Set(["WIND", "SOLAR", "HYDRO", "GEOTHERMAL"])
    
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
            add_single_thermal_generator!(sys, row, bus)
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
            add_single_renewable_generator!(sys, row, bus)
            renewable_count += 1
        catch e
            @warn "Failed to add renewable generator $(row.name): $e"
        end
    end
    
    @info "Added $thermal_count thermal and $renewable_count renewable generators from unified file"
end

function add_single_thermal_generator!(sys::System, row, bus::ACBus)
    # Create cost structure
    variable_cost = get(row, :variable, 50.0)
    variable_linear_curve = LinearCurve(variable_cost)
    variable_cost_curve = CostCurve(variable_linear_curve)
    
    operation_cost = ThermalGenerationCost(
        variable = variable_cost_curve,
        fixed = get(row, :fixed_cost, 0.0),
        start_up = get(row, :startup, 1000.0),
        shut_down = get(row, :shutdown, 0.0)
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

function add_single_renewable_generator!(sys::System, row, bus::ACBus)
    # Create cost structure
    variable_cost = get(row, :variable, 0.0)
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

function add_loads_from_csv!(sys::System, load_file::String)
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

function add_branches_from_csv!(sys::System, branch_file::String)
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

function add_dc_branches_from_csv!(sys::System, dc_branch_file::String)
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

function add_storage_from_csv!(sys::System, storage_file::String)
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
            # Create cost structure for storage
            variable_cost = get(row, :variable_cost, 0.0)
            variable_linear_curve = LinearCurve(variable_cost)
            variable_cost_curve = CostCurve(variable_linear_curve)
            operation_cost = StorageCost(variable_cost_curve)
            
            # Create EnergyReservoirStorage
            storage = EnergyReservoirStorage(
                name = string(row.name),
                available = get(row, :available, true),
                bus = bus,
                prime_mover_type = PrimeMovers.BA,  # Battery
                storage_technology_type = StorageTech.BATTERY,  # Battery technology
                storage_capacity = get(row, :energy_capacity, 100.0),  # Energy capacity in MWh
                initial_energy = get(row, :initial_energy, 50.0),
                state_of_charge_limits = (min = 0.1, max = 0.9),
                rating = get(row, :output_active_power_limits, 50.0),  # Power rating
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
            # Print more detailed error information
            @warn "Error details: $(sprint(showerror, e))"
        end
    end
    
    @info "Added $(nrow(df)) storage units"
end

function add_time_series_data!(sys::System, data_dir::String)
    @info "Checking for time series data..."
    
    # Check for time series metadata
    ts_metadata_file = joinpath(data_dir, "timeseries_metadata.json")
    if !isfile(ts_metadata_file)
        @info "No timeseries_metadata.json found - skipping time series"
        return
    end
    
    try
        metadata = JSON3.read(read(ts_metadata_file, String))
        @info "Found time series metadata with $(length(metadata)) entries"
        
        # Process time series (simplified - would need full implementation)
        @info "Time series processing would be implemented here"
        
    catch e
        @warn "Failed to process time series data: $e"
    end
end

function validate_system_components(sys::System)
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

function post_process_system!(sys::System)
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
    
    # Calculate total capacity
    total_thermal_cap = sum(get_active_power_limits(g).max for g in thermal_gens)
    total_renewable_cap = sum(get_rating(g) for g in renewable_gens)
    
    # Storage capacity (if any)
    total_storage_cap = 0.0
    try
        storage_units = collect(get_components(EnergyReservoirStorage, sys))
        total_storage_cap = sum(get_rating(s) for s in storage_units)
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

function run_simple_economic_dispatch(sys::System)
    @info "Running simple Economic Dispatch..."
    
    try
        template = ProblemTemplate(NetworkModel(StandardPTDFModel))
        set_device_model!(template, ThermalStandard, ThermalBasicUnitCommitment)
        set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
        set_device_model!(template, PowerLoad, StaticPowerLoad)
        set_device_model!(template, Line, StaticBranch)
        set_device_model!(template, TwoTerminalHVDCLine, HVDCTwoTerminalLossless)
        
        decision_model = DecisionModel(
            template,
            sys;
            name = "Simple_ED",
            optimizer = HiGHS.Optimizer,
            optimizer_solve_log_print = false
        )
        
        @info "✓ Model created"
        
        build!(decision_model; output_dir = pwd())
        solve!(decision_model)
        
        results = OptimizationProblemResults(decision_model)
        @info "✓ Economic Dispatch completed successfully"
        
        return results
        
    catch e
        @error "Economic Dispatch failed: $e"
        return nothing
    end
end

# Main function
function main()
    data_dir = expanduser("~/showcase/pypsa-rsa/networks/sienna/TEST/dispatch_2030/")
    
    @info "Custom Sienna System Builder (No Parser)"
    @info "Data directory: $data_dir"
    
    # Build system
    sys = build_sienna_system(data_dir)
    
    # Save system
    save_system(sys, "custom_pypsa_system.json")
    
    # Run dispatch
    results = run_simple_economic_dispatch(sys)
    
    if results !== nothing
        @info "✅ SUCCESS: System built and Economic Dispatch completed!"
    else
        @info "✅ System built successfully (dispatch failed but system is usable)"
    end
    
    return sys
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
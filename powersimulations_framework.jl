#!/usr/bin/env julia

"""
PowerSimulations.jl Framework v0.30.2
=====================================

Complete framework for PowerSimulations v0.30.2 + PowerSystems v4.6.2 + InfrastructureSystems v2.6.0
Pure modern API throughout - no legacy compatibility layers

Features:
- Modern DecisionModel API with explicit optimizer handling
- Updated device formulation names for v0.30.2
- Proper constraint handling for InfrastructureSystems v2.6.0
- Clean time series integration with PowerSystems v4.6.2
- Comprehensive error handling with actionable feedback
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
using TimeSeries

global_logger(ConsoleLogger(stderr, Logging.Info))

"""
    create_operations_problem_template(sys::System, formulation_type::String="economic_dispatch")

Create operations problem template using PowerSimulations v0.30.2 API.
"""
function create_operations_problem_template(sys::System, formulation_type::String="economic_dispatch")
    @info "Creating operations template for: $formulation_type"
    
    # Create template with modern API
    template = ProblemTemplate()
    device_counts = Dict{String, Int}()
    
    # Network Model - CopperPlate for reliability with v0.30.2
    set_network_model!(template, NetworkModel(PSI.CopperPlatePowerModel))
    @info "  ✓ Network: CopperPlatePowerModel"
    
    # === THERMAL GENERATORS ===
    thermal_gens = get_components(ThermalStandard, sys)
    if !isempty(thermal_gens)
        if formulation_type == "economic_dispatch"
            set_device_model!(template, ThermalStandard, PSI.ThermalBasicDispatch)
        elseif formulation_type == "unit_commitment"
            set_device_model!(template, ThermalStandard, PSI.ThermalBasicUnitCommitment)
        else
            error("Unknown formulation type: $formulation_type")
        end
        device_counts["ThermalStandard"] = length(thermal_gens)
        @info "  ✓ ThermalStandard: $(length(thermal_gens)) units"
    end
    
    # === RENEWABLE DISPATCH ===
    renewable_dispatch = get_components(RenewableDispatch, sys)
    if !isempty(renewable_dispatch)
        set_device_model!(template, RenewableDispatch, PSI.RenewableFullDispatch)
        device_counts["RenewableDispatch"] = length(renewable_dispatch)
        @info "  ✓ RenewableDispatch: $(length(renewable_dispatch)) units"
    end
    
    # === NON-DISPATCHABLE RENEWABLES ===
    renewable_nondispatch = get_components(RenewableNonDispatch, sys)
    if !isempty(renewable_nondispatch)
        set_device_model!(template, RenewableNonDispatch, PSI.FixedOutput)
        device_counts["RenewableNonDispatch"] = length(renewable_nondispatch)
        @info "  ✓ RenewableNonDispatch: $(length(renewable_nondispatch)) units"
    end
    
    # === HYDRO GENERATORS ===
    hydro_dispatch = get_components(HydroDispatch, sys)
    if !isempty(hydro_dispatch)
        set_device_model!(template, HydroDispatch, PSI.HydroDispatchRunOfRiver)
        device_counts["HydroDispatch"] = length(hydro_dispatch)
        @info "  ✓ HydroDispatch: $(length(hydro_dispatch)) units"
    end
    
    hydro_energy = get_components(HydroEnergyReservoir, sys)
    if !isempty(hydro_energy)
        set_device_model!(template, HydroEnergyReservoir, PSI.HydroDispatchReservoirBudget)
        device_counts["HydroEnergyReservoir"] = length(hydro_energy)
        @info "  ✓ HydroEnergyReservoir: $(length(hydro_energy)) units"
    end
    
    # === STORAGE ===
    # Check for different storage types in v4.6.2
    try
        battery_storage = get_components(GenericBattery, sys)
        if !isempty(battery_storage)
            set_device_model!(template, GenericBattery, PSI.BookKeeping)
            device_counts["GenericBattery"] = length(battery_storage)
            @info "  ✓ GenericBattery: $(length(battery_storage)) units"
        end
    catch
        # GenericBattery might not exist, try alternatives
        try
            energy_storage = get_components(EnergyReservoirStorage, sys)
            if !isempty(energy_storage)
                set_device_model!(template, EnergyReservoirStorage, PSI.BookKeeping)
                device_counts["EnergyReservoirStorage"] = length(energy_storage)
                @info "  ✓ EnergyReservoirStorage: $(length(energy_storage)) units"
            end
        catch
            # No storage found, continue
        end
    end
    
    # === LOADS ===
    loads = get_components(PowerLoad, sys)
    if !isempty(loads)
        set_device_model!(template, PowerLoad, PSI.StaticPowerLoad)
        device_counts["PowerLoad"] = length(loads)
        @info "  ✓ PowerLoad: $(length(loads)) units"
    end
    
    # === INTERRUPTIBLE LOADS ===
    try
        interruptible_loads = get_components(InterruptiblePowerLoad, sys)
        if !isempty(interruptible_loads)
            set_device_model!(template, InterruptiblePowerLoad, PSI.DispatchablePowerLoad)
            device_counts["InterruptiblePowerLoad"] = length(interruptible_loads)
            @info "  ✓ InterruptiblePowerLoad: $(length(interruptible_loads)) units"
        end
    catch
        # InterruptiblePowerLoad may not exist in system
    end
    
    # === AC BRANCHES ===
    lines = get_components(Line, sys)
    if !isempty(lines)
        set_device_model!(template, Line, PSI.StaticBranch)
        device_counts["Line"] = length(lines)
        @info "  ✓ Line: $(length(lines)) branches"
    end
    
    # === TRANSFORMERS ===
    transformers = get_components(Transformer2W, sys)
    if !isempty(transformers)
        set_device_model!(template, Transformer2W, PSI.StaticBranch)
        device_counts["Transformer2W"] = length(transformers)
        @info "  ✓ Transformer2W: $(length(transformers)) units"
    end
    
    # === TAP TRANSFORMERS ===
    try
        tap_transformers = get_components(TapTransformer, sys)
        if !isempty(tap_transformers)
            set_device_model!(template, TapTransformer, PSI.StaticBranch)
            device_counts["TapTransformer"] = length(tap_transformers)
            @info "  ✓ TapTransformer: $(length(tap_transformers)) units"
        end
    catch
        # TapTransformer may not exist
    end
    
    # === DC LINES (Fixed for v0.30.2) ===
    dc_lines = get_components(TwoTerminalHVDCLine, sys)
    if !isempty(dc_lines)
        try
            # Try modern v0.30.2 DC line formulation
            set_device_model!(template, TwoTerminalHVDCLine, PSI.HVDCTwoTerminalDispatch)
            device_counts["TwoTerminalHVDCLine"] = length(dc_lines)
            @info "  ✓ TwoTerminalHVDCLine: $(length(dc_lines)) lines"
        catch e1
            @warn "Modern DC formulation failed: $e1"
            @info "Trying alternative DC line formulations..."
            
            # Try alternative formulations for v0.30.2
            try
                # Alternative 1: Try HVDCTwoTerminalLossless
                set_device_model!(template, TwoTerminalHVDCLine, PSI.HVDCTwoTerminalLossless)
                device_counts["TwoTerminalHVDCLine"] = length(dc_lines)
                @info "  ✓ TwoTerminalHVDCLine: $(length(dc_lines)) lines (lossless formulation)"
            catch e2
                @warn "Lossless DC formulation failed: $e2"
                try
                    # Alternative 2: Try basic dispatch formulation
                    set_device_model!(template, TwoTerminalHVDCLine, PSI.HVDCTwoTerminalPowerFlow)
                    device_counts["TwoTerminalHVDCLine"] = length(dc_lines)
                    @info "  ✓ TwoTerminalHVDCLine: $(length(dc_lines)) lines (power flow formulation)"
                catch e3
                    @error "All DC line formulations failed: $e1, $e2, $e3"
                    @error "DC lines will be excluded from optimization"
                    @error "This will significantly impact South African power system results!"
                end
            end
        end
    end
    
    @info "✓ Template created with $(length(device_counts)) device types"
    return template, device_counts
end

"""
    setup_optimizer(solver_type::String="HiGHS"; time_limit::Float64=300.0, mip_gap::Float64=0.01, threads::Int=0)

Setup optimizer for PowerSimulations v0.30.2 with proper parameter handling.
"""
function setup_optimizer(solver_type::String="HiGHS"; 
                        time_limit::Float64=300.0, 
                        mip_gap::Float64=0.01,
                        threads::Int=0)
    
    # Validate parameters
    if time_limit <= 0
        @warn "Invalid time_limit ($time_limit), using default 300.0 seconds"
        time_limit = 300.0
    end
    
    if mip_gap <= 0 || mip_gap >= 1
        @warn "Invalid mip_gap ($mip_gap), using default 0.01"
        mip_gap = 0.01
    end
    
    if threads < 0
        @warn "Invalid threads ($threads), using auto-detection"
        threads = 0
    end
    
    if solver_type == "HiGHS"
        # Conservative thread count for stability
        thread_count = threads > 0 ? threads : min(Sys.CPU_THREADS, 4)
        
        # Modern HiGHS configuration for v0.30.2
        optimizer = optimizer_with_attributes(
            HiGHS.Optimizer,
            "time_limit" => time_limit,
            "mip_rel_gap" => mip_gap,
            "threads" => thread_count,
            "output_flag" => false,
            "presolve" => "on",
            "parallel" => "on"
        )
        
        @info "✓ HiGHS optimizer configured:"
        @info "  Time limit: $(time_limit) seconds"
        @info "  MIP gap: $(mip_gap)"
        @info "  Threads: $(thread_count)"
        
    else
        error("❌ Solver $solver_type not supported. Available: HiGHS")
    end
    
    return optimizer
end

"""
    run_operations_problem(sys::System, formulation_type::String; kwargs...)

Run operations problem using pure PowerSimulations v0.30.2 API.
"""
function run_operations_problem(sys::System, 
                               formulation_type::String="economic_dispatch";
                               horizon_hours::Int=24,
                               output_dir::String=pwd(),
                               solver_type::String="HiGHS",
                               save_results::Bool=true,
                               use_default_template::Bool=false,
                               time_limit::Float64=0.0,
                               mip_gap::Float64=0.0,
                               threads::Int=0,
                               kwargs...)
    
    @info "="^60
    @info "RUNNING $(uppercase(formulation_type)) OPERATIONS PROBLEM"
    @info "="^60
    @info "PowerSimulations: v0.30.2"
    @info "PowerSystems: v4.6.2"
    @info "InfrastructureSystems: v2.6.0"
    @info ""
    @info "System: $(get_name(sys))"
    @info "Devices: $(length(get_components(Device, sys)))"
    @info "Horizon: $horizon_hours hours"
    @info "Solver: $solver_type"
    
    # Validate system has forecast data
    validate_system_for_operations(sys, formulation_type)
    
    # Create output directory
    mkpath(output_dir)
    
    # Create or use template
    if use_default_template
        @info "Using default PSI template..."
        try
            if formulation_type == "economic_dispatch"
                template = PSI.template_economic_dispatch()
                @info "✓ Default economic dispatch template"
            elseif formulation_type == "unit_commitment"  
                template = PSI.template_unit_commitment()
                @info "✓ Default unit commitment template"
            else
                error("Unknown formulation type: $formulation_type")
            end
        catch e
            @warn "Default template failed: $e"
            @info "Falling back to custom template..."
            template, device_counts = create_operations_problem_template(sys, formulation_type)
        end
    else
        @info "Creating custom template..."
        template, device_counts = create_operations_problem_template(sys, formulation_type)
    end
    
    # Setup optimizer with intelligent defaults
    if time_limit > 0.0 && mip_gap > 0.0
        @info "Using provided solver settings"
        optimizer = setup_optimizer(solver_type, 
                                   time_limit=time_limit, 
                                   mip_gap=mip_gap, 
                                   threads=threads)
    else
        @info "Using intelligent defaults for $formulation_type"
        default_time_limit = formulation_type == "unit_commitment" ? 600.0 : 300.0
        default_mip_gap = formulation_type == "unit_commitment" ? 0.02 : 0.01
        optimizer = setup_optimizer(solver_type, 
                                   time_limit=default_time_limit, 
                                   mip_gap=default_mip_gap, 
                                   threads=threads)
    end
    
    # Create DecisionModel with modern v0.30.2 API
    @info "Creating DecisionModel..."
    problem = try
        DecisionModel(
            template, 
            sys; 
            optimizer = optimizer, 
            horizon = Hour(horizon_hours),
            name = "$(formulation_type)_problem"
        )
    catch e
        @error "Failed to create DecisionModel: $e"
        @error "This suggests a fundamental compatibility issue"
        rethrow(e)
    end
    @info "✓ DecisionModel created successfully"
    
    # Build model
    @info "Building optimization model..."
    build_time = @elapsed begin
        try
            build!(problem, output_dir = output_dir)
        catch e
            @error "Model build failed: $e"
            @error "Check system components and formulations"
            rethrow(e)
        end
    end
    @info "✓ Model built in $(round(build_time, digits=2)) seconds"
    
    # Solve with proper error handling
    @info "Solving optimization problem..."
    solve_time = @elapsed begin
        try
            solve!(problem)
        catch e
            @error "Solve failed: $e"
            @error "This may indicate:"
            @error "  1. Infeasible problem formulation"
            @error "  2. Missing or invalid time series data"
            @error "  3. Constraint encoding issues (InfrastructureSystems v2.6.0)"
            rethrow(e)
        end
    end
    
    # Get results using v0.30.2 API
    @info "Extracting results..."
    results = try
        OptimizationProblemResults(problem)
    catch e
        @error "Failed to extract results: $e"
        rethrow(e)
    end
    @info "✓ Results extracted successfully"
    
    @info "✓ Operations problem completed in $(round(solve_time, digits=2)) seconds"
    
    # Print summary
    print_operations_summary(results, sys, formulation_type)
    
    # Save results if requested
    if save_results
        results_dir = save_operations_results(results, sys, formulation_type, output_dir)
        @info "✓ Results saved to: $results_dir"
    end
    
    return results, problem
end

"""
    validate_system_for_operations(sys::System, formulation_type::String)

Validate system is ready for operations problems with v0.30.2.
"""
function validate_system_for_operations(sys::System, formulation_type::String)
    @info "Validating system for operations..."
    
    # Check basic components
    buses = get_components(ACBus, sys)
    if isempty(buses)
        error("❌ System has no buses")
    end
    
    loads = get_components(PowerLoad, sys)
    if isempty(loads)
        @warn "⚠️  System has no loads - unusual but may be valid"
    end
    
    # Check generators
    thermal_gens = get_components(ThermalStandard, sys)
    renewable_gens = get_components(RenewableDispatch, sys)
    total_gens = length(thermal_gens) + length(renewable_gens)
    
    if total_gens == 0
        error("❌ System has no generators")
    end
    
    # Check for unit commitment requirements
    if formulation_type == "unit_commitment" && isempty(thermal_gens)
        @warn "⚠️  Unit commitment requested but no thermal generators found"
    end
    
    # Check time series data - crucial for v0.30.2
    @info "Checking time series data..."
    
    # Check loads have time series
    load_ts_count = 0
    for load in loads
        ts_keys = get_time_series_keys(load)
        if !isempty(ts_keys)
            load_ts_count += 1
        end
    end
    
    # Check renewables have time series
    renewable_ts_count = 0
    for gen in renewable_gens
        ts_keys = get_time_series_keys(gen)
        if !isempty(ts_keys)
            renewable_ts_count += 1
        end
    end
    
    @info "Time series status:"
    @info "  Loads with data: $load_ts_count/$(length(loads))"
    @info "  Renewables with data: $renewable_ts_count/$(length(renewable_gens))"
    
    # For v0.30.2, we need some time series data
    if load_ts_count == 0 && renewable_ts_count == 0
        @warn "⚠️  No time series data found - operations may use static values only"
    end
    
    @info "✓ System validation completed"
end

"""
    print_operations_summary(results::OptimizationProblemResults, sys::System, formulation_type::String)

Print comprehensive operations results summary.
"""
function print_operations_summary(results::OptimizationProblemResults, sys::System, formulation_type::String)
    @info "\n" * "="^60
    @info "$(uppercase(formulation_type)) RESULTS SUMMARY"
    @info "="^60
    
    # Objective value
    try
        obj_value = get_objective_value(results)
        @info "Objective Value: \$$(round(obj_value, digits=2))"
    catch e
        @warn "Could not retrieve objective value: $e"
    end
    
    # Optimizer statistics
    try
        stats = get_optimizer_stats(results)
        if nrow(stats) > 0
            @info "Solve Time: $(round(stats.solve_time[1], digits=2)) seconds"
            @info "Termination Status: $(stats.termination_status[1])"
            @info "Primal Status: $(stats.primal_status[1])"
        end
    catch e
        @warn "Could not retrieve optimizer stats: $e"
    end
    
    # Generation analysis
    analyze_generation_dispatch(results, sys)
    
    # Unit commitment analysis (if applicable)
    if formulation_type in ["unit_commitment"]
        analyze_unit_commitment(results, sys)
    end
    
    @info "="^60
end

"""
    analyze_generation_dispatch(results::OptimizationProblemResults, sys::System)

Analyze generation dispatch from results.
"""
function analyze_generation_dispatch(results::OptimizationProblemResults, sys::System)
    @info "\nGeneration Dispatch Analysis:"
    
    try
        var_names = list_variable_names(results)
        total_generation = Dict{String, Float64}()
        
        # Thermal generation
        if "ActivePowerVariable__ThermalStandard" in var_names
            thermal_data = read_variable(results, "ActivePowerVariable__ThermalStandard")
            if ncol(thermal_data) > 1
                total_generation["Thermal"] = sum(sum(eachcol(thermal_data[:, 2:end])))
                max_thermal = maximum(sum(eachrow(thermal_data[:, 2:end])))
                @info "  Thermal: $(round(total_generation["Thermal"], digits=1)) MWh total, $(round(max_thermal, digits=1)) MW peak"
            end
        end
        
        # Renewable generation
        if "ActivePowerVariable__RenewableDispatch" in var_names
            renewable_data = read_variable(results, "ActivePowerVariable__RenewableDispatch")
            if ncol(renewable_data) > 1
                total_generation["Renewable"] = sum(sum(eachcol(renewable_data[:, 2:end])))
                max_renewable = maximum(sum(eachrow(renewable_data[:, 2:end])))
                @info "  Renewable: $(round(total_generation["Renewable"], digits=1)) MWh total, $(round(max_renewable, digits=1)) MW peak"
            end
        end
        
        # Hydro generation
        if "ActivePowerVariable__HydroDispatch" in var_names
            hydro_data = read_variable(results, "ActivePowerVariable__HydroDispatch")
            if ncol(hydro_data) > 1
                total_generation["Hydro"] = sum(sum(eachcol(hydro_data[:, 2:end])))
                max_hydro = maximum(sum(eachrow(hydro_data[:, 2:end])))
                @info "  Hydro: $(round(total_generation["Hydro"], digits=1)) MWh total, $(round(max_hydro, digits=1)) MW peak"
            end
        end
        
        # Generation mix
        if !isempty(total_generation)
            total = sum(values(total_generation))
            @info "  Total Generation: $(round(total, digits=1)) MWh"
            
            for (tech, gen) in total_generation
                if gen > 0
                    pct = (gen / total) * 100
                    @info "  $(tech) share: $(round(pct, digits=1))%"
                end
            end
        else
            @warn "  No generation variables found in results"
        end
        
    catch e
        @warn "Generation analysis failed: $e"
    end
end

"""
    analyze_unit_commitment(results::OptimizationProblemResults, sys::System)

Analyze unit commitment decisions.
"""
function analyze_unit_commitment(results::OptimizationProblemResults, sys::System)
    @info "\nUnit Commitment Analysis:"
    
    try
        var_names = list_variable_names(results)
        
        if "OnVariable__ThermalStandard" in var_names
            on_status = read_variable(results, "OnVariable__ThermalStandard")
            thermal_gens = get_components(ThermalStandard, sys)
            
            if ncol(on_status) > 1
                total_units = length(thermal_gens)
                hourly_committed = [sum(on_status[h, 2:end]) for h in 1:nrow(on_status)]
                
                avg_committed = mean(hourly_committed)
                max_committed = maximum(hourly_committed)
                min_committed = minimum(hourly_committed)
                
                @info "  Total thermal units: $total_units"
                @info "  Average committed: $(round(avg_committed, digits=1)) units ($(round(avg_committed/total_units*100, digits=1))%)"
                @info "  Peak commitment: $max_committed units ($(round(max_committed/total_units*100, digits=1))%)"
                @info "  Minimum commitment: $min_committed units ($(round(min_committed/total_units*100, digits=1))%)"
            end
        end
        
        # Start-up analysis
        if "StartVariable__ThermalStandard" in var_names
            start_data = read_variable(results, "StartVariable__ThermalStandard")
            if ncol(start_data) > 1
                total_starts = sum(sum(eachcol(start_data[:, 2:end])))
                @info "  Total start-ups: $(Int(total_starts))"
            end
        end
        
    catch e
        @warn "Unit commitment analysis failed: $e"
    end
end

"""
    save_operations_results(results::OptimizationProblemResults, sys::System, formulation_type::String, output_dir::String)

Save operations results with v0.30.2 compatible format.
"""
function save_operations_results(results::OptimizationProblemResults, 
                                sys::System, 
                                formulation_type::String, 
                                output_dir::String)
    
    # Create timestamped results directory
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    results_dir = joinpath(output_dir, "$(formulation_type)_$(timestamp)")
    mkpath(results_dir)
    
    # Save results summary
    summary_file = joinpath(results_dir, "results_summary.json")
    try
        results_summary = Dict(
            "formulation_type" => formulation_type,
            "timestamp" => timestamp,
            "objective_value" => get_objective_value(results),
            "optimizer_stats" => get_optimizer_stats(results),
            "system_info" => Dict(
                "name" => get_name(sys),
                "base_power" => get_base_power(sys),
                "device_count" => length(get_components(Device, sys))
            ),
            "variables_available" => list_variable_names(results),
            "parameters_available" => list_parameter_names(results),
            "powersimulations_version" => "v0.30.2",
            "powersystems_version" => "v4.6.2"
        )
        
        open(summary_file, "w") do f
            JSON3.pretty(f, results_summary)
        end
    catch e
        @warn "Failed to save results summary: $e"
    end
    
    # Save key variables as CSV
    try
        var_names = list_variable_names(results)
        if !isempty(var_names)
            variables_dir = joinpath(results_dir, "variables")
            mkpath(variables_dir)
            
            # Important variables for analysis
            important_vars = [
                "ActivePowerVariable__ThermalStandard",
                "ActivePowerVariable__RenewableDispatch", 
                "OnVariable__ThermalStandard",
                "StartVariable__ThermalStandard",
                "StopVariable__ThermalStandard",
                "ActivePowerVariable__HydroDispatch"
            ]
            
            for var_name in important_vars
                if var_name in var_names
                    try
                        var_data = read_variable(results, var_name)
                        var_file = joinpath(variables_dir, "$(var_name).csv")
                        CSV.write(var_file, var_data)
                    catch e
                        @warn "Failed to save variable $var_name: $e"
                    end
                end
            end
        end
    catch e
        @warn "Failed to save variables: $e"
    end
    
    return results_dir
end

"""
    run_operations_suite(sys::System; kwargs...)

Run comprehensive operations suite with v0.30.2.
"""
function run_operations_suite(sys::System; 
                             output_dir::String=pwd(),
                             horizon_hours::Int=24,
                             include_unit_commitment::Bool=true,
                             use_default_templates::Bool=false,
                             solver_type::String="HiGHS",
                             time_limit_ed::Float64=300.0,
                             time_limit_uc::Float64=600.0,
                             mip_gap_ed::Float64=0.01,
                             mip_gap_uc::Float64=0.02,
                             threads::Int=0)
    
    @info "="^60
    @info "RUNNING OPERATIONS SUITE (PSI v0.30.2)"
    @info "="^60
    @info "System: $(get_name(sys))"
    @info "Output: $output_dir"
    
    mkpath(output_dir)
    suite_results = Dict()
    
    # Economic Dispatch
    @info "\n>>> Economic Dispatch"
    try
        ed_results, ed_problem = run_operations_problem(
            sys, "economic_dispatch",
            horizon_hours=horizon_hours,
            output_dir=joinpath(output_dir, "economic_dispatch"),
            solver_type=solver_type,
            use_default_template=use_default_templates,
            time_limit=time_limit_ed,
            mip_gap=mip_gap_ed,
            threads=threads
        )
        suite_results["economic_dispatch"] = Dict(
            "status" => "SUCCESS",
            "objective_value" => get_objective_value(ed_results),
            "solve_time" => get_optimizer_stats(ed_results).solve_time[1]
        )
    catch e
        @error "Economic Dispatch failed: $e"
        suite_results["economic_dispatch"] = Dict("status" => "FAILED", "error" => string(e))
    end
    
    # Unit Commitment
    if include_unit_commitment
        @info "\n>>> Unit Commitment"
        try
            uc_results, uc_problem = run_operations_problem(
                sys, "unit_commitment",
                horizon_hours=horizon_hours,
                output_dir=joinpath(output_dir, "unit_commitment"),
                solver_type=solver_type,
                use_default_template=use_default_templates,
                time_limit=time_limit_uc,
                mip_gap=mip_gap_uc,
                threads=threads
            )
            suite_results["unit_commitment"] = Dict(
                "status" => "SUCCESS",
                "objective_value" => get_objective_value(uc_results),
                "solve_time" => get_optimizer_stats(uc_results).solve_time[1]
            )
        catch e
            @error "Unit Commitment failed: $e"
            suite_results["unit_commitment"] = Dict("status" => "FAILED", "error" => string(e))
        end
    end
    
    # Save suite summary
    suite_summary_file = joinpath(output_dir, "operations_suite_summary.json")
    suite_summary = Dict(
        "timestamp" => string(now()),
        "system_name" => get_name(sys),
        "horizon_hours" => horizon_hours,
        "powersimulations_version" => "v0.30.2",
        "solver_settings" => Dict(
            "solver_type" => solver_type,
            "time_limit_ed" => time_limit_ed,
            "time_limit_uc" => time_limit_uc,
            "mip_gap_ed" => mip_gap_ed,
            "mip_gap_uc" => mip_gap_uc,
            "threads" => threads
        ),
        "results" => suite_results
    )
    
    open(suite_summary_file, "w") do f
        JSON3.pretty(f, suite_summary)
    end
    
    @info "\n" * "="^60
    @info "OPERATIONS SUITE COMPLETE"
    for (problem_type, result) in suite_results
        status = result["status"]
        if status == "SUCCESS"
            obj_val = round(result["objective_value"], digits=2)
            solve_time = round(result["solve_time"], digits=2)
            @info "  ✓ $problem_type: \$obj_val in $(solve_time)s"
        else
            @info "  ✗ $problem_type: FAILED"
        end
    end
    @info "Suite summary: $suite_summary_file"
    @info "="^60
    
    return suite_results
end

"""
    load_system_and_run_operations(system_path::String, formulation_type::String; kwargs...)

Load system and run operations in one step.
"""
function load_system_and_run_operations(system_path::String, 
                                       formulation_type::String="economic_dispatch"; 
                                       kwargs...)
    @info "Loading system from: $system_path"
    
    if !isfile(system_path)
        error("❌ System file not found: $system_path")
    end
    
    # Load system
    sys = System(system_path)
    @info "✓ System loaded: $(get_name(sys))"
    
    # Run operations problem
    return run_operations_problem(sys, formulation_type; kwargs...)
end

"""
    run_simulation(sys_da::System, sys_rt::System; kwargs...)

Run multi-stage simulation with modern v0.30.2 API.
"""
function run_simulation(sys_da::System, 
                       sys_rt::System=sys_da;
                       output_dir::String=pwd(),
                       steps::Int=7,
                       simulation_name::String="simulation",
                       solver_type::String="HiGHS",
                       time_limit::Float64=600.0,
                       mip_gap::Float64=0.02,
                       threads::Int=0,
                       kwargs...)
    
    @info "="^60
    @info "RUNNING MULTI-STAGE SIMULATION (PSI v0.30.2)"
    @info "="^60
    @info "DA System: $(get_name(sys_da))"
    @info "RT System: $(get_name(sys_rt))"
    @info "Steps: $steps"
    
    # Setup optimizer
    optimizer = setup_optimizer(solver_type, 
                               time_limit=time_limit, 
                               mip_gap=mip_gap, 
                               threads=threads)
    
    # Create templates
    template_uc = PSI.template_unit_commitment()
    template_ed = PSI.template_economic_dispatch()
    
    # Create Decision Models with modern API
    decision_model_uc = DecisionModel(
        template_uc, 
        sys_da; 
        optimizer = optimizer, 
        name = "UC"
    )
    
    decision_model_ed = DecisionModel(
        template_ed, 
        sys_rt; 
        optimizer = optimizer, 
        name = "ED"
    )
    
    # Create SimulationModels
    models = SimulationModels(
        decision_models = [decision_model_uc, decision_model_ed]
    )
    
    # Define sequence with feedforward
    sequence = SimulationSequence(
        models = models,
        feedforwards = Dict(
            "ED" => [
                PSI.SemiContinuousFeedforward(
                    component_type = ThermalStandard,
                    source = PSI.OnVariable,
                    affected_values = [PSI.ActivePowerVariable],
                ),
            ],
        ),
    )
    
    # Create and run simulation
    sim = Simulation(
        name = simulation_name,
        steps = steps,
        models = models,
        sequence = sequence,
        simulation_folder = output_dir
    )
    
    @info "Building simulation..."
    build!(sim)
    
    @info "Executing simulation..."
    execute!(sim)
    
    @info "✓ Multi-stage simulation completed!"
    
    results = SimulationResults(sim)
    return results, sim
end

"""
    demo_workflow(system_path::String; kwargs...)

Complete demonstration workflow with v0.30.2.
"""
function demo_workflow(system_path::String; 
                      time_limit::Float64=300.0, 
                      mip_gap::Float64=0.01, 
                      threads::Int=0)
    @info "="^60
    @info "DEMO: POWERSIMULATIONS v0.30.2 WORKFLOW"
    @info "="^60
    
    # Load system
    sys = System(system_path)
    
    # Create output directory
    output_dir = joinpath(dirname(system_path), "powersimulations_v030_results")
    
    # Run operations suite
    operations_results = run_operations_suite(
        sys, 
        output_dir=output_dir,
        use_default_templates=true,
        time_limit_ed=time_limit,
        time_limit_uc=time_limit * 2,
        mip_gap_ed=mip_gap,
        mip_gap_uc=mip_gap * 2,
        threads=threads
    )
    
    @info "Demo complete! Check results in: $output_dir"
    return operations_results
end

"""
    run_operations_on_exported_system(export_info::Dict; kwargs...)

Run operations on system exported by build_sienna_system.jl.
"""
function run_operations_on_exported_system(export_info::Dict; 
                                          time_limit_ed::Float64=300.0,
                                          time_limit_uc::Float64=600.0,
                                          mip_gap_ed::Float64=0.01,
                                          mip_gap_uc::Float64=0.02,
                                          threads::Int=0,
                                          kwargs...)
    @info "Running operations on exported system: $(export_info["case_name"])"
    
    # Load system
    system_file = export_info["system_file"]
    if !isfile(system_file)
        error("❌ Exported system file not found: $system_file")
    end
    
    sys = System(system_file)
    
    # Create output directory
    output_dir = joinpath(dirname(system_file), "operations_results")
    
    # Run operations suite
    return run_operations_suite(
        sys, 
        output_dir=output_dir, 
        use_default_templates=true,
        time_limit_ed=time_limit_ed,
        time_limit_uc=time_limit_uc,
        mip_gap_ed=mip_gap_ed,
        mip_gap_uc=mip_gap_uc,
        threads=threads;
        kwargs...
    )
end

"""
    validate_powersimulations_compatibility()

Check PowerSimulations.jl version compatibility.
"""
function validate_powersimulations_compatibility()
    @info "Validating PowerSimulations.jl compatibility..."
    
    try
        # Check that we can create basic components
        template = ProblemTemplate()
        set_network_model!(template, NetworkModel(PSI.CopperPlatePowerModel))
        
        # Check modern formulations exist
        formulations_to_check = [
            PSI.ThermalBasicDispatch,
            PSI.ThermalBasicUnitCommitment,
            PSI.RenewableFullDispatch,
            PSI.StaticPowerLoad,
            PSI.StaticBranch
        ]
        
        for formulation in formulations_to_check
            try
                # Just check the formulation exists
                typeof(formulation)
            catch e
                @error "Missing formulation: $formulation"
                return false
            end
        end
        
        @info "✓ PowerSimulations v0.30.2 compatibility confirmed"
        return true
        
    catch e
        @error "PowerSimulations compatibility check failed: $e"
        return false
    end
end

"""
    create_minimal_test_problem(sys::System)

Create minimal test problem for debugging.
"""
function create_minimal_test_problem(sys::System)
    @info "Creating minimal test problem for debugging..."
    
    template = ProblemTemplate()
    set_network_model!(template, NetworkModel(PSI.CopperPlatePowerModel))
    
    # Only add essential components
    thermal_gens = get_components(ThermalStandard, sys)
    if !isempty(thermal_gens)
        set_device_model!(template, ThermalStandard, PSI.ThermalBasicDispatch)
        @info "  ✓ Added ThermalStandard (ED only)"
    end
    
    loads = get_components(PowerLoad, sys)
    if !isempty(loads)
        set_device_model!(template, PowerLoad, PSI.StaticPowerLoad)
        @info "  ✓ Added PowerLoad"
    end
    
    # Skip renewables and complex components for minimal test
    @info "✓ Minimal template created (thermal + loads only)"
    return template
end

# Export main functions
export run_operations_problem, run_operations_suite, load_system_and_run_operations,
       run_simulation, demo_workflow, run_operations_on_exported_system,
       setup_optimizer, create_operations_problem_template, 
       validate_powersimulations_compatibility, create_minimal_test_problem

# Main execution when run directly
if abspath(PROGRAM_FILE) == @__FILE__
    @info "PowerSimulations.jl Framework v0.30.2"
    @info "Compatible with:"
    @info "  • PowerSimulations v0.30.2"
    @info "  • PowerSystems v4.6.2" 
    @info "  • InfrastructureSystems v2.6.0"
    @info ""
    @info "Key Functions:"
    @info "  • run_operations_problem(sys, \"economic_dispatch\")"
    @info "  • run_operations_suite(sys, use_default_templates=true)"
    @info "  • load_system_and_run_operations(\"path/to/system.json\")"
    @info "  • validate_powersimulations_compatibility()"
end
#!/usr/bin/env julia

"""
Test what data format Deterministic actually expects
"""

using PowerSystems
using TimeSeries
using Dates
using DataStructures
using InfrastructureSystems


function test_data_formats()
    println("Testing different data formats for Deterministic...")
    
    # Create test data
    timestamps = DateTime(2024,1,1):Hour(1):DateTime(2024,1,1,5)
    values = [100.0, 110.0, 120.0, 115.0, 105.0, 95.0]
    
    # Format 1: Plain Dict{DateTime, Float64}
    println("\n1. Testing: Dict{DateTime, Float64}")
    try
        data_dict = Dict(zip(timestamps, values))
        forecast = Deterministic(
            name = "max_active_power",
            data = data_dict,
            resolution = Hour(1)
        )
        println("✓ SUCCESS! Dict{DateTime, Float64} works")
        return "dict_datetime_float"
    catch e
        println("✗ FAILED: $e")
    end
    
    # Format 2: SortedDict
    println("\n2. Testing: SortedDict{DateTime, Float64}")
    try

        data_dict = SortedDict(zip(timestamps, values))
        forecast = Deterministic(
            name = "max_active_power",
            data = data_dict,
            resolution = Hour(1)
        )
        println("✓ SUCCESS! SortedDict works")
        return "sorted_dict"
    catch e
        println("✗ FAILED: $e")
    end
    
    # Format 3: Vector of values with separate timestamps
    println("\n3. Testing: Just the values array")
    try
        forecast = Deterministic(
            name = "max_active_power",
            data = values,  # Just the values
            resolution = Hour(1)
        )
        println("✓ SUCCESS! Values array works")
        return "values_array"
    catch e
        println("✗ FAILED: $e")
    end
    
    # Format 4: Check what convert_data expects
    println("\n4. Checking convert_data function...")
    try
        methods(InfrastructureSystems.convert_data)
    catch e
        println("Can't access convert_data methods: $e")
    end
    
    # Format 5: Try the old-style constructor patterns
    println("\n5. Testing: Old-style constructors")
    
    # Check what RawTimeSeries expects  
    try
        # This might be what we need to create first
        println("Checking RawTimeSeries...")
        raw_ts = InfrastructureSystems.RawTimeSeries()
        println("RawTimeSeries type: $(typeof(raw_ts))")
    catch e
        println("RawTimeSeries failed: $e")
    end
    
    return "none_worked"
end

# Run the test
result = test_data_formats()
println("\nResult: $result")
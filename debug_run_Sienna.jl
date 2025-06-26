#!/usr/bin/env julia

# Minimal test to check if everything works

println("Testing Julia environment...")

# Test 1: Package loading
try  
    using PowerSystems, PowerSimulations, HiGHS
    println("✓ Packages loaded successfully")
catch e
    println("✗ Package loading failed: $e")
    exit(1)
end

# Test 2: File inclusion
try
    include("build_sienna_system.jl")
    println("✓ build_sienna_system.jl included successfully")
catch e
    println("✗ File inclusion failed: $e")
    exit(1)
end

# Test 3: Directory access
data_dir = expanduser("~/showcase/pypsa-rsa/networks/sienna/TEST/dispatch_2030/")
println("Checking directory: $data_dir")

if isdir(data_dir)
    println("✓ Directory exists")
    files = readdir(data_dir)
    println("Found $(length(files)) files:")
    for file in files[1:min(3, end)]
        println("  - $file")
    end
else
    println("✗ Directory not found")
    println("Current working directory: $(pwd())")
    exit(1)
end

# Test 4: Try building system
try
    println("Testing system build...")
    sys = build_sienna_system(data_dir; validate_system=false, load_timeseries=false)
    println("✓ System built successfully!")
    println("System has $(length(get_components(Bus, sys))) buses")
catch e
    println("✗ System build failed: $e")
end

println("Test complete!")
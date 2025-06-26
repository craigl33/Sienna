"""
Main Julia file for building and running PowerSystems.jl and PowerSimulations.jl

Example  showing how to use the Sienna system builder as a module
"""

# Include the system builder (like Python's import)
include("build_sienna_system.jl")

# Now you can use all the functions from build_sienna_system.jl

function run_my_analysis()
    # Set your data directory
    data_dir = "~/showcase/pypsa-rsa/networks/sienna/TEST/dispatch_2030/"
    
    # Build the system with custom parameters
    println("Building Sienna system...")
    sys = build_sienna_system(data_dir; 
                             validate_system=true, 
                             load_timeseries=true,
                             base_power=100.0)
    
    # Save the system
    save_system(sys, "my_pypsa_system.json")
    
    # Run custom analysis
    println("Running economic dispatch...")
    results = run_economic_dispatch_example(sys; horizon=168)  # One week
    
    if results !== nothing
        # Export results
        export_results(results, "my_results")
        
        # Do additional analysis
        analyze_dispatch_results(results, sys)
    end
    
    return sys, results
end

function analyze_dispatch_results(results, sys)
    """Custom analysis of dispatch results"""
    println("Analyzing results...")
    
    # Get generator dispatch
    try
        gen_results = read_realized_variables(results, names=[:ActivePowerVariable])
        println("Successfully extracted generator dispatch data")
        
        # You can add custom analysis here
        # - Plot generation by technology
        # - Calculate costs
        # - Analyze renewable curtailment
        # etc.
        
    catch e
        println("Could not extract detailed results: $e")
    end
end

# Run the analysis
if abspath(PROGRAM_FILE) == @__FILE__
    sys, results = run_my_analysis()
    println("Analysis complete!")
end
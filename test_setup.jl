using PowerSystems
using PowerSystemCaseBuilder
using PowerFlows

# Load a test system
sys = build_system(PSISystems, "c_sys5_pglib")
println("System loaded: $(get_name(sys))")
println("Buses: $(length(get_components(Bus, sys)))")
println("Loads: $(length(get_components(PSY.StaticLoad, sys)))")
println("Generators: $(length(get_components(PSY.PowerPlant, sys)))")

# Solve the power flow
result = run_acpf(sys, PowerFlowPowerModels)
println("Power flow solved: $(result.solution.status)")

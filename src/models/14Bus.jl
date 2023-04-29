using PowerSystems
using NLsolve
const PSY = PowerSystems
using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics

include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
include(joinpath(dirname(@__FILE__), "customControllers.jl"))
############### Data Network ########################
sys_dir = joinpath(dirname(@__FILE__), "14Bus.raw")
sys = System(sys_dir, runchecks = false)
set_units_base_system!(sys, "DEVICE_BASE")

############### Data Dynamic devices ########################
function inv_gfoll(static_device)
    return PSY.DynamicInverter(
        get_name(static_device),
        1.0, #ω_ref
        converter_low_power(), #converter
        outer_control_gfoll(), #outercontrol
        current_mode_inner(), #inner_control
        dc_source_lv(),
        reduced_pll(),
        filt_gfoll(),
    ) #pss
end

function dyn_gen_marconato_tg(generator)
    return PSY.DynamicGenerator(
        name = get_name(generator),
        ω_ref = 1.0, # ω_ref,
        machine = machine_marconato(), #machine
        shaft = shaft_no_damping(), #shaft
        avr = avr_type2(), #avr
        prime_mover = tg_type1(), #tg
        pss = pss_none(),
    )
end

# Add dynamic generators to the system (each gen is linked through a static one)
for g in get_components(Generator, sys)
    if get_number(get_bus(g)) in [1, 2, 8]
        case_gen = dyn_gen_marconato_tg(g)
        add_component!(sys, case_gen, g)
    elseif get_number(get_bus(g)) in [3, 6]
        case_gen=inv_gfoll(g)
        # case_gen=add_damping_controller(g)
        add_component!(sys, case_gen, g)
    end
end

IBR = collect(get_components(DynamicInverter, sys))
SG = collect(get_components(DynamicGenerator, sys))
loads = collect(get_components(StandardLoad, sys))

set_units_base_system!(sys, "DEVICE_BASE")

set_kp_pll!(IBR[1].freq_estimator, 0.20) # Change the proportional gain on the PLL 2.0
set_kp_pll!(IBR[2].freq_estimator, 0.20) # Change the proportional gain on the PLL 2.0
set_kp_pll!(IBR[3].freq_estimator, 0.20) # Change the proportional gain on the PLL 2.0


for l in get_components(StandardLoad, sys)
    set_impedance_active_power!(l, get_constant_active_power(l))
    set_impedance_reactive_power!(l, get_constant_reactive_power(l))
    set_constant_active_power!(l, 0.0)
    set_constant_reactive_power!(l, 0.0)
end


path = (joinpath(pwd(), "src/outputLog"))

for b in get_components(Line, sys)
    if get_name(b) != "Bus 02-Bus 05-i_1"
        dyn_branch = PowerSystems.DynamicBranch(b)
        add_component!(sys, dyn_branch)
    end
end


l_device = get_component(ElectricLoad, sys, "load131")
l_change = LoadChange(0.5, l_device, :P_ref_impedance, 0.2)
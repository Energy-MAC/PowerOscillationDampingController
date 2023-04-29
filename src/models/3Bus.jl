using PowerSystems
using NLsolve
const PSY = PowerSystems

include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
include(joinpath(dirname(@__FILE__), "customControllers.jl"))
############### Data Network ########################
sys_dir = joinpath(dirname(@__FILE__), "ThreeBusMultiLoad.raw")
sys = System(sys_dir, runchecks = false)
set_units_base_system!(sys, "DEVICE_BASE")


function add_battery(sys, battery_name, bus_name, capacity, P, Q)
    return GenericBattery(
        name = battery_name,
        bus = get_component(Bus, sys, bus_name),
        available = true,
        prime_mover = PrimeMovers.BA,
        active_power = P,
        reactive_power = Q,
        rating = 1.1,
        base_power = capacity,
        initial_energy = 50.0,
        state_of_charge_limits = (min = 5.0, max = 100.0),
        input_active_power_limits = (min = 0.0, max = 1.0),
        output_active_power_limits = (min = 0.0, max = 1.0),
        reactive_power_limits = (min = -1.0, max = 1.0),
        efficiency = (in = 0.80, out = 0.90),
    )
end



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

for g in get_components(Generator, sys)
    if get_number(get_bus(g)) == 101
        case_gen = dyn_gen_marconato_tg(g)
        add_component!(sys, case_gen, g)
    elseif get_number(get_bus(g)) == 102
        case_gen=inv_gfoll(g)
        add_component!(sys, case_gen, g)
    end
end

SG = get_component(Generator, sys, "generator-101-1")
Gf = get_component(Generator, sys, "generator-102-1")
set_units_base_system!(sys, "DEVICE_BASE")
SG.dynamic_injector.base_power=100
Gf.dynamic_injector.base_power=100

for l in get_components(StandardLoad, sys)
    set_impedance_active_power!(l, get_constant_active_power(l))
    set_impedance_reactive_power!(l, get_constant_reactive_power(l))
    set_constant_active_power!(l, 0.0)
    set_constant_reactive_power!(l, 0.0)
end

l_device = get_component(ElectricLoad, sys, "load1031")
l_change = LoadChange(0.5, l_device, :P_ref, 1.6)


set_kp_pll!(Gf.dynamic_injector.freq_estimator, 0.2)
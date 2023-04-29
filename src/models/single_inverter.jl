using PowerSystems
using NLsolve
const PSY = PowerSystems

include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
omib_file_dir = joinpath(dirname(@__FILE__), "OMIB.raw")
sys = System(omib_file_dir, runchecks = false)
add_source_to_ref(sys)
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

#Attach dynamic generator. Currently use PSS/e format based on bus #.
gen = [g for g in get_components(Generator, sys)][1]
case_gen=inv_gfoll(gen)
add_component!(sys, case_gen, gen)
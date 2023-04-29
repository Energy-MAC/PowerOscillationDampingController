#Define the complete struct
mutable struct pssReactivePower <: PSY.ReactivePowerControl
        "Proportional Gain"
        Kp_q::Float64
        "Integral Gain"
        Ki_q::Float64
        "filter frequency cutoff"
        ωf::Float64
        "stabilizing gain"
        K1::Float64
        "stabilizing gain"
        K2::Float64
        "stabilizing gain"
        K3::Float64
        "stabilizing gain"
        K4::Float64
        "stabilizing gain"
        K5::Float64
        "stabilizing gain"
        K6::Float64
        "stabilizing gain"
        K7::Float64
        "stabilizing gain"
        K8::Float64
        "Washout time constant"
        Tw::Float64
        "Lead-lag time constant 1"
        T1::Float64
        "Lead-lag time constant 2"
        T2::Float64
        "Lead-lag time constant 3"
        T3::Float64
        "Lead-lag time constant 4"
        T4::Float64
        "Voltage Set-point"
        V_ref::Float64
        "Reactive Power Set-point"
        Q_ref::Float64
        ext::Dict{String, Any}
        "The states of the pssReactivePower model are:
        σq_oc: Integrator state of the PI Controller,
        q_oc: Measured reactive power of the inverter model
        v1: Output of the washout filter
        v2: Output of first lead-lag controller
        v3: Output of second lead-lag controller"
        states::Vector{Symbol}
        "pssReactiveProbe has five states"
        n_states::Int
end

function pssReactivePower(Kp_q, Ki_q, ωf, K1, K2, K3, K4, K5, K6, K7, K8, Tw, T1, T2, T3, T4, V_ref=1.0, Q_ref=1.0, ext=Dict{String, Any}(), )
    pssReactivePower(Kp_q, Ki_q, ωf, K1, K2, K3, K4, K5, K6, K7, K8, Tw, T1, T2, T3, T4, V_ref, Q_ref, ext, [:σq_oc, :q_oc, :v1, :v2, :v3], 5, )
end

function pssReactivePower(; Kp_q, Ki_q, ωf, K1, K2, K3, K4, K5, K6, K7, K8, Tw, T1, T2, T3, T4, V_ref=1.0, Q_ref=1.0, ext=Dict{String, Any}(), states=[:σq_oc, :q_oc, :v1, :v2, :v3], n_states=5, )
    pssReactivePower(Kp_q, Ki_q, ωf, K1, K2, K3, K4, K5, K6, K7, K8, Tw, T1, T2, T3, T4, V_ref, Q_ref, ext, states, n_states, )
end

"""Get [`pssReactivePower`](@ref) `Kp_q`."""
get_Kp_q(value::pssReactivePower) = value.Kp_q
"""Get [`pssReactivePower`](@ref) `Ki_q`."""
get_Ki_q(value::pssReactivePower) = value.Ki_q
"""Get [`pssReactivePower`](@ref) `ωf`."""
get_ωf(value::pssReactivePower) = value.ωf
"""Get [`pssReactivePower`](@ref) `ext`."""
get_ext(value::pssReactivePower) = value.ext
"""Get [`pssReactivePower`](@ref) `states`."""
get_states(value::pssReactivePower) = value.states
"""Get [`pssReactivePower`](@ref) `n_states`."""
get_n_states(value::pssReactivePower) = value.n_states
"""Get [`pssReactivePower`](@ref) `K1`."""
get_K1(value::pssReactivePower) = value.K1
"""Get [`pssReactivePower`](@ref) `K2`."""
get_K2(value::pssReactivePower) = value.K2
"""Get [`pssReactivePower`](@ref) `K3`."""
get_K3(value::pssReactivePower) = value.K3
"""Get [`pssReactivePower`](@ref) `K4`."""
get_K4(value::pssReactivePower) = value.K4
"""Get [`pssReactivePower`](@ref) `K5`."""
get_K5(value::pssReactivePower) = value.K5
"""Get [`pssReactivePower`](@ref) `K6`."""
get_K6(value::pssReactivePower) = value.K6
"""Get [`pssReactivePower`](@ref) `K7`."""
get_K7(value::pssReactivePower) = value.K7
"""Get [`pssReactivePower`](@ref) `K8`."""
get_K8(value::pssReactivePower) = value.K8
"""Get [`pssReactivePower`](@ref) `Tw`."""
get_Tw(value::pssReactivePower) = value.Tw
"""Get [`pssReactivePower`](@ref) `T1`."""
get_T1(value::pssReactivePower) = value.T1
"""Get [`pssReactivePower`](@ref) `T2`."""
get_T2(value::pssReactivePower) = value.T2
"""Get [`pssReactivePower`](@ref) `T3`."""
get_T3(value::pssReactivePower) = value.T3
"""Get [`pssReactivePower`](@ref) `T4`."""
get_T4(value::pssReactivePower) = value.T4

PSY.get_Q_ref(value::pssReactivePower) = value.Q_ref

PSY.set_Q_ref!(value::pssReactivePower, val) = value.Q_ref = val
"""Get [`pssReactivePower`](@ref) `ωz`."""

PSY.get_V_ref(value::pssReactivePower) = value.V_ref

PSY.set_V_ref!(value::pssReactivePower, val) = value.V_ref = val
"""Get [`pssReactivePower`](@ref) `ωz`."""

set_K1(value::pssReactivePower, val) = value.K1 = val
set_K2(value::pssReactivePower, val) = value.K2 = val
set_K3(value::pssReactivePower, val) = value.K3 = val
set_K4(value::pssReactivePower, val) = value.K4 = val
set_K5(value::pssReactivePower, val) = value.K5 = val
set_K6(value::pssReactivePower, val) = value.K6 = val
set_K7(value::pssReactivePower, val) = value.K7 = val
set_K8(value::pssReactivePower, val) = value.K8 = val

set_Tw(value::pssReactivePower, val) = value.Tw = val
set_T1(value::pssReactivePower, val) = value.T1 = val
set_T2(value::pssReactivePower, val) = value.T2 = val
set_T3(value::pssReactivePower, val) = value.T3 = val
set_T4(value::pssReactivePower, val) = value.T4 = val

function PSID.mdl_outer_ode!(
    device_states,
    output_ode,
    inner_vars,
    _,
    dynamic_device::PSID.DynamicWrapper{PSY.DynamicInverter{
            C,
            PSY.OuterControl{ActivePowerPI, pssReactivePower},
            IC,
            DC,
            P,
            F,
        },
        },
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = PSID.get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{ActivePowerPI, pssReactivePower},
    )
    
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]
    Ir_filter = device_states[external_ix[3]]
    Ii_filter = device_states[external_ix[4]]

    γd_ic = device_states[8]
    γq_ic = device_states[9]
    vq_pll = device_states[10]
    ε_pll = device_states[11]

    #Obtain inner variables for component
    θ_pll = inner_vars[PSID.θ_freq_estimator_var]
    ω_pll = inner_vars[PSID.ω_freq_estimator_var]
    
    #Get Active Power Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    active_power_control = outer_control.active_power_control#PSY.get_active_power(outer_control)
    Kp_p = get_Kp_p(active_power_control) #Proportional Gain
    Ki_p = get_Ki_p(active_power_control) #Integral Gain
    ωz = get_ωz(active_power_control) #Frequency cutoff frequency


    #Get Reactive Power Controller parameters
    reactive_power_control = outer_control.reactive_power_control#PSY.get_reactive_power(outer_control)
    Kp_q = get_Kp_q(reactive_power_control) #Proportional Gain
    Ki_q = get_Ki_q(reactive_power_control) #Integral Gain
    ωf = get_ωf(reactive_power_control) #Reactive power filter cutoff frequency
    K1 = get_K1(reactive_power_control) #
    K2 = get_K2(reactive_power_control) #
    K3 = get_K3(reactive_power_control) #
    K4 = get_K4(reactive_power_control) #
    K5 = get_K5(reactive_power_control) #
    K6 = get_K6(reactive_power_control) #
    K7 = get_K7(reactive_power_control) #
    K8 = get_K8(reactive_power_control) #
    Tw = get_Tw(reactive_power_control) #
    T1 = get_T1(reactive_power_control) #
    T2 = get_T2(reactive_power_control) #
    T3 = get_T3(reactive_power_control) #
    T4 = get_T4(reactive_power_control) #

    #Obtain external parameters
    p_ref = PSID.get_P_ref(dynamic_device)
    q_star = PSID.get_Q_ref(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = PSID.get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{ActivePowerPI, pssReactivePower},
    )


    #Define internal states for outer control
    internal_states = @view device_states[local_ix]
    σp_oc = internal_states[1]
    p_oc = internal_states[2]
    σq_oc = internal_states[3]
    q_oc = internal_states[4]
    v1 = internal_states[5]
    v2 = internal_states[6]
    v3 = internal_states[7]

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    pss_in = K1*σp_oc + K2*p_oc + K3*σq_oc  + K4*q_oc  + K5*γd_ic + K6*γq_ic  + K7*vq_pll + K8*ε_pll  

    Δq_ref = v3 + (T3/T4)*(v2 + (T1/T2)*(pss_in+v1))
    q_ref = q_star + Δq_ref
    

    #Compute 7 states ODEs
    output_ode[local_ix[1]] = p_ref - p_oc
    output_ode[local_ix[2]] = ωz * (p_elec_out - p_oc)
    output_ode[local_ix[3]] = q_ref - q_oc
    output_ode[local_ix[4]] = ωf * (q_elec_out - q_oc)
    output_ode[local_ix[5]] = (-pss_in -v1)/Tw
    output_ode[local_ix[6]] = ((1-(T1/T2))*(pss_in + v1) -v2)/T2
    output_ode[local_ix[7]] = ( (1- (T3/T4) ) * (v2+( (T1/T2) * (pss_in + v1) ) )-v3)/T4

    #Update inner vars
    inner_vars[PSID.θ_oc_var] = θ_pll
    inner_vars[PSID.ω_oc_var] = ω_pll
    inner_vars[PSID.Iq_oc_var] = Kp_p * (p_ref - p_oc) + Ki_p * σp_oc
    inner_vars[PSID.Id_oc_var] = Kp_q * (q_ref - q_oc) + Ki_q * σq_oc
end

function PSID.initialize_outer!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSID.DynamicWrapper{PSY.DynamicInverter{
            C,
            PSY.OuterControl{PSY.ActivePowerPI, pssReactivePower},
            IC,
            DC,
            P,
            F,
        },
        },
    inner_vars,
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = PSID.get_input_port_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerPI, pssReactivePower},
    )
    Vr_filter = device_states[external_ix[1]]
    Vi_filter = device_states[external_ix[2]]
    Ir_filter = device_states[external_ix[3]]
    Ii_filter = device_states[external_ix[4]]
    Ir_cnv = device_states[external_ix[5]]
    Ii_cnv = device_states[external_ix[6]]

    #Obtain additional expressions
    θ0_oc = inner_vars[PSID.θ_freq_estimator_var]
    I_dq_cnv = PSID.ri_dq(θ0_oc + pi / 2) * [Ir_cnv; Ii_cnv]
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Get Outer Controller parameters
    outer_control = PSY.get_outer_control(dynamic_device)
    # active_power_control = PSY.get_active_power(outer_control)
    Ki_p = get_Ki_p(outer_control.active_power_control) #Integral Gain
    # reactive_power_control = PSY.get_reactive_power(outer_control)
    Ki_q = get_Ki_q(outer_control.reactive_power_control) #Integral Gain

    #Update inner_vars
    inner_vars[PSID.P_ES_var] = p_elec_out
    #Update states
    outer_ix = PSID.get_local_state_ix(
        dynamic_device,
        PSY.OuterControl{PSY.ActivePowerPI, pssReactivePower},
    )

    outer_states = @view device_states[outer_ix]
    outer_states[1] = I_dq_cnv[PSID.q] / Ki_p #σp_oc
    outer_states[2] = p_elec_out #p_oc
    outer_states[3] = I_dq_cnv[PSID.d] / Ki_q #σq_oc
    outer_states[4] = q_elec_out #q_oc
    outer_states[5] = 0.0
    outer_states[6] = 0.0
    outer_states[7] = 0.0

    #Update inner vars
    inner_vars[PSID.θ_oc_var] = θ0_oc
    inner_vars[PSID.ω_oc_var] = PSID.get_ω_ref(dynamic_device)
    inner_vars[PSID.Id_oc_var] = I_dq_cnv[PSID.d]
    inner_vars[PSID.Iq_oc_var] = I_dq_cnv[PSID.q]
    #Update Q_ref. Initialization assumes q_ref = q_elec_out from PF solution
    PSID.set_P_ref(dynamic_device, p_elec_out)
    # PSY.set_P_ref!(PSY.get_active_power(PSY.get_outer_control(dynamic_device)), p_elec_out)
    PSY.set_P_ref!(outer_control.active_power_control, p_elec_out)
    
    PSID.set_Q_ref(dynamic_device, q_elec_out)
    # PSY.set_Q_ref!(PSY.get_reactive_power(PSY.get_outer_control(dynamic_device)),q_elec_out)
    PSY.set_Q_ref!(outer_control.reactive_power_control, q_elec_out)
end

function outer_damping_control()
    function activePI()
        return ActivePowerPI(Kp_p = 2.0, Ki_p = 30.0, ωz = 0.132 * 2 * pi * 50)
    end
    function reactivePSS()
        return pssReactivePower(Kp_q = 2.0, Ki_q = 30.0, ωf = 0.132 * 2 * pi * 50.0, K1=0.0, K2=0.0, K3=0.0, K4=0.0, K5=0.0, K6=0.0, K7=0.0, K8=0.0, Tw=1.0, T1=1.0, T2=1.0, T3=1.0, T4=1.0)
    end
    return OuterControl(activePI(), reactivePSS())
end

function add_damping_controller(static_device)
    return PSY.DynamicInverter(
        get_name(static_device),
        1.0, #ω_ref
        converter_low_power(), #converter
        outer_damping_control(), #outercontrol
        current_mode_inner(), #inner_control
        dc_source_lv(),
        reduced_pll(),
        filt_gfoll(),
    ) #pss
end






























# #Define the complete struct
# mutable struct ssidActiveProbe <: PSY.ActivePowerControl
#     "Proportional Gain"
#     Kp_p::Float64
#     "Integral Gain"
#     Ki_p::Float64
#     "filter frequency cutoff"
#     ωz::Float64
#     "Active Loop Status"
#     c1::Float64
#     "Probing Status"
#     c2::Float64
#     "Probe input"
#     u::Float64
#     "Reference Power Set-point"
#     P_ref::Float64
#     ext::Dict{String, Any}
#     "The states of the ActivePowerPI model are:
#     σp_oc: Integrator state of the PI Controller,
#     p_oc: Measured active power of the inverter model"
#     states::Vector{Symbol}
#     "ssidActiveProbe has two states"
#     n_states::Int
# end

# function ssidActiveProbe(Kp_p, Ki_p, ωz, c1, c2, u, P_ref=1.0, ext=Dict{String, Any}(), )
#     ssidActiveProbe(Kp_p, Ki_p, ωz, c1, c2, u, P_ref, ext, [:σp_oc, :p_oc], 2, )
# end

# function ssidActiveProbe(; Kp_p, Ki_p, ωz, c1, c2, u, P_ref=1.0, ext=Dict{String, Any}(), states=[:σp_oc, :p_oc], n_states=2, )
#     ssidActiveProbe(Kp_p, Ki_p, ωz, c1, c2, u, P_ref, ext, [:σp_oc, :p_oc], 2, )
# end


# mutable struct ssidReactiveProbe <: PSY.ReactivePowerControl
#     "Proportional Gain"
#     Kp_q::Float64
#     "Integral Gain"
#     Ki_q::Float64
#     "filter frequency cutoff"
#     ωf::Float64
#     "Active Loop Status"
#     c1::Float64
#     "Probing Status"
#     c2::Float64
#     "Voltage Set-point"
#     V_ref::Float64
#     "Reactive Power Set-point"
#     Q_ref::Float64
#     ext::Dict{String, Any}
#     "The states of the ssidReactiveProbe model are:
# 	σq_oc: Integrator state of the PI Controller,
# 	q_oc: Measured reactive power of the inverter model"
#     states::Vector{Symbol}
#     "ssidReactiveProbe has two states"
#     n_states::Int
# end

# function ssidReactiveProbe(Kp_q, Ki_q, ωf, c1, c2, V_ref=1.0, Q_ref=1.0, ext=Dict{String, Any}(), )
#     ssidReactiveProbe(Kp_q, Ki_q, ωf, c1, c2, V_ref, Q_ref, ext, [:σq_oc, :q_oc], 2, )
# end

# function ssidReactiveProbe(; Kp_q, Ki_q, ωf, c1, c2, V_ref=1.0, Q_ref=1.0, ext=Dict{String, Any}(), states=[:σq_oc, :q_oc], n_states=2, )
#     ssidReactiveProbe(Kp_q, Ki_q, ωf, c1, c2, V_ref, Q_ref, ext, states, n_states, )
# end



# """Get [`ssidActiveProbe`](@ref) `Kp_p`."""
# get_Kp_p(value::ssidActiveProbe) = value.Kp_p
# """Get [`ssidActiveProbe`](@ref) `Ki_p`."""
# get_Ki_p(value::ssidActiveProbe) = value.Ki_p
# """Get [`ssidActiveProbe`](@ref) `ωz`."""
# get_ωz(value::ssidActiveProbe) = value.ωz
# """Get [`ssidActiveProbe`](@ref) `ωz`."""
# get_c1(value::ssidActiveProbe) = value.c1
# """Get [`ssidActiveProbe`](@ref) `ωz`."""
# get_c2(value::ssidActiveProbe) = value.c2
# """Get [`ssidActiveProbe`](@ref) `ωz`."""
# get_u(value::ssidActiveProbe) = value.u
# """Get [`ssidActiveProbe`](@ref) `ωr`."""
# PSY.get_P_ref(value::ssidActiveProbe) = value.P_ref

# PSY.set_P_ref!(value::ssidActiveProbe, val) = value.P_ref = val
# """Get [`ssidActiveProbe`](@ref) `ωz`."""
# set_u(value::ssidActiveProbe, val) = value.u = val

# set_c1(value::ssidActiveProbe, val) = value.c1 = val

# set_c2(value::ssidActiveProbe, val) = value.c2 = val


# """Get [`ssidReactiveProbe`](@ref) `Kp_q`."""
# get_Kp_q(value::ssidReactiveProbe) = value.Kp_q
# """Get [`ssidReactiveProbe`](@ref) `Ki_q`."""
# get_Ki_q(value::ssidReactiveProbe) = value.Ki_q
# """Get [`ssidReactiveProbe`](@ref) `ωf`."""
# get_ωf(value::ssidReactiveProbe) = value.ωf
# """Get [`ssidReactiveProbe`](@ref) `V_ref`."""
# get_V_ref(value::ssidReactiveProbe) = value.V_ref
# """Get [`ssidReactiveProbe`](@ref) `Q_ref`."""
# get_Q_ref(value::ssidReactiveProbe) = value.Q_ref
# """Get [`ssidReactiveProbe`](@ref) `ext`."""
# get_ext(value::ssidReactiveProbe) = value.ext
# """Get [`ssidReactiveProbe`](@ref) `states`."""
# get_states(value::ssidReactiveProbe) = value.states
# """Get [`ssidReactiveProbe`](@ref) `n_states`."""
# get_n_states(value::ssidReactiveProbe) = value.n_states
# """Get [`ssidActiveProbe`](@ref) `ωz`."""
# get_c1(value::ssidReactiveProbe) = value.c1
# """Get [`ssidActiveProbe`](@ref) `ωz`."""
# get_c2(value::ssidReactiveProbe) = value.c2

# PSY.get_Q_ref(value::ssidReactiveProbe) = value.Q_ref

# PSY.set_Q_ref!(value::ssidReactiveProbe, val) = value.Q_ref = val
# """Get [`ssidActiveProbe`](@ref) `ωz`."""

# PSY.get_V_ref(value::ssidReactiveProbe) = value.V_ref

# PSY.set_V_ref!(value::ssidReactiveProbe, val) = value.V_ref = val
# """Get [`ssidActiveProbe`](@ref) `ωz`."""

# set_c1(value::ssidReactiveProbe, val) = value.c1 = val

# set_c2(value::ssidReactiveProbe, val) = value.c2 = val



# function PSID.mdl_outer_ode!(
#     device_states,
#     output_ode,
#     inner_vars,
#     _,
#     dynamic_device::PSID.DynamicWrapper{PSY.DynamicInverter{
#             C,
#             PSY.OuterControl{ssidActiveProbe, ssidReactiveProbe},
#             IC,
#             DC,
#             P,
#             F,
#         },
#         },
# ) where {
#     C <: PSY.Converter,
#     IC <: PSY.InnerControl,
#     DC <: PSY.DCSource,
#     P <: PSY.FrequencyEstimator,
#     F <: PSY.Filter,
# }

#     #Obtain external states inputs for component
#     external_ix = PSID.get_input_port_ix(
#         dynamic_device,
#         PSY.OuterControl{ssidActiveProbe, ssidReactiveProbe},
#     )
    
#     Vr_filter = device_states[external_ix[1]]
#     Vi_filter = device_states[external_ix[2]]
#     Ir_filter = device_states[external_ix[3]]
#     Ii_filter = device_states[external_ix[4]]

#     #Obtain inner variables for component
#     θ_pll = inner_vars[PSID.θ_freq_estimator_var]
#     ω_pll = inner_vars[PSID.ω_freq_estimator_var]
    
#     #Get Active Power Controller parameters
#     outer_control = PSY.get_outer_control(dynamic_device)
#     active_power_control = PSY.get_active_power(outer_control)
#     Kp_p = get_Kp_p(active_power_control) #Proportional Gain
#     Ki_p = get_Ki_p(active_power_control) #Integral Gain
#     ωz = get_ωz(active_power_control) #Frequency cutoff frequency
#     p1 = get_c1(active_power_control) #Attack Frequency
#     p2 = get_c2(active_power_control) #Attack Frequency
#     u = get_u(active_power_control) #Attack Frequency

#     #Get Reactive Power Controller parameters
#     reactive_power_control = PSY.get_reactive_power(outer_control)
#     Kp_q = get_Kp_q(reactive_power_control) #Proportional Gain
#     Ki_q = get_Ki_q(reactive_power_control) #Integral Gain
#     ωf = get_ωf(reactive_power_control) #Reactive power filter cutoff frequency
#     q1 = get_c1(reactive_power_control) #Attack Frequency
#     q2 = get_c2(reactive_power_control) #Attack Frequency

#     #Obtain external parameters
#     p_ref = PSID.get_P_ref(dynamic_device)
#     q_ref = PSID.get_Q_ref(dynamic_device)

#     #Obtain indices for component w/r to device
#     local_ix = PSID.get_local_state_ix(
#         dynamic_device,
#         PSY.OuterControl{ssidActiveProbe, ssidReactiveProbe},
#     )


#     #Define internal states for outer control
#     internal_states = @view device_states[local_ix]
#     σp_oc = internal_states[1]
#     p_oc = internal_states[2]
#     σq_oc = internal_states[3]
#     q_oc = internal_states[4]

#     #Obtain additional expressions
#     p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
#     q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter
        
#     #Compute 4 states ODEs
#     output_ode[local_ix[1]] = p1*(p_ref - p_oc)
#     output_ode[local_ix[2]] = ωz * (p_elec_out - p_oc)
#     output_ode[local_ix[3]] = q1*(q_ref - q_oc)
#     output_ode[local_ix[4]] = ωf * (q_elec_out - q_oc)

#     #Update inner vars
#     inner_vars[PSID.θ_oc_var] = θ_pll
#     inner_vars[PSID.ω_oc_var] = ω_pll
#     inner_vars[PSID.Iq_oc_var] = p1*(Kp_p * (p_ref - p_oc) + Ki_p * σp_oc) + p2*(Ki_p * σp_oc + u)
#     inner_vars[PSID.Id_oc_var] = q1*(Kp_q * (q_ref - q_oc) + Ki_q * σq_oc) + q2*(Ki_q * σq_oc)
# end

# function PSID.initialize_outer!(
#     device_states,
#     static::PSY.StaticInjection,
#     dynamic_device::PSID.DynamicWrapper{PSY.DynamicInverter{
#             C,
#             PSY.OuterControl{ssidActiveProbe, ssidReactiveProbe},
#             IC,
#             DC,
#             P,
#             F,
#         },
#         },
#     inner_vars,
# ) where {
#     C <: PSY.Converter,
#     IC <: PSY.InnerControl,
#     DC <: PSY.DCSource,
#     P <: PSY.FrequencyEstimator,
#     F <: PSY.Filter,
# }

#     #Obtain external states inputs for component
#     external_ix = PSID.get_input_port_ix(
#         dynamic_device,
#         PSY.OuterControl{ssidActiveProbe, ssidReactiveProbe},
#     )
#     Vr_filter = device_states[external_ix[1]]
#     Vi_filter = device_states[external_ix[2]]
#     Ir_filter = device_states[external_ix[3]]
#     Ii_filter = device_states[external_ix[4]]
#     Ir_cnv = device_states[external_ix[5]]
#     Ii_cnv = device_states[external_ix[6]]

#     #Obtain additional expressions
#     θ0_oc = inner_vars[PSID.θ_freq_estimator_var]
#     I_dq_cnv = PSID.ri_dq(θ0_oc + pi / 2) * [Ir_cnv; Ii_cnv]
#     p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
#     q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

#     #Get Outer Controller parameters
#     outer_control = PSY.get_outer_control(dynamic_device)
#     active_power_control = PSY.get_active_power(outer_control)
#     Ki_p = get_Ki_p(active_power_control) #Integral Gain
#     reactive_power_control = PSY.get_reactive_power(outer_control)
#     Ki_q = get_Ki_q(reactive_power_control) #Integral Gain

#     #Update inner_vars
#     inner_vars[PSID.P_ES_var] = p_elec_out
#     #Update states
#     outer_ix = PSID.get_local_state_ix(
#         dynamic_device,
#         PSY.OuterControl{ssidActiveProbe, ssidReactiveProbe},
#     )
#     outer_states = @view device_states[outer_ix]
#     outer_states[1] = I_dq_cnv[PSID.q] / Ki_p #σp_oc
#     outer_states[2] = p_elec_out #p_oc
#     outer_states[3] = I_dq_cnv[PSID.d] / Ki_q #σq_oc
#     outer_states[4] = q_elec_out #q_oc

#     #Update inner vars
#     inner_vars[PSID.θ_oc_var] = θ0_oc
#     inner_vars[PSID.ω_oc_var] = PSID.get_ω_ref(dynamic_device)
#     inner_vars[PSID.Id_oc_var] = I_dq_cnv[PSID.d]
#     inner_vars[PSID.Iq_oc_var] = I_dq_cnv[PSID.q]
#     #Update Q_ref. Initialization assumes q_ref = q_elec_out from PF solution
#     PSID.set_P_ref(dynamic_device, p_elec_out)
#     PSY.set_P_ref!(PSY.get_active_power(PSY.get_outer_control(dynamic_device)), p_elec_out)
#     PSID.set_Q_ref(dynamic_device, q_elec_out)
#     PSY.set_Q_ref!(PSY.get_reactive_power(PSY.get_outer_control(dynamic_device)),q_elec_out)
# end

# function outer_damping_control()
#     function activePI()
#         return ActivePowerPI(Kp_p = 2.0, Ki_p = 30.0, ωz = 0.132 * 2 * pi * 50)
#     end
#     function reactivePSS()
#         return pssReactivePower(Kp_q = 2.0, Ki_q = 30.0, ωf = 0.132 * 2 * pi * 50.0, K1=0.0, K2=0.0, K3=0.0, K4=0.0, K5=0.0, K6=0.0, K7=0.0, K8=0.0, Tw=1.0, T1=1.0, T2=1.0, T3=1.0, T4=1.0)
#     end
#     return OuterControl(activePI(), reactivePSS())
# end

# function add_damping_controller(static_device)
#     return PSY.DynamicInverter(
#         get_name(static_device),
#         1.0, #ω_ref
#         converter_low_power(), #converter
#         outer_damping_control(), #outercontrol
#         current_mode_inner(), #inner_control
#         dc_source_lv(),
#         reduced_pll(),
#         filt_gfoll(),
#     ) #pss
# end
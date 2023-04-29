function add_damping_controller!(sys, Gf)
    remove_component!(sys, get_dynamic_injector(Gf))
    case_gen=add_damping_controller(Gf)
    add_component!(sys, case_gen, Gf)
    set_kp_pll!(Gf.dynamic_injector.freq_estimator, 0.2)
    Gf.dynamic_injector.base_power=100
end

function build_sim_pss(sys, Gf, optParams, tspan, perturbation)
    IBR_params=zeros(11)
    IBR_params[9]=IBR_params[10]=IBR_params[11]=1.0
    if length(sys.bus_numbers)==3
        change_pss_gains!(Gf, IBR_params)
    end

    sim = Simulation!(ResidualModel, sys, path, tspan, all_lines_dynamic = true)
    x0=sim.x0_init

    if length(sys.bus_numbers)==3
        change_pss_gains!(Gf, optParams)
    end

    sim = Simulation!(ResidualModel, sys, path, tspan, initialize_simulation = false, initial_conditions = x0) 
    execute!(sim, Sundials.IDA())
    results = read_results(sim)
    x0=deepcopy(results.solution.u[end])


    sim = Simulation!(ResidualModel, sys, path, tspan, perturbation, initialize_simulation = false, initial_conditions = x0, all_lines_dynamic = true)

    return sim 
end

function get_exp_data(ibr, exp_time, tsave, res)
    
    gf_states=collect(values(res.global_index[ibr]))[1:8]
    exp_window=exp_time:tsave:exp_time+20
    y=Array(res.solution(exp_window,idxs=gf_states))
    y=y .- mean(y, dims=2)

    return y
end 

function build_rom(data, n, u, tsave)
    
    if length(data) < 100
        data_mean = sum(data)/length(data)
    else
        data_mean=data
    end

    d = iddata(data_mean, u, tsave)
    est_sys = subspaceid(d, n, verbose=true, W=:MOESP, r=100, s1=10, s2=10, zeroD=true, stable = false)
    rom = d2c(est_sys.sys, :zoh)
    
    svc=zeros(8)
    for i = 1:8
        svc[i] = est_sys.s.S[i+1]^2 + (i*size(data_mean)[1] + i*(size(data_mean)[1]+1) + size(data_mean)[1])*log(length(u))/length(u) 
    end

    return rom, svc
end

function tune_pss(rom)

    minParams=[0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, -0.5, 0.2, 0.01, 0.01]
    maxParams=[0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 0.5]

    nParticles=300
    iters=200
    c1 = 1.496 
    c2 = 1.496 
    w = 0.7298

    optParams, minValue=optimizeGains_V2(rom, minParams, maxParams, nParticles, c1, c2, w, iters)

    return optParams, minValue
end 

function build_Acl_V2(rom, paramVec)
    K1, K2, K3, K4, K5, K6, K7, K8, Tw, T1, T2= paramVec
    
    Apss = [-1/Tw 0 0 
            (T2-T1)/T2^2 -1/T2 0
            ((T2-T1)/T2^2)*(T1/T2) ((T2-T1)/T2^2)*T2 -1/T2]

    Bpss = [K1 K2 K3 K4 K5 K6 K7 K8
            K1 K2 K3 K4 K5 K6 K7 K8
            K1 K2 K3 K4 K5 K6 K7 K8]

    Bpss[1,:] = (-1/Tw)*Bpss[1,:]
    Bpss[2,:] = ((T2-T1)/T2^2)*Bpss[2,:]
    Bpss[3,:] = ((T2-T1)/T2^2)*(T1/T2)*Bpss[3,:]

    Cpss = [(T1*T1)/(T2*T2) T1/T2 1]

    Dpss = (T1*T1)/(T2*T2)*[K1 K2 K3 K4 K5 K6 K7 K8]
    
    A_cl = [rom.A + rom.B*Dpss*rom.C      rom.B*Cpss
            Bpss*rom.C                        Apss]

    return A_cl
end

function calLoss(X)
    A_cl = build_Acl(X)
    eigs=eigvals(A_cl)
    damping=-real(eigs)./sqrt.(real(eigs).^2 .+ imag(eigs).^2)
    
    damping_diff = 0.7 .-damping
    
    R_eigVec = eigvecs(A_cl)
    L_eigVec = inv(R_eigVec)
    par_fact = zeros(size(L_eigVec))
    for (ix, eigs) in enumerate(eigs)
        den = sum(abs.(L_eigVec[:, ix]) .* abs.(R_eigVec[ix, :]))
        par_fact[ix, :] = abs.(L_eigVec[:, ix]) .* abs.(R_eigVec[ix, :]) ./ den
    end

    sum_part = vec(sum(par_fact[1:2, :], dims=1))

    loss = 10*dot(Int.(damping_diff .> 0.0), damping_diff.^2) + sum(sum_part.*real(eigs)) + 1e10*sum(real(eigs) .> 0)
    return loss
end

function calLoss_V2(rom, X)
    A_cl = build_Acl_V2(rom, X)
    eigs=eigvals(A_cl)
    damping=-real(eigs)./sqrt.(real(eigs).^2 .+ imag(eigs).^2)

    damping_diff = 0.7 .-damping
    
    R_eigVec = eigvecs(A_cl)
    L_eigVec = inv(R_eigVec)
    par_fact = zeros(size(L_eigVec))
    for (ix, eigs) in enumerate(eigs)
        den = sum(abs.(L_eigVec[:, ix]) .* abs.(R_eigVec[ix, :]))
        par_fact[ix, :] = abs.(L_eigVec[:, ix]) .* abs.(R_eigVec[ix, :]) ./ den
    end

    sum_part = vec(sum(par_fact[1:2, :], dims=1))

    loss = 10*dot(Int.(damping_diff .> 0.0), damping_diff) + sum(sum_part.*real(eigs)) + 1e10*sum(real(eigs) .> 0)
    return loss
end

function optimizeGains(minParams, maxParams, nParticles, c1, c2, w, iters)
    nParams=length(minParams)
    paramRange=maxParams-minParams
    X = minParams .+ rand(nParams,nParticles) .* paramRange
    V = rand(nParams,nParticles) .* 0.1
    pbest=X
    pbest_obj=vec(mapslices(calLoss, X,  dims = [1]))
    gbest=X[:, argmin(pbest_obj)]
    gbest_obj = minimum(pbest_obj)
    
    for i = 1:iters
        r1, r2 = rand(2)
        V = w * V + c1*r1*(pbest - X) + c2*r2*(gbest .-X)
        X = X + V
        X = transpose(reduce(hcat, [clamp.(X[i, :], minParams[i], maxParams[i]) for i in eachindex(minParams)]))
        obj = vec(mapslices(calLoss, X,  dims = [1]))
        pbest[:, (pbest_obj .> obj)] = X[:, (pbest_obj .> obj)]
        pbest_obj=vec(minimum(hcat(pbest_obj, obj), dims=2))
        gbest=X[:, argmin(pbest_obj)]
        gbest_obj = minimum(pbest_obj)
    end
    return gbest, gbest_obj
end

function optimizeGains_V2(rom, minParams, maxParams, nParticles, c1, c2, w, iters)
    nParams=length(minParams)
    paramRange=maxParams-minParams
    X = minParams .+ rand(nParams,nParticles) .* paramRange
    V = rand(nParams,nParticles) .* 0.1
    pbest=X
    pbest_obj=vec(calLoss_V2.(Ref(rom), eachcol(X)))
    gbest=X[:, argmin(pbest_obj)]
    gbest_obj = minimum(pbest_obj)
    
    for i = 1:iters
        r1, r2 = rand(2)
        V = w * V + c1*r1*(pbest - X) + c2*r2*(gbest .-X)
        X = X + V
        X = transpose(reduce(hcat, [clamp.(X[i, :], minParams[i], maxParams[i]) for i in eachindex(minParams)]))
        obj = vec(calLoss_V2.(Ref(rom), eachcol(X)))
        pbest[:, (pbest_obj .> obj)] = X[:, (pbest_obj .> obj)]
        pbest_obj=vec(minimum(hcat(pbest_obj, obj), dims=2))
        gbest=X[:, argmin(pbest_obj)]
        gbest_obj = minimum(pbest_obj)
    end
    return gbest, gbest_obj
end


function change_pss_gains!(Gf, Params)
    set_K1(Gf.dynamic_injector.outer_control.reactive_power_control, Params[1])
    set_K2(Gf.dynamic_injector.outer_control.reactive_power_control, Params[2])
    set_K3(Gf.dynamic_injector.outer_control.reactive_power_control, Params[3])
    set_K4(Gf.dynamic_injector.outer_control.reactive_power_control, Params[4])
    set_K5(Gf.dynamic_injector.outer_control.reactive_power_control, Params[5])
    set_K6(Gf.dynamic_injector.outer_control.reactive_power_control, Params[6])
    set_K7(Gf.dynamic_injector.outer_control.reactive_power_control, Params[7])
    set_K8(Gf.dynamic_injector.outer_control.reactive_power_control, Params[8])
    set_Tw(Gf.dynamic_injector.outer_control.reactive_power_control, Params[9])
    set_T1(Gf.dynamic_injector.outer_control.reactive_power_control, Params[10])
    set_T3(Gf.dynamic_injector.outer_control.reactive_power_control, Params[10])
    set_T2(Gf.dynamic_injector.outer_control.reactive_power_control, Params[11])
    set_T4(Gf.dynamic_injector.outer_control.reactive_power_control, Params[11])
    
end

function generate_probing_signal(probe_amp, tsave, tstop, ω1, ω2)
    t_control = 0.0:tsave:tstop
    u = probe_amp*sign.(sin.(ω1*t_control .+ ((ω2-ω1)*(t_control) .^ 2)/(2*tstop))) 
    u[end]=0

    return u
end

function append_perturb_vector!(u, perturb_vector, ibr, t_control)

    c=findall(x->x>0, abs.(diff(u))).+1
    Q0=get_Q_ref(ibr.outer_control.reactive_power_control)

    for i in eachindex(c)
        push!(perturb_vector, ControlReferenceChange(t_control[c[i]], ibr, :Q_ref, Q0+u[c[i]]))
    end

    return perturb_vector

end

function generate_perturb_vector(u, ibr, t_control, loads)

    c=findall(x->x>0, abs.(diff(u))).+1
    Q0=get_Q_ref(ibr.outer_control.reactive_power_control)

    n4sid_perturb =  ControlReferenceChange(t_control[c[1]], ibr, :Q_ref, Q0+u[c[1]])

    load_step_times = round.(34*rand(20), digits=3)
    load_changed = rand(loads)
    l_change = LoadChange(load_step_times[1], load_changed, :P_ref, get_impedance_active_power(load_changed)+0.005*randn(1)[1])

    perturb_vector=[n4sid_perturb, l_change]
    for i in 2:length(c)
        push!(perturb_vector, ControlReferenceChange(t_control[c[i]], ibr, :Q_ref, Q0+u[c[i]]))
    end

    for i in 2:length(load_step_times)
        load_changed = rand(loads)
        push!(perturb_vector, LoadChange(load_step_times[i], load_changed, :P_ref, get_impedance_active_power(load_changed)+0.005*randn(1)[1]))
    end

    return perturb_vector

end
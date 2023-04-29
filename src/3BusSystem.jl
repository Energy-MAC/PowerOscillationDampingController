using Pkg
Pkg.activate(".")
using PowerSimulationsDynamics
using DifferentialEquations
using Sundials
using Plots
using PowerSystems 
using DataFrames
using CSV
using Logging
configure_logging(console_level = Logging.Warn)
using ControlSystems
using LinearAlgebra
using ControlSystemIdentification
const PSY = PowerSystems
const PSID = PowerSimulationsDynamics
using Random 
using LaTeXStrings
using Statistics
using Distributions


##################################################
############### LOAD DATA ########################
##################################################
include(joinpath(pwd(), "src/ssid_functions.jl"))
include(joinpath(pwd(), "src/models/3Bus.jl"))

tstop=5.0
tspan=(0.0, tstop)
tsave=0.005
path = (joinpath(pwd(), "src/outputLog"))

sim = Simulation!(ResidualModel, sys, path, (0.0, 3.0), l_change, all_lines_dynamic = true)

small_sig = small_signal_analysis(sim)
full_eigs = small_sig.eigenvalues
part_factors = small_sig.participation_factors;
part_states=Gf.dynamic_injector.states
sum_part = sum(part_factors["generator-102-1"][i] for i in part_states)

execute!(sim, Sundials.IDA())
results = read_results(sim)
apower = PSID.get_activepower_series(results, Gf.name, dt=tsave)
rpower = PSID.get_reactivepower_series(results, Gf.name, dt=tsave)

base_case_load=DataFrame(Time=apower[1], P=apower[2], Q=rpower[2])


# #Chirp
tstop=20.0
probe_amp=0.005
t_control = 0.0:tsave:tstop
ω1=10
ω2=15
u = probe_amp*sign.(sin.(ω1*t_control .+ ((ω2-ω1)*(t_control) .^ 2)/(2*tstop))) # sample a control input for identification

c=findall(x->x>0, abs.(diff(u))).+1

Q0=get_Q_ref(Gf.dynamic_injector.outer_control.reactive_power_control)


n4sid_perturb =  ControlReferenceChange(t_control[c[1]], Gf.dynamic_injector, :Q_ref, Q0+u[c[1]])

perturb_vector=[n4sid_perturb]
for i in 2:length(c)
    push!(perturb_vector, ControlReferenceChange(t_control[c[i]], Gf.dynamic_injector, :Q_ref, Q0+u[c[i]]))
end

sim = Simulation!(ResidualModel, sys, path, (0.0, tstop), perturb_vector, all_lines_dynamic = true)

execute!(sim, Sundials.IDA(), saveat=tsave)

res = read_results(sim)

# Get data for SSID and estimate 

y=get_exp_data(Gf.name, 0.0, tsave, res)
sampling_period=tsave
d = iddata(y, u, sampling_period)

plot(0:tsave:tstop, y[2,:], label="Δp [p.u.]", xlabel="Time [s]", xlim=(0.0, 5.0), ylim=[-0.015, 0.016])
plot!(0:tsave:tstop, y[4,:], label="Δq [p.u.]",)
plot!(0:tsave:tstop, u, label="u=Δq⋆")
savefig("src/results/figs/3Bus/open_loop_probing.pdf")
CSV.write("src/results/data/3Bus/open_loop_probing.csv", DataFrame(hcat(Vector(0:tsave:tstop), u, y[2:4,:]'), :auto), header = false)

est_sys = subspaceid(d, :2, verbose=true, r=100, s1=10, s2=10, W=:MOESP, zeroD=true, scaleU=true, stable=false)
rom = d2c(est_sys.sys, :zoh)
eigs_est=eigvals(rom.A)


scatter(real(full_eigs[28:end]), imag(full_eigs[28:end]), marker_z=sum_part[28:end],  markersize=8, color =:viridis,  xlabel=L"\operatorname{Re}(\lambda)", ylabel=L"\operatorname{Im} (\lambda)", label="FOM",  colorbar_title="Sum of CIG state participation")
scatter!(real(eigs_est), imag(eigs_est), label="ROM", markersize=5, markershape=:star5, color=:red, legend=:topright)
savefig("src/results/figs/3Bus/eig_plot.pdf")
CSV.write("src/results/data/3Bus/full_eigs.csv", DataFrame(hcat(real(full_eigs[28:end]), imag(full_eigs[28:end]), sum_part[28:end]), :auto), header = false)
CSV.write("src/results/data/3Bus/rom_eigs.csv", DataFrame(hcat(real(eigs_est), imag(eigs_est)), :auto), header = false)

max_order=6
svc=zeros(max_order)
for i = 1:max_order
    svc[i] = est_sys.s.S[i+1]^2 + (i*size(y)[1] + i*(size(y)[1]+1) + size(y)[1])*log(length(u))/length(u) 
end


scatter(1:max_order, svc, label="", xlabel="Reduced model order n", ylabel=L"SVC(n)", markersize=5)
plot!(1:max_order, svc, label="")
savefig("src/results/figs/3Bus/ROM_order.pdf")
CSV.write("src/results/data/3Bus/rom_order.csv", DataFrame(hcat(Vector(1:max_order), svc), :auto), header = false)

opt_gains, minVal=tune_pss(rom)

add_damping_controller!(sys, Gf)

sim = build_sim_pss(sys, Gf, opt_gains, tspan, l_change)

small_sig = small_signal_analysis(sim)
new_eigs = small_sig.eigenvalues

scatter(real(new_eigs), imag(new_eigs),  color=:lightgrey, xlim=[-20, 0], ylim=[-100, 100], xlabel=L"\operatorname{Re}(\lambda)", ylabel=L"\operatorname{Im} (\lambda)", label="", markersize=9, colorbar_title="Sum of CIG state participation")
scatter!(real(new_eigs[30:31]), imag(new_eigs[30:31]), color=:salmon, markersize=9, label="With propsoed controller")
scatter!(real(full_eigs[31:32]), imag(full_eigs[31:32]), color=:lightskyblue, markersize=9, label="Base case", legend=:bottomright)
savefig("src/results/figs/3Bus/eig_comp.pdf")

CSV.write("src/results/data/3Bus/pss_eigs.csv", DataFrame(hcat(real(new_eigs[30:end]), imag(new_eigs[30:end])), :auto), header = false)
CSV.write("src/results/data/3Bus/moved_eigs.csv", DataFrame(hcat(real(full_eigs[31:32]), imag(full_eigs[31:32]), real(new_eigs[30:31]), imag(new_eigs[30:31])), :auto), header = false)


execute!(sim, Sundials.IDA())
results = read_results(sim)
apower = PSID.get_activepower_series(results, Gf.name, dt=tsave)
rpower = PSID.get_reactivepower_series(results, Gf.name, dt=tsave)

pod_load=DataFrame(Time=apower[1], P=apower[2], Q=rpower[2])


plot(base_case_load.Time, base_case_load.P, xlim=[0, 2], label="Base case", xlabel="Time [s]", ylabel="Gf Active Power [p.u.]")
plot!(pod_load.Time, pod_load.P, label="With damping controller", legend=:bottomright)
savefig("src/results/figs/3Bus/GF_active.pdf")

CSV.write("src/results/data/3Bus/base_case.csv", base_case_load, header = false)
CSV.write("src/results/data/3Bus/pod.csv", pod_load, header = false)


plot(base_case_load.Time, base_case_load.Q, xlim=[0, 2], label="Base case", xlabel="Time [s]", ylabel="Gf Reactive Power [p.u.]")
plot!(pod_load.Time, pod_load.Q, label="With damping controller", legend=:bottomright)
savefig("src/results/figs/3Bus/GF_reactive.pdf")


Q0=get_Q_ref(Gf.dynamic_injector.outer_control.reactive_power_control)


n4sid_perturb =  ControlReferenceChange(t_control[c[1]], Gf.dynamic_injector, :Q_ref, Q0+u[c[1]])

perturb_vector=[n4sid_perturb]
for i in 2:length(c)
    push!(perturb_vector, ControlReferenceChange(t_control[c[i]], Gf.dynamic_injector, :Q_ref, Q0+u[c[i]]))
end

tspan = (0.0, 20.0)
sim = build_sim_pss(sys, Gf, opt_gains, tspan, perturb_vector)
execute!(sim, Sundials.IDA(), saveat=tsave)

res = read_results(sim)

y_pss=get_exp_data(Gf.name, 0.0, tsave, res)

plot(0:tsave:20, y_pss[2,:], label="Δp [p.u.]", xlabel="Time [s]", xlim=(0.0, 3.0), ylim=[-0.015, 0.016])
plot!(0:tsave:20, y_pss[4,:], label="Δq [p.u.]", xlim=(0.0, 3.0))

CSV.write("src/results/data/3Bus/closed_loop_probing.csv", DataFrame(hcat(Vector(0:tsave:tstop), u, y_pss[2:4,:]'), :auto), header = false)

# IEEE-39 model parameters

# Moved to testing files because used in LOAD_CHANGE_EVENTS_test
#const SYS_BASE_MVA = 100.0     #System base[MVA]
const F_BASE       = 60.0      #Nominal frequency[Hz]
const Ω_BASE_RAD_S = 2π * F_BASE     #Nominal angular frequency[rad/s]

const TF_SIM             = 10.0          #Simulation time[s]
const EXPORT_DT          = 0.01          #CSV/plot export step[s]
const RESAMPLE_FOR_CSV   = true      #Resample exported CSV data
const RESAMPLE_FOR_PLOTS = true       #Resample plotted data


const TG_TS  = 0.5    #Governor servo time constant[s]
const TG_EPS = 0.001   # other internal turbine stages are treated as very fast time (user can change)

const PSS_KS = 0.0    #PSS gain PSSFixed(V_pss = 0.0) made fixed 
#Solver settings

const SOLVER_ABSTOL    = 1e-5      #IDA absolute tolerance
const SOLVER_RELTOL    = 1e-5      #IDA relative tolerance
const SOLVER_MAX_ORDER = 5      #IDA maximum order
const SOLVER_INIT_DT   = 1e-3     #Initial solver step[s]
const SOLVER_DTMAX     = 1e-3   #Maximum solver step[s] #  IDA maximum step[s]   ( here  An adaptive timestepBDF DAE solver)


#SG/GFM/GFL bus groups


SG_BUSES  = SG_BUSES_test      #Synchronous-generator buses
GFM_BUSES = GFM_BUSES_test       #Grid-forming inverter buses
GFL_BUSES = GFL_BUSES_test         #Grid-following inverter buses
PLOT_BUSES_EXT = vcat(SG_BUSES, GFM_BUSES, GFL_BUSES)               #Buses shown in plots

const R_SG_PU  = 0.05              #SG speed droop[p.u.]
const R_GFM_PU = 0.08      #GFM speed droop[p.u.]


#Note to user: default parameters  maps are included  for all buses to make the script flexible when changing bus  assets assignments. 
#This prevents missing-parameter errors after changing bus assignments, and the defaults can be replaced with detailed values later.

#PLL and inverter droop defaults


const GFL_PLL_KP = GFL_PLL_KP_test    #PLL Kp
const GFL_PLL_KI = GFL_PLL_KI_test         #PLL Ki


#Important Note to User  regarding GFL  kω_droop::Float64 = GFL_KW_DEFAULT
# NOTE: I modified my local PSID source  to incoprate a FFR for my GFL outer power loop: eg 
# C:\Users\<my_username>\.julia\dev\PowerSimulationsDynamics\src\models\inverter_models\outer_control_models.jl
# The GFL ActivePowerPI/ReactivePowerPI mdl_outer_ode! block now includes FFR droop:
#   p_ref_eff = p_ref - Kω * (ω_pll - 1.0)
# and uses (p_ref_eff - p_oc) instead of the original fixed-P_ref error (p_ref - p_oc). User can go to the source file to make similkar chnages Else
# Users keeping the original PSID source file  should not pass Kω/ext=Dict("Kω"=>...);
# use the default ActivePowerPI constructor   



const GFL_KW_DEFAULT = 75.0         #GFL frequency-droop gain
const GFL_KW_map_ext = Dict{Int,Float64}(
    30 => GFL_KW_DEFAULT,              #Bus 30 GFL frequency-droop gain
    31 => GFL_KW_DEFAULT,                    #Bus 31 GFL frequency-droop gain
)

const GFM_TA_DEFAULT = 2.0 * 4.2#GFM virtual-inertia time constant[s]
const GFM_TA_map_ext = Dict{Int,Float64}(
    31 => 2 * 3.030,#Bus 31 virtual-inertia time constant[s]
    32 => 2 * 3.580,#Bus 32 virtual-inertia time constant[s]
    35 => 3.480,#Bus 35 virtual-inertia time constant[s]
    34 => 2 * 4.33,#Bus 34 virtual-inertia time constant[s]
)

const GFM_KW_DEFAULT = 12.5               #GFM frequency-droop gain
const GFM_KW_map_ext = Dict{Int,Float64}(
    31 => GFM_KW_DEFAULT,                # Bus 31 frequency-droop gain
    32 => GFM_KW_DEFAULT,         #Bus 32 frequency-droop gain
    35 => GFM_KW_DEFAULT,              #Bus 35 frequency-droop gain
    38 => GFM_KW_DEFAULT,             #Bus 38 frequency-droop gain
)

const GFM_KD_DEFAULT = 12.5             #GFM damping gain
const GFM_KD_map_ext = Dict{Int,Float64}(
    31 => GFM_KD_DEFAULT,            #Bus 31 damping gain
    32 => GFM_KD_DEFAULT,         #Bus 32 damping gain
    35 => GFM_KD_DEFAULT,          #Bus 35 damping gain
    38 => GFM_KD_DEFAULT,                    #Bus 38 damping gain
)


#Plot and export settings


const ABS_VOLTAGE_YLIM   = (0.94, 1.08)           #Voltage plot y-limits[p.u.]
const ABS_VOLTAGE_YTICKS = 0.94:0.02:1.08                #Voltage plot y-ticks[p.u.]

const BUS_PLOT_COLORS_EXT = Dict{Int,String}(
    30 => "#1f77b4",    #Bus 30 plot color
    31 => "#9467bd",      #Bus 31 plot color
    32 => "#ff7f0e",  # Bus 32 plot color
    33 => "#8c564b",   #Bus 33 plot color
    34 => "#d62728",   #Bus 34 plot color
    35 => "#e377c2",         #Bus 35 plot color
    36 => "#7f7f7f",      #Bus 36 plot color
    37 => "#bcbd22",      #Bus 37 plot color
    38 => "#17becf",   #Bus 38 plot color
    39 => "#2ca02c",     #Bus 39 plot color
)

bus_plot_color(bus::Int) = get(BUS_PLOT_COLORS_EXT, bus, "#333333")


# Network P/Q export settings


const SAVE_DIAGNOSTIC_NETWORK_PQ = true  # this bSave network P/Q diagnostic CSVs
const NETWORK_PQ_ALL_BUSES       = collect(1:39)         #Buses used in network P/Q reconstruction
const NETWORK_PQ_EXPORT_BUSES    = collect(30:39)   #Generator-side buses exported by default

const IEEE39_LOAD_BUSES_ALL = [1, 3, 4, 7, 8, 9, 12, 15, 16, 18, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 35]#IEEE-39 load buses
const IEEE39_EXTRA_LOAD_BUSES = [b for b in IEEE39_LOAD_BUSES_ALL if !(b in collect(30:39))]#Load buses outside 30..39
const EXPORT_EXTRA_LOAD_BUS_CSVS = true#Save load-bus voltage/PQ CSVs


# Inertia and damping


const H_ALL_map_ext = Dict{Int,Float64}(
    30 => 4.200,#Bus 30 inertia constant[s]
    31 => 3.030,#Bus 31 inertia constant[s]
    32 => 3.580,   #Bus 32 inertia constant[s]
    33 => 2.860,   #Bus 33 inertia constant[s]
    34 => 4.333,   #Bus 34 inertia constant[s]
    35 => 3.480,#Bus 35 inertia constant[s]
    36 => 3.480,#Bus 36 inertia constant[s]
    37 => 2.430,#Bus 37 inertia constant[s]
    38 => 3.450,#Bus 38 inertia constant[s]
    39 => 5.000,#Bus 39 inertia constant[s]
)

const D_SCALE      = 1.0    #Damping scale factor
const D_SG_DEFAULT = 40.0 / D_SCALE        #Default SG damping coefficient

const D_SG_map_ext = Dict{Int,Float64}(
    33 => 40.0 / D_SCALE,   #     Bus 33 damping coefficient
    36 => 40.0 / D_SCALE,    #Bus 36 damping coefficient
    37 => 40.0 / D_SCALE,  #Bus 37 damping coefficient
    38 => 40.0 / D_SCALE,   #Bus 38 damping coefficient
    39 => 48.0 / D_SCALE,     #Bus 39 damping coefficient
)


#bus helpers

is_sg_bus(bus::Int)  = bus in SG_BUSES
is_gfm_bus(bus::Int) = bus in GFM_BUSES
is_gfl_bus(bus::Int) = bus in GFL_BUSES

function model_label(bus::Int)
    if is_sg_bus(bus)
        return "SG"
    elseif is_gfm_bus(bus)
        return "GFM"
    elseif is_gfl_bus(bus)
        return "GFL"
    else
        return "Unknown"
    end
end

function getH(bus::Int)
    haskey(H_ALL_map_ext, bus) && return H_ALL_map_ext[bus]
    @warn "No inertia H found for bus; using H = 4.0" bus = bus model = model_label(bus)
    return 4.0
end

function getDsg(bus::Int)
    haskey(D_SG_map_ext, bus) && return D_SG_map_ext[bus]
    @warn "No SG damping D found; using D_SG_DEFAULT" bus = bus model = model_label(bus)
    return D_SG_DEFAULT
end

function getGfmTa(bus::Int)
    haskey(GFM_TA_map_ext, bus) && return GFM_TA_map_ext[bus]
    @warn "No GFM Ta found; using GFM_TA_DEFAULT" bus = bus model = model_label(bus)
    return GFM_TA_DEFAULT
end

function getGfmKw(bus::Int)
    haskey(GFM_KW_map_ext, bus) && return GFM_KW_map_ext[bus]
    @warn "No GFM kω found; using GFM_KW_DEFAULT" bus = bus model = model_label(bus)
    return GFM_KW_DEFAULT
end

function getGfmKd(bus::Int)
    haskey(GFM_KD_map_ext, bus) && return GFM_KD_map_ext[bus]
    @warn "No GFM kd found; using GFM_KD_DEFAULT" bus = bus model = model_label(bus)
    return GFM_KD_DEFAULT
end

function getGflKw(bus::Int)
    haskey(GFL_KW_map_ext, bus) && return GFL_KW_map_ext[bus]
    @warn "No GFL kω found; using GFL_KW_DEFAULT" bus = bus model = model_label(bus)
    return GFL_KW_DEFAULT
end


# Resource-model construction defaults


const DEVICE_W_REF = 1.0    #Device speed reference[p.u.]

const CONVERTER_RATED_VOLTAGE = 1.0    #Converter rated voltage[p.u.]
const CONVERTER_RATED_CURRENT = 1.0    #Converter rated current[p.u.]
const DC_SOURCE_VOLTAGE       = 1200.0       #DC link voltage[V]

const VIRTUAL_INERTIA_DEFAULTS = (
    Ta    = 0.1,    #Virtual-inertia time constant[s]
    kd    = 0.0,      #Virtual damping gain
    kω    = 0.0,       #Frequency-droop gain
    P_ref = 1.0,           #Active-power reference[p.u.]
)

const SG_MACHINE_SATURATION_SE = (0.0, 0.0)     #Machine saturation coefficients
const SG_AVR_VA_LIM            = (-10.0, 10.0)      #AVR output limits
const SG_AVR_AE                = 0.001       #Exciter saturation A
const SG_AVR_BE                = 1.0      #Exciter saturation B
const SG_TG_VALVE_POSITION_LIMITS = (min = -5.0, max = 5.0)        #Valve position limits
const SG_PSS_FIXED_VPSS        = 0.0#Fixed PSS output

const GFM_REACTIVE_POWER_DROOP_WF = 10.0        #GFM reactive-power droop filter

const GFM_DEFAULT_PARAMS = (
    kq   = 0.05,         #Reactive-power droop gain
    Ta   = GFM_TA_DEFAULT,#Virtual-inertia time constant[s]
    kd   = GFM_KD_DEFAULT,#Virtual damping gain
    kω   = GFM_KW_DEFAULT,#Frequency-droop gain
    kpv  = 1.0,#Voltage controller proportional gain
    kiv  = 1000.0,#Voltage controller integral gain
    kpc  = 3.0,   # Current controller proportional gain
    kic  = 600.0,#Current controller integral gain
    kffv = 0.0,    #Voltage feedforward gain
    kffi = 0.0,#Current feedforward gain
    rv   = 0.0,#Virtual resistance
    lv   = 0.0,#Virtual inductance
    ωad  = 2000.0,     #Active damping cutoff frequency[rad/s]
    kad  = 10.0,#Active damping gain
    lf   = 0.08,#LCL inverter-side inductance
    rf   = 0.05,#LCL inverter-side resistance
    cf   = 0.074,#LCL shunt capacitance
    lg   = 0.2,#LCL grid-side inductance
    rg   = 0.04,#LCL grid-side resistance
)

const GFL_KW_DROOP_CAP = 35.0#GFL droop gain cap

const GFL_DEFAULT_PARAMS = (
    kp_p     = 0.8,#Active-power PI proportional gain
    ki_p     = 8.0,#Active-power PI integral gain
    ωz_p     = 80.0,#Active-power PI zero shaping frequency
    kp_q     = 0.8,#Reactive-power PI proportional gain
    ki_q     = 12.0,#Reactive-power PI integral gain
    ωf_q     = 120.0,#Reactive loop filter cutoff frequency
    kω_droop = min(GFL_KW_DEFAULT, GFL_KW_DROOP_CAP),#Frequency-droop gain
    kpc      = 1.2,#Current controller proportional gain
    kic      = 120.0,#Current controller integral gain
    kffv     = 0.0,#Voltage feedforward gain
    pll_ωlp  = 80.0,#PLL low-pass filter bandwidth
    pll_kp   = GFL_PLL_KP,#PLL proportional gain
    pll_ki   = GFL_PLL_KI,#PLL integral gain
    lf       = 0.08,#LCL inverter-side inductance
    rf       = 0.05,#LCL inverter-side resistance
    cf       = 0.074,#LCL shunt capacitance
    lg       = 0.2,#LCL grid-side inductance
    rg       = 0.04,#LCL grid-side resistance
)

const GFL_MAIN_CALL_OVERRIDES = (
    kp_q = 2.0,#Reactive-power PI proportional gain
    ki_q = 30.0,#Reactive-power PI integral gain
    ωf_q = 500.0,#Reactive loop filter cutoff frequency
    kic  = 300.0,#Current controller integral gain
)


#Generator base and SG parameter tables


const GEN_MBASE_DEFAULT_MVA = 1000.0#Default generator base[MVA]
const GEN_MBASE_SPECIAL_MVA = Dict{Int,Float64}(
    39 => 10000.0,#Bus 39 generator base[MVA]
    34 => 600.0,#Bus 34 generator base[MVA]
)

function generator_mbase_map()
    out = Dict{Int,Float64}()
    for b in 30:39
        out[b] = GEN_MBASE_DEFAULT_MVA
    end
    for (b, mbase) in GEN_MBASE_SPECIAL_MVA
        out[b] = mbase
    end
    return out
end

function default_sg_machine_params()
    return Dict{Symbol,Float64}(
        :Rs     => 0.0025,#Stator resistance[p.u.]
        :Xl     => 0.20,#Leakage reactance[p.u.]
        :Xd     => 2.0,#d-axis synchronous reactance[p.u.]
        :Xq     => 1.9,#q-axis synchronous reactance[p.u.]
        :Xd_p   => 0.40,#d-axis transient reactance[p.u.]
        :Xq_p   => 0.60,#q-axis transient reactance[p.u.]
        :Xd_pp  => 0.30,#d-axis subtransient reactance[p.u.]
        :Tdo_p  => 7.0,#d-axis transient open-circuit time constant[s]
        :Tqo_p  => 0.70,#q-axis transient open-circuit time constant[s]
        :Tdo_pp => 0.035,#d-axis subtransient open-circuit time constant[s]
        :Tqo_pp => 0.070,#q-axis subtransient open-circuit time constant[s]
    )
end

function sg_machine_params_by_bus()
    return Dict{Int,Dict{Symbol,Float64}}(
        b => default_sg_machine_params() for b in 30:39
    )
end

function stable_avr_params()
    return Dict{Symbol,Float64}(
        :Ka => 15.0,#AVR gain
        :Ke => 1.0,#Exciter gain
        :Kf => 0.08,#Stabilizer gain
        :Ta => 0.03,#AVR time constant[s]
        :Te => 0.8,#Exciter time constant[s]
        :Tf => 1.0,#Stabilizer time constant[s]
        :Tr => 0.08,#Transducer time constant[s]
    )
end

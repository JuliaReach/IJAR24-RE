# ----------------------------------------------------
#  Case 1: small uncertainty Wpos = Wvel = 0.1
#  Case 2: intermediate uncertainty Wpos = Wvel = 0.4
# ----------------------------------------------------

function _initial_states_quadrotor(; Wpos=0.8, Wvel=0.8, NSteps=NSTEPS)
    X0c = zeros(12)

    # fix the discretization steps
    PBA.setSteps(1)
    x1 = U(X0c[1] - Wpos, X0c[1] + Wpos)
    x2 = U(X0c[2] - Wpos, X0c[2] + Wpos)

    x1n = normalize(x1)
    x2n = normalize(x2)

    PBA.setSteps(NSteps)
    x3 = U(X0c[3] - Wpos, X0c[3] + Wpos)
    x3n = normalize(x3)

    PBA.setSteps(1)
    x4 = U(X0c[4] - Wvel, X0c[4] + Wvel)
    x5 = U(X0c[5] - Wvel, X0c[5] + Wvel)
    x6 = U(X0c[6] - Wvel, X0c[6] + Wvel)

    x4n = normalize(x4)
    x5n = normalize(x5)
    x6n = normalize(x6)

    aux = U(0, 0.1)
    x7 = aux
    x8 = aux
    x9 = aux
    x10 = aux
    x11 = aux
    x12 = aux

    x7n = normalize(x7)
    x8n = normalize(x8)
    x9n = normalize(x9)
    x10n = normalize(x10)
    x11n = normalize(x11)
    x12n = normalize(x12)

    C = IndepCopula(12)
    U0 = C([x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12])
    U0n = C([x1n, x2n, x3n, x4n, x5n, x6n, x7n, x8n, x9n, x10n, x11n, x12n])

    return U0, U0n
end

function _initial_states_quadrotor_samesteps(; Wpos=0.8, Wvel=0.8, Nsteps=1000)
    X0c = zeros(12)

    # fix the discretization steps
    PBA.setSteps(Nsteps)

    x1 = U(X0c[1] - Wpos, X0c[1] + Wpos)
    x2 = U(X0c[2] - Wpos, X0c[2] + Wpos)

    x1n = normalize(x1)
    x2n = normalize(x2)

    x3 = U(X0c[3] - Wpos, X0c[3] + Wpos)
    x3n = normalize(x3)

    x4 = U(X0c[4] - Wvel, X0c[4] + Wvel)
    x5 = U(X0c[5] - Wvel, X0c[5] + Wvel)
    x6 = U(X0c[6] - Wvel, X0c[6] + Wvel)

    x4n = normalize(x4)
    x5n = normalize(x5)
    x6n = normalize(x6)

    aux = U(0, 0.1)
    x7 = aux
    x8 = aux
    x9 = aux
    x10 = aux
    x11 = aux
    x12 = aux

    x7n = normalize(x7)
    x8n = normalize(x8)
    x9n = normalize(x9)
    x10n = normalize(x10)
    x11n = normalize(x11)
    x12n = normalize(x12)

    C = IndepCopula(12)
    U0 = C([x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12])
    U0n = C([x1n, x2n, x3n, x4n, x5n, x6n, x7n, x8n, x9n, x10n, x11n, x12n])

    return U0, U0n
end

if SAMESTEPS
    U0_case1, U0n_case1 = _initial_states_quadrotor_samesteps(; Wpos=0.4, Wvel=0.4)
    U0_case2, U0n_case2 = _initial_states_quadrotor_samesteps(; Wpos=0.8, Wvel=0.8)
else
    U0_case1, U0n_case1 = _initial_states_quadrotor(; Wpos=0.4, Wvel=0.4)
    U0_case2, U0n_case2 = _initial_states_quadrotor(; Wpos=0.8, Wvel=0.8)
end
alg = TMJets20(; abstol=1e-7, orderT=5, orderQ=1, adaptive=false);
prob_case1 = @ivp(x' = quadrotor!(x), dim:12, x(0) ∈ U0_case1);
prob_case2 = @ivp(x' = quadrotor!(x), dim:12, x(0) ∈ U0_case2);

@time sol_case1 = solve(prob_case1; tspan=Tspan, alg=alg, ensemble=false, trajectories=10,
                        initial_pbox_normalized=U0n_case1);
solz_case1 = overapproximate(sol_case1, Zonotope);

@time sol_case2 = solve(prob_case2; tspan=Tspan, alg=alg, ensemble=false, trajectories=10,
                        initial_pbox_normalized=U0n_case2);
solz_case2 = overapproximate(sol_case2, Zonotope);

# ------------------------------------------------------------------------------
# Property 1: b1 = (x[3] < 1.4) for all times
unsafe1 = HalfSpace(-Vector(v3), -1.4) # unsafe: x3 >= 1.4
# check if the flowpipe is disjoint with the unsafe states
sat_cond1_case1 = all(is_intersection_empty(unsafe1, set(R)) for R in solz_case1) # true
@assert sat_cond1_case1
sat_cond1_case2 = all(is_intersection_empty(unsafe1, set(R)) for R in solz_case2) # false
@assert !sat_cond1_case2

pfail_cond1_case1 = failure_probability(sol_case1, 0 .. 5, unsafe1, U0n_case1; normalize=false,
                                        cons=false) # [0, 0]
pfail_cond1_case2 = failure_probability(sol_case2, 0 .. 5, unsafe1, U0n_case2; normalize=false,
                                        cons=false) # [0.0299999, 0.250001]

# ------------------------------------------------------------------------------
# Property 2: x[3] > 0.9 for t ≥ 1.0
unsafe2 = HalfSpace(v3, 0.9) # unsafe: x3 <= 0.9
sat_cond2_case1 = all(is_intersection_empty(unsafe2, set(R)) for R in solz_case1(1.0 .. 5.0)) # true
@assert sat_cond2_case1
sat_cond2_case2 = all(is_intersection_empty(unsafe2, set(R)) for R in solz_case2(1.0 .. 5.0)) # false
@assert !sat_cond2_case2

pfail_cond2_case1 = failure_probability(sol_case1, 1 .. 5, unsafe2, U0n_case1; normalize=false,
                                        cons=false) # [0, 0]
pfail_cond2_case2 = failure_probability(sol_case2, 1 .. 5, unsafe2, U0n_case2; normalize=false,
                                        cons=false) # [0.129999, 0.220001]

# ------------------------------------------------------------------------------
# Property 3: x[3] ⊆ Interval(0.98, 1.02) for t ≥ 5.0
sat_cond3_case1 = set(project(solz_case1[end]; vars=(3))) ⊆ LazySets.Interval(0.98, 1.02) # true
@assert sat_cond3_case1
sat_cond3_case2 = set(project(solz_case2[end]; vars=(3))) ⊆ LazySets.Interval(0.98, 1.02) # true
@assert sat_cond3_case2

unsafe3a = HalfSpace(v3, 0.98) # x[3] <= 0.98
unsafe3b = HalfSpace(-v3, -1.02) # x[3] >= 1.02
pfail_cond3_case1 = failure_probability(sol_case1, 5.0, unsafe3a ∪ unsafe3b, U0n_case1;
                                        normalize=false, cons=false) # [0, 0]
pfail_cond3_case2 = failure_probability(sol_case2, 5.0, unsafe3a ∪ unsafe3b, U0n_case2;
                                        normalize=false, cons=false) # [0, 0]

# ------------------------------------------------------------------------------
# Properties 1-3
ts = [0 .. 5, 1 .. 5, 5.0]
Us = [unsafe1, unsafe2, unsafe3a ∪ unsafe3b]
pfail_cond13_case1 = failure_probability(sol_case1, ts, Us, U0n_case1; normalize=false, cons=false) # [0, 0]
pfail_cond13_case2 = failure_probability(sol_case2, ts, Us, U0n_case2; normalize=false, cons=false) # [0, 0]

# ------------------------------------------------------------------------------
# Compute confidence
pconf = BENCHMARK_MODE ? [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001] :
        [1.0, 0.25, 0.5, 0.75, 0.001]
PBA.setSteps(NSTEPS) # reset step size
conf_case1 = confidence(sol_case1, pconf; cons=false, pbinit=U0n_case1, normalize=false)
conf_case2 = confidence(sol_case2, pconf; cons=false, pbinit=U0n_case2, normalize=false)

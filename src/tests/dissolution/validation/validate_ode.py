"""
Monolithic vs Operator-Split ODE Validation for FoxFlow Dissolution Model
=========================================================================

Solves the 1D pressure diffusion equation on a linear chain of 11 nodes:
    FluidSource -- Fluid -- Fluid -- ... -- Fluid -- FluidSink

Three solvers:
  1. Monolithic: Full system solved as one ODE  (ground truth)
  2. Operator-split (geocell): Domain decomposed into geocells, each solved
     independently per SSA step, mimicking DGGML's approximate SSA.
  3. Analytical steady-state: P_i = 1 - i/N

Compares time evolution and error.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os

# ═══════════════════════════════════════════════════════════════════════════════
# Model Parameters (from dissolution model)
# ═══════════════════════════════════════════════════════════════════════════════

N = 11              # total nodes: 1 source + 9 interior + 1 sink
ALPHA = 11.11       # diffusion coefficient = k / (mu * c * dx^2) ≈ 1/0.09
DX = 0.3            # spacing between nodes
P_SOURCE = 1.0      # fixed source pressure (node 0)
P_SINK = 0.0        # fixed sink pressure (node N-1)
DELTA = 0.1         # SSA macro time step (ODE integration window)
T_TOTAL = 10.0      # total simulation time

# Initial conditions: source=1, all others=0
P0 = np.zeros(N)
P0[0] = P_SOURCE


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Monolithic Solver (Ground Truth)
# ═══════════════════════════════════════════════════════════════════════════════

def rhs_monolithic(t, P):
    """Full system RHS: dP_i/dt = alpha * sum_neighbors(P_j - P_i)
    with fixed boundary conditions at node 0 (source) and node N-1 (sink)."""
    dPdt = np.zeros_like(P)
    for i in range(1, N - 1):  # interior nodes only
        # left neighbor
        dPdt[i] += ALPHA * (P[i - 1] - P[i])
        # right neighbor
        dPdt[i] += ALPHA * (P[i + 1] - P[i])
    # Boundary nodes: dP/dt = 0 (fixed)
    dPdt[0] = 0.0
    dPdt[N - 1] = 0.0
    return dPdt


print("=" * 70)
print("Solving MONOLITHIC system (ground truth)...")
print("=" * 70)

t_eval = np.linspace(0, T_TOTAL, 201)
sol_mono = solve_ivp(rhs_monolithic, [0, T_TOTAL], P0, t_eval=t_eval,
                     method='RK45', rtol=1e-10, atol=1e-12)
print(f"  Steps: {sol_mono.nfev} function evaluations")
print(f"  Success: {sol_mono.success}")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Operator-Split Solver (Mimics DGGML geocell decomposition)
# ═══════════════════════════════════════════════════════════════════════════════

def assign_geocells(N, n_cells):
    """Assign each node to a geocell based on spatial position.
    Simple 1D partitioning: divide the Z range into n_cells equal intervals."""
    z_positions = np.array([i * DX for i in range(N)])
    z_min, z_max = z_positions[0], z_positions[-1]
    cell_width = (z_max - z_min + 1e-12) / n_cells
    cells = np.minimum(((z_positions - z_min) / cell_width).astype(int), n_cells - 1)
    return cells, z_positions


def solve_operator_split(P_init, n_geocells, delta, t_total, sub_method='RK45'):
    """
    Operator-split solver mimicking DGGML's approximate SSA:

    For each macro step of size DELTA:
      For each geocell k:
        1. Gather all rule matches assigned to geocell k
        2. Build local ODE system from those matches
        3. Solve for duration DELTA
        4. Write results back to global state (copy_back)

    Rule assignment: An edge (i,j) match with anchor=i is assigned to
    the geocell containing node i. The reverse match (j,i) with anchor=j
    goes to node j's geocell.
    """
    geocells, z_pos = assign_geocells(N, n_geocells)
    n_steps = int(t_total / delta)

    # Storage for time history
    times = [0.0]
    states = [P_init.copy()]
    P = P_init.copy()

    for step in range(n_steps):
        t_start = step * delta
        t_end = (step + 1) * delta

        # Process each geocell independently
        for cell_k in range(n_geocells):
            # Find all edges assigned to this geocell.
            # Edge (i, i+1): anchor=i → assigned to geocell of node i
            # Edge (i+1, i): anchor=i+1 → assigned to geocell of node i+1
            local_edges = []  # list of (anchor, neighbor) pairs

            for i in range(N - 1):
                # Forward edge: anchor=i, neighbor=i+1
                if geocells[i] == cell_k:
                    local_edges.append((i, i + 1))
                # Reverse edge: anchor=i+1, neighbor=i
                if geocells[i + 1] == cell_k:
                    local_edges.append((i + 1, i))

            if not local_edges:
                continue

            # Determine which nodes are in the local ODE system
            local_nodes = sorted(set(n for edge in local_edges for n in edge))
            node_to_idx = {n: idx for idx, n in enumerate(local_nodes)}
            n_local = len(local_nodes)

            # Initial conditions from global state
            y0_local = np.array([P[n] for n in local_nodes])

            # Determine which are boundary nodes (source=0, sink=N-1)
            # AND which nodes are "owned" by this geocell (their anchor matches
            # appear here) vs "borrowed" (appear only as neighbors in other
            # geocell's matches)
            fixed = set()
            if 0 in local_nodes:
                fixed.add(0)
            if (N - 1) in local_nodes:
                fixed.add(N - 1)

            # Build local RHS from assigned edges
            def make_local_rhs(edges, nodes, n2i, fixed_nodes):
                def local_rhs(t, y_local):
                    dydt = np.zeros(len(nodes))
                    for anchor, neighbor in edges:
                        a_idx = n2i[anchor]
                        n_idx = n2i[neighbor]

                        # Rule: dP_anchor/dt += alpha * (P_neighbor - P_anchor)
                        # BUT boundary nodes are fixed
                        if anchor not in fixed_nodes:
                            dydt[a_idx] += ALPHA * (y_local[n_idx] - y_local[a_idx])

                    # Fixed boundary nodes
                    for fn in fixed_nodes:
                        if fn in n2i:
                            dydt[n2i[fn]] = 0.0
                    return dydt
                return local_rhs

            local_rhs_func = make_local_rhs(local_edges, local_nodes,
                                            node_to_idx, fixed)

            # Solve local ODE
            sol_local = solve_ivp(local_rhs_func, [t_start, t_end], y0_local,
                                  method=sub_method, rtol=1e-8, atol=1e-10)

            if not sol_local.success:
                print(f"  WARNING: geocell {cell_k} failed at step {step}")

            # Copy back: write final values to global state
            # (Only update nodes that are NOT fixed boundary)
            for n in local_nodes:
                if n not in fixed:
                    P[n] = sol_local.y[node_to_idx[n], -1]

        times.append(t_end)
        states.append(P.copy())

    return np.array(times), np.array(states)


# Solve with different numbers of geocells
print("\n" + "=" * 70)
print("Solving OPERATOR-SPLIT systems...")
print("=" * 70)

# 1 geocell = monolithic (no splitting error)
# 4 geocells matches the DGGML settings (CELL_NZ=4 for Z-direction)
# More geocells = more splitting, more error
geocell_configs = [1, 2, 4, 6, 11]
split_solutions = {}

for n_gc in geocell_configs:
    print(f"\n  Geocells: {n_gc}, DELTA: {DELTA}")
    times_split, states_split = solve_operator_split(P0, n_gc, DELTA, T_TOTAL)
    split_solutions[n_gc] = (times_split, states_split)

    # Report geocell assignment
    cells, _ = assign_geocells(N, n_gc)
    for k in range(n_gc):
        nodes_in_cell = [i for i in range(N) if cells[i] == k]
        print(f"    Cell {k}: nodes {nodes_in_cell}")


# Also solve with smaller DELTA for the 4-geocell case
delta_configs = [0.1, 0.05, 0.02, 0.01, 0.005]
delta_solutions = {}

print("\n" + "=" * 70)
print("Solving 4-GEOCELL with varying DELTA...")
print("=" * 70)

for dt in delta_configs:
    print(f"  DELTA: {dt}")
    times_dt, states_dt = solve_operator_split(P0, 4, dt, T_TOTAL)
    delta_solutions[dt] = (times_dt, states_dt)


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Analytical Steady State
# ═══════════════════════════════════════════════════════════════════════════════

P_steady = np.array([P_SOURCE * (1.0 - i / (N - 1)) for i in range(N)])


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Visualization
# ═══════════════════════════════════════════════════════════════════════════════

output_dir = os.path.dirname(os.path.abspath(__file__))

# --- Figure 1: Pressure profiles at snapshots ---
fig1 = plt.figure(figsize=(16, 12))
gs = GridSpec(2, 2, figure=fig1, hspace=0.3, wspace=0.3)

snapshots = [0.05, 0.2, 1.0, 10.0]
# Only plot interior nodes (exclude source=0 and sink=N-1)
interior = slice(1, N - 1)
x_nodes = np.arange(1, N - 1)

for idx, t_snap in enumerate(snapshots):
    ax = fig1.add_subplot(gs[idx // 2, idx % 2])

    # Monolithic (ground truth)
    i_mono = np.argmin(np.abs(sol_mono.t - t_snap))
    ax.plot(x_nodes, sol_mono.y[interior, i_mono], 'k-o', linewidth=2, markersize=6,
            label='Monolithic (truth)', zorder=10)

    # Operator-split with different geocells
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']
    for ci, n_gc in enumerate(geocell_configs):
        times_s, states_s = split_solutions[n_gc]
        i_split = np.argmin(np.abs(times_s - t_snap))
        ls = '--' if n_gc > 1 else '-'
        ax.plot(x_nodes, states_s[i_split, interior], ls, color=colors[ci], linewidth=1.5,
                marker='s', markersize=4, alpha=0.8,
                label=f'Split {n_gc} geocells')

    # Steady state
    ax.plot(x_nodes, P_steady[interior], ':', color='gray', linewidth=1, label='Steady state')

    ax.set_title(f't = {t_snap:.2f} s', fontsize=13, fontweight='bold')
    ax.set_xlabel('Interior node index')
    ax.set_ylabel('Pressure')
    ax.set_ylim(-0.05, 1.1)
    ax.grid(True, alpha=0.3)
    if idx == 0:
        ax.legend(fontsize=8, loc='upper right')

fig1.suptitle('Pressure Profiles: Monolithic vs Operator-Split\n'
              f'(α={ALPHA}, N={N} nodes, DELTA={DELTA}s)',
              fontsize=14, fontweight='bold')

fig1_path = os.path.join(output_dir, 'pressure_comparison.png')
fig1.savefig(fig1_path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {fig1_path}")


# --- Figure 2: Error vs time for different geocell counts ---
fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(14, 6))

# Left: max |error| vs time for different geocell counts
for ci, n_gc in enumerate(geocell_configs):
    if n_gc == 1:
        continue  # skip 1 geocell (no split)
    times_s, states_s = split_solutions[n_gc]
    errors = []
    for i, t_s in enumerate(times_s):
        i_mono = np.argmin(np.abs(sol_mono.t - t_s))
        err = np.max(np.abs(states_s[i, 1:-1] - sol_mono.y[1:-1, i_mono]))
        errors.append(err)
    ax2a.plot(times_s, errors, '-', linewidth=1.5,
              label=f'{n_gc} geocells', color=colors[geocell_configs.index(n_gc)])

ax2a.set_xlabel('Time (s)')
ax2a.set_ylabel('Max |P_split - P_mono|')
ax2a.set_title('Splitting Error vs Time\n(varying geocell count, DELTA=0.1s)')
ax2a.legend()
ax2a.grid(True, alpha=0.3)
ax2a.set_yscale('log')

# Right: max |error| vs time for different DELTA (4 geocells)
delta_colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(delta_configs)))
for di, dt in enumerate(delta_configs):
    times_dt, states_dt = delta_solutions[dt]
    errors = []
    for i, t_d in enumerate(times_dt):
        i_mono = np.argmin(np.abs(sol_mono.t - t_d))
        err = np.max(np.abs(states_dt[i, 1:-1] - sol_mono.y[1:-1, i_mono]))
        errors.append(err)
    ax2b.plot(times_dt, errors, '-', linewidth=1.5,
              label=f'Δ={dt}', color=delta_colors[di])

ax2b.set_xlabel('Time (s)')
ax2b.set_ylabel('Max |P_split - P_mono|')
ax2b.set_title('Splitting Error vs Time\n(4 geocells, varying DELTA)')
ax2b.legend()
ax2b.grid(True, alpha=0.3)
ax2b.set_yscale('log')

fig2.suptitle('Operator-Splitting Error Analysis', fontsize=14, fontweight='bold')
fig2.tight_layout()

fig2_path = os.path.join(output_dir, 'splitting_error.png')
fig2.savefig(fig2_path, dpi=150, bbox_inches='tight')
print(f"Saved: {fig2_path}")


# --- Figure 3: Space-time heatmap comparison ---
fig3, axes3 = plt.subplots(1, 3, figsize=(18, 6))

# Monolithic (interior nodes only)
im0 = axes3[0].imshow(sol_mono.y[1:-1, :], aspect='auto', origin='lower',
                       extent=[0, T_TOTAL, 1, N - 2], vmin=0, vmax=1,
                       cmap='hot')
axes3[0].set_title('Monolithic (Truth)')
axes3[0].set_xlabel('Time (s)')
axes3[0].set_ylabel('Interior node index')
plt.colorbar(im0, ax=axes3[0], label='Pressure')

# Split (4 geocells, DELTA=0.1) — interior nodes only
times_4, states_4 = split_solutions[4]
im1 = axes3[1].imshow(states_4[:, 1:-1].T, aspect='auto', origin='lower',
                       extent=[0, T_TOTAL, 1, N - 2], vmin=0, vmax=1,
                       cmap='hot')
axes3[1].set_title(f'Split (4 geocells, Δ={DELTA})')
axes3[1].set_xlabel('Time (s)')
axes3[1].set_ylabel('Interior node index')
plt.colorbar(im1, ax=axes3[1], label='Pressure')

# Error heatmap (interior nodes only)
from scipy.interpolate import interp1d
mono_interp = interp1d(sol_mono.t, sol_mono.y[1:-1, :], axis=1, kind='linear')
mono_at_split_times = mono_interp(np.clip(times_4, sol_mono.t[0], sol_mono.t[-1]))
error_matrix = np.abs(states_4[:, 1:-1].T - mono_at_split_times)

im2 = axes3[2].imshow(error_matrix, aspect='auto', origin='lower',
                       extent=[0, T_TOTAL, 1, N - 2],
                       cmap='RdYlBu_r')
axes3[2].set_title('|Error| (Split - Monolithic)')
axes3[2].set_xlabel('Time (s)')
axes3[2].set_ylabel('Interior node index')
plt.colorbar(im2, ax=axes3[2], label='|ΔP|')

fig3.suptitle('Space-Time Pressure Evolution', fontsize=14, fontweight='bold')
fig3.tight_layout()

fig3_path = os.path.join(output_dir, 'spacetime_comparison.png')
fig3.savefig(fig3_path, dpi=150, bbox_inches='tight')
print(f"Saved: {fig3_path}")


# --- Print summary table ---
print("\n" + "=" * 70)
print("ERROR SUMMARY")
print("=" * 70)

print(f"\n{'Config':<30} {'Max Error @ t=0.1':<18} {'Max Error @ t=1.0':<18} {'Max Error @ t=10':<18}")
print("-" * 84)

for n_gc in geocell_configs:
    times_s, states_s = split_solutions[n_gc]
    errs = []
    for t_check in [0.1, 1.0, 10.0]:
        i_s = np.argmin(np.abs(times_s - t_check))
        i_m = np.argmin(np.abs(sol_mono.t - t_check))
        errs.append(np.max(np.abs(states_s[i_s, 1:-1] - sol_mono.y[1:-1, i_m])))
    print(f"  {n_gc} geocells, Δ={DELTA:<13} {errs[0]:<18.6e} {errs[1]:<18.6e} {errs[2]:<18.6e}")

print()
for dt in delta_configs:
    times_dt, states_dt = delta_solutions[dt]
    errs = []
    for t_check in [0.1, 1.0, 10.0]:
        i_s = np.argmin(np.abs(times_dt - t_check))
        i_m = np.argmin(np.abs(sol_mono.t - t_check))
        errs.append(np.max(np.abs(states_dt[i_s, 1:-1] - sol_mono.y[1:-1, i_m])))
    print(f"  4 geocells, Δ={dt:<13} {errs[0]:<18.6e} {errs[1]:<18.6e} {errs[2]:<18.6e}")

# Steady-state error (interior only)
i_m_final = np.argmin(np.abs(sol_mono.t - T_TOTAL))
ss_err = np.max(np.abs(sol_mono.y[1:-1, i_m_final] - P_steady[1:-1]))
print(f"\n  Monolithic vs analytical steady state at t={T_TOTAL}: {ss_err:.6e}")

print("\nDone!")

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE / "files" / "rzut_ukosny"
OUT.mkdir(parents=True, exist_ok=True)

G = 9.81
V0 = 33.0
ANGLE_DEG = 40.0
MASS = 0.145
RHO = 1.225
CD = 0.47
DIAMETER = 0.070
AREA = np.pi * (DIAMETER / 2) ** 2
DRAG_COEFFICIENT = 0.5 * RHO * CD * AREA
DT = 0.002

BLUE = "#2463A6"
ORANGE = "#D47722"
GREEN = "#168578"
INK = "#233044"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelcolor": INK,
    "text.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.2,
    "lines.linewidth": 2.2,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def derivative(state):
    x, y, vx, vy = state
    speed = np.hypot(vx, vy)
    factor = DRAG_COEFFICIENT * speed / MASS
    return np.array([vx, vy, -factor * vx, -G - factor * vy])


def rk4_step(state, dt):
    k1 = derivative(state)
    k2 = derivative(state + 0.5 * dt * k1)
    k3 = derivative(state + 0.5 * dt * k2)
    k4 = derivative(state + dt * k3)
    return state + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)


def simulate(angle_deg, dt=DT):
    angle = np.deg2rad(angle_deg)
    state = np.array([0.0, 0.0, V0 * np.cos(angle), V0 * np.sin(angle)])
    rows = [state.copy()]
    times = [0.0]

    while state[1] >= 0.0 and len(rows) < 20000:
        next_state = rk4_step(state, dt)

        if next_state[1] < 0.0:
            fraction = state[1] / (state[1] - next_state[1])
            next_state = state + fraction * (next_state - state)
            next_state[1] = 0.0
            rows.append(next_state)
            times.append(times[-1] + fraction * dt)
            break

        rows.append(next_state)
        times.append(times[-1] + dt)
        state = next_state

    return np.asarray(rows), np.asarray(times)


def save(fig, filename):
    fig.savefig(OUT / filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    angle = np.deg2rad(ANGLE_DEG)
    vx0 = V0 * np.cos(angle)
    vy0 = V0 * np.sin(angle)
    flight_time_ideal = 2 * vy0 / G
    range_ideal = vx0 * flight_time_ideal
    max_height_ideal = vy0**2 / (2 * G)

    states, time = simulate(ANGLE_DEG)
    x, y, vx, vy = states.T
    speed = np.hypot(vx, vy)
    energy_kinetic = 0.5 * MASS * speed**2
    energy_potential = MASS * G * y

    print(f"v0 = {V0:.2f} m/s; alpha = {ANGLE_DEG:.1f} deg")
    print(f"Ideal: flight time = {flight_time_ideal:.3f} s, range = {range_ideal:.2f} m, max height = {max_height_ideal:.2f} m")
    print(f"With drag: flight time = {time[-1]:.3f} s, range = {x[-1]:.2f} m, max height = {y.max():.2f} m")
    print(f"Drag coefficient c = {DRAG_COEFFICIENT:.6f} kg/m; figures: {OUT}")

    ideal_time = np.linspace(0, flight_time_ideal, 600)
    ideal_x = vx0 * ideal_time
    ideal_y = vy0 * ideal_time - 0.5 * G * ideal_time**2

    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    ax.plot(ideal_x, ideal_y, color=BLUE, label="bez oporu")
    ax.plot(x, y, color=ORANGE, label="opór kwadratowy")
    ax.scatter([0], [0], color=INK, zorder=4, label="punkt wyrzutu")
    ax.scatter([range_ideal, x[-1]], [0, 0], color=[BLUE, ORANGE], zorder=4)
    ax.set(title="Tor ruchu: model idealny i opór powietrza",
           xlabel="x [m]", ylabel="y [m]")
    ax.legend(frameon=False)
    save(fig, "tory_ruchu.png")

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout="constrained")
    axes[0].plot(time, x, color=BLUE)
    axes[0].set(ylabel="x [m]", title="Współrzędne w funkcji czasu (opór kwadratowy)")
    axes[1].plot(time, y, color=GREEN)
    axes[1].set(xlabel="t [s]", ylabel="y [m]")
    save(fig, "wspolrzedne_czas.png")

    fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
    ax.plot(time, vx, color=BLUE, label="$v_x$")
    ax.plot(time, vy, color=ORANGE, label="$v_y$")
    ax.axhline(0, color=INK, lw=1)
    ax.set(title="Składowe prędkości w czasie (opór kwadratowy)",
           xlabel="t [s]", ylabel="prędkość [m/s]")
    ax.legend(frameon=False)
    save(fig, "predkosci_czas.png")

    fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
    ax.plot(time, energy_kinetic, color=BLUE, label="$E_k$")
    ax.plot(time, energy_potential, color=GREEN, label="$E_p$")
    ax.plot(time, energy_kinetic + energy_potential, color=ORANGE,
            label="$E_m=E_k+E_p$")
    ax.set(title="Energia mechaniczna maleje wskutek oporu",
           xlabel="t [s]", ylabel="energia [J]")
    ax.legend(frameon=False)
    save(fig, "energia.png")

    angles = np.linspace(10, 80, 36)
    ranges_drag = np.array([simulate(value)[0][-1, 0] for value in angles])
    ranges_ideal = V0**2 * np.sin(2 * np.deg2rad(angles)) / G
    fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
    ax.plot(angles, ranges_ideal, color=BLUE, label="bez oporu")
    ax.plot(angles, ranges_drag, color=ORANGE, label="opór kwadratowy")
    ax.axvline(ANGLE_DEG, color=INK, ls="--", lw=1.2,
               label=fr"kąt w raporcie: {ANGLE_DEG:.0f}°")
    ax.set(title="Zasięg a kąt wyrzutu przy stałej prędkości początkowej",
           xlabel="kąt wyrzutu [°]", ylabel="zasięg [m]")
    ax.legend(frameon=False)
    save(fig, "zasieg_kat.png")


if __name__ == "__main__":
    main()
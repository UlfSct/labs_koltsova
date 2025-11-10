import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.signal import find_peaks

def get_neighbors(i, j, size):
    neighbors = []
    for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        ni, nj = (i + di) % size, (j + dj) % size
        neighbors.append((ni, nj))
    return neighbors

def get_all_neighbor_pairs(size):
    pairs = set()
    for i in range(size):
        for j in range(size):
            neighbors = get_neighbors(i, j, size)
            for ni, nj in neighbors:
                pair = tuple(sorted([(i, j), (ni, nj)]))
                pairs.add(pair)
    return list(pairs)

def calculate_rates(lattice, size, P_CO):
    k01, k_1, k02, k3, k4, k5, k6, k7, k8 = 1e4, 3e2, 2.5e3, 2.5e4, 0.11, 6.5e-3, 1e5, 1e5, 1e5
    k1 = P_CO * k01
    events, rates = [], []

    for i in range(size):
        for j in range(size):
            if lattice[i, j] == 0:
                events.append(('adsorb_CO', i, j))
                rates.append(k1)
            if lattice[i, j] == 1:
                events.append(('desorb_CO', i, j))
                rates.append(k_1)
            if lattice[i, j] == 2:
                events.append(('create_vacancy', i, j))
                rates.append(k4)

    neighbor_pairs = get_all_neighbor_pairs(size)
    for (i1, j1), (i2, j2) in neighbor_pairs:
        state1, state2 = lattice[i1, j1], lattice[i2, j2]
        if state1 == 0 and state2 == 0:
            events.append(('adsorb_O2', i1, j1, i2, j2))
            rates.append(k02 / 2)
        if (state1 == 1 and state2 == 2) or (state1 == 2 and state2 == 1):
            events.append(('react_CO_O', i1, j1, i2, j2))
            rates.append(k3 / 2)
        if (state1 == 1 and state2 == 3) or (state1 == 3 and state2 == 1):
            events.append(('react_CO_OV', i1, j1, i2, j2))
            rates.append(k5 / 2)
        if state1 == 1 and state2 == 0:
            events.append(('migrate_CO', i1, j1, i2, j2))
            rates.append(k6)
        elif state1 == 0 and state2 == 1:
            events.append(('migrate_CO', i2, j2, i1, j1))
            rates.append(k6)
        if state1 == 2 and state2 == 0:
            events.append(('migrate_O', i1, j1, i2, j2))
            rates.append(k7)
        elif state1 == 0 and state2 == 2:
            events.append(('migrate_O', i2, j2, i1, j1))
            rates.append(k7)
        if state1 == 3 and state2 == 0:
            events.append(('migrate_OV', i1, j1, i2, j2))
            rates.append(k8)
        elif state1 == 0 and state2 == 3:
            events.append(('migrate_OV', i2, j2, i1, j1))
            rates.append(k8)
    return events, rates

def execute_event(lattice, event):
    event_type = event[0]
    if event_type == 'adsorb_CO':
        _, i, j = event
        lattice[i, j] = 1
    elif event_type == 'desorb_CO':
        _, i, j = event
        lattice[i, j] = 0
    elif event_type == 'adsorb_O2':
        _, i1, j1, i2, j2 = event
        lattice[i1, j1] = 2
        lattice[i2, j2] = 2
    elif event_type == 'create_vacancy':
        _, i, j = event
        lattice[i, j] = 3
    elif event_type == 'react_CO_O':
        _, i1, j1, i2, j2 = event
        lattice[i1, j1] = 0
        lattice[i2, j2] = 0
    elif event_type == 'react_CO_OV':
        _, i1, j1, i2, j2 = event
        lattice[i1, j1] = 0
        lattice[i2, j2] = 0
    elif event_type == 'migrate_CO':
        _, from_i, from_j, to_i, to_j = event
        lattice[to_i, to_j] = 1
        lattice[from_i, from_j] = 0
    elif event_type == 'migrate_O':
        _, from_i, from_j, to_i, to_j = event
        lattice[to_i, to_j] = 2
        lattice[from_i, from_j] = 0
    elif event_type == 'migrate_OV':
        _, from_i, from_j, to_i, to_j = event
        lattice[to_i, to_j] = 3
        lattice[from_i, from_j] = 0

def get_coverage(lattice):
    return np.sum(lattice == 1) / lattice.size

def simulate_co_oxidation(size=20, initial_state='empty', P_CO=0.3, total_time=10.0, save_interval=0.1):
    lattice = np.zeros((size, size), dtype=int)
    if initial_state == 'oxygen':
        total_sites = size * size
        oxygen_sites = int(0.4 * total_sites)
        vacancy_sites = int(0.4 * total_sites)
        positions = random.sample([(i, j) for i in range(size) for j in range(size)], oxygen_sites)
        for i, j in positions:
            lattice[i, j] = 2
        remaining_positions = [(i, j) for i in range(size) for j in range(size) if lattice[i, j] == 0]
        vacancy_positions = random.sample(remaining_positions, vacancy_sites)
        for i, j in vacancy_positions:
            lattice[i, j] = 3

    time_val, step_count = 0.0, 0
    time_history, coverage_history = [0.0], [get_coverage(lattice)]
    next_save_time = save_interval

    while time_val < total_time:
        print(time_val)
        events, rates = calculate_rates(lattice, size, P_CO)
        if not events:
            break

        total_rate = sum(rates)
        if total_rate == 0:
            break

        r = random.random() * total_rate
        cumulative_rate, selected_event = 0, None
        for event, rate in zip(events, rates):
            cumulative_rate += rate
            if cumulative_rate >= r:
                selected_event = event
                break

        if selected_event:
            execute_event(lattice, selected_event)
            dt = -np.log(random.random()) / total_rate
            time_val += dt
            step_count += 1

            if time_val >= next_save_time:
                coverage_history.append(get_coverage(lattice))
                time_history.append(time_val)
                next_save_time += save_interval

    coverage_history.append(get_coverage(lattice))
    time_history.append(time_val)
    return time_history, coverage_history, lattice

def analyze_oscillations(time_history, coverage_history):
    if len(coverage_history) < 10:
        return None, None

    coverage = np.array(coverage_history)
    time_arr = np.array(time_history)

    window_size = min(5, len(coverage) // 10)
    if window_size > 1:
        coverage_smooth = np.convolve(coverage, np.ones(window_size)/window_size, mode='valid')
        time_smooth = time_arr[window_size-1:]
    else:
        coverage_smooth = coverage
        time_smooth = time_arr

    peaks, _ = find_peaks(coverage_smooth, height=0.1, distance=5)
    troughs, _ = find_peaks(-coverage_smooth, height=-0.9, distance=5)

    if len(peaks) > 1 and len(troughs) > 1:
        peak_times = time_smooth[peaks]
        periods = np.diff(peak_times)
        avg_period = np.mean(periods)
        amplitudes = coverage_smooth[peaks] - coverage_smooth[troughs[:len(peaks)]]
        avg_amplitude = np.mean(amplitudes)
        return avg_period, avg_amplitude

    return None, None

def plot_results(time_history, coverage_history, size, P_CO):
    plt.figure(figsize=(12, 6))
    plt.plot(time_history, coverage_history, 'b-', linewidth=2)
    plt.xlabel('Время, с')
    plt.ylabel('Доля адсорбированного CO')
    plt.title(f'Динамика адсорбции CO (решетка {size}x{size}, P_CO = {P_CO} Торр)')
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()

def study_size_effect():
    sizes = [10, 20, 40]
    results = {}
    plt.figure(figsize=(12, 8))

    for size in sizes:
        time_history, coverage_history, _ = simulate_co_oxidation(size=size, total_time=10.0)
        results[size] = {'time': time_history, 'coverage': coverage_history}
        plt.plot(time_history, coverage_history, label=f'{size}x{size}', linewidth=2)
        period, amplitude = analyze_oscillations(time_history, coverage_history)
        if period:
            print(f"Размер {size}x{size}: период {period:.3f} с, амплитуда {amplitude:.3f}")

    plt.xlabel('Время, с')
    plt.ylabel('Доля адсорбированного CO')
    plt.title('Влияние размера решетки на динамику системы')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1)
    plt.show()
    return results

def study_pressure_effect():
    pressures = [0.05, 0.1, 0.2, 0.3, 0.4]
    results = {}
    plt.figure(figsize=(12, 8))

    for P_CO in pressures:
        time_history, coverage_history, _ = simulate_co_oxidation(P_CO=P_CO, total_time=10.0)
        results[P_CO] = {'time': time_history, 'coverage': coverage_history}
        plt.plot(time_history, coverage_history, label=f'P_CO = {P_CO} Торр', linewidth=2)
        period, amplitude = analyze_oscillations(time_history, coverage_history)
        if period:
            print(f"Давление {P_CO} Торр: период {period:.3f} с, амплитуда {amplitude:.3f}")

    plt.xlabel('Время, с')
    plt.ylabel('Доля адсорбированного CO')
    plt.title('Влияние давления CO на динамику системы')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1)
    plt.show()
    return results

if __name__ == "__main__":
    time_history, coverage_history, _ = simulate_co_oxidation()
    plot_results(time_history, coverage_history, 20, 0.3)
    analyze_oscillations(time_history, coverage_history)
    study_size_effect()
    study_pressure_effect()
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import random


def get_neighbors(i, j, X, Y):
    neighbors = []
    for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        ni, nj = (i + di) % X, (j + dj) % Y
        neighbors.append((ni, nj))
    return neighbors


def calculate_rates(lattice, X, Y, W, Q):
    events = []
    rates = []
    total_rate = 0.0
    for i in range(X):
        for j in range(Y):
            if lattice[i, j] == 1:
                neighbors = get_neighbors(i, j, X, Y)
                for ni, nj in neighbors:
                    if lattice[ni, nj] == 0:
                        events.append(('move', i, j, ni, nj))
                        rates.append(W)
                        total_rate += W
                    elif lattice[ni, nj] == 1:
                        events.append(('react', i, j, ni, nj))
                        rates.append(Q)
                        total_rate += Q
    return events, rates, total_rate


def execute_event(lattice, event):
    event_type, i, j, ni, nj = event
    if event_type == 'move':
        lattice[ni, nj] = 1
        lattice[i, j] = 0
    elif event_type == 'react':
        lattice[i, j] = 2
        lattice[ni, nj] = 0


def simulate_step(lattice, time, step_count, X, Y, W, Q):
    events, rates, total_rate = calculate_rates(lattice, X, Y, W, Q)
    r = random.random() * total_rate
    sum_rate = 0
    selected_event = None
    for event, rate in zip(events, rates):
        sum_rate += rate
        if sum_rate >= r:
            selected_event = event
            break
    execute_event(lattice, selected_event)
    dt = -np.log(random.random()) / total_rate
    time += dt
    step_count += 1
    return time, step_count


def save_state(history, time, step_count, lattice):
    state = {
        'time': time,
        'step': step_count,
        'lattice': lattice.copy(),
        'count_A': np.sum(lattice == 1),
        'count_B': np.sum(lattice == 2),
        'count_empty': np.sum(lattice == 0)
    }
    history.append(state)


def visualize(history):
    for el in history:
        fig, ax1 = plt.subplots(1, 1, figsize=(15, 6))
        cmap = mcolors.ListedColormap(['white', 'blue', 'red'])
        bounds = [0, 1, 2, 3]
        norm = mcolors.BoundaryNorm(bounds, cmap.N)
        im = ax1.imshow(el['lattice'], cmap=cmap, norm=norm,
                        interpolation='nearest')
        ax1.set_title(f"t = {el['time']:.4f}; N = {el['step']:.1f}")
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        cbar = plt.colorbar(im, ax=ax1, ticks=[0.5, 1.5, 2.5])
        cbar.ax.set_yticklabels(['Пусто', 'A', 'B'])
        plt.tight_layout()
        plt.show()

    fig, ax2 = plt.subplots(1, 1, figsize=(12, 5))
    times = [state['time'] for state in history]
    counts_A = [state['count_A'] for state in history]
    counts_B = [state['count_B'] for state in history]
    ax2.plot(times, counts_A, 'bo-', label='Частицы A', markersize=4)
    ax2.plot(times, counts_B, 'ro-', label='Частицы B', markersize=4)
    ax2.set_xlabel('Время')
    ax2.set_ylabel('Количество частиц')
    ax2.set_title('Динамика реакции A + A → B')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    for i in range(1, len(history)):
        if i % (len(history) // 10) == 0 or i == len(history) - 1:
            ax2.axvline(x=times[i], color='gray', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()


def main():
    X, Y = 8, 12
    N_A = 7
    W = 100
    Q = 0.1
    T_total = 100.0

    # Инициализация решетки
    lattice = np.zeros((X, Y), dtype=int)
    positions = random.sample([(i, j) for i in range(X) for j in range(Y)], N_A)
    for i, j in positions:
        lattice[i, j] = 1

    time = 0.0
    step_count = 0
    history = []

    # Симуляция
    output_interval = T_total / 10
    next_output_time = output_interval
    print(f"[N = 0]: t = {time:.4f}")
    save_state(history, time, step_count, lattice)

    while time < T_total:
        time, step_count = simulate_step(lattice, time, step_count, X, Y, W, Q)
        if time >= next_output_time:
            print(f"[N = {step_count}]: t = {time:.4f}")
            save_state(history, time, step_count, lattice)
            next_output_time += output_interval

    visualize(history)


if __name__ == "__main__":
    main()
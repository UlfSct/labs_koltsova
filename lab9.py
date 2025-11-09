import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import random


class LatticeGasReaction:
    def __init__(self, X, Y, N_A, W, Q, T):
        self.X = X
        self.Y = Y
        self.N_A = N_A
        self.W = W
        self.Q = Q
        self.T_total = T
        self.lattice = np.zeros((X, Y), dtype=int)

        positions = random.sample([(i, j) for i in range(X) for j in range(Y)], N_A)
        for i, j in positions:
            self.lattice[i, j] = 1

        self.time = 0.0
        self.step_count = 0
        self.history = []

    def get_neighbors(self, i, j):
        neighbors = []
        for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            ni, nj = (i + di) % self.X, (j + dj) % self.Y
            neighbors.append((ni, nj))
        return neighbors

    def calculate_rates(self):
        events = []
        rates = []
        total_rate = 0.0
        for i in range(self.X):
            for j in range(self.Y):
                if self.lattice[i, j] == 1:
                    neighbors = self.get_neighbors(i, j)
                    for ni, nj in neighbors:
                        if self.lattice[ni, nj] == 0:
                            events.append(('move', i, j, ni, nj))
                            rates.append(self.W)
                            total_rate += self.W
                        elif self.lattice[ni, nj] == 1:
                            events.append(('react', i, j, ni, nj))
                            rates.append(self.Q)
                            total_rate += self.Q
        return events, rates, total_rate

    def execute_event(self, event):
        event_type, i, j, ni, nj = event
        if event_type == 'move':
            self.lattice[ni, nj] = 1
            self.lattice[i, j] = 0
        elif event_type == 'react':
            self.lattice[i, j] = 2
            self.lattice[ni, nj] = 0

    def step(self):
        events, rates, total_rate = self.calculate_rates()
        r = random.random() * total_rate
        sum_rate = 0
        selected_event = None
        for event, rate in zip(events, rates):
            sum_rate += rate
            if sum_rate >= r:
                selected_event = event
                break
        self.execute_event(selected_event)
        dt = -np.log(random.random()) / total_rate
        self.time += dt
        self.step_count += 1
        return True

    def simulate(self):
        output_interval = self.T_total / 10
        next_output_time = output_interval
        print(f"[N = 0]: t = {self.time:.4f}")
        self.save_state()
        while self.time < self.T_total:
            self.step()
            if self.time >= next_output_time:
                print(f"[N = {self.step_count}]: t = {self.time:.4f}")
                self.save_state()
                next_output_time += output_interval

    def save_state(self):
        state = {
            'time': self.time,
            'step': self.step_count,
            'lattice': self.lattice.copy(),
            'count_A': np.sum(self.lattice == 1),
            'count_B': np.sum(self.lattice == 2),
            'count_empty': np.sum(self.lattice == 0)
        }
        self.history.append(state)

    def visualize(self):
        for el in self.history:
            fig, (ax1) = plt.subplots(1, 1, figsize=(15, 6))
            cmap = mcolors.ListedColormap(['white', 'blue', 'red'])
            bounds = [0, 1, 2, 3]
            norm = mcolors.BoundaryNorm(bounds, cmap.N)
            im = ax1.imshow(el['lattice'], cmap=cmap, norm=norm,
                            interpolation='nearest')
            ax1.set_title(f"t = {el['time']:.4f}")
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            cbar = plt.colorbar(im, ax=ax1, ticks=[0.5, 1.5, 2.5])
            cbar.ax.set_yticklabels(['Пусто', 'A', 'B'])
            plt.tight_layout()
            plt.show()

        fig, (ax2) = plt.subplots(1, 1, figsize=(12, 5))
        times = [state['time'] for state in self.history]
        counts_A = [state['count_A'] for state in self.history]
        counts_B = [state['count_B'] for state in self.history]
        ax2.plot(times, counts_A, 'bo-', label='Частицы A', markersize=4)
        ax2.plot(times, counts_B, 'ro-', label='Частицы B', markersize=4)
        ax2.set_xlabel('Время')
        ax2.set_ylabel('Количество частиц')
        ax2.set_title('Динамика реакции A + A → B')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        for i in range(1, len(self.history)):
            if i % (len(self.history) // 10) == 0 or i == len(self.history) - 1:
                ax2.axvline(x=times[i], color='gray', linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    X, Y = 8, 12
    N_A = 7
    W = 100
    Q = 0.1
    T = 100.0

    system = LatticeGasReaction(X, Y, N_A, W, Q, T)
    system.simulate()
    system.visualize()

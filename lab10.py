import numpy as np
import matplotlib.pyplot as plt
import random
from collections import defaultdict
import time

class COOxidationCatalyst:
    def __init__(self, size=20, initial_state='empty', P_CO=0.3):
        """
        Модель каталитического окисления СО

        Parameters:
        size - размер решетки (size x size)
        initial_state - начальное состояние: 'empty' или 'oxygen'
        P_CO - давление CO в Торр
        """
        self.size = size
        self.P_CO = P_CO

        # Константы скорости из условия
        self.k01 = 1e4  # c⁻¹ Торр⁻¹
        self.k_1 = 3e2  # c⁻¹
        self.k02 = 2.5e3  # c⁻¹ (уже с учетом PO2)
        self.k3 = 2.5e4  # c⁻¹
        self.k4 = 0.11  # c⁻¹
        self.k5 = 6.5e-3  # c⁻¹
        self.k6 = self.k7 = self.k8 = 1e5  # c⁻¹

        # Рассчитываем k1 с учетом давления CO
        self.k1 = self.P_CO * self.k01

        # 0 - (*) 1 - [CO] 2 - [O] 3 - [O]v
        self.lattice = np.zeros((size, size), dtype=int)

        # Инициализация начального состояния
        if initial_state == 'oxygen':
            total_sites = size * size
            oxygen_sites = int(0.4 * total_sites)
            vacancy_sites = int(0.4 * total_sites)

            positions = random.sample([(i, j) for i in range(size) for j in range(size)], oxygen_sites)
            for i, j in positions:
                self.lattice[i, j] = 2

            remaining_positions = [(i, j) for i in range(size) for j in range(size) if self.lattice[i, j] == 0]
            vacancy_positions = random.sample(remaining_positions, vacancy_sites)
            for i, j in vacancy_positions:
                self.lattice[i, j] = 3

        self.time = 0.0
        self.step_count = 0
        self.history = []

        # Статистика
        self.coverage_history = []
        self.time_history = []

    def get_neighbors(self, i, j):
        """Получить соседей с периодическими граничными условиями"""
        neighbors = []
        for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            ni, nj = (i + di) % self.size, (j + dj) % self.size
            neighbors.append((ni, nj))
        return neighbors

    def get_all_neighbor_pairs(self):
        """Получить все уникальные пары соседей (без дубликатов)"""
        pairs = set()
        for i in range(self.size):
            for j in range(self.size):
                neighbors = self.get_neighbors(i, j)
                for ni, nj in neighbors:
                    # Создаем уникальный идентификатор пары
                    pair = tuple(sorted([(i, j), (ni, nj)]))
                    pairs.add(pair)
        return list(pairs)

    def calculate_rates(self):
        """Рассчитать все возможные события и их скорости"""
        events = []
        rates = []

        for i in range(self.size):
            for j in range(self.size):
                # Адсорбция CO на свободный центр
                if self.lattice[i, j] == 0:
                    events.append(('adsorb_CO', i, j))
                    rates.append(self.k1)
                # Десорбция CO
                if self.lattice[i, j] == 1:  # [CO]
                    events.append(('desorb_CO', i, j))
                    rates.append(self.k_1)
                # Образование кислородных вакансий
                if self.lattice[i, j] == 2:  # [O]
                    events.append(('create_vacancy', i, j))
                    rates.append(self.k4)

        neighbor_pairs = self.get_all_neighbor_pairs()
        for (i1, j1), (i2, j2) in neighbor_pairs:
            # Адсорбция O2 на два соседних свободных центра
            state1 = self.lattice[i1, j1]
            state2 = self.lattice[i2, j2]
            if state1 == 0 and state2 == 0:
                events.append(('adsorb_O2', i1, j1, i2, j2))
                rates.append(self.k02 / 2)
            # Реакция [CO] + [O] -> CO2 + 2*
            if (state1 == 1 and state2 == 2) or (state1 == 2 and state2 == 1):
                events.append(('react_CO_O', i1, j1, i2, j2))
                rates.append(self.k3 / 2)

            # Реакция [CO] + [O]v -> CO2 + 2*
            if (state1 == 1 and state2 == 3) or (state1 == 3 and state2 == 1):
                events.append(('react_CO_OV', i1, j1, i2, j2))
                rates.append(self.k5 / 2)

            # Миграция [CO]
            if state1 == 1 and state2 == 0:
                events.append(('migrate_CO', i1, j1, i2, j2))
                rates.append(self.k6)
            elif state1 == 0 and state2 == 1:
                events.append(('migrate_CO', i2, j2, i1, j1))
                rates.append(self.k6)

            # Миграция [O]
            if state1 == 2 and state2 == 0:
                events.append(('migrate_O', i1, j1, i2, j2))
                rates.append(self.k7)
            elif state1 == 0 and state2 == 2:
                events.append(('migrate_O', i2, j2, i1, j1))
                rates.append(self.k7)

            # Миграция [O]v
            if state1 == 3 and state2 == 0:
                events.append(('migrate_OV', i1, j1, i2, j2))
                rates.append(self.k8)
            elif state1 == 0 and state2 == 3:
                events.append(('migrate_OV', i2, j2, i1, j1))
                rates.append(self.k8)
        return events, rates

    def execute_event(self, event):
        """Выполнить выбранное событие"""
        event_type = event[0]

        if event_type == 'adsorb_CO':
            _, i, j = event
            self.lattice[i, j] = 1  # [CO]

        elif event_type == 'desorb_CO':
            _, i, j = event
            self.lattice[i, j] = 0  # Свободно

        elif event_type == 'adsorb_O2':
            _, i1, j1, i2, j2 = event
            self.lattice[i1, j1] = 2  # [O]
            self.lattice[i2, j2] = 2  # [O]

        elif event_type == 'create_vacancy':
            _, i, j = event
            self.lattice[i, j] = 3  # [O]v

        elif event_type == 'react_CO_O':
            _, i1, j1, i2, j2 = event
            self.lattice[i1, j1] = 0  # Свободно
            self.lattice[i2, j2] = 0  # Свободно

        elif event_type == 'react_CO_OV':
            _, i1, j1, i2, j2 = event
            self.lattice[i1, j1] = 0  # Свободно
            self.lattice[i2, j2] = 0  # Свободно

        elif event_type == 'migrate_CO':
            _, from_i, from_j, to_i, to_j = event
            self.lattice[to_i, to_j] = 1  # [CO]
            self.lattice[from_i, from_j] = 0  # Свободно

        elif event_type == 'migrate_O':
            _, from_i, from_j, to_i, to_j = event
            self.lattice[to_i, to_j] = 2  # [O]
            self.lattice[from_i, from_j] = 0  # Свободно

        elif event_type == 'migrate_OV':
            _, from_i, from_j, to_i, to_j = event
            self.lattice[to_i, to_j] = 3  # [O]v
            self.lattice[from_i, from_j] = 0  # Свободно

    def step(self):
        """Выполнить один шаг алгоритма BKL"""
        events, rates = self.calculate_rates()

        if len(events) == 0:
            return False  # Нет возможных событий

        total_rate = sum(rates)
        if total_rate == 0:
            return False

        # Выбор события с вероятностью пропорциональной скорости
        r = random.random() * total_rate
        cumulative_rate = 0
        selected_event = None

        for event, rate in zip(events, rates):
            cumulative_rate += rate
            if cumulative_rate >= r:
                selected_event = event
                break

        # Выполнение события
        if selected_event:
            self.execute_event(selected_event)

            # Обновление времени
            dt = -np.log(random.random()) / total_rate
            self.time += dt
            self.step_count += 1

            return True

        return False

    def get_coverage(self):
        """Получить долю адсорбированного CO"""
        total_sites = self.size * self.size
        co_sites = np.sum(self.lattice == 1)
        return co_sites / total_sites

    def get_composition(self):
        """Получить состав поверхности"""
        total_sites = self.size * self.size
        return {
            'empty': np.sum(self.lattice == 0) / total_sites,
            'CO': np.sum(self.lattice == 1) / total_sites,
            'O': np.sum(self.lattice == 2) / total_sites,
            'O_v': np.sum(self.lattice == 3) / total_sites
        }

    def simulate(self, total_time=10.0, save_interval=0.1):
        """Провести моделирование"""
        print(f"Запуск моделирования на {total_time} секунд...")
        print(f"Размер решетки: {self.size}x{self.size}")
        print(f"Давление CO: {self.P_CO} Торр")

        next_save_time = 0.0

        # Начальное состояние
        coverage = self.get_coverage()
        self.coverage_history.append(coverage)
        self.time_history.append(self.time)
        print(f"Время: {self.time:.3f} с, Доля CO: {coverage:.4f}")

        while self.time < total_time:
            print(self.time)
            if not self.step():
                print("Нет возможных событий. Моделирование остановлено.")
                break

            # Сохраняем состояние через заданные интервалы
            if self.time >= next_save_time:
                coverage = self.get_coverage()
                self.coverage_history.append(coverage)
                self.time_history.append(self.time)

                if len(self.coverage_history) % 100 == 0:
                    print(f"Время: {self.time:.3f} с, Доля CO: {coverage:.4f}")

                next_save_time += save_interval

        # Финальное состояние
        coverage = self.get_coverage()
        self.coverage_history.append(coverage)
        self.time_history.append(self.time)
        print(f"Финальное время: {self.time:.3f} с, Доля CO: {coverage:.4f}")

        return self.time_history, self.coverage_history

    def plot_results(self, show_composition=True):
        """Построить графики результатов"""
        fig, axes = plt.subplots(2 if show_composition else 1, 1,
                               figsize=(12, 8 if show_composition else 6))

        if show_composition:
            ax1, ax2 = axes
        else:
            ax1 = axes
            ax2 = None

        # График доли CO
        ax1.plot(self.time_history, self.coverage_history, 'b-', linewidth=2)
        ax1.set_xlabel('Время, с')
        ax1.set_ylabel('Доля адсорбированного CO')
        ax1.set_title(f'Динамика адсорбции CO (решетка {self.size}x{self.size}, P_CO = {self.P_CO} Торр)')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 1)

        # График состава поверхности (если нужно)
        if show_composition and hasattr(self, 'composition_history'):
            time_comp = [t for t, _ in self.composition_history]
            empty = [comp['empty'] for _, comp in self.composition_history]
            co = [comp['CO'] for _, comp in self.composition_history]
            oxygen = [comp['O'] for _, comp in self.composition_history]
            oxygen_v = [comp['O_v'] for _, comp in self.composition_history]

            ax2.plot(time_comp, empty, 'k-', label='Свободно', linewidth=2)
            ax2.plot(time_comp, co, 'b-', label='CO', linewidth=2)
            ax2.plot(time_comp, oxygen, 'r-', label='O', linewidth=2)
            ax2.plot(time_comp, oxygen_v, 'g-', label='O_v', linewidth=2)
            ax2.set_xlabel('Время, с')
            ax2.set_ylabel('Доля поверхности')
            ax2.set_title('Состав поверхности')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            ax2.set_ylim(0, 1)

        plt.tight_layout()
        plt.show()

    def analyze_oscillations(self):
        """Анализ автоколебаний"""
        if len(self.coverage_history) < 10:
            print("Недостаточно данных для анализа колебаний")
            return

        # Поиск максимумов и минимумов
        coverage = np.array(self.coverage_history)
        time_arr = np.array(self.time_history)

        # Усреднение для сглаживания шума
        window_size = min(5, len(coverage) // 10)
        if window_size > 1:
            coverage_smooth = np.convolve(coverage, np.ones(window_size)/window_size, mode='valid')
            time_smooth = time_arr[window_size-1:]
        else:
            coverage_smooth = coverage
            time_smooth = time_arr

        # Поиск экстремумов
        from scipy.signal import find_peaks
        peaks, _ = find_peaks(coverage_smooth, height=0.1, distance=5)
        troughs, _ = find_peaks(-coverage_smooth, height=-0.9, distance=5)

        if len(peaks) > 1 and len(troughs) > 1:
            # Период колебаний
            peak_times = time_smooth[peaks]
            periods = np.diff(peak_times)
            avg_period = np.mean(periods)

            # Амплитуда
            amplitudes = coverage_smooth[peaks] - coverage_smooth[troughs[:len(peaks)]]
            avg_amplitude = np.mean(amplitudes)

            print(f"Анализ автоколебаний:")
            print(f"  Количество периодов: {len(periods)}")
            print(f"  Средний период: {avg_period:.3f} с")
            print(f"  Средняя амплитуда: {avg_amplitude:.3f}")

            return avg_period, avg_amplitude
        else:
            print("Явные автоколебания не обнаружены")
            return None, None


def study_size_effect():
    """Исследование влияния размера решетки"""
    sizes = [10, 20, 40]
    results = {}

    for size in sizes:
        print(f"\n{'='*50}")
        print(f"Моделирование для решетки {size}x{size}")
        print(f"{'='*50}")

        model = COOxidationCatalyst(size=size, initial_state='empty', P_CO=0.3)
        time_history, coverage_history = model.simulate(total_time=10.0)

        results[size] = {
            'time': time_history,
            'coverage': coverage_history,
            'model': model
        }

        # Анализ колебаний
        model.analyze_oscillations()

    # Сравнительный график
    plt.figure(figsize=(12, 8))
    for size, data in results.items():
        plt.plot(data['time'], data['coverage'], label=f'{size}x{size}', linewidth=2)

    plt.xlabel('Время, с')
    plt.ylabel('Доля адсорбированного CO')
    plt.title('Влияние размера решетки на динамику системы')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1)
    plt.show()

    return results


def study_pressure_effect():
    """Исследование влияния давления CO"""
    pressures = [0.05, 0.1, 0.2, 0.3, 0.4]
    results = {}

    for P_CO in pressures:
        print(f"\n{'='*50}")
        print(f"Моделирование для P_CO = {P_CO} Торр")
        print(f"{'='*50}")

        model = COOxidationCatalyst(size=20, initial_state='empty', P_CO=P_CO)
        time_history, coverage_history = model.simulate(total_time=10.0)

        results[P_CO] = {
            'time': time_history,
            'coverage': coverage_history,
            'model': model
        }

        # Анализ колебаний
        period, amplitude = model.analyze_oscillations()
        if period is not None:
            print(f"  Период колебаний: {period:.3f} с")
            print(f"  Амплитуда: {amplitude:.3f}")

    # Сравнительный график
    plt.figure(figsize=(12, 8))
    for P_CO, data in results.items():
        plt.plot(data['time'], data['coverage'], label=f'P_CO = {P_CO} Торр', linewidth=2)

    plt.xlabel('Время, с')
    plt.ylabel('Доля адсорбированного CO')
    plt.title('Влияние давления CO на динамику системы')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1)
    plt.show()

    return results


if __name__ == "__main__":
    model = COOxidationCatalyst(size=20, initial_state='empty', P_CO=0.3)
    model.simulate(total_time=10.0)
    model.plot_results(show_composition=True)
    model.analyze_oscillations()
    size_results = study_size_effect()
    pressure_results = study_pressure_effect()

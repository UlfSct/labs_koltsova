import numpy as np
import matplotlib.pyplot as plt


def read_data_file(filename):
    """Чтение файла данных"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Получаем x_values из первой строки (пропускаем "t\x")
    first_line = lines[0].split()
    x_values = np.array(list(map(float, first_line[1:])))

    # Читаем остальные данные
    data_lines = lines[1:]
    data = []
    y_values = []

    for line in data_lines:
        values = list(map(float, line.split()))
        y_values.append(values[0])  # первое значение - y
        data.append(values[1:])  # остальные - значения функции

    return np.array(data), np.array(y_values), x_values


def plot_selected_lines(filename, num_lines=10):
    """Строит 10 равномерно распределенных линий на одном графике с легендой"""

    # Читаем данные
    data, y_values, x_values = read_data_file(filename)

    # Выбираем индексы для 10 равномерно распределенных линий
    total_lines = len(y_values)
    indices = np.linspace(0, total_lines - 1, num_lines, dtype=int)

    # Создаем график
    plt.figure(figsize=(12, 8))

    # Строим выбранные линии
    for idx in indices:
        y_val = y_values[idx]
        row_data = data[idx]
        plt.plot(x_values, row_data, linewidth=2, label=f't = {y_val:.2f}')

    # Настраиваем график
    plt.xlabel('X', fontsize=12)
    plt.ylabel('Z', fontsize=12)
    plt.title(f'График из файла: {filename}', fontsize=14)
    plt.grid(True, alpha=0.3)

    # Размещаем легенду
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.tight_layout()
    plt.show()

    # Выводим информацию о выбранных линиях
    print(f"Всего строк в файле: {total_lines}")
    print("Отображены строки с y-values:")
    for idx in indices:
        print(f"  y = {y_values[idx]:.2f}")


plot_selected_lines("./results/lab5/n_1001_j_.txt")
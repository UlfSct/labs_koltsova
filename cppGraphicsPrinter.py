import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob


def read_data_file(filename):
    """Чтение одного файла данных"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Пропускаем первую строку с заголовками
    data_lines = lines[1:]

    # Извлекаем данные
    data = []
    y_values = []

    for line in data_lines:
        values = list(map(float, line.split()))
        y_values.append(values[0])  # первое значение - y
        data.append(values[1:])  # остальные - значения функции

    # Преобразуем в numpy массивы
    data_array = np.array(data)
    y_array = np.array(y_values)

    return data_array, y_array


def create_simple_3d_plot(file_signature):
    """Создание одного 3D графика из первого найденного файла"""

    # Находим файлы по сигнатуре
    file_pattern = f"{file_signature}"
    files = sorted(glob.glob(file_pattern))

    if not files:
        print(f"Файлы по шаблону '{file_pattern}' не найдены!")
        return

    # Берем первый файл
    filename = files[0]
    print(f"Используем файл: {filename}")

    # Читаем данные
    data, y_values = read_data_file(filename)

    # Получаем x_values из первой строки файла
    with open(filename, 'r') as f:
        first_line = f.readline().split()
        x_values = np.array(list(map(float, first_line[1:])))  # пропускаем "t\x"

    # Создаем сетку
    X, Y = np.meshgrid(x_values, y_values)

    # Создаем 3D график
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Рисуем поверхность
    surf = ax.plot_surface(X, Y, data, cmap='viridis', alpha=0.9)

    # Настройки графика
    ax.set_xlabel('L')
    ax.set_ylabel('T')
    ax.set_zlabel('F')
    ax.set_title(f'3D график: {filename}')

    # Цветовая шкала
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

    plt.show()


# Использование
create_simple_3d_plot("./results/lab1/ksi_n_10.txt")
create_simple_3d_plot("./results/lab1/psi_n_10.txt")
create_simple_3d_plot("./results/lab1/vx_n_10.txt")
create_simple_3d_plot("./results/lab1/vy_n_10.txt")
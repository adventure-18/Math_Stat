import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma, gammainc


def LowerIncGamma(x, alpha):
    return gamma(alpha) * gammainc(alpha, x)


# Создаем массивы параметров
a_values = np.array([2, 1, 2, 5, 7])
d_values = np.array([0.5, 1, 1, 1, 1])
p_values = np.array([0.5, 0.5, 2, 5, 7])

# Создаем массив значений x от 0.0001 до 7 с шагом 0.01
x_values = np.arange(0.0001, 7, 0.01)

# Создаем пустые списки для хранения результатов
pdf_values = []
cdf_values = []
params_values = []

# Вычисляем значения PDF и CDF для всех комбинаций параметров и значений x
for a, d, p in zip(a_values, d_values, p_values):
    params = f'a={a}, d={d}, p={p}'
    params_values.extend([params] * len(x_values))

    pdf = (p / a ** d) * x_values ** (d - 1) * np.exp(-(x_values / a) ** p) / gamma(d / p)
    pdf_values.extend(pdf)

    cdf = gammainc(d / p, (x_values / a) ** p)
    cdf_values.extend(cdf)

# Выводим результаты
for params, x, pdf, cdf in zip(params_values, x_values, pdf_values, cdf_values):
    print(f"Params: {params}, x: {x}, PDF: {pdf}, CDF: {cdf}")

# Построение графиков
for i in range(len(a_values)):
    plt.figure(figsize=(10, 4))
    plt.plot(x_values, pdf_values[i * len(x_values):(i + 1) * len(x_values)], label='PDF')
    plt.plot(x_values, cdf_values[i * len(x_values):(i + 1) * len(x_values)], label='CDF')
    plt.title(f'a={a_values[i]}, d={d_values[i]}, p={p_values[i]}')
    plt.xlabel('x')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True)
    plt.show()

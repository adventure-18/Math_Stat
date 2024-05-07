import csv
import random
import math
import matplotlib.pyplot as plt
import scipy.special
import scipy.special as sc


def read_data():
    with open("7. Обобщенное_гамма-распределение_var_7.csv", encoding='utf-8') as r_file:
        # Создаем объект reader, указываем символ-разделитель ","
        file_reader = csv.reader(r_file, delimiter=",")
        nums = []
        count = 0
        for row in file_reader:
            nums.append(float(row[0]))
            count += 1
        return nums


def sum_func(nums):
    sum = 0
    for i in range(0, len(nums)):
        sum += nums[i]
    # print(f"Sum {sum}")
    return sum


def srednee_func(nums):
    srednee = sum_func(nums) / len(nums)
    # print(f"Srednee {srednee}")
    return srednee


def median(arr):
    sorted_arr = sorted(arr)
    n = len(arr)
    if n % 2 == 0:
        # Если количество элементов четное, медиана - среднее двух центральных элементов
        return (sorted_arr[n // 2 - 1] + sorted_arr[n // 2]) / 2
    else:
        # Если количество элементов нечетное, медиана - центральный элемент
        return sorted_arr[n // 2]


def mode(arr):
    # Создаем словарь для подсчета количества вхождений каждого элемента
    counts = {}
    for num in arr:
        if num in counts:
            counts[num] += 1
        else:
            counts[num] = 1

    # Находим максимальное количество вхождений
    max_count = max(counts.values())

    # Создаем список, содержащий все элементы, которые встречаются максимальное количество раз
    modes = [num for num, count in counts.items() if count == max_count]

    return modes[0]


def range_func(data):
    return max(data) - min(data)


def shifted_disp(data, xmed):
    disp_sqr = 0
    count = 0
    for xi in data:
        count += 1
        disp_sqr += (xi - xmed) ** 2
    return disp_sqr / count


def non_shifted_disp(data, xmed):
    disp_sqr = 0
    count = 0
    for xi in data:
        count += 1
        disp_sqr += (xi - xmed) ** 2
    return disp_sqr / (count - 1)


def start_moment(data, k):
    res = 0
    count = 0
    for xi in data:
        count += 1
        res += xi ** k
    return res / count


def mid_moment(data, xmed, k):
    res = 0
    count = 0
    for xi in data:
        count += 1
        res += (xi - xmed) ** k
    return res / count


def empirical_func(data, elem):
    n = len(data)
    data = sorted(data)
    # elem = data
    ecdf = []
    for point in elem:
        count = sum(1 for d in data if d <= point)
        ecdf.append(count / n)
    return ecdf


def draw_graph(ecdf, x, k):
    plt.step(x, ecdf, where='post', label='Эмпирическая функция распределения')
    plt.xlabel('X')
    plt.ylabel('F(X)')
    plt.title(f'Эмпирическая функция распределения для выборки {k}')
    plt.legend()
    plt.grid(True)
    plt.show()


def draw_gisto(data, bins_param):
    plt.hist(data, bins=int(bins_param * 0.75), density=True, alpha=0.75, edgecolor='black', label='Гистограмма данных')
    plt.xlabel('X')
    plt.ylabel(f"Плотность вероятности для подвыборки")
    plt.title(f'Гистограмма данных для подвыборки из {bins_param} элементов')
    plt.legend()
    plt.grid(True)
    plt.show()

from scipy import integrate

def integrand(t, s):
    return t**(s-1) * math.exp(-t)

def gamma_low(s, x):
    result, _ = integrate.quad(integrand, 0, x, args=(s,))
    return result
from scipy import integrate
import math

def integrand_g(t, s):
    return t**(s-1) * math.exp(-t)

def gamma(s, x):
    result, _ = integrate.quad(integrand_g, x, float('inf'), args=(s,))
    return result

# Пример использования
s = 0
x = 2
result = gamma(s, x)
print("Γ(", s, ",", x, ") =", result)




# Пример использования
s = 2
x = 5
result = gamma_low(s, x)
print("γ(", s, ",", x, ") =", result)
print(f"Г встроенная {scipy.special.gamma(s)}")
def graphics_and_gistos():
    sorted_data = sorted(data)

    sub_data = random.sample(data, k=10)
    # empirical func for 10
    emp_res = empirical_func(sub_data, sorted_data)
    draw_graph(emp_res, sorted_data, 10)
    draw_gisto(sub_data, 10)

    sub_data = random.sample(data, k=100)
    # empirical func for 100
    emp_res = empirical_func(sub_data, sorted_data)
    draw_graph(emp_res, sorted_data, 100)
    draw_gisto(sub_data, 100)

    sub_data = random.sample(data, k=200)
    # empirical func for 200
    emp_res = empirical_func(sub_data, sorted_data)
    draw_graph(emp_res, sorted_data, 200)
    draw_gisto(sub_data, 200)

    # empirical func for all
    emp_res = empirical_func(data, sorted_data)
    draw_graph(emp_res, sorted_data, len(data))
    # draw_gisto(emp_res)


def my_gamma(z):
    return sc.gamma(z)


# Пример использования
z = 0.6
print("Гамма-функция для z =", z, ":", my_gamma(z))


def incomplete_gamma_function(s, x):
    return sc.gammainc(s, x)


z = 0.6
print("Гамма-функция для z =", z, ":", incomplete_gamma_function(8, 3))


def generalized_gamma_distrib(t, θ, k, δ):
    pass


data = read_data()

if __name__ == '__main__':
    summa = sum_func(data)
    srednee = srednee_func(data)
    mediana = median(data)
    moda = mode(data)
    _range = range_func(data)
    disp_sqr = shifted_disp(data, srednee)
    non_disp_sqr = non_shifted_disp(data, srednee)
    start_mom = start_moment(data, 3)  # k
    mid_mom = mid_moment(data, srednee, 3)  # k
    # graphics_and_gistos()

    print(summa)
    print(srednee)
    print(mediana)
    print(moda)
    print(_range)
    print(disp_sqr)
    print(non_disp_sqr)
    print(start_mom)
    print(mid_mom)

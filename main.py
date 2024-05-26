import csv
import random
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import scipy.special as sc
import scipy.stats as stats
import mpmath


# uerieirireire

def read_data():
    with open("7. Обобщенное_гамма-распределение_var_7.csv", encoding='utf-8') as r_file:
        # Создаем объект reader, указываем символ-разделитель ","
        file_reader = csv.reader(r_file, delimiter=",")
        nums = []
        count = 0
        for row in file_reader:
            nums.append(float(row[0]))
            count += 1
        print('we')
        return nums


def sum_func(nums):
    sum = 0
    for i in range(0, len(nums)):
        sum += nums[i]
    return sum


def average_func(nums):
    avg = sum_func(nums) / len(nums)
    return avg


def median(arr):
    sorted_arr = sorted(arr)
    n = len(arr)
    if n % 2 == 0:
        return (sorted_arr[n // 2 - 1] + sorted_arr[n // 2]) / 2
    else:
        return sorted_arr[n // 2]


def mode(arr):
    counts = {}
    for num in arr:
        if num in counts:
            counts[num] += 1
        else:
            counts[num] = 1

    max_count = max(counts.values())

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
    return t ** (s - 1) * math.exp(-t)


def gamma_low(s, x):
    result, _ = integrate.quad(integrand, 0, x, args=(s,))
    return result


from scipy import integrate
import math


def integrand_g(t, s):
    return t ** (s - 1) * math.exp(-t)


def gamma(s, x):
    result, _ = integrate.quad(integrand_g, x, float('inf'), args=(s,))
    return result


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


def LowerIncGamma(x, alpha):
    return my_gamma(alpha) * sc.gammainc(alpha, x)


def incomplete_gamma_function(s, x):
    return sc.gammainc(s, x)


# print("Гамма-функция для z =", z, ":", incomplete_gamma_function(8, 3))


def generalized_gamma_distrib(d, p, a, x):
    return incomplete_gamma_function(d / p, (x / a) ** p) / my_gamma(d / p)


def gen2(d, p, a, x):
    return mpmath.gammainc((x / a) ** p, d / p) / my_gamma(d / p)


data = read_data()

from scipy.stats import gengamma
from scipy.stats import gamma

import sympy as sp


def pdf_gener(x, b, c, theta):
    return (abs(c) * ((x / theta) ** (b - 1)) / (theta * my_gamma(b / c))) * sp.exp(-1 * ((x / theta) ** c))


def gener_graph(x_values, y_values, a, d, p):
    plt.plot(x_values, y_values)
    plt.xlabel('Аргументы')
    plt.ylabel('Значения функции')
    plt.title(f'График функции\na={a}, d={d:.2}, p={p:.2}')
    plt.grid(True)
    plt.show()


def pdf_grapher():
    res_mas = [pdf_gener(i, 5.0, 1.0, 5) for i in range(50)]
    gener_graph([i for i in range(50)], res_mas, 2.0, 0.5, 0.5)


def integral_lower_gamma(t, c, b):
    x = sp.symbols('x')

    f = (x ** (b - 1)) * sp.exp(-1 * (x ** c))

    integral = sp.integrate(f, (x, 0, t))

    return integral


def new_distr(t, b, c, a):
    tt = t / a

    integral = integral_lower_gamma(tt, c, b)

    res = abs(c) * integral / my_gamma(b / c)

    return res


def sum1_func(a, p):
    s = 0
    for i in range(len(data)):
        s += (abs(data[i]) / a) ** p
    return s


def sum2_func(a):
    s = 0
    for i in range(len(data)):
        s += math.log(abs(data[i]) / a)
    return s


def sum3_func(a, p):
    s = 0
    for i in range(len(data)):
        s += ((abs(data[i]) / a) ** p) * math.log(abs(data[i]) / a)
    return s


def truefunc(d, p, a):
    prod = 1
    for x in data:
        xi = abs(x)
        res = (abs(p) * (xi / a) ** (d - 1)) * math.exp(-1 * (xi / a) ** p) / (a * my_gamma(d / p))
        prod = prod * res
    print('res1', prod)
    res1 = prod
    d2 = 3
    p2 = 1
    a2 = 1
    prod = 1
    for x in data:
        xi = abs(x)
        res33 = (abs(p2) * (xi / a2) ** (d2 - 1)) * math.exp(-1 * (xi / a2) ** p2) / (a2 * my_gamma(d2 / p2))
        prod = prod * res33
    print('res2', prod)
    res2 = prod
    y = 2 * math.log(res1 / res2)
    print(y)


if __name__ == '__main__':
    n = len(data)
    for a in range(1, 5):
        for d in range(1, 100):
            p = 1
            eq1 = (d - (p / n) * sum1_func(a, p))
            eq2 = sc.digamma(d / p) - (p / n) * sum2_func(a)
            # eq3 = 1 / p + (d / (p ** 2)) * sc.digamma(d / p) - (1 / n) * sum3_func(a, p)
            if abs(eq1) < 1 and abs(eq2) < 1:
                print(f' d = {d}, p = {p}, a = {a} res = {abs(eq1) + abs(eq2)}')
    # d = 1, p = 2, a = 4
    truefunc(36, 9, 3)
    # PPP d = 36 p = 9 a = 3
    # PPP d = 2 p = 1 a = 2
    exit(0)
    # d = 3, p = 1, a = 1
    # pass
    # print(data)
    # summa = sum_func(data)
    # srednee = average_func(data)
    # mediana = median(data)
    # moda = mode(data)
    # _range = range_func(data)
    # disp_sqr = shifted_disp(data, srednee)
    # non_disp_sqr = non_shifted_disp(data, srednee)
    # start_mom = start_moment(data, 3)  # k
    # mid_mom = mid_moment(data, srednee, 3)  # k
    # graphics_and_gistos()
    #
    # print(summa)
    # print(srednee)
    # print(mediana)
    # print(moda)
    # print(_range)
    # print(disp_sqr)
    # print(non_disp_sqr)
    # print(start_mom)
    # print(mid_mom)
    # print('-----------------------')
    # vrand = random.randint(1, 15) / random.randint(1, 2) + 1
    # krand = random.randint(0, 10) / random.randint(5, 23)
    # sigmarand = random.randint(0, 10) / random.randint(5, 23)
    #
    # x = sp.symbols('x')
    #
    # #pdf_grapher()

    # Определяем функцию
    # f = x ** 2

    # Вычисляем определенный интеграл от 0 до 1
    # integral = sp.integrate(f, (x, 0, 1))

    # sdf = integral * 3

    a = 2
    d = 1.0
    p = 2.0
    cpdata = data
    cpdata = sorted(cpdata)
    sub = random.sample(cpdata, k=200)
    print(sub)
    shape_parameter = d
    rv = gengamma(d, p)
    cdf_value = gengamma.cdf(4, d, p)
    # cdf_value = gamma.cdf(4, a=1, scale=1)
    th = my_gamma(4)
    # res = generalized_gamma_distrib(d, p, a, 4)
    res = sc.gammainc(d / p, (4 / a) ** p)
    # scipy.gengamma()
    buffer = []
    xs = []
    x_values = np.arange(0.0001, 7, 0.01)
    for ind in sorted(data):
        res = generalized_gamma_distrib(d, p, a, abs(ind))
        xs.append(res)
    gener_graph([abs(it) for it in sorted(data)], xs, a, d, p)
    exit(2)
    # print(vrand, krand, sigmarand)
    for i in range(0, 10):
        # res = gen2(d, p, a, i)
        # res = sc.gammainc(d / p, (i / a) ** p)
        # res = new_distr(i, d, p, a)
        res = generalized_gamma_distrib(d, p, a, i)
        print(f'Iter={i}: x = {sub[i]}, res = {res}\n')
        buffer.append(res)
    # gener_graph(sub, buffer, a, d, p)
    gener_graph([i for i in range(0, 10)], buffer, a, d, p)

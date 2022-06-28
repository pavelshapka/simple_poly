# -*- coding: utf-8 -*-
import numpy as np


def search_inv(a, P):  # Поиск обратного элемента в поле характеристики P
    for i in range(1, P):
        if (a * i) % P == 1:
            return i
    return 0


P = int(input("Ввдети характеристику поля: "))
st = input("Введите строку коэффициентов: ")
lst = list(map(int, list(st)))
p_0 = np.poly1d(np.array(lst) % P)  # Исходный многочлен
n = np.arange(0, P, 1)
deg = p_0.order  # Степень исходного многочлена

der = np.poly1d(p_0.deriv().c % P)  # Производная исходного многочлена
left = p_0
right = der
r0 = left

if der == np.poly1d([0]):  # Если производная 0 - то она делится на сам многочлен - сразу выходим
    print("False")
    exit()

elif p_0.order <= 1:  # Если степень 1 - сразу выходим
    print("True")
    exit()

elif p_0.order <= 3:  # Если степень <= 3 и нет корней - выходим
    koef = list(reversed(p_0.c))
    for i in range(P):
        tmp_sum = 0
        for j in range(len(koef)):
            tmp_sum += i ** j * koef[j]
        if tmp_sum % P == 0:
            print("False")
            exit()
    print("True")
    exit()
else:  # Проверяем на взаимную простоту
    while r0.order >= right.order and r0.order != 0:  # Проверка на простоту со своей производной
        # Ищем значение, на которе надо домножить
        first = left.c[0]
        tmp = right
        tmp *= np.poly1d([1, 0]) ** (left.order - right.order)
        k = first * search_inv(tmp.c[0], P) % P
        r0 = np.poly1d((left - k * tmp).c % P)
        left = r0
        right = right
    if r0.c[0] == 0:  # Если не взаимно просты, сразу выходим
        print('False')
        exit()

pattern = np.poly1d(-p_0.c[1:] * search_inv(p_0.c[0], P) % P)
lst_obr = []

# Базис
for i in range(deg):
    koef_tmp = np.zeros(i + 1)
    koef_tmp[0] = 1
    p_tmp = np.poly1d(koef_tmp)
    # Ищем образы:
    p_obr = np.poly1d((p_tmp ** P - p_tmp).c % P)
    while p_obr.order >= deg:
        fir = p_obr.c[0]  # первый коэффициент
        k = fir * search_inv(p_0.c[0], P) % P
        tmp_pattern = pattern
        tmp_pattern *= np.poly1d([1, 0]) ** (p_obr.order - deg)  # Домножаем шаблон нужное кол-во раз на x
        p_obr += np.poly1d(k * tmp_pattern)  # Прибавляем шаблон, умноженный на первый коэффицинет
        p_obr = np.poly1d(p_obr.c[1:] % P)  # Убираем старшую степень
    koef = list(reversed(np.array(p_obr.c) % P))
    while len(koef) < deg:
        koef.append(0)  # (x^2 + x + 1)^2 + x^2 + x + 1 = x^4 + x^2 + 1 + x^2 + x + 1 = x^4 + x = x^2 + x + 1
    lst_obr.append(np.array(koef))

# Приведение к диагональному виду
for i in range(1, deg):
    if lst_obr[i][i] == 0:
        for j in range(deg):
            if lst_obr[j][i] != 0:
                lst_obr[i], lst_obr[j] = lst_obr[j], lst_obr[i]
    if lst_obr[i][i] != 0:
        for j in range(deg):
            if i != j:
                lst_obr[j] = (lst_obr[j] - lst_obr[i] * search_inv(lst_obr[i][i], P) * lst_obr[j][i]) % P

matrix = np.asmatrix(lst_obr).transpose()
rank = np.linalg.matrix_rank(matrix)
print(matrix)
if rank == deg - 1:
    print("True")
else:
    print("False")

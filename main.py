import numpy as np
import itertools as it

def search_inv(a, P): # Поиск обратного элемента в поле характеристики P
    for i in range(2, P):
        if (a * i) % P == 1:
            return i
    return 1

P = int(input("Ввдети характеристику поля: "))
st = input("Введите строку коэффициентов: ")
lst = list(map(int, list(st)))
pol = np.poly1d(np.array(lst) % P) # Исходный многочлен
n = np.arange(0, P, 1)
deg = pol.order # Степень исходного многочлена
data = list(it.product(n, repeat=deg)) # Создаем все вычеты по модулю многочлена pol

der = np.poly1d(pol.deriv().c % P) # Производная исходного многочлена
left = pol
right = der
r0 = left
while r0.order >= right.order and r0.order != 0: # Проверка на простоту со своей производной
    # Ищем значение, на которе надо домножить
    first = left.c[0]
    tmp = right
    for i in range(left.order - right.order):
        tmp *= np.poly1d([1, 0])
    q = search_inv(tmp.c[0], P)
    k = first * search_inv(tmp.c[0], P) % P
    r0 = np.poly1d((left - k*tmp).c % P)
    left = r0
    right = right


if r0.c[0] == 0: # Если не взаимно просты, сразу выходим
    print('False')
    exit()

pattern = np.poly1d(-pol.c[1:] * search_inv(pol.c[0], P) % P) # Ищем правило выражения в поле вида x^n = x^(n-1) + x^(n-2) + ... + x + a
count = 0
for i in data: # Проверка на количество корней уравнения h^P - h = 0 (Проверяем все элементы поля)
    tmp = np.poly1d(i)
    m = np.poly1d(i)
    for q in range(P - 1): # Возводим в степень
        tmp *= m
    tmp = np.poly1d(tmp.c % P)
    while tmp.order >= deg:
        fir = tmp.c[0] # первый коэффициент
        tmp_pattern = pattern
        for j in range(deg, tmp.order):
            tmp_pattern *= np.poly1d([1, 0]) # Домножаем шаблон нужное кол-во раз на x
        tmp += np.poly1d(fir * tmp_pattern) # Прибавляем шаблон, умноженный на первый коэффицинет
        tmp = np.poly1d(tmp.c[1:] % P) # Убираем старшую степень
    if (m - tmp) == np.poly1d(0): # Увеличиваем счетчик, если данный многочлен - решение
        count += 1


if count != P:
    print('False')
else:
    print('True')



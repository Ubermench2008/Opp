#!/bin/bash

# Функция вывода справки
usage() {
    echo "Usage: $0 -v1|-v2 mode [matrix_size]"
    echo "  -v1          Использовать версию 1 (lab1_omp_var1.cpp -> Omp_variant1)"
    echo "  -v2          Использовать версию 2 (lab1_omp_var2.cpp -> Omp_var2)"
    echo "  mode         Режим работы: 1 или 2"
    echo "  matrix_size  Размер матрицы (только для режима 2)"
    exit 1
}

# Проверка количества аргументов
if [ $# -lt 2 ]; then
    usage
fi

# Инициализация переменных
version_flag=""
mode=""
matrix_size=""

# Разбор аргументов
while [[ $# -gt 0 ]]; do
    case $1 in
        -v1)
            version_flag="v1"
            shift
            ;;
        -v2)
            version_flag="v2"
            shift
            ;;
        1|2)
            mode=$1
            shift
            ;;
        *)
            if [ -z "$matrix_size" ]; then
                matrix_size=$1
                shift
            else
                usage
            fi
            ;;
    esac
done

# Проверка корректности аргументов
if [ -z "$version_flag" ] || [ -z "$mode" ]; then
    usage
fi

if [ "$mode" -eq 2 ] && [ -z "$matrix_size" ]; then
    echo "Ошибка: для режима 2 необходимо указать размер матрицы."
    usage
fi

# Опционально: задаём число потоков OpenMP (можно изменить или убрать)
export OMP_NUM_THREADS=2

# Компиляция и запуск в зависимости от выбранной версии
if [ "$version_flag" == "v1" ]; then
    echo "Компиляция версии 1 (lab1_omp_var1.cpp)..."
    g++ -fopenmp -O2 lab1_omp_var1.cpp matrix_generators.cpp -o Omp_variant1
    if [ "$mode" -eq 1 ]; then
        echo "Запуск Omp_variant1 с режимом 1"
        ./Omp_variant1 1
    else
        echo "Запуск Omp_variant1 с режимом 2 и размером матрицы $matrix_size"
        ./Omp_variant1 2 $matrix_size
    fi
elif [ "$version_flag" == "v2" ]; then
    echo "Компиляция версии 2 (lab1_omp_var2.cpp)..."
    g++ -fopenmp -O2 lab1_omp_var2.cpp matrix_generators.cpp -o Omp_var2
    if [ "$mode" -eq 1 ]; then
        echo "Запуск Omp_var2 с режимом 1"
        ./Omp_var2 1
    else
        echo "Запуск Omp_var2 с режимом 2 и размером матрицы $matrix_size"
        ./Omp_var2 2 $matrix_size
    fi
else
    usage
fi

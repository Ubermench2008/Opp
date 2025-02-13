#include <iostream>                     // Для стандартного ввода-вывода
#include <vector>                       // Для динамических массивов (векторов)
#include <cmath>                        // Для математических функций (sqrt, sin и т.д.)
#include <chrono>                       // Для измерения времени выполнения
#include <iomanip>                      // Для манипуляции форматированием вывода
#include <cstdlib>                      // Для стандартных функций (например, atoi)
#include <string>                       // Для работы со строками
#include "matrix_generators.h"          // Пользовательский заголовочный файл с функциями генерации матриц/векторов

using namespace std;
using namespace std::chrono;

// Функция для вывода подсказки по использованию программы
void printUsage(const char* progName) {                
    cerr << "Использование:\n"
         << progName << " 1 # Ручной ввод\n"  
         << progName << " 2 # Генерация матрицы/вектора\n";
}

// Функция для запроса целочисленного ввода с сообщением
int promptInt(const string &message) {                
    int value;
    cout << message;
    cin >> value;
    return value;
}

// Функция для запроса вещественного ввода с сообщением (не используется в данном варианте)
double promptDouble(const string &message) {          
    double value;
    cout << message;
    cin >> value;
    return value;
}

// Функция для ручного ввода матрицы и вектора (режим 1)
void inputMatrixAndVector(int &N, vector<vector<double>> &A, vector<double> &b) {  
    N = promptInt("Введите размер матрицы [N]: ");  
    A.assign(N, vector<double>(N, 0.0));             
    b.assign(N, 0.0);                                

    cout << "Введите матрицу A (" << N << "x" << N << "):\n";  
    for (int i = 0; i < N; ++i) {                    
        for (int j = 0; j < N; ++j) {                
            cin >> A[i][j];                          
        }
    }

    cout << "Введите вектор b (" << N << " элементов):\n";  
    for (int i = 0; i < N; ++i) {                    
        cin >> b[i];                               
    }
}

// Функция для получения значения N из аргументов командной строки или установки значения по умолчанию
int getNFromArgsOrDefault(int argc, char* argv[], int defaultN) {  
    int N = defaultN;                              
    if (argc >= 3) {                               
        N = atoi(argv[2]);                         
        if (N <= 0) {                              
            cerr << "Некорректное N, используем N=" << defaultN << " по умолчанию\n";  
            N = defaultN;                          
        }
    } else {                                       
        cerr << "N не задано, используем N=" << defaultN << "\n";  
    }
    return N;                                      
}

// Функция для вычисления нормы (евклидовой) вектора
double computeNorm(const vector<double>& vec) {    
    double sum = 0.0;                              
    for (double x : vec) {                         
        sum += x * x;                              
    }
    return sqrt(sum);                              
}

// Функция для вычисления невязки r = A*x - b
vector<double> computeResidual(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {  
    int N = x.size();                              
    vector<double> r(N, 0.0);                        
    for (int i = 0; i < N; ++i) {                    
        double sum = 0.0;                          
        for (int j = 0; j < N; ++j) {                
            sum += A[i][j] * x[j];                  
        }
        r[i] = sum - b[i];                         
    }
    return r;                                      
}

int main(int argc, char* argv[]) {                 
    if (argc < 2) {                                
        printUsage(argv[0]);                       
        return 1;                                  
    }

    int mode = atoi(argv[1]);                      
    int N = 0;                                     
    vector<vector<double>> A;                      
    vector<double> b;                              

    if (mode == 1) {                               
        inputMatrixAndVector(N, A, b);             
    } else if (mode == 2) {                        
        N = getNFromArgsOrDefault(argc, argv, 1000);  
        generateMatrix(N, A, b);                   
        cout << "Сгенерирована матрица " << N << "x" << N 
             << "\nСгенерирован вектор b со значениями: " << (N + 1) << "\n";
    } else {                                       
        printUsage(argv[0]);                       
        return 1;                                  
    }

    // Используем значение tau, предлагаемое функцией suggestTau
    double tau = suggestTau(N, 0.95);
    cout << "\nИспользуем tau = " << tau << "\n";
    
    const double epsilon = 0.00001;                
    cout << "Точность epsilon = " << epsilon << endl;
    int maxIterations = promptInt("Введите maxIterations: ");

    vector<double> x(N, 0.0);
    double normB = computeNorm(b);
    if (normB == 0.0) {                            
        normB = 1.0;
    }

    auto start = high_resolution_clock::now();
    int iteration = 0;
    while (iteration < maxIterations) {
        vector<double> r = computeResidual(A, x, b);
        double normR = computeNorm(r);

        if (normR / normB < epsilon) {             
            break;
        }

        for (int i = 0; i < N; ++i) {
            x[i] -= tau * r[i];
        }
        ++iteration;
    }

    auto end = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(end - start);

    double elapsedSeconds = static_cast<double>(duration.count()) / 1000000.0;

    cout << "\nМетод завершён за " << iteration << " итераций.\n";
    cout << "Время выполнения: " << fixed << setprecision(6) << elapsedSeconds << " секунд.\n";

    vector<double> rFinal = computeResidual(A, x, b);
    double normFinal = computeNorm(rFinal);
    cout << "\n||A*x - b|| = " << fixed << setprecision(6) << normFinal << "\n\n";

    cout << fixed << setprecision(6);
    int nPrint = (N < 50 ? N : 50);
    cout << "Первые " << nPrint << " компонент x:\n";
    for (int i = 0; i < nPrint; ++i) {
        cout << "x[" << i << "] = " << x[i] << "\n";
    }

    return 0;
}

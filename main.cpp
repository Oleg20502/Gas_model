#include <iostream>
#include <fstream>
#include <string>

#include "molecule.h"
#include "vector_functions.h"
#include "LJP.h"
#include "space.h"


int CountLinesInFile(std::string filename)
{
  std::ifstream F(filename, std::ios::in);
  if (!F) {
    return -1;
  }
  int count = 0;
  char buffer[1000];
  while (!F.eof()){
    ++count;
    F.getline(buffer, 1000);
  }
  F.close();
  return count;
}

int GetStringsFromFileC(std::string filename, char*** _lines = nullptr)
{
  // 1. Дополнительные переменные
  char** lines;
  int n = CountLinesInFile(filename); // получить количество строк в файле

  // 2. Проверка, известно ли количество строк
  if (n == -1) return -1;

  // 3. Объявить файловую переменную
  std::ifstream F(filename); // открыть файл для чтения

  // 4. Проверка, открылся ли файл
  if (!F) return -1;

  // 5. Попытка выделить память для n строк
  try
  {
    lines = new char* [n];
  }
  catch (std::bad_alloc e)
  {
    // если невозможно выделить память, то выход
    std::cout << e.what() << std::endl; // вывести сообщение об ошибке
    F.close(); // закрыть файл
    return -1;
  }

  // 6. Цикл чтения строк и выделения памяти для них
  int len; // длина одной строки
  char buffer[1000]; // память, куда записывается одна строка из файла

  for (int i = 0; i < n; i++)
  {
    // 6.1. Считать строку из файла
    F.getline(buffer, 1000);

    // 6.2. Определить длину прочитанной строки
    for (len = 0; buffer[len] != '\0'; len++);

    // 6.3. Выделить память для len+1 символов
    lines[i] = new char[len + 1];

    // 6.4. Скопировать данные из buffer в lines[i]
    for (int j = 0; j < len; j++)
      lines[i][j] = buffer[j];
    lines[i][len] = '\0'; // добавить символ конца строки
  }

  // 7. Закрыть файл
  F.close();

  // 8. Записать результат
  *_lines = lines;
  return n;
}

bool ChangeStringInFileC(std::string filename, int position, std::string str)
{
  // 1. Получить строки файла в виде списка
  char** lines; // список строк файла
  int count; // количество строк файла
  count = GetStringsFromFileC(filename, &lines); // получить список lines

  // 2. Проверка, корректно ли прочитан файл
  if (count < 0) return false;

  // 3. Проверка, корректна ли позиция 0 <= position < count
  if ((position < 0) || (position >= count)) return false;

  // 4. Запись строк lines в файл до позиции position
  std::ofstream F(filename); // открыть файл для записи

  // 5. Проверка, открылся ли файл корректно - функция is_open()
  if (!F.is_open()) return false;

  if (position < count-1){
      for (int i = 0; i < position; i++)
        F << lines[i] << std::endl; // вывести строку в файл

      // 6. Запись строки с позиции position
      F << str.c_str() << std::endl; // здесь пишется строка str

      // 7. Запись строк после позиции position
      for (int i = position + 1; i < count - 1; i++)
        F << lines[i] << std::endl;

      // 8. Записать последнюю строку без символа '\n'
      F << lines[count - 1];
  }
  else{
      for (int i = 0; i < count - 1; i++)
        F << lines[i] << std::endl;

      // 8. Записать последнюю строку без символа '\n'
      F << str.c_str();
  }

  // 9. Закрыть файл
  F.close();

  // 10. Освободить память, выделенную для строк
  for (int i = 0; i < count; i++)
    delete lines[i];

  // 11. Освободить указатели на строки
  delete[] lines;
  return true;
}


template<typename type>
void process(std::string path)
{
    unsigned int N;
    int a;
    type T, dt, tau, x_size, y_size, z_size, eps, sigma, k, rmin, v_mean, v_std;
    bool save;
    std::ifstream out(path);
    std::string str1, str2;
    if (out.is_open()){
        getline(out, str2);
        getline(out, str1);
        N = stoi(str1);

        getline(out, str2);
        getline(out, str1);
        x_size = stod(str1);

        getline(out, str2);
        getline(out, str1);
        y_size = stod(str1);

        getline(out, str2);
        getline(out, str1);
        z_size = stod(str1);

        getline(out, str2);
        getline(out, str1);
        T = stod(str1);

        getline(out, str2);
        getline(out, str1);
        dt = stod(str1);

        getline(out, str2);
        getline(out, str1);
        tau = stod(str1);

        getline(out, str2);
        getline(out, str1);
        eps = stod(str1);

        getline(out, str2);
        getline(out, str1);
        sigma = stod(str1);

        getline(out, str2);
        getline(out, str1);
        k = stod(str1);

        getline(out, str2);
        getline(out, str1);
        rmin = stod(str1);

        getline(out, str2);
        getline(out, str1);
        v_mean = stod(str1);

        getline(out, str2);
        getline(out, str1);
        v_std = stod(str1);

        getline(out, str2);
        getline(out, str1);
        save = stoi(str1);

        getline(out, str2);
        getline(out, str1);
        a = stoi(str1);
    }
    out.close();

    std::cout << N << " particles" << '\n';
    std::cout << T << " model time" << '\n';

    std::ofstream file("Data/Parameters_"+std::to_string(a)+".txt");
    if(file.is_open()){
        file << N <<' '<< x_size <<' '<< y_size <<' '<< z_size <<' '<<
            T <<' '<< dt <<' '<< tau <<' '<< eps <<' '<< sigma <<' '<<
            k <<' '<< rmin <<' '<< v_mean <<' '<< v_std <<'\n';
    }
    file.close();

    Space<double> s(N, x_size, y_size, z_size, eps, sigma, k);

    //s.set_random_points(rmin);
    s.set_crystal_cell();
    s.set_random_speed(v_mean, v_std);

    //s.load_points("Data/Points_data.txt");
    //s.load_speed("Data/Speed_data.txt");

    s.run(T, dt, tau, save, std::to_string(a));

    ChangeStringInFileC(path, CountLinesInFile(path)-1, std::to_string(a+1));
}


int main()
{
    std::string path = "Params.txt";
    process<double>(path);
    std::cout << "Finished\n";
    return 0;
}

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
  char** lines;
  int n = CountLinesInFile(filename); // получить количество строк в файле
  if (n == -1) return -1;
  std::ifstream F(filename); // открыть файл для чтения
  if (!F) return -1;
  try
  {
    lines = new char* [n];
  }
  catch (std::bad_alloc e)
  {
    std::cout << e.what() << std::endl; // вывести сообщение об ошибке
    F.close(); // закрыть файл
    return -1;
  }
  int len; // длина одной строки
  char buffer[1000]; // память, куда записывается одна строка из файла
  for (int i = 0; i < n; i++)
  {
    F.getline(buffer, 1000);
    for (len = 0; buffer[len] != '\0'; len++);
    lines[i] = new char[len + 1];
    for (int j = 0; j < len; j++)
      lines[i][j] = buffer[j];
    lines[i][len] = '\0'; // добавить символ конца строки
  }
  F.close();
  *_lines = lines;
  return n;
}

bool ChangeStringInFileC(std::string filename, int position, std::string str)
{
  char** lines; // список строк файла
  int count; // количество строк файла
  count = GetStringsFromFileC(filename, &lines); // получить список lines
  if (count < 0) return false;
  if ((position < 0) || (position >= count)) return false;
  std::ofstream F(filename); // открыть файл для записи
  if (!F.is_open()) return false;
  if (position < count-1){
      for (int i = 0; i < position; i++)
        F << lines[i] << std::endl; // вывести строку в файл
      F << str.c_str() << std::endl; // здесь пишется строка str
      for (int i = position + 1; i < count - 1; i++)
        F << lines[i] << std::endl;
      F << lines[count - 1];
  }
  else{
      for (int i = 0; i < count - 1; i++)
        F << lines[i] << std::endl;
      F << str.c_str();
  }
  F.close();
  for (int i = 0; i < count; i++)
    delete lines[i];
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

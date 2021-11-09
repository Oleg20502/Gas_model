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
  // 1. �������������� ����������
  char** lines;
  int n = CountLinesInFile(filename); // �������� ���������� ����� � �����

  // 2. ��������, �������� �� ���������� �����
  if (n == -1) return -1;

  // 3. �������� �������� ����������
  std::ifstream F(filename); // ������� ���� ��� ������

  // 4. ��������, �������� �� ����
  if (!F) return -1;

  // 5. ������� �������� ������ ��� n �����
  try
  {
    lines = new char* [n];
  }
  catch (std::bad_alloc e)
  {
    // ���� ���������� �������� ������, �� �����
    std::cout << e.what() << std::endl; // ������� ��������� �� ������
    F.close(); // ������� ����
    return -1;
  }

  // 6. ���� ������ ����� � ��������� ������ ��� ���
  int len; // ����� ����� ������
  char buffer[1000]; // ������, ���� ������������ ���� ������ �� �����

  for (int i = 0; i < n; i++)
  {
    // 6.1. ������� ������ �� �����
    F.getline(buffer, 1000);

    // 6.2. ���������� ����� ����������� ������
    for (len = 0; buffer[len] != '\0'; len++);

    // 6.3. �������� ������ ��� len+1 ��������
    lines[i] = new char[len + 1];

    // 6.4. ����������� ������ �� buffer � lines[i]
    for (int j = 0; j < len; j++)
      lines[i][j] = buffer[j];
    lines[i][len] = '\0'; // �������� ������ ����� ������
  }

  // 7. ������� ����
  F.close();

  // 8. �������� ���������
  *_lines = lines;
  return n;
}

bool ChangeStringInFileC(std::string filename, int position, std::string str)
{
  // 1. �������� ������ ����� � ���� ������
  char** lines; // ������ ����� �����
  int count; // ���������� ����� �����
  count = GetStringsFromFileC(filename, &lines); // �������� ������ lines

  // 2. ��������, ��������� �� �������� ����
  if (count < 0) return false;

  // 3. ��������, ��������� �� ������� 0 <= position < count
  if ((position < 0) || (position >= count)) return false;

  // 4. ������ ����� lines � ���� �� ������� position
  std::ofstream F(filename); // ������� ���� ��� ������

  // 5. ��������, �������� �� ���� ��������� - ������� is_open()
  if (!F.is_open()) return false;

  if (position < count-1){
      for (int i = 0; i < position; i++)
        F << lines[i] << std::endl; // ������� ������ � ����

      // 6. ������ ������ � ������� position
      F << str.c_str() << std::endl; // ����� ������� ������ str

      // 7. ������ ����� ����� ������� position
      for (int i = position + 1; i < count - 1; i++)
        F << lines[i] << std::endl;

      // 8. �������� ��������� ������ ��� ������� '\n'
      F << lines[count - 1];
  }
  else{
      for (int i = 0; i < count - 1; i++)
        F << lines[i] << std::endl;

      // 8. �������� ��������� ������ ��� ������� '\n'
      F << str.c_str();
  }

  // 9. ������� ����
  F.close();

  // 10. ���������� ������, ���������� ��� �����
  for (int i = 0; i < count; i++)
    delete lines[i];

  // 11. ���������� ��������� �� ������
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

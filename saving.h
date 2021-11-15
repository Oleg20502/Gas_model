#ifndef saving_h
#define saving_h

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
  int n = CountLinesInFile(filename); // �������� ���������� ����� � �����
  if (n == -1) return -1;
  std::ifstream F(filename); // ������� ���� ��� ������
  if (!F) return -1;
  try
  {
    lines = new char* [n];
  }
  catch (std::bad_alloc e)
  {
    std::cout << e.what() << std::endl; // ������� ��������� �� ������
    F.close(); // ������� ����
    return -1;
  }
  int len; // ����� ����� ������
  char buffer[1000]; // ������, ���� ������������ ���� ������ �� �����
  for (int i = 0; i < n; i++)
  {
    F.getline(buffer, 1000);
    for (len = 0; buffer[len] != '\0'; len++);
    lines[i] = new char[len + 1];
    for (int j = 0; j < len; j++)
      lines[i][j] = buffer[j];
    lines[i][len] = '\0'; // �������� ������ ����� ������
  }
  F.close();
  *_lines = lines;
  return n;
}

bool ChangeStringInFileC(std::string filename, int position, std::string str)
{
  char** lines; // ������ ����� �����
  int count; // ���������� ����� �����
  count = GetStringsFromFileC(filename, &lines); // �������� ������ lines
  if (count < 0) return false;
  if ((position < 0) || (position >= count)) return false;
  std::ofstream F(filename); // ������� ���� ��� ������
  if (!F.is_open()) return false;
  if (position < count-1){
      for (int i = 0; i < position; i++)
        F << lines[i] << std::endl; // ������� ������ � ����
      F << str.c_str() << std::endl; // ����� ������� ������ str
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

#endif // saving_h

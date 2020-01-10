using namespace std;

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

int bad(char c)
{
    return c != 'A' && c != 'C' && c != 'G' && c != 'T';
}

int main(int argc, char **argv)
{
  ifstream input(argv[1]);
  string cur;
  ostringstream out("");
  cout << "Reading reference genome" << endl;
  while (getline(input, cur))
  {
      if(cur[0] != '>')
      {
          for(int i = 0; i<cur.length(); i++)
          {
              if(cur[i] >= 'a' && cur[i] <= 'z') cur[i] += 'A' - 'a';
              if(!bad(cur[i]))
                 out << cur[i];
          }
      }
  }
  string reference = out.str();
  size_t n = reference.length();
  cout << n << endl;

  FILE *refoutfile = fopen ((string(argv[1]) + string(".ref")).c_str(), "wb");
  fwrite(&reference[0], sizeof(char), n, refoutfile);
}

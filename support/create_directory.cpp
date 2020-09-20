#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>

using namespace std;

struct stat st = {0};

main()
{
  string dir = "some1";
  if (stat(dir.c_str(), &st) == -1) {
      mkdir(dir.c_str(), 0700);
  }
}

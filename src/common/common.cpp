#include "global.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <errno.h>
#include <cstring>

using namespace std;

#ifdef _WIN32
#include <windows.h>
#else
#define MAX_PATH 4096
#endif

double RandomDouble(double min, double max)
{
    return min + double(rand())/(double(RAND_MAX/(max - min)));
}

void RenameConsole(const string &title)
{
//	LPCWSTR t = title.c_str();
//    SetConsoleTitleA(title.c_str());
}

string CreateUniqueFileName(const string &filename, const string &ext)
{
    string name = filename + ext;

    for (int i = 1; ifstream(name).is_open(); ++i)
    {
        name = filename + '(' + to_string(i) + ')' + ext;
    }

    return name;
}

string CreateFolder(string &name)
{
#ifdef _WIN32
    char curDir[MAX_PATH] = "";
    char newDir[MAX_PATH];

    if (!GetCurrentDirectoryA(MAX_PATH, curDir))
    {
        cerr << "Error getting current directory: #" << GetLastError();
    }

    strcat(curDir, "\\");
    strcpy(newDir, curDir);
    strcat(newDir, name.c_str());

    char dirN[MAX_PATH];
    strcpy(dirN, newDir);
    string num = "";

    for (int i = 1; !CreateDirectoryA(dirN, NULL); ++i)
    {
        num = "(" + to_string(i) + ")";
        string dirName = newDir + num;
        strcpy(dirN, dirName.c_str());
    }

    name += num;
    return curDir;
#else
    // Extract basename from path (e.g. "results/run_123" -> "run_123")
    string dirPath = name;
    string baseName = name;
    size_t slash = name.rfind('/');
    if (slash != string::npos)
        baseName = name.substr(slash + 1);

    // Create directory, append (N) if already exists
    string num = "";
    int rc = mkdir(dirPath.c_str(), 0755);
    if (rc != 0 && errno == EEXIST)
    {
        for (int i = 1; ; ++i)
        {
            num = "(" + to_string(i) + ")";
            dirPath = name + num;
            rc = mkdir(dirPath.c_str(), 0755);
            if (rc == 0 || errno != EEXIST) break;
        }
    }
    if (rc != 0)
        cerr << "WARNING: could not create directory '" << dirPath
             << "': " << strerror(errno) << endl;
    name = baseName + num;
    return dirPath + "/";
#endif
}

string CreateDir(const string &name)
{
    string dirName;
#ifdef _WIN32
    char cdir[MAX_PATH];
    cdir[0] = '\0';
    strcat(cdir, name.c_str());

    char dirN[MAX_PATH];
    strcpy(dirN, cdir);

    for (int i = 1; !CreateDirectoryA(cdir, NULL); ++i)
    {
        string dirName = dirN;
        dirName += "(" + to_string(i) + ")";
        strcpy(cdir, dirName.c_str());
    }

    strcat(cdir, "\\");
    dirName = cdir;
#else
    dirName = name;
#endif
    return dirName;
}

void EraseConsoleLine(int lenght)
{
    cout << '\r';

    for (int i = 0; i < lenght; ++i)
    {
        cout << ' ';
    }

    cout << '\r';
}

double DegToRad(double deg)
{
    return (deg*M_PI)/180;
}

double RadToDeg(double rad)
{
    return (rad*180)/M_PI;
}

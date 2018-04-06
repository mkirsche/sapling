#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

std::string replaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
    return str;
}

string rep(string s)
{
    string res = s;
    res = replaceAll(res, "N", "");
    res = replaceAll(res, "M", "");
    res = replaceAll(res, "R", "");
    res = replaceAll(res, "W", "");
    res = replaceAll(res, "X", "");
    res = replaceAll(res, "Y", "");
    res = replaceAll(res, "Z", "");
    res = replaceAll(res, "J", "");
    res = replaceAll(res, "K", "");
    res = replaceAll(res, "L", "");
    return res;
}

#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <string>

#include "zstr/zstr.hpp"

char* GetOpt(char **begin, char **end, const std::string &option) {
  char **it = std::find(begin, end, option);
  return ((it != end && ++it != end) ? *it : 0);
}

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

std::unique_ptr<std::istream> OpenInstream(char **begin, char **end, const std::string &option) {
  const std::string &filename = GetOpt(begin, end, option);
  return std::unique_ptr<zstr::ifstream>(new zstr::ifstream(filename));
}

#endif

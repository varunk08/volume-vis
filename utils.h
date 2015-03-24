#ifndef UTILS_H
#define UTILS_H
namespace utils
{
template<typename T>
T min(T a, T b)
{
  return a < b?a : b;
}

template <typename T>
T max(T a, T b)
{
  return a > b? a:b;
}
}
#endif

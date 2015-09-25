#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <stdexcept>

//------------------------------------------------------------------------------
//
// d2s<T> - convert data of type T to string
//
//------------------------------------------------------------------------------
/**
 * Convert the data of any type which has oveloaded operator '<<' to string
 * @param data - the data
 * @param scientific - use scientific (exponential) format or not
 * @param precision - if scientific format is used, we can change the precision
 * @return data in string format
 */
template <typename T>
inline std::string d2s(T data,
                       bool scientific = false,
                       int precision = 6,
                       bool noperiod = false)
{
  std::ostringstream o;
  if (scientific)
  {
    o.setf(std::ios::scientific);
    o.precision(precision);
  }

  if (!(o << data))
    throw std::runtime_error("Bad conversion of data to string!");

  if (noperiod) // eliminate a period in case of 'double' numbers
  {
    std::string res = o.str(); // resulting string
    std::string::size_type pos = res.find('.');
    if (pos != std::string::npos)
      res.erase(pos, 1);
    return res;
  }

  return o.str();
}

#endif // UTILITIES_HPP

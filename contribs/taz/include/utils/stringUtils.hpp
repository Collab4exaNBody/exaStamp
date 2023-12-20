/// @file 
/// @brief Some elements for string manipulation (used in input)

#ifndef __STRING_UTILS_HPP_INCLUDED
#define __STRING_UTILS_HPP_INCLUDED


#include <algorithm>
#include <iostream>
#include <string>

#include "utils/array/array.hpp"
#include "utils/vec3/vec3.hpp"


/// @brief Compare two strings
/// @param [in] str1 First string
/// @param [in] str2 Second string
inline bool is_equal(const std::string& str1, const std::string& str2) {
  return str1.compare(str2)==0;
}


/// @brief Remove spaces in front and at the end of a given string
/// @param [in,out] line Line to trim
/// @param [in] otherSpaces Other characters counting as spaces concatenated in a stream, default=none
inline void trim_spaces(std::string& line, const std::string otherSpaces="") {
  auto position = line.find_first_not_of(" "+otherSpaces);
  if (position!=std::string::npos) line.erase(0, position);
  position = line.find_last_not_of(" "+otherSpaces);
  if (position!=std::string::npos) line.erase(position+1, line.size());
}


/// @brief Remove string content when a given character comment is found
/// @param [in] line Line to trim
/// @param [in] comment Character indicating the begining of the comments
inline void trim_comment(std::string& line, const char& comment) {
  auto position = line.find_first_of(comment);
  if (position!=std::string::npos) line.erase(position, line.size()-position);
}


/// @brief Get a line without spaces at the beginning and end from the stream
/// @param [in,out] in Source stream
/// @param [in,out] line Extracted line
/// @param [in] otherSpaces Other characters counting as spaces concatenated in a stream, default=tabulation
inline void getTrimedLine(std::istream& in, std::string& line, const std::string otherSpaces="\t") {
	getline(in,line);
	trim_spaces(line,otherSpaces);
}


/// @brief Split a line containing the character '=' into strings which 
/// contains the field to be assigned and its value
/// @param [in] line The line
/// @param [out] field Name of the field
/// @param [in] separator Character '=' or another separation character
/// @param [out] value Value to be assigned
inline void extract_content(const std::string& line, std::string& field, char separator, std::string& value) {

  auto position = line.find(separator);

  field = line.substr(0, position);
  trim_spaces(field, "\t");

  value = line.substr(position+1, line.size());
  trim_spaces(value, "\t");

}


/// @brief Convert a string value into type
/// @tparam T Type to convert to
/// @param [in] field Name of the field where the converted value will go (used for error print)
/// @param [in] str String to convert
/// @return Converted value
template <class T> T from_string(const std::string& field,const std::string& str);


/// @version Specialization for a bool
template <> inline bool from_string(const std::string& field,const std::string& str) {
  if      (is_equal(str, "true"))  return true;
  else if (is_equal(str, "false")) return false;
  else if (is_equal(str, "yes"))   return true;
  else if (is_equal(str, "no"))    return false;
  else if (is_equal(str, "1"))     return true;
  else if (is_equal(str, "0"))     return false;
  else return false;
}


/// @brief Convert a formated string into an array of data
/// @tparam T Type of the element of the array
/// @param [in] field Name of the field where the converted values will go (used for error print)
/// @param [in] str_ String to convert
/// @param [in] delimiter_left Left delimiter for an array
/// @param [in] delimiter_right Right delimiter for an array
/// @param [in] separator Element separator
/// @return Converted array
template <class T> Array<T> from_string_array(const std::string& field, const std::string& str_, char delimiter_left, char delimiter_right, char separator) {

  std::string str = str_;

  // Reduce string content to what is inside the delimiters
  auto position = str.find_first_of(delimiter_left);
  if (position!=std::string::npos) str.erase(0, 1+position);
  position = str.find_last_of(delimiter_right);
  if (position!=std::string::npos) str.erase(position, str.size());

  // Count number of elements in array
  int size = 1+std::count(str.begin(), str.end(), separator);

  // Exit if the array is empty
  if (size==1 && str.size()==0) return Array<T>(0);

  // Alloc array
  Array<T> array(size);
  // For each element
  for (int i=0; i<size; ++i) {
  	// Get the string until next separator
    position = str.find_first_of(separator);
    // Convert to the right type
    array[i] = from_string<T>(field, str.substr(0, position));
    // Delete this part
    str.erase(0, position+1);
    trim_spaces(str);
  }

  return array;

}


/// @brief Convert a formated string into a 3D vector of data
/// @tparam T Type of the element of the array
/// @param [in] field Name of the field where the converted values will go (used for error print)
/// @param [in] str_ String to convert
/// @param [in] delimiter_left Left delimiter for an array
/// @param [in] delimiter_right Right delimiter for an array
/// @param [in] separator Element separator
/// @return Converted array
template <class T>
vec3<T> from_string_vec3(const std::string& field, const std::string& str_, char delimiter_left, char delimiter_right, char separator) {

	// Get an array from the string
  Array<T> array = from_string_array<T>(field, str_, delimiter_left, delimiter_right, separator);
  // Check that the size of the array is 3 (if not abort)
  if (array.size()!=VEC3_NDIMS) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'from_string_vec3(const std::string& field, const std::string&, char, char, char)' : wrong size for a vec3<>. STOP."
	     << std::endl;
    exit(-1);
  }

  // Put the values into a 3D vector
  return vec3<T>(array[0], array[1], array[2]);

}


/// @brief Interpret a string into a floating point number
/// @tparam T Type of the floating point number
/// @param [in] str String to interpret
/// @param [out] idx Pointer to the number of character read, default=null
/// @return Interpreted value
template <class T> T floating_point_from_string(const std::string& str, size_t* idx=0);


/// @version Specialization for a double
template <> inline double floating_point_from_string(const std::string& str, size_t* idx) {
  return std::stod(str, idx);
}


/// @version Specialization for a float
template <> inline float floating_point_from_string(const std::string& str, size_t* idx) {
  return std::stof(str, idx);
}


/// @version Specialization for a long double
template <> inline long double floating_point_from_string(const std::string& str, size_t* idx) {
  return std::stold(str, idx);
}


/// @brief Interpret a string into an integer number
/// @tparam T Type of the integer number
/// @param [in] str String to interpret
/// @param [out] idx Pointer to the number of character read, default=null
/// @param [in] base Base of the number, default=10
/// @return Interpreted value
template <class T> T integer_from_string(const std::string& str, size_t* idx=0, int base=10);


/// @version Specialization for an int
template <> inline int integer_from_string(const std::string& str, size_t* idx, int base) {
  return std::stoi(str, idx, base);
}


/// @version Specialization for a long
template <> inline long integer_from_string(const std::string& str, size_t* idx, int base) {
  return std::stol(str, idx, base);
}


/// @version Specialization for a long long
template <> inline long long integer_from_string(const std::string& str, size_t* idx, int base) {
  return std::stoll(str, idx, base);
}


/// @version Specialization for an unsigned long
template <> inline unsigned long integer_from_string(const std::string& str, size_t* idx, int base) {
  return std::stoul(str, idx, base);
}


/// @version Specialization for an unsigned long long
template <> inline unsigned long long integer_from_string(const std::string& str, size_t* idx, int base) {
  return std::stoull(str, idx, base);
}


/// @version Specialization for a double
template <> inline double from_string(const std::string& field, const std::string& str) {
  return floating_point_from_string<double>(str);
}


/// @version Specialization for a float
template <> inline float from_string(const std::string& field, const std::string& str) {
  return floating_point_from_string<float>(str);
}


/// @version Specialization for a long double
template <> inline long double from_string(const std::string& field, const std::string& str) {
  return floating_point_from_string<long double>(str);
}


/// @version Specialization for an integer
template <> inline int from_string(const std::string& field, const std::string& str) {
  return integer_from_string<int>(str);
}


/// @version Specialization for a long integer
template <> inline long from_string(const std::string& field, const std::string& str) {
  return integer_from_string<long>(str);
}


/// @version Specialization for a long long integer
template <> inline long long from_string(const std::string& field, const std::string& str) {
  return integer_from_string<long long>(str);
}


/// @version Specialization for a unsigned long integer
template <> inline unsigned long from_string(const std::string& field, const std::string& str) {
  return integer_from_string<unsigned long>(str);
}


/// @version Specialization for a unsigned long long integer
template <> inline unsigned long long from_string(const std::string& field, const std::string& str) {
  return integer_from_string<unsigned long long>(str);
}


/// @version Specialization for a unsigned string
template <> inline std::string from_string(const std::string& field, const std::string& str) {
  std::string str_ = str;
  trim_spaces(str_);
  return str_;
}

#endif // __STRING_UTILS_HPP_INCLUDED

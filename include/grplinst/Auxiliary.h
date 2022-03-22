#ifndef GRPLINST_AUXILIARY_H
#define GRPLINST_AUXILIARY_H

#include <cassert>
#include <cstdarg>
#include <string>
#include <vector>




/**
 */
inline const std::string formatString(const char* fmt, ...) {
	// arguments
	va_list listOfArgs;
	va_start(listOfArgs, fmt);

	// dummy copy of args to avoid future memory problems
	va_list listOfArgs2;
    va_copy(listOfArgs2, listOfArgs);
	int length = std::vsnprintf(nullptr, 0, fmt, listOfArgs2);
    va_end(listOfArgs2);
	
	std::vector<char> s(length + 1);
	std::vsnprintf(s.data(), s.size(), fmt, listOfArgs);
	va_end(listOfArgs);

	std::string formattedString = std::string(s.data(), length);

	return formattedString;
}


#endif
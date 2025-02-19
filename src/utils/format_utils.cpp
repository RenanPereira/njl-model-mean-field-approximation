#include "utils/format_utils.h"

using namespace std;


string trim0ToDot0(const double value) 
{
    string result = to_string(value);

    // Remove unnecessary trailing zeroes
    result.erase(result.find_last_not_of('0') + 1);

    // Ensure at least one digit remains after the decimal point (e.g., 0. -> 0.0)
    if (result.back() == '.') 
    {
        result += '0';
    }

    return result;
}

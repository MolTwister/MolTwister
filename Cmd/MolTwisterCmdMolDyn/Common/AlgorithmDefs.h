#pragma once

#include <algorithm>
#include <vector>
template<class InputIt, class OutputIt, class UnaryOperation>
OutputIt mttransform(InputIt first1, InputIt last1, OutputIt d_first, UnaryOperation unary_op )
{
    return std::transform(first1, last1, d_first, unary_op);
}
template<class T> using hostvector = std::vector<T>;
template<class T> using devicevector = std::vector<T>;

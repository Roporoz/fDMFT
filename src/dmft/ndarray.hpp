#pragma once

#include <vector>

// Replaces legacy headers.new.h's new_double3/new_complex3/new_int3,
// new_int4/new_complex4 and their matching delete_*() calls. Those were
// hand-rolled T***/T**** allocators, each requiring the caller to remember
// the exact (n1, n2, n3[, n4]) triple/quad again at the matching delete_*()
// call — get it wrong (or forget the call, as the original does for
// new_int4, which has no delete_int4 counterpart at all) and it's a crash or
// a leak. Plain nested std::vector gives the same "array of arrays" shape
// with automatic, correct cleanup, so no delete_* equivalents are needed
// here at all.
namespace dmft {

template <typename T>
using Array3 = std::vector<std::vector<std::vector<T>>>;

template <typename T>
using Array4 = std::vector<std::vector<std::vector<std::vector<T>>>>;

template <typename T>
Array3<T> make_array3(int n1, int n2, int n3, T fill = T{}) {
    return Array3<T>(n1, std::vector<std::vector<T>>(n2, std::vector<T>(n3, fill)));
}

template <typename T>
Array4<T> make_array4(int n1, int n2, int n3, int n4, T fill = T{}) {
    return Array4<T>(
        n1, std::vector<std::vector<std::vector<T>>>(
                n2, std::vector<std::vector<T>>(n3, std::vector<T>(n4, fill))));
}

}  // namespace dmft

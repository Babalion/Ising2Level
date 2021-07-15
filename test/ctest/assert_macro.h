//
// Created by chris on 26.04.21.
//

#ifndef NUMERISCHEMETHODENSTATISTISCHENPHYSIK_ASSERT_MACRO_H
#define NUMERISCHEMETHODENSTATISTISCHENPHYSIK_ASSERT_MACRO_H
#include <iostream>
#include <sstream>


#define assertEqual( ... )               \
do {                                            \
    if( !( __VA_ARGS__ ) ) {                     \
        std::cerr << "Unit test assert [ " \
        << ( #__VA_ARGS__ )             \
        << " ] failed in line [ "       \
        << __LINE__                     \
        << " ] file [ "                 \
        << __FILE__ << " ]"             \
        << std::endl;                     \
        err_code = 1;                           \
    }                                            \
} while( false )
#endif //NUMERISCHEMETHODENSTATISTISCHENPHYSIK_ASSERT_MACRO_H

// hex_index_zorder.h
// Z-Order Curve Indexing for Apertures 3, 4, 7
// Based on DGGRID's DgZOrderStringRF.cpp

#pragma once

#include <string>

namespace hexify {
namespace zorder {

// Z-Order for Aperture 3
std::string encode_ap3(long long i, long long j, int resolution);
void decode_ap3(const std::string& z_str, int resolution, 
                long long& i, long long& j);

// Z-Order for Aperture 4
std::string encode_ap4(long long i, long long j, int resolution);
void decode_ap4(const std::string& z_str, long long& i, long long& j);

// Z-Order for Aperture 7
std::string encode_ap7(long long i, long long j, int resolution);
void decode_ap7(const std::string& z_str, int resolution, 
                long long& i, long long& j);

} // namespace zorder
} // namespace hexify

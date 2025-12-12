#pragma once
// cell_index.h - Unified space-filling curve indexing for all apertures

#include <string>
#include <vector>
#include <cstdint>

namespace hexify {

enum class IndexType {
  AUTO,
  ZORDER,
  Z3,
  Z7
};

std::string cell_to_index(int face, long long i, long long j, 
                          int resolution, int aperture,
                          IndexType index_type = IndexType::AUTO);

void index_to_cell(const std::string& index, int aperture,
                   IndexType index_type,
                   int& face, long long& i, long long& j, int& resolution);

std::string cell_to_index_ap34(int face, long long i, long long j,
                               const std::vector<int>& ap_seq);

void index_to_cell_ap34(const std::string& index,
                        const std::vector<int>& ap_seq,
                        int& face, long long& i, long long& j);

uint64_t index_to_uint64(const std::string& index, int aperture,
                         IndexType index_type);

std::string uint64_to_index(uint64_t value, int resolution, int aperture,
                            IndexType index_type);

std::string get_parent_index(const std::string& index, int aperture,
                             IndexType index_type);

std::vector<std::string> get_children_indices(const std::string& index, 
                                              int aperture,
                                              IndexType index_type);

int compare_indices(const std::string& idx1, const std::string& idx2);

int get_index_resolution(const std::string& index, int aperture,
                         IndexType index_type);

bool is_valid_index_type(int aperture, IndexType index_type);

IndexType get_default_index_type(int aperture);

} // namespace hexify

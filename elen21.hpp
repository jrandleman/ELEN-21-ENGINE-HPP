// Author: Jordan Randleman - elen21.hpp
// About: SOP-POS Function/Kmap/TruthTable Generator.

#ifndef ELEN_21_HPP_
#define ELEN_21_HPP_

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#define K_ERR_HEADER "\x1b[1m" + std::string(__FILE__)    + ":" +\
                                 std::string(__func__)    + ":" +\
                                 std::to_string(__LINE__) + ":" +\
                                 "\x1b[31mERROR:\x1b[0m\x1b[1m INVALID "

using Bit = unsigned long long;
using Bits = std::vector<Bit>;
using Bit_Matrix = std::vector<Bits>;

const Bit DONT_CARE_BIT = 2;
const Bit BLANK_BIT = 3;

/******************************************************
* OPERATION * MINIMUM SOP-POS NOTATION * C++ NOTATION *
*******************************************************
*   NOT     *            A'            *     !A       *
*           *                          *              *
*   AND     *          (A.B)           *   (A & B)    *
*   OR      *          (A+B)           *   (A | B)    *
*   XOR     *          (A^B)           *   (A ^ B)    *
*           *                          *              *
*   NAND    *          (A.B)'          *   !(A & B)   *
*   NOR     *          (A+B)'          *   !(A | B)   *
*   XNOR    *          (A'^B)          *   !(A ^ B)   *
******************************************************/

/*
 * 16 FUNCTIONS FOR USERS:
 *     1) SOP_str: returns SOP string for fcn, given string vector w/ var names 
 *                 & "Bits" of fcn results as if from a truth table
 *     2) SOP_str: returns SOP string for fcn, given string vector w/ var names 
 *                 & a fcn to apply to a generated truth table
 *     3) POS_str: returns POS string for fcn, given string vector w/ var names 
 *                 & "Bits" of fcn results as if from a truth table
 *     4) POS_str: returns POS string for fcn, given string vector w/ var names 
 *                 & a fcn to apply to a generated truth table
 *     => NOTE: STRING "VARIABLE NAME" VECTOR ARGS FOR THE ABOVE 4 FCNS SHALL
 *                 BE LISTED IN DESCENDING BIT SIGNIFICANCE
 *
 *     5) xor_replace:    returns string w/ fcn's expanded XOR/XNOR's abbreviated
 *     6) compl_literals: returns string w/ fcn's literals complemented
 *     7) dual_fcn:       returns string w/ fcn's dual
 *     8) compl_fcn:      composes "compl_literals" & "dual_fcn"
 *
 *     9) BCD_decoder: convert vector of bits into decimal # (helps w/ dont cares)
 *
 *    10) Kmap_print: given bit#
 *    11) Kmap_print: given "Bits" of fcn results (as if from a truth table)
 *    12) Kmap_print: given bit# & fcn to apply across generated truthtable
 *
 *    13) TruthTable_print: given bit#
 *    14) TruthTable_print: given bit# & vector of fcns to apply to table
 *    15) TruthTable_print: given vector of "Bits" (results of truth table fcns)
 *
 *    16) TruthTable_fcn: given bit# & a fcn, returns vector of fcn results 
 *                        from truth table
 */

/*
 * GENERAL CUSTOM FCN LAYOUT TO BE PASSED TO SOP/POS & TRUTHTABLE/KMAP FCNS:
 *   => !!! NOTE the mandatory "Bit" return & "Bits" arg types !!!
 *   => the "Bits" arg shall be a vector sorted by descending bit significance,
 *      contining a "row" of truth table bits as its arg
 *
 * Bit your_func_name_here(Bits bit_vector) { 
 *    Bit dont_care_row = 3; // if we don't care about row 3 (rows starting from 0)
 *    
 *    // BCD_decoder returns decimal value of given bit vector
 *    if(BCD_decoder(bit_vector) == dont_care_row) 
 *      return DONT_CARE_BIT; // use predefined "DONT_CARE_BIT" to denote don't cares
 *
 *    return (bit_vector[0] & !bit_vector[1]) ^ bit_vector[3]; // any operation
 * }
 */


// These 3 function templates allow us to use sorted, unqiue-elt 
// "std::vector"s: avoids overhead of std::set what w/ its typical red-black
// tree implementation, & this program's use of iteration FAR more frequently
// than insertion/erasure means we ought to simulate std::set w/ std::vector.

// Check out this PDF for more info on std::set vs std::vector:
// http://lafstern.org/matt/col1.pdf

template<typename T>
static void sorted_insert(std::vector<T>&v, const T elt) {
  // "lower_bound" performs binary search
  auto position = std::lower_bound(v.begin(), v.end(), elt);
  if(position == v.end() || elt < *position)
    v.insert(position, elt);
}

template<typename T>
static void sorted_erase(std::vector<T>&v, const T elt) {
  // "lower_bound" performs binary search
  auto position = std::lower_bound(v.begin(), v.end(), elt);
  if(position != v.end() && *position == elt)
    v.erase(position);
}

template<typename T>
static bool sorted_find(const std::vector<T>&v, const T elt) {
  // "lower_bound" performs binary search
  auto position = std::lower_bound(v.begin(), v.end(), elt);
  return (position != v.end() && *position == elt);
}

/******************************************************************************
* ARITHMETIC HELPER FUNCTIONS
******************************************************************************/

// returns if a is a perfect power of b => UNSIGNED INTS ONLY (ie no fractionals)
constexpr inline bool a_isPowOf_b(const Bit a, const Bit b) {
  return (a==b||a==1) ? true : ((a%b)||(a==0)) ? false : a_isPowOf_b(a/b, b);
}


// returns log2(n):
//   @err: returns max unsigned (ie -1 overflowed) if n != perfect power of 2
static const Bit BAD_LOG2_ARG = -1;
constexpr inline Bit log2(const Bit n) {
  if(!a_isPowOf_b(n, 2)) return BAD_LOG2_ARG;
  Bit counter = 1, result = 0;
  while(counter < n) counter <<= 1, ++result;
  return result;
}

/******************************************************************************
* HYPERCUBE ADJACENCY MATRIX FUNCTIONS
******************************************************************************/

/******************************************************************************
 * Labeled Points on Hypercube of Dimension D = 3 & Corresponding 3-Var K-Map *
 *                                                                            *
 *              2 ---------- 3                                                *
 *             /|           /|             ----- ----- ----- -----            *
 *            / |          / |            |  0  |  0  |  0  |  0  |           *
 *           0 ---------- 1  |            |00)  |01)  |03)  |02)  |           *
 *           |  |         |  |             ----- ----- ----- -----            *
 *           |  6 --------|- 7            |  0  |  0  |  0  |  0  |           *
 *           | /          | /             |04)  |05)  |07)  |06)  |           *
 *           |/           |/               ----- ----- ----- -----            *
 *           4 ---------- 5                                                   *
 *                                                                            *
 *****************************************************************************/

// Find all points connected to point P in hypercube of dimension D.
//   => Dimension D repns K-map of D vars, and adjacent points on hypercube
//      repn logically adjacent cells on K-map (thus P == K-map cell idx).
// How: recursively half dimensions to find the point connected to P
//      on the opposite side of the current dimension.
// Ie:  find 4th dim projection of P, then 3rd, then 2nd, and finally 1st
// Time & Space Complexity: O(D)
static void find_adjacent_hypercube_points(const Bit D, const Bit original_D,
                                           const Bit P, const Bit original_P,
                                           Bits &points_connected_to_P){
  if(!D || P > (1<<original_D)) return;
  const Bit half_of_dimension = (1<<(D-1));
  if(P >= half_of_dimension) {
    sorted_insert(points_connected_to_P, original_P-half_of_dimension);
    find_adjacent_hypercube_points(D-1, original_D, P-half_of_dimension, original_P, points_connected_to_P);
  } else {
    sorted_insert(points_connected_to_P, original_P+half_of_dimension);
    find_adjacent_hypercube_points(D-1, original_D, P, original_P, points_connected_to_P);
  }
}

/******************************************************************************
* NESTED HYPERCUBE CONTAINMENT FUNCTION
******************************************************************************/

// Confirms set "points" all fall w/in dimension D's range
// Time  Complexity: O(points.size())
// Space Complexity: O(1)
static bool points_are_within_D(const Bit D, const Bits &points){
  const Bit total_D_points = (1<<D);
  if(points.empty() || points.size()>total_D_points) return false;
  for(Bit pt : points) if(pt >= total_D_points) return false;
  return true;
}

/******************************************************************************
* NESTED POINT-OVERLAPPING HYPERCUBES HELPER FUNCTIONS
******************************************************************************/

// Returns a set of indices where "function_bit_values" has a "sought" value
// => "sought_bit" = 1 for SOP (minterm) && 0 for POS (maxterm)
static void convert_fcn_to_point_set(const Bits &function_bit_values, 
  const Bit sought_bit,Bits &need_to_cover, Bits &points) {
  const Bit n = function_bit_values.size();
  for(Bit i = 0; i < n; ++i)
    if(function_bit_values[i] == sought_bit)
      sorted_insert(points, i), sorted_insert(need_to_cover, i);
    else if(function_bit_values[i] == DONT_CARE_BIT)
      sorted_insert(points, i);
}


// Returns greatest hypercube possible w/ "total_vertices" # of vertices
static Bit max_formable_hypercube(Bit max_inner_hypercube_vertices, 
const Bit total_vertices) {
  while(max_inner_hypercube_vertices > total_vertices) 
    max_inner_hypercube_vertices >>= 1;
  return log2(max_inner_hypercube_vertices);
}

/******************************************************************************
* EPI PARSING FROM PI MATRIX (FIND FEWEST LARGEST PI'S TO COVER ALL NEEDED PTS)
******************************************************************************/

// Return whether NOT all rows in "PI_matrix" are empty
static bool PI_matrix_is_not_empty(const std::vector<Bit_Matrix>&PI_matrix){
  for(auto PI_sets : PI_matrix)
    if(!PI_sets.empty()) return true;
  return false;
}


// Find idx of smallest "values" list (ie the set of PI's containing the least 
// frequent desired elt)
static Bit idx_of_row_with_least_values_in_PI_matrix(
const std::vector<Bit_Matrix>&PI_matrix) {
  bool first_non_empty_row = true;
  Bit shortest_row_idx = 0, count = 0;
  for(auto PI_set : PI_matrix) {
    if(!PI_set.empty() && 
    (first_non_empty_row || PI_set.size()<PI_matrix[shortest_row_idx].size())){
      shortest_row_idx = count;
      first_non_empty_row = false;
    }
    ++count;
  }
  return shortest_row_idx;
}


// Rm excess PI's w/ points in the current "epi" from the "PI_matrix"
static void rm_redundant_PIs_from_PI_matrix(const Bits &epi,
std::vector<Bit_Matrix> &PI_matrix,
const Bits &desired_keys) {
  // erase rows for epi's points
  const Bit total_desired = desired_keys.size();
  for(auto desired_point_key : epi) {
    Bit i = 0;
    while(i < total_desired && desired_keys[i] != desired_point_key) ++i;
    if(i < total_desired) PI_matrix[i].clear();
  }
  
  // erase redundant rows for epi's points
  const Bit n = PI_matrix.size();
  Bit_Matrix redundant_PI;
  for(Bit i = 0; i < n; ++i) {
    // no redundant sets if only 1 or 0 sets remaining in a list of sets 
    // containing a desired point outside of the current "epi"
    if(PI_matrix[i].size() < 2) continue;

    Bit count = PI_matrix[i].size();
    redundant_PI.clear();

    for(auto PI : PI_matrix[i]) {
      if(count < 2) break; 
      for(auto point : epi)
        if(sorted_find(PI, point)) { 
          sorted_insert(redundant_PI, PI), --count;
          break;
        }
    }

    for(auto PI : redundant_PI) sorted_erase(PI_matrix[i], PI);
  }
}


// Rm all non-EPIs from "cube_matrix" of cube-PI's
static void rm_non_EPIs_from_PI_cube_matrix(Bit_Matrix &cube_matrix,
const Bits &need_to_cover) {
  Bits desired_keys;

  // Put intersection of "need_to_cover" & "cube_matrix" in "desired_keys".
  // These are the elts the EPI set must contain.
  for(auto cube : cube_matrix)
    for(auto point : cube)
      if(sorted_find(need_to_cover, point))
        sorted_insert(desired_keys, point);


  const Bit total_desired = desired_keys.size();

  // Idxs = "keys", where "key" = elt in "desired_keys".
  //   => IE: if desired_keys[0] = 7, 7 = key for PI_matrix[0]
  // Matrices of sets (matrices of PI's) = "values" list, where each set 
  //   (PI) in "values" contains its "key" as an elt.
  std::vector<Bit_Matrix>PI_matrix(total_desired,Bit_Matrix{});

  // Initialize PI_matrix according to the parameters above its definition
  for(Bit i = 0; i < total_desired; ++i)
    for(auto cube : cube_matrix)
      if(sorted_find(cube, desired_keys[i]))
        sorted_insert(PI_matrix[i], cube);

  cube_matrix.clear(); // will hold the set of EPI's

  // Parse EPIs from w/in PI_matrix
  Bits epi;
  while(PI_matrix_is_not_empty(PI_matrix)) {

    // Find idx of smallest "values" list (ie the set of PI's 
    // containing the least frequent desired elt)
    Bit shortest_row_idx = idx_of_row_with_least_values_in_PI_matrix(PI_matrix);

    // Selected an epi (could be any in "PI_matrix[shortest_row_idx]" set)
    epi = PI_matrix[shortest_row_idx][0];

    // Eliminate sets of PI's containing any points from the current "epi" 
    //   (given they're now covered), SO LONG as such a set is not the ONLY set
    //   left in a row's "values" list (otherwise would lose all sets
    //   containing the row's desired elt "key").
    // Furthermore, also eliminate all rows w/ points in "epi" as their key:
    //   since "epi" already covers the point, no need for a list of other sets 
    //   containing it.
    rm_redundant_PIs_from_PI_matrix(epi,PI_matrix,desired_keys);
    sorted_insert(cube_matrix, epi);
  }
}

/******************************************************************************
* PI DERIVATION (FINDING INNER HYPERCUBES COVERING NEEDED POINTS)
******************************************************************************/

// Clears row at idx "rm" & erases all elt instances of "rm" from "adj_matrix"
static void erase_row_and_all_elt_instances_of(const Bit rm,
Bit_Matrix &adj_matrix) {
  adj_matrix[rm].clear(); // erase row
  for(Bit i = 0, n = adj_matrix.size(); i < n; ++i) // erase elts
    sorted_erase(adj_matrix[i], rm);
}


// Removes all references & rows of points NOT in "need_to_cover" set from "adj_matrix"
static void rm_dont_needToCovers_from_adj_matrix(const Bit MAX_VERTICES,
const Bits &need_to_cover,Bit_Matrix &adj_matrix){
  // Get which bits we don't need to try covering
  std::vector<bool> dont_need_to_cover_hash_table(MAX_VERTICES, true);
  for(auto v : need_to_cover) dont_need_to_cover_hash_table[v] = false;
  // rm these "don't care bits" from the adjacency list matrix cpy
  // once complete, we have an adj. list matrix ONLY w/ "need_to_cover" bits
  for(Bit i = 0; i < MAX_VERTICES; ++i)
    if(dont_need_to_cover_hash_table[i])
      erase_row_and_all_elt_instances_of(i, adj_matrix);
}


// Returns whether "adj_matrix" has enough points to form an "nth" sub-dimension w/in
static bool has_dimension_n(const Bit_Matrix &adj_matrix,const Bit n){
  const Bit minimum_points_in_dimension_n = (1<<n);
  const Bit adj_matrix_size = adj_matrix.size();
  Bit total_points = 0;
  for(Bit i = 0; i < adj_matrix_size; ++i) {
    if(adj_matrix[i].size() >= n)
      ++total_points;
    if(total_points >= minimum_points_in_dimension_n) 
      return true;
  }
  return false;
}


// Returns set matrix of combinations size "n" w/in iterator range [begin,end)
static void get_n_combinations_in_range(const Bit n, 
const Bits::iterator begin, const Bits::iterator end, 
Bits seen, Bit_Matrix &combo_matrix) {
  Bits::iterator temp;
  for(Bits::iterator it = begin; it != end;) {
    sorted_insert(seen, *it);
    temp = it++;
    switch(n) {
      case 1:  sorted_insert(combo_matrix, seen); break;
      default: get_n_combinations_in_range(n-1, it, end, seen, combo_matrix);
    } 
    sorted_erase(seen, *temp);
  }
}


// Remove all elts in "seen" set from set-matrix "val_list"
static void rm_seen_elts_from_valList(const Bits &seen,Bit_Matrix &val_list){
  for(auto point : seen)
    for(Bits &adj_list : val_list)
      sorted_erase(adj_list, point);
}


// Refills "val_list" w/ adjacency lists from "adj_matrix" & push points from
// "pt_rows" to "seen" set
static void refill_valList_and_push_to_seen(Bits &seen,Bit_Matrix &val_list,
const Bits &pt_rows, const Bit_Matrix &adj_matrix) {
  for(auto point : pt_rows) {
    val_list.push_back(adj_matrix[point]);
    sorted_insert(seen, point);
  }
}


// If all rows in val_list share "edge" distinct elts w/ "n - edge" other rows,
// return set of these shared elts - else, return empty set
static void shared_edges_set_forming_sub_dimension(const Bit n,const Bit edge,
Bits &pt_rows,const Bit_Matrix &val_list) {
  pt_rows.clear();

  Bit this_row_idx = 0;
  for(auto this_row : val_list) {

    Bit elt_count = 0;
    for(auto elt : this_row) {

      Bit row_count = 0, row_idx = 0;
      for(auto row : val_list) {
        if(sorted_find(row, elt) && row_idx != this_row_idx)
          ++row_count;
        ++row_idx;
      }

      if(row_count == (n - edge)) {
        sorted_insert(pt_rows, elt);
        ++elt_count;
      }
    }

    if(elt_count != edge) {
      pt_rows.clear();
      return;
    }

    ++this_row_idx;
  }
}


// Get "nth" sub-dimensions w/ point "P" in "adj_matrix"
static void get_nth_sub_cubes_containing_point_P(const Bit P,const Bit n,
Bit_Matrix &cube_matrix, Bit_Matrix &adj_matrix){
  // Get "n"-length combinations of P's values to try and form cubes w/
  Bit_Matrix combos;
  get_n_combinations_in_range(n,adj_matrix[P].begin(),adj_matrix[P].end(),{},combos);

  Bits seen;
  Bit_Matrix val_list;

  // get all sub-dimensions of point P connected to different combinations of "edge" points in its adjacency list
  for(auto pt_rows : combos) { // "pt_rows" = the "idx" keys for "val_list"'s rows
    bool is_cube = true;
    seen = {P};
    val_list.clear();
    
    // get initial series of point adj lists & save the to set of "seen" points
    refill_valList_and_push_to_seen(seen,val_list,pt_rows,adj_matrix);

    // check whether point & "edge" points combo in its adjacency list forms an actual sub-dimension or not
    for(Bit edge = n-1; edge > 0; --edge) {
      rm_seen_elts_from_valList(seen,val_list);

      // if all rows in val_list share "edge" distinct elts w/ "n - edge" other rows, 
      // put set of these shared elts into "pt_rows"
      shared_edges_set_forming_sub_dimension(n, edge, pt_rows, val_list);

      if(!pt_rows.empty()) {
        val_list.clear();
        refill_valList_and_push_to_seen(seen,val_list,pt_rows,adj_matrix);
      } else {
        is_cube = false;
        break;
      }
    }

    if(is_cube) sorted_insert(cube_matrix, seen);
  }
}


// Confirms cube = subset of any higher-dimensional cubes in "cube_matrix"
static bool cube_subset_of_larger_cubes(const Bits &cube,
const Bit_Matrix &cube_matrix){
  for(auto hcube : cube_matrix) {
    if(hcube.size() <= cube.size()) 
      continue;
    bool hcube_contains_cube = true;
    for(auto point : cube)
      if(!sorted_find(hcube, point)) {
        hcube_contains_cube = false;
        break;
      }
    if(hcube_contains_cube) 
      return true; // cube IS a subset of higher dimension in "cube_matrix"
  }
  return false; // cube is NOT a subset of higher dimension in "cube_matrix"
}


// FIRST parses out PI's for dimension (sub-dimensions spanning "need_to_cover"
// points), THEN rms any non-essential PI's (thus only leaving EPIs)
void get_EPI_sub_dimension_matrix_for_D(const Bit D, const Bit max_inner_hypercube,
Bits &need_to_cover, const Bit_Matrix &full_adj_matrix,
Bit_Matrix &cube_matrix,const Bits &points) {
  // rm all references & rows of points NOT in "need_to_cover" set from "adj_matrix"
  Bit_Matrix adj_matrix(full_adj_matrix);
  rm_dont_needToCovers_from_adj_matrix((1<<D), points, adj_matrix);

  Bit_Matrix local_adj_matrix;
  Bit_Matrix local_cube_matrix;
  Bit_Matrix subset_cubes;

  // For each possible sub-dimension in "adj_matrix" of "need_to_cover" points
  for(Bit n = max_inner_hypercube; n >= 0 && !need_to_cover.empty(); --n) {

    // if "full_adj_matrix" has a sub-dimension size n
    if(has_dimension_n(adj_matrix,n)) {
      local_adj_matrix = adj_matrix;
      local_cube_matrix.clear();

      // Get sub-cubes for the current "nth" dimension in "adj_matrix"
      for(Bit i = 0; i < local_adj_matrix.size(); ++i)
        if(local_adj_matrix[i].size() >= n) {
          // Add nth dimensional cubes that can be formed including point "i"
          get_nth_sub_cubes_containing_point_P(i,n,local_cube_matrix,local_adj_matrix);
          // Rm instances of point "i" in "local_adj_matrix", given
          // we already have all the sub-hypercubes that can be formed w/ it
          erase_row_and_all_elt_instances_of(i, local_adj_matrix);
        }

      // Rm all local cubes that are subsets of a pre-existing larger cube
      subset_cubes.clear();
      for(auto cube : local_cube_matrix)
        if(cube_subset_of_larger_cubes(cube, cube_matrix)) 
          sorted_insert(subset_cubes, cube);
      for(auto cube : subset_cubes)
        sorted_erase(local_cube_matrix, cube);

      // Rm all non-essential PI's from "local_cube_matrix" 
      //   => HENCEFORTH "local_cube_matrix" SHALL ONLY CONTAIN "EPI"'S !!!
      rm_non_EPIs_from_PI_cube_matrix(local_cube_matrix,need_to_cover); 

      // push "local_cube_matrix" cubes to "cube_matrix", 
      // & rm their pts from "need_to_cover"
      for(auto cube : local_cube_matrix) {
        sorted_insert(cube_matrix, cube);
        for(auto point : cube) 
          sorted_erase(need_to_cover, point);
      }
    }

    // End of iteration: manual break since 'Bit' type is unsigned
    if(n == 0) break;
  }

  // insert the remaining points left over in "need_to_cover" set, which
  // (having no connections to other points in the adjacency list matrix)
  // are their own 0th dimensional cubes
  for(auto point : need_to_cover)
    sorted_insert(cube_matrix, Bits{point});
}

/******************************************************************************
* NESTED POINT-OVERLAPPING HYPERCUBES MAIN FUNCTION
******************************************************************************/

// Returns hypercube overlap matrix of vertices in set "points" w/in dimension D
// => IE the sets of vertices w/in "points" that are K-map logically adjacent!
Bit_Matrix max_hypercube_overlap_matrix(const Bit D, 
const Bits function_bit_values, const Bit sought_bit) {
  Bits need_to_cover, points;
  convert_fcn_to_point_set(function_bit_values,sought_bit,need_to_cover,points);
  if(D == 0 || !points_are_within_D(D,points)) 
    return Bit_Matrix{};

  // find greatest formable hypercube from set "points"
  Bit max_inner_hypercube = max_formable_hypercube((1<<D), points.size());
  if(max_inner_hypercube > D) return Bit_Matrix{};

  // adajcency matrix for points in dimension D 
  //   => ie: logical adjacency matrix for cells in K-map of D vars
  const Bit TOTAL_VERTICES = (1 << D);
  Bit_Matrix adj_matrix(TOTAL_VERTICES, Bits{});
  for(Bit P = 0; P < TOTAL_VERTICES; ++P) // generates adjacency matrix per vertex
    find_adjacent_hypercube_points(D,D,P,P, adj_matrix[P]); 

  // Each inner set = set of points forming a sub-dimension/inner-hypercube
  Bit_Matrix EPI_matrix;
  get_EPI_sub_dimension_matrix_for_D(D, max_inner_hypercube, need_to_cover, adj_matrix, EPI_matrix, points);

  return EPI_matrix;
}

/******************************************************************************
* SOP-POS FUNCTION GENERATION HELPER FUNCTIONS
******************************************************************************/

// Confirms valid SOP / POS variable names given
static void validate_sop_pos_vars(const std::vector<std::string>var_names){
  const std::string reserved("()+.^'");
  for(auto name : var_names)
    for(auto ch : reserved)
      if(name.find(ch) != std::string::npos) {
        std::cerr << K_ERR_HEADER << " SOP/POS VARIABLE NAME \"" << name 
                  << "\" CONTAINS RESERVED CHAR '" << ch << "'\x1b[0m\n"
                  << "=> NAMES SHALL NOT CONTAIN THE FOLLOWING: \"()+.^'\"\n" 
                  << "Terminating program.\n";
        exit(EXIT_FAILURE);
      }
}


// Check whether fcn is only 0's or 1's.
// Returns 0 for false, else the uniform bit+1 (ie returning 1 == only 0's)
static Bit uniform_fcn_result(const Bits fcn_bit_values){
  if(fcn_bit_values.empty()) return 0;
  Bit val1 = 2;
  for(auto elt : fcn_bit_values) 
    if(elt != 2) { val1 = elt; break; }
  if(val1 == 2) return 1; // returns uniform 0 for all "Don't Cares"
  for(auto elt : fcn_bit_values) 
    if(elt != val1 && elt != DONT_CARE_BIT) 
      return 0;
  return val1+1;
}


// Applies a function across a generated Truth Table of n vars & 
// returns a vector w/ the results. 
// => NOTE: THE FUNCTION SHALL BE OF TYPE "Bit" & ACCEPT AN ARG OF TYPE "Bits"
template<typename Fcn>
static Bits apply_user_TT_func(const Bit n, Fcn user_func){
  Bits fcn_result;
  if(n == 0) return fcn_result;
  for(Bit counter = 0, rows = 1<<n; counter < rows; ++counter) {
    Bits row_result;
    for(Bit shift = n-1; shift >= 0; --shift) {
      row_result.push_back((counter >> shift) & 1);
      if(shift == 0) break;
    }
    fcn_result.push_back(user_func(row_result));
  }
  return fcn_result;
}


// Returns whether var v is w/in domain of var #N for dimension D
static bool v_in_var_N_domain(const Bit D, const Bit N, const Bit v) {
  if(N >= D) return false;
  if(N == 0) return (v&1); // odds
  if(N == 1) return (!(v&1) && (v&3)) || (!((v-1)&1) && ((v-1)&3)); // %2 & !%4
  const Bit start = (1<<N), end = (1<<(N+1)), max_k = (1<<(D-N-1));
  Bit k = 0;
  while(k < max_k) {
    if(v >= start+(end*k) && v < end+(end*k)) return true;
    ++k;
  }
  return false;
}


// Returns whether all "overlap_cube" points are in domain of variable N (1)
//         or whether none of the points are in the domain              (-1)
//         or some are in the domain but others aren't                   (0)
static int total_or_empty_set_intersection(const Bit D, const Bit N, 
const Bits overlap_cube) {
  Bit total_overlap = 0, no_overlap = 0;
  for(auto v : overlap_cube) { // for each point in cube
    if(v_in_var_N_domain(D, N, v))
      ++total_overlap;
    else
      ++no_overlap;
    if(no_overlap && total_overlap) return 0;
  }
  return total_overlap ? 1 : -1;
}


// Given variable names and fcn bit values, returns string of fcn in SOP/POS form
// => NOTE: variable names in "var_names" shall be listed in descending bit 
//          significance (ie var_names[0] == MSB && var_names[n-1] == LSB).
// => NOTE: bits in "fcn_bit_values" shall be listed in top-down order as if
//          read from a truth table.
// => "sought_bit" = 1 for SOP (minterm) && 0 for POS (maxterm)
static std::string get_SOP_or_POS(std::vector<std::string> var_names, 
const Bits fcn_bit_values, const Bit sought_bit) {
  validate_sop_pos_vars(var_names);
  const Bit D = log2(fcn_bit_values.size());
  if(D == -1 || var_names.size() != D || fcn_bit_values.size() > (1 << D)) 
    return std::string{};
  if(Bit uniform = uniform_fcn_result(fcn_bit_values)) 
    return std::to_string(uniform-1);

  std::reverse(var_names.begin(), var_names.end()); // var-pt domain fcns work backwards
  int intersection = 0;
  std::string SOP_POS_FORM;
  char inner_operation = sought_bit ? '.' : '+'; // SOP or POS
  char outer_operation = sought_bit ? '+' : '.'; // SOP or POS
  auto overlapping_cubes = max_hypercube_overlap_matrix(D, 
                           fcn_bit_values, sought_bit);
  
  Bits cube_vars, cube_compl_vars;
  for(auto overlap_cube : overlapping_cubes) { // for each cube
    cube_vars.clear(), cube_compl_vars.clear();
    
    for(Bit N = 0; N < D; ++N) { // for each var
      intersection = total_or_empty_set_intersection(D, N, overlap_cube);
      if(intersection == 1)
        cube_vars.push_back(N);
      else if(intersection == -1)
        cube_compl_vars.push_back(N);
    }

    if(!SOP_POS_FORM.empty()) SOP_POS_FORM += outer_operation;
    SOP_POS_FORM += '(';

    const Bit total_vars       = cube_vars.size();
    const Bit total_compl_vars = cube_compl_vars.size();

    for(Bit i = 0; i < total_vars; ++i) {
      SOP_POS_FORM += var_names[cube_vars[i]];
      if(!sought_bit) SOP_POS_FORM += '\''; // negate POS here
      if(i < total_vars-1) SOP_POS_FORM += inner_operation;
    }
    if(!cube_compl_vars.empty() && 
      SOP_POS_FORM[SOP_POS_FORM.size()-1] != inner_operation && 
      SOP_POS_FORM[SOP_POS_FORM.size()-1] != '(')
      SOP_POS_FORM += inner_operation;
    for(Bit i = 0; i < total_compl_vars; ++i) {
      SOP_POS_FORM += var_names[cube_compl_vars[i]];
      if(sought_bit) SOP_POS_FORM += '\''; // negate SOP here
      if(i < total_compl_vars-1) SOP_POS_FORM += inner_operation;
    }
    SOP_POS_FORM += ')';
  }

  return SOP_POS_FORM;
}

/******************************************************************************
* FUNCTION MANIPULATION HELPER FUNCTION
******************************************************************************/

// Returns fcn w/o double paren-based negations ( ie "(var')'" => "(var)" )
static std::string rmv_double_compl(std::string fcn) {
  if(fcn.empty()) return "";
  const std::string reserved_chars("()+.^'");
  std::string cleaned_fcn;

  for(auto ch = fcn.begin(); ch != fcn.end(); ++ch) {
    // find next parens
    while(ch != fcn.end() && *ch != '(') cleaned_fcn += *(ch++);
    if(ch == fcn.end())  return cleaned_fcn;
    
    // determine if parens only contain 1 literal
    auto scout = ch+1;
    while(scout!=fcn.end() && reserved_chars.find(*scout)==std::string::npos)
      ++scout;
    while(ch != scout) cleaned_fcn += *(ch++);
    
    // skip double negation if present, else cpy string as is
    if(scout != fcn.end() && scout+1 != fcn.end() && scout+2 != fcn.end() && 
      *scout == '\'' && *(scout+1) == ')' && *(scout+2) == '\'')
      cleaned_fcn += ')', ch += 2;
    else if(ch != fcn.end())
      cleaned_fcn += *ch;
  }
  return cleaned_fcn;
}

/******************************************************************************
* TRUTH TABLE PRINTING HELPER FUNCTION
******************************************************************************/

// Verifies matrix of truth table fcn results = same size & a power of 2
static inline void verify_valid_truth_table_fcn_result_matrix(
const Bit_Matrix fcn_results) {
  const std::string clearNewline = "\x1b[0m\n";
  if(fcn_results.empty()) {
    std::cerr << K_ERR_HEADER 
              << "(empty) FUNCTION RESULTS MATRIX GIVEN!"
              << clearNewline << "Terminating Program.\n";
    exit(EXIT_FAILURE);
  }

  const Bit n = fcn_results[0].size();
  for(auto res : fcn_results)
    if(res.size() != n) {
      std::cerr << K_ERR_HEADER 
                << "FUNCTION RESULTS MATRIX GIVEN!"
                << clearNewline
                << "Result vectors != same length,"
                << " hence cannot derive truth table.\n"
                << "Terminating Program.\n";
      exit(EXIT_FAILURE);
    }

  if(log2(n) == BAD_LOG2_ARG) {
    std::cerr << K_ERR_HEADER 
              << "FUNCTION RESULTS MATRIX GIVEN, RESULT SIZES != POWER OF 2!\x1b[0m\n"
              << "Terminating Program.\n";
    exit(EXIT_FAILURE);
  }
}

/******************************************************************************
* K-MAP PRINTING HELPER FUNCTION
******************************************************************************/

// KMAP-Printer function error handler
static inline void verify_vector_is_kmappable(const Bit TOTAL_RESULTS, 
const Bit LOG_VAL) {
  const std::string clearNewline = "\x1b[0m\n";
  if(TOTAL_RESULTS == 0) {
    std::cerr << K_ERR_HEADER 
              << "(EMPTY) VECTOR OF FUNCTION VALUES != PRINTABLE KMAP!"
              << clearNewline;
  } else if(LOG_VAL == BAD_LOG2_ARG) {
    std::cerr << K_ERR_HEADER 
              << "VECTOR GIVEN FOR KMAP: LENGTH != POWER OF 2!"
              << clearNewline
              << " => HENCE \x1b[1mNOT\x1b[0m A TRUTHTABLE RESULT "
              << "(AS REQUIRED)!\n";
  }
  if(TOTAL_RESULTS == 0 || LOG_VAL == BAD_LOG2_ARG) {
    std::cerr << "Terminating Program.\n";
    exit(EXIT_FAILURE);
  }
}

/******************************************************************************
* XOR-SUBSTITUTION FOR FCN-STRS MAIN FUNCTION
******************************************************************************/

std::string xor_replace(std::string fcn) {
  if(fcn.empty()) return "";
  const std::string AND   = "\\.", OR  = "\\+", ONE = "\\1", TWO = "\\2";
  const std::string BEGIN = "\\(", END = "\\)";
  const std::string VAR   = R"(([^\(\)\+\.\^']+))";
  const std::string SOP_CONJUNCT = R"(\)(([^'|.]*)\+(.*))\()";
  const std::string POS_CONJUNCT = R"(\)(([^'|.]*)\.(.*))\()";
  #define NOT(NOT_VAR) NOT_VAR+"'"
  
  const Bit TOTAL_XORS   = 4; 
  const Bit XOR_RANGE    = (TOTAL_XORS/2)-1;
  const Bit BTWN_XOR_STR = 3, XOR_LHS = 1, XOR_RHS = 2;

  /*
   * Parser's XOR-XNOR Interpretation:
   * xor_SOP  = A.B'  + B.A'
   * xor_POS  = A+B   . A'+B'
   * xnor_SOP = A'.B' + A.B
   * xnor_POS = A'+B  . B'+A
   */
  
  const std::regex xors[TOTAL_XORS] = { // xor_SOP, xor_POS, xnor_SOP, xnor_POS
    std::regex(BEGIN+VAR+AND+NOT(VAR)+SOP_CONJUNCT+TWO+AND+NOT(ONE)+END),
    std::regex(BEGIN+VAR+OR+VAR+POS_CONJUNCT+NOT(ONE)+OR+NOT(TWO)+END),
    std::regex(BEGIN+NOT(VAR)+AND+NOT(VAR)+SOP_CONJUNCT+ONE+AND+TWO+END),
    std::regex(BEGIN+NOT(VAR)+OR+VAR+POS_CONJUNCT+NOT(TWO)+OR+ONE+END)
  };
  #undef NOT

  for(Bit i = 0; i < TOTAL_XORS; ++i)
    for(std::sregex_iterator it(fcn.begin(), fcn.end(), xors[i]); it != std::sregex_iterator{}; ++it)
      if(it->suffix().str()[0] != '\'')
        fcn = it->prefix().str() + it->str(BTWN_XOR_STR).substr(1) + 
              "(" + it->str(XOR_LHS) +
              ((i>XOR_RANGE) ? "'" : "") +
              "^" + it->str(XOR_RHS) + 
              ")" + it->suffix().str();
  return fcn;
}

/******************************************************************************
* FUNCTION MANIPULATION MAIN FUNCTIONS
******************************************************************************/

// Returns fcn w/ its literals complemented
std::string compl_literals(std::string fcn) {
  if(fcn.empty()) return "";
  const std::string reserved_chars("()+.^");
  auto ch = fcn.rbegin();
  while(ch != fcn.rend()) {
    // if at a literal, rmv or add complement
    if(reserved_chars.find(*ch)==std::string::npos &&
      ((*ch == '\'' && *(ch+1) != ')') || *ch != '\'')) {
      (*ch != '\'') ? fcn.insert(ch.base(), '\'') : fcn.erase(((ch++)+1).base());
      while(ch!=fcn.rend() && reserved_chars.find(*ch)==std::string::npos)++ch;
    } else ++ch;
  }
  return rmv_double_compl(fcn);
}


// Returns dual of the fcn (swap '+' & '.', negate '^' (as XOR dual == XNOR))
std::string dual_fcn(std::string fcn) {
  if(fcn.empty()) return "";
  // swap '+' & '.', negate '^' to become XNOR (or XNOR to XOR)
  for(auto ch = fcn.begin(); ch != fcn.end(); ++ch) {
    switch(*ch) {
      case '+': *ch = '.'; break;
      case '.': *ch = '+'; break;
      case '^': 
        switch(*(ch-1)) {
          case '\'': fcn.erase(ch-1); break;
          default:   fcn.insert(ch++, '\'');
        }
        break;
    }
  }
  return rmv_double_compl(fcn);
}


// Returns complement of the fcn (dual + negated literals)
inline std::string compl_fcn(std::string fcn) {
  return compl_literals(dual_fcn(fcn));
}

/******************************************************************************
* SOP-POS FUNCTION GENERATION MAIN FUNCTIONS
******************************************************************************/

// => NOTE: Both "SOP_str" & "POS_str" have 2 versions, the 1st arg for all of 
//          which is a string vector of variable names, and the 2nd arg a:
//           - fcn to invoke on a generated truth table for the given variables
//               => NOTE: THE FCN SHALL BE TYPE "Bit" & ACCEPT ARG TYPE "Bits"
//           - unsigned vector of fcn results (ie a truth table's fcn "column")

// Returns string of given fcn in SOP form using variable names from vector
// => NOTE: variable names in "var_names" shall be listed in descending bit 
//          significance (ie var_names[0] == MSB && var_names[n-1] == LSB).
// => NOTE: bits in "fcn_bit_values" shall be listed in top-down order as if
//          read from a truth table.
std::string SOP_str(std::vector<std::string>var_names,
const Bits fcn_bit_values){
  return get_SOP_or_POS(var_names, fcn_bit_values, 1);
}
template<typename Fcn>
std::string SOP_str(std::vector<std::string> var_names, Fcn bit_row_fcn){
  return SOP_str(var_names,apply_user_TT_func(var_names.size(),bit_row_fcn));
}


// Returns string of given fcn in POS form using variable names from vector
// => NOTE: variable names in "var_names" shall be listed in descending bit 
//          significance (ie var_names[0] == MSB && var_names[n-1] == LSB).
// => NOTE: bits in "fcn_bit_values" shall be listed in top-down order as if
//          read from a truth table.
std::string POS_str(std::vector<std::string> var_names, 
const Bits fcn_bit_values) {
  return get_SOP_or_POS(var_names, fcn_bit_values, 0);
}
template<typename Fcn>
std::string POS_str(std::vector<std::string> var_names, Fcn bit_row_fcn){
  return POS_str(var_names,apply_user_TT_func(var_names.size(),bit_row_fcn));
}

/******************************************************************************
* TRUTH TABLE FUNCTION APPLICATION MAIN FUNCTION
******************************************************************************/

// Applies given fcn across an n-bit truth table, returning vector of results.
template<typename Fcn>
Bits TruthTable_fcn(const Bit n, Fcn func) {
  Bits results, row;
  if(n == 0) return results;
  for(Bit counter = 0, rows = (1<<n); counter < rows; ++counter){
    row.clear();
    for(Bit shift = n-1; shift >= 0; --shift) {
      row.push_back((counter >> shift) & 1);
      if(shift == 0) break;
    }
    results.push_back(func(row));
  }
  return results;
}

/******************************************************************************
* TRUTH TABLE PRINTING MAIN FUNCTIONS
******************************************************************************/

// Prints Truth Table for given a matrix of truth table function results.
void TruthTable_print(const Bit_Matrix fcn_results) {
  verify_valid_truth_table_fcn_result_matrix(fcn_results);
  const Bit bit_No = log2(fcn_results[0].size());
  if(bit_No == 0) return;
  const Bit ROWS = (1 << bit_No);
  const bool has_fcns = !fcn_results.empty();
  const Bit fcn_header_length = (has_fcns?(4+2*fcn_results.size()):0);
  const std::string table_header_sperator(2*bit_No+4+fcn_header_length, '-');

  std::cout << "Bit:";
  for(Bit i = 0; i < bit_No; ++i) std::cout << " " << i;
  if(has_fcns) {
    std::cout << " | F";
    for(Bit i = 0; i < fcn_results.size(); ++i) std::cout << " " << i;
  }
  std::cout << std::endl << table_header_sperator << std::endl;

  for(Bit counter = 0; counter < ROWS; ++counter) {
    std::cout << std::setw(3) << counter << ": ";
    for(Bit shift = bit_No-1; shift >= 0; --shift) {
      std::cout << ((counter >> shift) & 1) << " ";
      if(shift == 0) break;
    }
    if(has_fcns) {
      std::cout << "|   ";
      for(auto res : fcn_results) 
        std::cout << ((res[counter] == DONT_CARE_BIT)
                      ?"X":std::to_string(res[counter])) 
                  << " ";
    }
    std::cout << std::endl;
  }
}

// Prints Truth Table for given 1 vector of truth table function results.
void TruthTable_print(Bits fcn_results) {
  TruthTable_print(Bit_Matrix{fcn_results});
}



// Prints Truth Table for "bit_No" bits & fcn results from a vector of fcns
//   to apply to the table.
void TruthTable_print(const Bit bit_No, 
const std::vector<std::function<Bit(Bits)>> fcns_vec){
  if(bit_No == 0) return;
  const Bit total_fcns = fcns_vec.size();
  const Bit ROWS = (1 << bit_No);
  Bit_Matrix fcn_results(total_fcns, Bits{});

  Bits row;
  for(Bit counter = 0; counter < ROWS; ++counter) {
    row.clear();
    // Get row of truth table bits
    for(Bit shift = bit_No-1; shift >= 0; --shift) {
      row.push_back((counter >> shift) & 1);
      if(shift == 0) break;
    }
    // Get result from each function acting on the row
    for(Bit i = 0; i < total_fcns; ++i)
      fcn_results[i].push_back(fcns_vec[i](row));
  }

  TruthTable_print(fcn_results);
}

// Prints Truth Table for "bit_No" bits & fcn results from 1 fcn to apply 
//   to the table.
template<typename Fcn>
void TruthTable_print(const Bit bit_No, Fcn func){
  TruthTable_print(bit_No, {func});
}



// Prints Truth Table for "bit_No" bits.
void TruthTable_print(const Bit bit_No) {
  if(bit_No == 0) return;
  const Bit ROWS = (1 << bit_No);
  const std::string table_header_sperator(2*bit_No+4, '-');

  std::cout << "Bit:";
  for(Bit i = 0; i < bit_No; ++i) std::cout << " " << i;
  std::cout << std::endl << table_header_sperator << std::endl;

  for(Bit counter = 0; counter < ROWS; ++counter) {
    std::cout << std::setw(3) << counter << ": ";
    for(Bit shift = bit_No-1; shift >= 0; --shift) {
      std::cout << ((counter >> shift) & 1) << " ";
      if(shift == 0) break;
    }
    std::cout << std::endl;
  }
}

/******************************************************************************
* K-MAP PRINTING MAIN FUNCTIONS
******************************************************************************/

// Given std::vector<unsigned> of truth table fcn values, prints
//   its kmap layout.
// SAMPLE KMAP LAYOUT INTERPRETATION: 
// => B(n), where lower n == more significant bit
//    0  0 B2 B2
// 0             0
// 0             B1
// B0            B1
// B0            0
//    0 B3 B3  0
void Kmap_print(const Bits function_bit_values, bool style_kmap_var_sections = true) {
  using std::cout;
  // Cell ASCII Layout
  const std::string cell_top(" -----"), 
                    cell_topleft("|  "), cell_topright("  "),   
                    cell_botleft("|"),   cell_botright("  "), 
                    new_cell_row("|\n");
  // KMAP Constants
  const Bit TOTAL_RESULTS = function_bit_values.size(),
            bit_No        = log2(TOTAL_RESULTS),
            ROWS          = (bit_No<4)?(bit_No<3)?1:2:4,
            COLS          = (bit_No==1)?2:4,
            FCN_VALS_PER_4VAR_KMAP = 16,
            kmap_truthTable_idx[4][4] = {{0,1,3,2},{4,5,7,6},{12,13,15,14},{8,9,11,10}};
  verify_vector_is_kmappable(TOTAL_RESULTS, bit_No);

  // Output Formatting to Discern Variable Map Sections
  Bit i = 0, j = 0;
  const std::string MSB3("\x1b[1m"),  // MOST-SIGNIFICANT-BIT[3] BOLD
                    MSB2("\x1b[4m"),  // MOST-SIGNIFICANT-BIT[2] UNDERLINED 
                    MSB1("\x1b[31m"), // MOST-SIGNIFICANT-BIT[1] RED FONT 
                    MSB0("\x1b[7m");  // MOST-SIGNIFICANT-BIT[0] REVERSE TEXT & BACKGROUND COLORS
  const std::string kmap_section_styles[FCN_VALS_PER_4VAR_KMAP] = 
  {
    "",        MSB0,           MSB0+MSB1,           MSB1, 
    MSB2,      MSB0+MSB2,      MSB0+MSB1+MSB2,      MSB1+MSB2, 
    MSB2+MSB3, MSB0+MSB2+MSB3, MSB0+MSB1+MSB2+MSB3, MSB1+MSB2+MSB3, 
    MSB3,      MSB0+MSB3,      MSB0+MSB1+MSB3,      MSB1+MSB3
  };

  // Formatting Lambdas
  auto format = [&](){
    if(i==ROWS)i=ROWS-1;
    return(style_kmap_var_sections)?kmap_section_styles[i*4+j]:"";
  };
  auto clearFormat = [&](){return(style_kmap_var_sections)?"\x1b[0m":"";};

  // Formatting Explanation:
  if(style_kmap_var_sections) {
    cout << std::endl << R"(>> "MSB[n]" == "M(ost)S(gfnt)B(it)[descending]")";
    if(bit_No >= 4) cout << "\n\x1b[1mMSB[3]\x1b[0m <BOLD>";
    if(bit_No >= 3) cout << "\n\x1b[4mMSB[2]\x1b[0m <UNDERLINE>";
    if(bit_No >= 2) cout << "\n\x1b[31mMSB[1]\x1b[0m <RED_FONT>";
    if(bit_No >= 1) cout << "\n\x1b[7mMSB[0]\x1b[0m <REVERSE_TEXT_BACKGROUND>";
    cout << std::endl;
  }

  for(Bit kmap4var_instanceNo = 0; kmap4var_instanceNo < TOTAL_RESULTS; 
  kmap4var_instanceNo+=FCN_VALS_PER_4VAR_KMAP) {
    Bit extension = std::to_string(kmap4var_instanceNo+FCN_VALS_PER_4VAR_KMAP-1).size();
    extension = (extension > 2) ? (extension - 2) : 0;
    std::string top_extend(extension, '-'), side_extend(extension, ' ');
    // for each 4-var kmap composing total kmap 
    // (ie 1 if 4-var fcn, 2 if 5-var, & 4 if 6-var)
    for(i = 0; i < ROWS; ++i) {
      // print kmap row delimiter
      for(j = 0; j < COLS; ++j) cout << cell_top+top_extend;
      cout << "\n";
      
      // print kmap function values
      for(j = 0; j < COLS; ++j) {
        const Bit func_bit = function_bit_values[kmap_truthTable_idx[i][j]+kmap4var_instanceNo];
        cout << cell_topleft << format() 
             << ((func_bit==DONT_CARE_BIT)
                 ? "X" : ((func_bit==BLANK_BIT)
                           ? " " : std::to_string(func_bit)))
             << clearFormat() << cell_topright+side_extend;
      }
      cout << new_cell_row;
      
      // print kmap idxs
      for(j = 0; j < COLS; ++j) {
        const Bit idx_length = std::to_string(kmap_truthTable_idx[i][j]+kmap4var_instanceNo).size();
        cout << cell_botleft << format()
             << ((kmap_truthTable_idx[i][j]+kmap4var_instanceNo<10)?"0":"") 
             << kmap_truthTable_idx[i][j]+kmap4var_instanceNo 
             << ")" << clearFormat() << cell_botright+((idx_length>1 && idx_length-2<extension)?" ":"");
      }
      cout << new_cell_row;
    }

    // print kmap footer delimiter
    for(j = 0; j < COLS; ++j) cout << cell_top+top_extend;
    cout << "\n\n";

  } // end "kmap4var_instanceNo" for loop
} // end function


// Given kmap's # of bits (n) & a fcn to apply across a truth table of n-bits, 
//   prints the corresponding kmap layout of the fcn.
template <typename Fcn>
void Kmap_print(const Bit n, Fcn bit_row_fcn, 
bool style_kmap_var_sections = true) {
  Bits fcn_result, row_result;
  if(n == 0) return;
  for(Bit counter = 0, rows = 1<<n; counter < rows; ++counter) {
    row_result.clear();
    for(Bit shift = n-1; shift >= 0; --shift) {
      row_result.push_back((counter >> shift) & 1);
      if(shift == 0) break;
    }
    fcn_result.push_back(bit_row_fcn(row_result));
  }
  Kmap_print(fcn_result, style_kmap_var_sections);
}


// Print a blank Kmap of n bits
void Kmap_print(const Bit n, bool style_kmap_var_sections = true) {
  Kmap_print(n, [](auto v){return 3;}, style_kmap_var_sections);
}

/******************************************************************************
* BCD DECODER FUNCTION FOR IDENTIFYING CURRENT TRUTH TABLE ROW
******************************************************************************/

// Returns decimal value for bit vector ordered in descending bit significance.
// Intended to help users identify which row they're currently on when making
//    custom fcns to operate on truth table rows (ie whether in a "Don't Care")
inline Bit BCD_decoder(const Bits bit_vec){
  Bit decimal_value = 0;
  for(Bit b : bit_vec)
    decimal_value = (decimal_value << 1) | b;
  return decimal_value;
}

#endif

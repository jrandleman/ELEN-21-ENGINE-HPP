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
#include <set>
#include <string>
#include <vector>

#define K_ERR_HEADER "\x1b[1m" + std::string(__FILE__)    + ":" +\
                                 std::string(__func__)    + ":" +\
                                 std::to_string(__LINE__) + ":" +\
                                 "\x1b[31mERROR:\x1b[0m\x1b[1m INVALID "

using Bit = unsigned long long;
using Bits = std::vector<Bit>;

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
                                           std::set<Bit>&points_connected_to_P){
  if(!D || P > (1 << original_D)) return;
  const Bit half_of_dimension = (1 << (D-1));
  if(P >= half_of_dimension) {
    points_connected_to_P.insert(original_P-half_of_dimension);
    find_adjacent_hypercube_points(D-1, original_D, P-half_of_dimension, original_P, points_connected_to_P);
  } else {
    points_connected_to_P.insert(original_P+half_of_dimension);
    find_adjacent_hypercube_points(D-1, original_D, P, original_P, points_connected_to_P);
  }
}


// Returns set of points connected to point P in hypercube of dimension D
// Time & Space Complexity: O(D)
std::set<Bit> adjacent_hypercube_points(const Bit D, const Bit P) {
  std::set<Bit> connected_points;
  find_adjacent_hypercube_points(D, D, P, P, connected_points);
  return connected_points;
}


// Returns adjacency list for points in dimension D via MxN vector-set matrix.
//   => M == parent & N == parent's corresponding adjacent points
// Time & Space Complexity: O(D*2^D)
std::vector<std::set<Bit>> dimension_adjacency_matrix(const Bit D) {
  std::vector<std::set<Bit>> adj_matrix;
  for(Bit P = 0, total_P = (1 << D); P < total_P; ++P)
    adj_matrix.push_back(adjacent_hypercube_points(D,P));
  return adj_matrix;
}

/******************************************************************************
* NESTED HYPERCUBE CONTAINMENT FUNCTIONS
******************************************************************************/

// Confirms set "points" all fall w/in dimension D's range
// Time  Complexity: O(points.size())
// Space Complexity: O(1)
static bool points_are_within_D(const Bit D, const std::set<Bit>points){
  const Bit total_D_points = (1 << D);
  if(points.empty() || points.size()>total_D_points) return false;
  for(Bit pt : points) if(pt >= total_D_points) return false;
  return true;
}


// Confirms set "points" forms a hypercube w/in another hypercube dimension's 
//   adjacency matrix.
//   => HYPERCUBE DEFN: 2^n points with n edges btwn one another.
// Time  Complexity: O(points.size()^2)
// Space Complexity: O(1)
bool adj_matrix_contains_hypercube(const std::vector<std::set<Bit>>&adj_matrix, 
const std::set<Bit>points) {
  const Bit edges_btwn_pts = log2(points.size()); // sought subdimension
  Bit edge_counter;
  for(Bit pt : points) {
    edge_counter = 0;
    for(Bit pt2 : points)
      if(adj_matrix[pt].find(pt2) != adj_matrix[pt].end())
        ++edge_counter;
    if(edge_counter != edges_btwn_pts) return false;
  }
  return true;
}

/******************************************************************************
* COMBINATIONS HELPER FUNCTIONS
******************************************************************************/

// Recursive fcn to mk a combinations set matrix for "combinations" fcn below
template <typename T>
static void find_combinations(const Bit r, 
const typename T::iterator n_begin, const typename T::iterator n_end, 
T seen, std::set<T>&combo_matrix) {
  typename T::iterator temp;
  for(auto it = n_begin; it != n_end;) {
    seen.insert(*it);
    temp = it++;
    switch(r) {
      case 1:  combo_matrix.insert(seen); break;
      default: find_combinations(r-1, it, n_end, seen, combo_matrix);
    } 
    seen.erase(*temp);
  }
}


// Returns a set matrix of possible combinations length r w/in container 
// group_n of type "T", where T is a std::set of any type of elt.
// => Hence a set of ints are just as easy combined as a set of sets.
// NOTE: r == group size w/in set "group_n" we are finding combos 
template <typename T>
std::set<T> combinations(const Bit r,const T group_n){
  std::set<T> combo_matrix;
  if(r > group_n.size())
    return combo_matrix;
  find_combinations(r, group_n.begin(), group_n.end(), {}, combo_matrix);
  return combo_matrix;
}

/******************************************************************************
* NESTED POINT-OVERLAPPING HYPERCUBES HELPER FUNCTIONS
******************************************************************************/

// Returns a set of indices where "function_bit_values" has a "sought" value
// => "sought_bit" = 1 for SOP (minterm) && 0 for POS (maxterm)
static void convert_fcn_to_point_set(
const Bits function_bit_values, const Bit sought_bit,std::set<Bit>&need_to_cover, 
std::set<Bit>&points){
  for(Bit i = 0; i < function_bit_values.size(); ++i)
    if(function_bit_values[i] == sought_bit)
      points.insert(i), need_to_cover.insert(i);
    else if(function_bit_values[i] == DONT_CARE_BIT)
      points.insert(i);
}


// Confirms given sets are disjoint
static inline bool sets_are_disjoint(
const std::set<Bit>&set1,const std::set<Bit>&set2){
  for(auto elt : set1)
    if(set2.find(elt) != set2.end())
      return false;
  return true;
}


// Returns greatest hypercube possible w/ "total_vertices" # of vertices
static Bit max_formable_hypercube(Bit max_inner_hypercube_vertices, 
const Bit total_vertices) {
  while(max_inner_hypercube_vertices > total_vertices) 
    max_inner_hypercube_vertices >>= 1;
  return log2(max_inner_hypercube_vertices);
}


// Returns vector of set matrices w/ possible hypercube coordinate arrangements
// for the given "points" set
static std::vector<std::set<std::set<Bit>>> possible_hypercube_combinations(
const Bit max_inner_hypercube, const std::set<Bit> points) {
  std::vector<std::set<std::set<Bit>>> hypercube_combos;
  for(Bit dim = max_inner_hypercube; dim >= 0; --dim) {
    hypercube_combos.push_back(combinations((1 << dim), points));
    if(dim == 0) return hypercube_combos;
  }
  return hypercube_combos;
}

/******************************************************************************
* ESSENTIAL PRIME IMPLICANTS FROM PI'S PARSING FUNCTIONS
******************************************************************************/

// Returns set of all unique elts within the dimensionSet matrix
std::set<Bit> get_set_of_elts_in_dimension(const std::set<std::set<Bit>>dimensionSet) {
  std::set<Bit> uniqueElts;
  for(auto dim : dimensionSet)
    for(auto point : dim)
      uniqueElts.insert(point);
  return uniqueElts;
}


// Return matrix of all possible combinations of start-to-end groups of
// PI's w/in dimensionSet
static std::set<std::set<std::set<std::set<Bit>>>> get_all_EPI_candidate_combos(
const Bit start, const Bit end, const std::set<std::set<Bit>> dimensionSet) {
  std::set<std::set<std::set<std::set<Bit>>>> possible_EPI_3D_matrix; // Big Oof
  for(Bit i = start; i <= end; ++i)
    possible_EPI_3D_matrix.insert(combinations(i, dimensionSet));
  return possible_EPI_3D_matrix;
}


// confirm whether point matrix (or epi set) contains all of needToCoverElts
static bool set_contains_all_elts(const std::set<std::set<Bit>>epiSet, 
const std::set<Bit> needToCoverElts) {
  for(auto elt : needToCoverElts) {
    bool found = false;
    for(auto epi : epiSet) {
      for(auto point : epi)
        if(point == elt) {
          found = true;
          break;
        }
      if(found) break;
    }
    if(!found) return false;
  }
  return true;
}


// Return EPIs from PIs given in dimensionSet 
// (EPI = Essential Prime Implicant, PI = Prime Implicant)
static std::set<std::set<Bit>> EPIs_within_dimensionSet(
const std::set<std::set<Bit>> dimensionSet, const std::set<Bit> desiredCoveredBits) {
  if(dimensionSet.empty()) return {};
  const Bit m = dimensionSet.size();    // number of sub dimensions
  const Bit n = (*dimensionSet.begin()).size(); // number of points per sub dimension
  auto elts_in_dimension = get_set_of_elts_in_dimension(dimensionSet);

  // Unique elts are those we seek to find EPI's covering -- the intersection
  // of desired bits we ultimately seek to completely cover on the Kmap with 
  // the derived set of bits that this current subdimension DOES cover (those 
  // remaining the "desiredCoveredBits" will be covered by a lower dimension)
  std::set<Bit> uniqueElts;
  for(auto elt : elts_in_dimension)
    if(desiredCoveredBits.find(elt) != desiredCoveredBits.end())
      uniqueElts.insert(elt);

  // minimum sets needed to cover all unique elts
  const Bit minimum_sets_to_cover_points = 1+((uniqueElts.size()-1)/n);  // cieling(uniqueElts.size()/n)

  // All potential EPI combos, sorted in ascending length (hence the BEST 
  // (fewer dependancies) epi set will be closer to the front)
  auto epi_3D_set_combinations = get_all_EPI_candidate_combos(minimum_sets_to_cover_points, m, dimensionSet); 

  for(auto epi_set_combination : epi_3D_set_combinations) // per 3D matrix in 4D matrix
    for(auto epi_set : epi_set_combination)               // per 2D matrix in 3D matrix
      if(set_contains_all_elts(epi_set, uniqueElts))
        return epi_set;

  return {};
}

/******************************************************************************
* NESTED POINT-OVERLAPPING HYPERCUBES MAIN FUNCTION
******************************************************************************/

// Returns hypercube overlap matrix of vertices in set "points" w/in dimension D
// => IE the sets of vertices w/in "points" that are K-map logically adjacent!
std::set<std::set<Bit>> max_hypercube_overlap_matrix(const Bit D, 
const Bits function_bit_values, const Bit sought_bit) {
  std::set<Bit> need_to_cover, points;
  convert_fcn_to_point_set(function_bit_values,sought_bit,need_to_cover,points);
  if(D == 0 || !points_are_within_D(D,points)) 
    return std::set<std::set<Bit>>{};

  // find greatest formable hypercube from set "points"
  Bit max_inner_hypercube = max_formable_hypercube((1<<D), points.size());
  if(max_inner_hypercube > D) return std::set<std::set<Bit>>{};

  // adajcency matrix for points in dimension D 
  //   => ie: logical adjacency matrix for cells in K-map of D vars
  const std::vector<std::set<Bit>>adj_matrix=dimension_adjacency_matrix(D);

  // max overlaps (prime implicants), points (cells) not yet covered, & all 
  //   possible hypercube combos made from "points"
  std::set<std::set<Bit>> max_hypercube_overlaps; 
  auto hypercube_combos = possible_hypercube_combinations(max_inner_hypercube, points);

  // get formable hypercubes from "points" (prime implicants)
  for(Bit i = 0; i <= max_inner_hypercube && !need_to_cover.empty(); ++i) {

    // track all hypercubes covering the locally critical points 
    // (points not previously covered by any higher dimension)
    // and parse out EPI's within from a seperate function 
    // (this fcn only serves to ID Prime Implicants not the ESSENTIALS)
    auto LOCAL_DIM_NEED_TO_COVER(need_to_cover);
    std::set<std::set<Bit>> local_dim_prime_implicants;

    for(auto cube_set : hypercube_combos[i]) {
      if(!sets_are_disjoint(cube_set, LOCAL_DIM_NEED_TO_COVER) && adj_matrix_contains_hypercube(adj_matrix, cube_set)) {
        // erase found cube's pts from set "need_to_cover"
        for(auto pt : cube_set) 
          need_to_cover.erase(pt);
        // save dimension/hypercube formed
        local_dim_prime_implicants.insert(cube_set);
        // erase subset dimensions of the found dimension
        for(Bit j = i+1; j <= max_inner_hypercube; ++j)
          for(auto sub_cube_set = hypercube_combos[j].begin(); 
              sub_cube_set != hypercube_combos[j].end();) {
            if(sets_are_disjoint(*sub_cube_set, need_to_cover))
              hypercube_combos[j].erase(*(sub_cube_set++));
            else
              ++sub_cube_set;
          }
      } 
    }

    // Parse & save essential prime implicants from the current dimension's PI's
    auto epi_set = EPIs_within_dimensionSet(local_dim_prime_implicants, LOCAL_DIM_NEED_TO_COVER);
    for(auto epi : epi_set)
      max_hypercube_overlaps.insert(epi);
  }

  return max_hypercube_overlaps;
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
const std::set<Bit> overlap_cube) {
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
  
  for(auto overlap_cube : overlapping_cubes) { // for each cube
    Bits cube_vars, cube_compl_vars;
    
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
const std::vector<Bits> fcn_results) {
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
  Bits results;
  if(n == 0) return results;
  for(Bit counter = 0, rows = (1<<n); counter < rows; ++counter){
    Bits row;
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
void TruthTable_print(const std::vector<Bits> fcn_results) {
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
  TruthTable_print(std::vector<Bits>{fcn_results});
}



// Prints Truth Table for "bit_No" bits & fcn results from a vector of fcns
//   to apply to the table.
void TruthTable_print(const Bit bit_No, 
const std::vector<std::function<Bit(Bits)>> fcns_vec){
  if(bit_No == 0) return;
  const Bit total_fcns = fcns_vec.size();
  const Bit ROWS = (1 << bit_No);
  std::vector<Bits> fcn_results(total_fcns, Bits{});

  for(Bit counter = 0; counter < ROWS; ++counter) {
    Bits row;
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
  Bits fcn_result;
  if(n == 0) return;
  for(Bit counter = 0, rows = 1<<n; counter < rows; ++counter) {
    Bits row_result;
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

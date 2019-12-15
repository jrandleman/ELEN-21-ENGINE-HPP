// Author: Jordan Randleman - elen21_SampleExec.cpp – elen21.hpp Demo File

#include <iostream>
#include <string>

#include "elen21.hpp"

/*
 * 15 FUNCTIONS FOR USERS:
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
 *    10) Kmap_print: given "Bits" of fcn results (as if from a truth table)
 *    11) Kmap_print: given bit# & fcn to apply across generated truthtable
 *
 *    12) TruthTable_print: given bit#
 *    13) TruthTable_print: given bit# & vector of fcns to apply to table
 *    14) TruthTable_print: given vector of "Bits" (results of truth table fcns)
 *
 *    15) TruthTable_fcn: given bit# & a fcn, returns vector of fcn results 
 *                        from truth table
 */

Bit xnor01(Bits TT_row) {return !TT_row[0] ^ TT_row[1];}
Bit xnor02(Bits TT_row) {return !TT_row[0] ^ TT_row[2];}
Bit xnor03(Bits TT_row) {return !TT_row[0] ^ TT_row[3];}

Bit func_1(Bits v) {return (v[0] & v[1]) ^ v[2];}

int main() {

  using namespace std;

  // Generating POS & SOP Forms
  cout << "\nGenerating POS & SOP Forms:\n";
  string sop = SOP_str({"W", "X", "Y", "Z"}, {0,0,1,0,1,1,1,1,1,1,1,0,0,0,1,1});
  string pos = POS_str({"W", "X", "Y", "Z"}, {0,0,1,0,1,1,1,1,1,1,1,0,0,0,1,1});
  cout << " => SOP = " << sop << endl 
       << " => POS = " << pos << "\n\n";



  // Use a function to generate it's own minimum SOP/POS form:
  cout << "Use a function to generate its own minimum SOP/POS form:\n => ";
  cout << SOP_str({"W", "X", "Y", "Z"}, func_1) << "\n => ";
  cout << POS_str({"W", "X", "Y", "Z"}, func_1) << "\n\n";



  // Using Don't Cares:
  // => NOTE: USE A "2" TO DECONTE A DON'T CARE BIT!
  cout << "Using Don't Cares (USE \"2\" TO DENOTE A \"DON'T CARE\" BIT!):"
       << "\n => Minterms = 0,3,4:                 " << SOP_str({"A", "B", "C"}, {1,0,0,1,1,0,0,0})
       << "\n => Minterms = 0,4 & Don't Cares = 3: " << SOP_str({"A", "B", "C"}, {1,0,0,2,1,0,0,0})
       << "\n\n";



  // Complementing a fcn, its literals, taking the dual, and replacing xors:
  cout << "Complementing a fcn, its literals, taking the dual, & replacing XORs:\n";
  string temp("(X.W')+(A'.Y')+(W.X')+(A.Y)");
  cout << " => Demo Function String  = " << temp << endl;
  cout << " => Xor's Replaced        = " << xor_replace(temp) << endl;
  cout << " => Complemented Fcn      = " << compl_fcn(temp) << endl;
  cout << " => Dual Fcn              = " << dual_fcn(temp) << endl;
  cout << " => Complemented Literals = " << compl_literals(temp) << "\n\n";



  // Getting the decimal value of a binary vector 
  // => helps ID don't care rows in custom fcns
  cout << "Decimal for Binary of \"0101\" = " 
       << BCD_decoder({0,1,0,1}) << "\n\n";



  // Printing Truth Tables:
  cout << "Printing a Truth Table of 4 Bits:\n";
  TruthTable_print(4);
  
  cout << "\nPrinting a Truth Table of 4 Bits,\n"
       << "  XORing the most & 2nd least significant bits:\n";
  TruthTable_print(4, [](auto v){return v[0] ^ v[2];});

  cout << "\nPrinting a Truth Table of 4 Bits\n"
       << "  and applying a vector of 3 fcns:\n";
  TruthTable_print(4, {xnor01, xnor02, xnor03});

  cout << "\nPrinting a Truth Table of 2 bits given a \"Bits\" vector\n"
       << "  of function results:\n";
  TruthTable_print({{0,1,1,0}, {1,0,0,1}});
  cout << "\n\n";


  // Printing K-Maps:
  cout << "Printing a K-Map from given \"Bits\" of a fcn's\n"
       << "  truth table results:\n";
  Kmap_print({1,1,0,0,1,1,0,0});

  cout << "\nPrinting a K-Map for the given # of variables\n"
       << "  + a fcn to apply to their truth table:\n";
  Kmap_print(4, xnor01);



  cout << "\nBye!\n\n";
  return 0;
}

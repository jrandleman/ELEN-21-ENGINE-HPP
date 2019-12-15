# ELEN21-HPP
## _SOP/POS Minimization For N bits + Kmap &amp; Truth Table Generation!_
### _Compile >= C++17!_
----------------------

## The Library's Main 4 Functions:
* ***NOTE:** `SOP_str` & `POS_str`'s variable name string vector arg shall be sorted in __DESCENDING__ bit significance!*
1) `SOP_str`: _Returns string with a function in minimum SOP form!_</br>
              => Given either a `vector<string>` with variable names & `vector<unsigned>` of fcn results (as if from a truth table)</br>
              => OR a `vector<string>` with variable names & a fcn to apply to a generated truth table
2) `POS_str`: _Returns string with a function in minimum POS form!_</br>
              => Given either a `vector<string>` with variable names & `vector<unsigned>` of fcn results (as if from a truth table)</br>
              => OR a `vector<string>` with variable names & a fcn to apply to a generated truth table
3) `Kmap_print`: _Prints a K-Map to `stdout`!_</br>
              => Given either a `vector<unsigned>` of fcn results (as if from a truth table)</br>
              => OR a # of bits (# of K-Map vars) along with a fcn to apply across a generated truthtable
4) `TruthTable_print`: _Prints a Truth Table to `stdout`!_</br>
              => Given either a # of bits (# of variables in Truth Table)</br>
              => OR a # of bits and a `vector` or fcns to apply to the table</br>
              => OR a `vector<vector<unsigned>>` matrix holding fcn results (as if from a truth table)

## The Library's 4 "Function String" Helper Functions:
* ***NOTE:** Used to further manipulate fcns returned by `SOP_str` or `POS_str`!*
1) `xor_replace`:    _Returns `string` with fcn's expanded XOR/XNOR's abbreviated!_</br>
2) `compl_literals`: _Returns `string` with fcn's literals complemented!_</br>
3) `dual_fcn`:       _Returns `string` with fcn's dual!_</br>
4) `compl_fcn`:      _Composes "`compl_literals`" & "`dual_fcn`"!_

## The Library's 2 Additional Functions & Universal Constant:
1) `BCD_decoder`: _Convert `vector<unsigned>` of bits into a decimal #!_</br>
   => Helps with implementing "don't cares" in user-defined custom fcns
2) `TruthTable_fcn`: _Given a # of bits & a fcn, returns `vector` of fcn results across truth table!_
3) `DONT_CARE_BIT`: _Universal predefined constant - to be returned by user fcns if at a "don't care" row!_

-------------
## Sample User-Defined Function:
```c++
// NOTE the mandatory predefined "Bit" return & "Bits" arg types!
// NOTE "Bits" = vector sorted by descending bit significance, 
//      contining a "row" of truth table bits as its arg!
 
Bit your_func_name_here(Bits bit_vector) { 
   Bit dont_care_row = 3; // if we don't care about row 3 (rows starting from 0)
   
   // BCD_decoder returns decimal value of given bit vector
   if(BCD_decoder(bit_vector) == dont_care_row) 
     return DONT_CARE_BIT; // use predefined "DONT_CARE_BIT" to denote don't cares
     
   return (bit_vector[0] & !bit_vector[1]) ^ bit_vector[3] // any operation
}
```

-----------
## Operators Guide:
```c++
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
```

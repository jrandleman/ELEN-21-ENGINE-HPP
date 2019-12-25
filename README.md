# ELEN21-HPP
## _SOP/POS Minimization For N bits + Kmap &amp; Truth Table Generation!_
### _Compile >= C++17!_
----------------------

## General Information For [`elen21.hpp`](https://github.com/jrandleman/ELEN-21-ENGINE-HPP/blob/master/elen21.hpp):
### Library's 2 Predefined Types:
1) `Bit`: `unsigned long long`</br>
2) `Bits`: `vector<Bit>`</br>
=> *Use to denote bit sequences & perform operations in `fcn`s passed to the library!*</br>

### README Notes:
* *Refer to [`elen21_SampleExec.cpp`](https://github.com/jrandleman/ELEN-21-ENGINE-HPP/blob/master/elen21_SampleExec.cpp) for a demo!*</br>
* *`#` = any number of type `Bit`*</br>
* *"`name vector`" = `vector<string>` of variable names, ordered by **DESCENDING** bit significance!*</br>
* *"`fcn`" = user-defined function: returns `Bit` & takes arg `Bits` ([more on this below](#sample-user-defined-function))*</br>

----------------------

## 4 Main Functions:

1) `SOP_str`: _Returns SOP `string`!_</br>
              => Given `name vector` & `Bits` of fcn results (as if from a truth table)</br>
              => Given `name vector` & `fcn` to apply on their truth table
2) `POS_str`: _Returns POS `string`!_</br>
              => Given `name vector` & `Bits` of fcn results (as if from a truth table)</br>
              => Given `name vector` & `fcn` to apply on their truth table
3) `Kmap_print`: _Print K-Map to `stdout`!_</br>
              => Given `#` of bits (# K-Map vars) for blank K-Map</br>
              => Given `#` of bits (# K-Map vars) & `fcn` to apply on their truth table<br>
              => Given `Bits` of fcn results (as if from a truth table)
4) `TruthTable_print`: _Print Truth Table to `stdout`!_</br>
              => Given `#` of bits (# Truth Table vars) for blank Truth Table</br>
              => Given `#` of bits & a single (or `vector` of) `fcn` to apply on their truth table</br>
              => Given a `Bits` (or `vector<Bits>` matrix) of fcn results (as if from a truth table)

----------------------

## 4 Helper Functions:
* ***NOTE:** Further manipulates fcns returned by `SOP_str` & `POS_str`!*
1) `xor_replace`:    _Returns `string` w/ fcn's XOR/XNOR's using the `^` operator!_</br>
2) `compl_literals`: _Returns `string` w/ fcn's literals complemented!_</br>
3) `dual_fcn`:       _Returns `string` w/ fcn's dual!_</br>
4) `compl_fcn`:      _Composes "`compl_literals`" & "`dual_fcn`"!_

----------------------

## 2 Bonus Functions & DONT-CARE Universal Constant:
1) `BCD_decoder`: _Convert `Bits` to a decimal #!_</br>
   => Helps ID the current row in `fcn` applied to a truth table!</br>
   => Used for identifying "_Don't Care_" rows as needed!
2) `TruthTable_fcn`: _Given `#` of bits & `fcn`, returns `Bits` of fcn results across truth table!_
3) `DONT_CARE_BIT`: _Universal **constant**: return from `fcn`s at "_Don't Care_" rows!_

----------------------

## Sample User-Defined Function:
```c++
// Suppose variable names: A,B,C
// NOTE TYPE CONSTRAINTS: "Bit" return & "Bits" arg!
 
Bit some_func_name(Bits v) {
  return (v[0] ^ v[1]) & !v[2]; // A xor B and not C
}


// Same fcn if we don't care beyond row 3:
Bit some_func_name(Bits v) { 
  if(BCD_decoder(v) > 3)
    return DONT_CARE_BIT;
  return (v[0] ^ v[1]) & !v[2]; // A xor B and not C
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

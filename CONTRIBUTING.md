Contributing to MGLET-base
==========================

We encourages users to contribute to the development of the code. We are very
grateful for any contributions, from minor typos and bug fixes to new models
and features.

***Note that no contribution may be accepted unless the contributor agrees
with the license terms found in the LICENSE file in the repository's
top source directory.***


Quick Fortran style guide
-------------------------

This is a quick summary of the Fortran style that are used in MGLET-base:

### Indentation
Use 4 space characters for indentation. Consider an extra level of indentation
for continued lines to emphasize the continuation. Do not use tabs for any
purpose, not even within comments (tabs are not allowed in Fortran source
files).

### Capitalization
Use lowercase for all names constructed by the user (variables, subroutines,
functions, types, modules etc). Fortran keywords (such as `INTEGER`, `REAL`,
`SUBROUTINE`, `MODULE`, `CONTAINS`, `IF` etc.) are preferred in upper
case.

### Line lengths
Code lines are not allowed to exceed 80 characters. You are encouraged to limit
the line lengths of comments to the same limit, but certain exceptions
are be allowed, such as a comment with a URL reference that exceed 80
characters.

### `END`'s
Always use `END <construct>`, for example `END DO` and not `END` or `ENDDO`.

### Whitespace
* Always use a single whitespace after a comma, for example `bp(k, j, i)`.
* There shall be no trailing whitespace. Empty lines are not allowed to have
whitespace.
* Surround relational operators with a whitespace, e.g. `IF (myid == 0)`.
* In mathematical expressions: if operators with different priorities are
used, consider adding whitespace around the operators with the lowest
priority(ies).

### Naming conventions
* Modules end with the suffix `_mod`
* Derived types shall have the suffix `_t`.
* `kk, jj, ii` are used for the grid shape, do not use `kmx, jmx, imx`.
* Avoid variable names that are the same as a Fortran keyword, for
instance `if`, despite that this is allowed in the Fortran language.

### Order of variables in argument lists
It is nearly impossible to deduce a general, hard set of rules for the order
of variables in argument lists. However, please try to follow the following
guidelines:

1. Array shapes **must** appear before the array declarations themselves
(this is a hard rule), for example use `SUBROUTINE myroutine(kk, jj, ii, bp)`
where `bp` is an array of shape `(kk, jj, ii)`.

2. Indices appear in the order z-y-x, e.g. `kk, jj, ii`.

3. Physical coordinates and quantities appear in the order x-y-z, e.g.
`(dx, dy, dz)` and `(u, v, w)`.

4. Return values and results, declared as `INTENT(out)` or `INTENT(inout)`,
appear *before* pure inputs that are `INTENT(in)`, e.g.
`CALL get_mgdims(kk, jj, ii, igrid)` because `kk, jj, ii` are the results
of the call and hence appear first, and `igrid` is an input and appear
afterwards.

Where Fortran language rules specify certain rules that conflict these, the
standard will of course be the deciding factor. For example, optional arguments
have to appear at the end, independent of the `INTENT`.

### Relational operators
Always use the "modern" operators `==`, `/=`, `>`, `<`, `>=`, `<=` and not the
legacy `.EQ.`, `.NE.`, `.GT.`, `.LT.`, `.GE.`, `.LE.`. Relational operators
shall always be surrounded by a whitespace.

### Implicit typing and module visibility scopes
Always have `IMPLICIT NONE (type, external)` in every module to enforce
all interfaces and variables to be explicit and set the default module
visibility to `PRIVATE` except in special cases.

Always declare the `kind` of `REAL` and `INTEGER` whenever possible.

### `INTENT`'s
Always declare the `INTENT` of an argument whenever possible.

### When in doubt
Look at existing code and use your own judgment.


Timers numbering scheme
-----------------------

MGLET has a simple internal performance monitoring system using timers. To
avoid collisions and conflict of interest between various regions in the
code using the same indexes for the timers, this is the rule of assigning
numbers:

* 1-199: Used from main program (`mglet.F90`) and main time integration
(`timeloop_mod.F90`) and core functionality in core_mod
* 200-299: IB methods and blocking
* 300-399: Flow solver
* 400-499: Scalar solver
* 500-799: Reserved for future use
* 800-999: Plugins

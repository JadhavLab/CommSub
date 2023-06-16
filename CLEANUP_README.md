CLEANUP README

# Phase 1 : Section Collapse
 - Collapsing each major section (events, windowing) into function
 - Return simple struct from the sections instead of a slew of variables

# Phase 2 : Section cleanup

## Raw data structures

Too many lines and too unstructured. We can create a constructor method to
reduce the overhead.

## Result save section

- Some of the tables have confusing naming conventions. 
- Hide table munging and hashed saves into methods



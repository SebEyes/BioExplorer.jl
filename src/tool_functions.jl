"""
This function returns true if x is:
 - missing;
 - an empty string;
 - a string containing only spaces.
"""

_blank(x)=ismissing(x)||isempty(strip(string(x)))
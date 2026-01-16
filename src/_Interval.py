#!/usr/bin/python3
# -*- coding: utf-8 -*-

# NATURAL_INTERVAL:   [0, inf)
# POSITIVE_INTERVAL:  [1, inf)

class Interval:
    def __init__(self, _lower:int=1, _upper:int=float('Inf')):
        self.lower = _lower
        self.upper = _upper

    def __add__(self, other):
        if isinstance(other, Interval):
            return Interval(self.lower + other.lower, self.upper + other.upper)
        elif isinstance(other, int):
            return Interval(self.lower + other, self.upper + other)        
        else:
            return NotImplemented

    def __sub__(self, other:int):
        if isinstance(other, Interval):
            return Interval(self.lower - other.lower, self.upper - other.upper)
        elif isinstance(other, int):
            return Interval(self.lower - other, self.upper - other)             
        else:
            return NotImplemented

    # Join |
    def __or__(self, other):
        if isinstance(other, Interval):
            return Interval(min(other.lower, self.lower), max(other.upper, self.upper))
        else:
            return NotImplemented

    # Merge &
    def __and__(self, other):
        if isinstance(other, Interval):
            return Interval(max(other.lower, self.lower), min(other.upper, self.upper))
        else:
            return NotImplemented

    # obj * other
    def __mul__(self, other):
        if isinstance(other, Interval):
            return Interval(self.lower * other.lower, self.upper * other.upper)
        elif isinstance(other, int):
            return Interval(self.lower * other, self.upper * other)
        else:
            return NotImplemented
    
    # in
    def __contains__(self, other):
        if isinstance(other, Interval):
            return self.lower <= other.lower and other.upper <= self.upper
        elif isinstance(other, int):
            return self.lower <= other <= self.upper
        else:
            return NotImplemented

    def __str__(self):
        return f"[{self.lower}, {self.upper}]"

    def iter(self):
        yield from range(self.lower , self.upper + 1)
#!/usr/bin/env python3


start = 0
end   = 99
divisor=7
print("Print out all numbers from",start,"to",end, "that are not divisible by",divisor)
for solution in range(0, 99):
  if solution % 7 != 0:
    print(solution)

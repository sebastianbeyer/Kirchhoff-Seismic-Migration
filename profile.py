#!/usr/bin/env python


import pstats, cProfile


import sebitest

cProfile.runctx("sebitest.Testit()", globals(), locals(), "Profile.prof")

s = pstats.Stats("Profile.prof")
s.sort_stats("cumulative").print_stats()

These are reminder when using CMake.

1. For those who want to use CMake, please check whether the path to GLUT's include and lib are found. If not, please manually specify them.

2. Do not mix CMake with autoconf. Otherwise, since the generated makefiles might have different names (Makefile vs. makefile), it can cause unexpected issues.

Teng-Yok Lee
07/01/2011




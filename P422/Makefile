##############################################################
#               CMake Project Wrapper Makefile               #
############################################################## 

SHELL := /bin/bash
RM    := rm -rf

./build/Makefile:
	@ (cd build >/dev/null 2>&1 && cmake ..)

$(MAKECMDGOALS): ./build/Makefile
	@ $(MAKE) --no-print-directory -C build $(MAKECMDGOALS)

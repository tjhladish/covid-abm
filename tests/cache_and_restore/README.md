# How to run cache and restore test suite
---
This test suite will run the provided model to test whether the caching and restoration functionality is operating properly.

Copy the main.cpp file of choice (as well as create links to any other necessary file dependencies that main.cpp file requires).
However, please update the following in the main file to ensure proper test operation:
* Change #inlcude to use "cache_and_restore_simulator.h" instead of "simulator.h"
* Change main.cpp:FORECAST_DURATION to desired run length
* Change cache_and_restore_simulator.h:CACHE_POINT to desired day when cache will be made
* Turn off simulation database generation
* Turn off simvis figure generation

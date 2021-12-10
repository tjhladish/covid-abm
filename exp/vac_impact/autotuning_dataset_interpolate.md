# How to use the autotuning dataset interpolation system
---
After an autotuning run, the simulator will generate a file called "autotuning_dataset.csv".
If this filename is passed to the parameter `par->autotuning_dataset` and `par->behavioral_autotuning` is `false`, these tuned values will be read in to the simulator.

If the user desires to add interpolated values beyond the tuned dataset, follow the below steps:

1. Below the last row of tuned data, add the following line: "interpolate,,,"
2. Add dates and values to use as interpolation anchor points (you do not need to add values for columns 2 or 3)

Example:

2022-01-16,,,0.0011160  
2022-01-17,,,0.0022321  
2022-01-18,,,0.0033482  
interpolate,,,  
2022-01-19,,,0  
2022-03-23,,,0.5  
2022-03-29,,,0.5  
2022-04-07,,,0.54  
2022-04-11,,,0.54  
2022-04-20,,,0.5  
2022-04-28,,,0.5  
2022-05-31,,,0.0727539  

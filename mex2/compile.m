clear all;
clc;
mex "MATLAB_interface.cpp" -l"cudart" -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.1\lib\x64" -I"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.1\include"
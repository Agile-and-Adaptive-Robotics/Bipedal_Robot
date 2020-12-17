clc                 %Command Line Clear
clear               %Clear the workspace of stored variables
close all           %Close all open figures

ka = 2.7*260^(-0.265)

Se = ka*1*0.85*0.89*0.702*100

sig = 15000/(1.5*0.25)

dH = 0.5/2

Kf = 1 + 0.95*(2.44-1)

sigMax = Kf * sig

a = (0.76*260)^2/25

b = -1/3*log10(0.76*260/25)

N = (sigMax/a)^(1/b)

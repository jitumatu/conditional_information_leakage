# conditional_information_leakage

This repository gives source codes for computing conditional information leakage.
This topic was discussed in:

[1] U.Michiwaki and Y. Jitsumatsu, "A Definition of Information Leakage Given a Wiretapped signal," the 41st Sympo. on Inform. Theory and Its Applications (SITA2018) (in Japanese) 

[2] Y. Jitsumatsu, U.Michiwaki, and Y.Oohama, "Conditional Information Leakage Given Eavesdropper's Received Signals in Wiretap Cahnnels," submitted to IEICE Transactions on Fundamentals. 

#########

naive_comp_CIL_OMP.c :  the navie computation of the conditional information leakage.  
fast_comp_CIL_OMP.c :  a fast computation of the conditional information leakage.

naive_comp_CIL_OMP-B.c : 16 cases of (n,m) are computed by the navie method.  
                         n is the code length, m is the lenght of the message . 
fast_comp_CIL_OMP-B.c :  16 cases of (n,m) are computed by the fast method. 

I checked that these files can be compiled with gcc version 7.4.0 on Ubuntu 16.04.

 

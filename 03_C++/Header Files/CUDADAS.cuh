#include        <iostream>
#include        <cmath>
#include        "Header Files/DAS.hpp"

void            indexgpuwrapper(Dataset data2, IndexMatrix indmat, double* reconptr);
void            mmultgpuwrapper(Dataset data2, IndexMatrix indmat, float* cpureconptr);
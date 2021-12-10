#!/bin/bash

ASPERA_SCP_PASS=U2kSkXcp /gpfs/home/aft19qdu/.aspera/connect/bin/ascp -P33001  -O33001 -QT -l300M -L- ./S0904-S09-L1* ega-box-175@fasp.ega.ebi.ac.uk:/.

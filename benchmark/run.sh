#!/bin/bash
model="CVM_1d.mdl"
period="period.txt"
/home/yma/rbh/PROGRAMS.330/bin/sprep96 -M $model -L -R -PARR $period -NMOD 1
#/home/yma/rbh/PROGRAMS.330/bin/sprep96 -M $model -R -PARR $period -NMOD 1
/home/yma/rbh/PROGRAMS.330/bin/sdisp96 -v
# output the partial derivative
/home/yma/rbh/PROGRAMS.330/bin/slegn96 -DER
/home/yma/rbh/PROGRAMS.330/bin/sregn96 -DER
/home/yma/rbh/PROGRAMS.330/bin/sdpder96 -L -TXT
/home/yma/rbh/PROGRAMS.330/bin/sdpder96 -R -TXT
/home/yma/rbh/PROGRAMS.330/bin/slegn96 
/home/yma/rbh/PROGRAMS.330/bin/sregn96
/home/yma/rbh/PROGRAMS.330/bin/sdpegn96 -L -U -C -ASC
/home/yma/rbh/PROGRAMS.330/bin/sdpegn96 -R -U -C -ASC

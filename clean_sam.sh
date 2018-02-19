#!/usr/bin/env bash
set –euo pipefail
for sam in *.sam
do
    cut -f3,4 $sam | sort -k2,2n > $sam.clean
done


#!/bin/bash
for i in {1..60}; do
    tail -l '1221array_$i.log'
    echo "dsa"
done
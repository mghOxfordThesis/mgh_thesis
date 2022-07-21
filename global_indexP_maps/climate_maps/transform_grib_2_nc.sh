#!/bin/bash

for i in *.grib; do
    cdo -f nc copy "$i" "$i.nc"
done

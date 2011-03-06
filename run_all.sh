#!/bin/bash
CURDIR=`pwd`

$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 0
$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 1
$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 2
$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 3
$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 4
$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 5
$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 6
$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 7
$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 8
$CURDIR/riciandenoise3d -m 60 -n 60 -p 60 -i f.bin -o denoise.out -b 9

#!/bin/bash
rm -rf result/*_0.h
path=result
files=$(ls $path)
for filename in $files
do
  txtfile_0=`echo $filename | sed 's/.h$/\_0.h/'`
  echo txtfile_0
  ../contrib/utilities/wrapcomments.py  $path/$filename >> $path/$txtfile_0
done
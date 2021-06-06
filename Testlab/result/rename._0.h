#!/bin/bash
for filename in *.h; do
        newfile=`echo $filename | sed 's/\_pre//'`
        mv $filename $newfile
done

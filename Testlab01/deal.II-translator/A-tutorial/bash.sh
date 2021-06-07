#!/bin/bash
for filename in *.temp; do
        newfile=`echo $filename | sed 's/.temp/.h/'`
        mv $filename $newfile
done

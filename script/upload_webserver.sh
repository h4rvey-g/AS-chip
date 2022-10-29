#!/bin/bash
# get username
username=$(whoami)
mkdir -p /data0/shared/XuLab/webserver/$username
webpath=/data0/shared/XuLab/webserver/$username
# create a for loop to create symbolic links for input files to webpath, rename as current time
for i in $#; do
    ln -s $i $webpath/$(basename $i)_$(date +%Y%m%d%H%M)
    # attach the link after http://10.181.6.195/, and echo the link to the screen
    echo "Filename: $i"
    # echo $i in blue color
    echo -e

    echo http://10.181.6.195/webserver/$username/$(basename $i)_$(date +%Y%m%d%H%M)
done

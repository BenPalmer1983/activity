#!/bin/bash 
cd /data/activity/
tar cvzf activityData.tar.gz * --directory="/data/activity" --exclude="/data/activity/tar.sh"  \
 --exclude="/data/activity/activityData.tar.gz"
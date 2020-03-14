#!/bin/bash 
tar cvzf activity.tar.gz * \
--exclude="*/.git" \
--exclude="make.sh" \
--exclude="*.docx" \
--exclude="*activity.tar.gz" \
--exclude="tar.sh" \
--exclude="git.sh" \
--exclude="untar.sh" \
--exclude="bin/*.x"

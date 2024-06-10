#!/bin/bash
#
# USE:
# 
#    startscanjobs.sh  [ /path/to/scan ]  [ /path/to/scanjobs.sh ]
#
# Attempts to log into neon and start scanjobs.sh (or other specified script,
# defaulting to ~/scanjobs.sh), which scans the specified directory on the 
# remote host, (default: ~/jobs)
#  
#

# Christopher Kovach 2016


loginnode="neon.hpc.uiowa.edu"
default_scandir='~/jobs';  # Remote directory to scan
scan_exec=${2:-'~/scanjobs.sh'}; # Remote location of the scanjobsd script.

# Command to be executed
#command="nohup $scan_exec ${1:-$default_scandir}  &"

printf "\nEnter neon login (usually same as hawkID): " 
read loginid

ssh -p 40 -T $loginid@$loginnode << ENDHERE
  nohup $scan_exec ${1:-$default_scandir}  &
  exit
ENDHERE



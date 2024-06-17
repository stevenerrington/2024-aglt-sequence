#!/bin/bash

# Copy files between an xneon submit folder and some other freely accessible 
# folder so that requests can be submitted to neon without having to log in.
#
# This script should run as a background process on some server with access
# to both folders. It uses the unison file synchronization utility, which 
# must be installed on the sever.


submit_from=${1:-~/hbrl2/HBRL_Upload/xneon_submit}
submit_to=${2:-~/neon/xneon_submit/}

run_every=5; #Run at this interval in seconds

unison $submit_from $submit_to -batch -times -force newer -ignore "Path .bridge"

touch $submit_from/.bridge

sleep $run_every && exec $0




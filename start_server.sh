#! /bin/bash

kill $(pgrep -f auto_clone_server)
rm -rf ./output_folders/*
nohup python auto_clone_server.py 2>&1 1> server_log.txt &
disown
 

#! /bin/bash

nohup python auto_clone_server.py 2>&1 1> server_log.txt &
disown
 

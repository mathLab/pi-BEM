#!/bin/bash

exec 2> /dev/null
EI=0

# If you want to simulate what will be run, just pass echo on the command
# line. 
PREPEND_COMMANDS="$1"

# The script runs all commands in the before_install and script sections 
# of the .travis.yml file, one time for each env variable set

while ENVSET=`cat .travis.yml | shyaml get-value env.$EI`; do 
    echo $ENVSET
    eval "$ENVSET"
    ((EI++))
    for section in before_install script; do 
	echo ============================================================
	echo ============================================================
	echo "Current phase: $section"
	echo ============================================================
	echo ============================================================
	BI=0
	while CMD=`cat .travis.yml | shyaml get-value $section.$BI`; do
	    ((BI++))
	    echo "Executing: $CMD"
	    echo ============================================================
	    if ! eval "$PREPEND_COMMANDS $CMD"; then
		echo ************************************************************
		echo ************************************************************
		echo "FAILED"
		echo ************************************************************
		echo ************************************************************
		exit 1
	    else 
		echo ============================================================
		echo "SUCCESS"
		echo ============================================================
	    fi
	done
    done
done
    

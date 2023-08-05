#!/bin/sh
./scripts/indent
git diff
git diff-files --quiet 

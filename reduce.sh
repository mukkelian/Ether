#!/bin/bash
for f in *\ *.dat ; do mv "$f" "${f// /}" ; done

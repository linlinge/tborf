#!/bin/sh
mv Result/ ..
mv build/ ..
git add .
git commit -m "update"
git push origin master
mv ../Result .
mv ../build .

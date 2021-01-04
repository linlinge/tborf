#!/bin/sh
mv Result/ ..
mv build/ ..
git add .
git commit -m "update"
git push origin
mv ../Result .
mv ../build .

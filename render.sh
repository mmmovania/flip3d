#!/bin/sh
i=$1
while [ $i -le $2 ]; do
f=render/sunflow/"$i"_scene.sc
if [ -s "$f" ]; then
echo "Starting scene $i..."
java -Xmx2G -jar sunflow/sunflow.jar -nogui -o render/sunflow/"$i"_scene.png render/sunflow/"$i"_scene.sc
i=`expr $i + 1`
else
sleep 10
fi
done

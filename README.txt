Sunflow How To Render:

Please download sunflow-bin-v0.07.2.zip from

http://sunflow.sourceforge.net/

Then, expand here.
To render sequential file, run:

for i in {1..400}; 
do 
	java -Xmx2G -jar sunflow/sunflow.jar -nogui -o render/sunflow/"$i"_scene.png render/sunflow/"$i"_scene.sc; 
done
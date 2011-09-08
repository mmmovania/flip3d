How To Compile And Run:

Just do:

> make run

That's All.

How To Render With Sunflow:

Please download sunflow-bin-v0.07.2.zip by running

> wget http://prdownloads.sourceforge.net/sunflow/sunflow-bin-v0.07.2.zip

Then, expand here.
To render sequential file from frame 1 to e.g. 40, run:

> ./render.sh 1 40

How To Increase Domain Size:

Edit "#define N 32" in flip3d.cpp to the number whatever you want.
( The number doesn't have to be a power of two, e.g. 150 is OK )
Watch for the file size and compile codes appropriately otherwise files won't be
successfully written. Don't forget to increase the memory size in the variable in 
render.sh if you plan to render with Sunflow renderer.

Enjoy !
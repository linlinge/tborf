str_cmd=strcat('"/home/llg/workspace_cmake/code-of-MODM/matlab_scripts/build/MGS" -i /home/llg/dataset/3Dlib/cat_colmap.ply -o /home/llg/workspace_cmake/code-of-MODM/matlab_scripts/1.csv -M OGEM')
system(str_cmd)
og=csvread('1.csv')

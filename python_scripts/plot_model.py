import open3d as o3d 
ply=o3d.io.read_point_cloud("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")
o3d.visualization.draw_geometries([ply])
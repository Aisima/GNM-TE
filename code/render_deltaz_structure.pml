load cftr_9mxl_deltaz.pdb, cftr
hide everything
bg_color white
set ray_opaque_background, 1
set antialias, 2
set cartoon_fancy_helices, 1
set specular, 0.12
set depth_cue, 0
set ray_shadows, 0
select nodz, cftr and b > 90
show cartoon, cftr and resi 390-650+1170-1480
color grey80, nodz
spectrum b, orange_white_marine, cftr and not nodz, -3, 3
set cartoon_transparency, 0.30
select kv, cftr and resi 455+549+551+560+1244+1251+1303+1349 and name CA
show spheres, kv
set sphere_scale, 1.5, kv
orient cftr and resi 390-650+1170-1480
turn x, 90
ray 1500, 1250
png s_nbd.png, dpi=300
# identical view, spheres only, unique flat colours
hide everything
show spheres, kv
set sphere_scale, 1.5, kv
set ray_trace_mode, 0
set specular, 0
set ambient, 1.0
set direct, 0.0
set reflect, 0.0
set_color c455, [1.0,0.0,0.0]
set_color c549, [0.0,1.0,0.0]
set_color c551, [0.0,0.0,1.0]
set_color c560, [1.0,1.0,0.0]
set_color c1244, [1.0,0.0,1.0]
set_color c1251, [0.0,1.0,1.0]
set_color c1303, [0.5,0.0,0.0]
set_color c1349, [0.0,0.5,0.0]
color c455, resi 455 and name CA
color c549, resi 549 and name CA
color c551, resi 551 and name CA
color c560, resi 560 and name CA
color c1244, resi 1244 and name CA
color c1251, resi 1251 and name CA
color c1303, resi 1303 and name CA
color c1349, resi 1349 and name CA
ray 1500, 1250
png s_key.png, dpi=300

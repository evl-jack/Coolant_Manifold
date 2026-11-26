from pathlib import Path
import json
import numpy as np
from shapely.geometry import shape, LineString
from skimage.measure import marching_cubes
R=Path(__file__).resolve().parent; out=R/'assets/viewer';out.mkdir(exist_ok=True)
v=R/'evidence/verification/input';n=R/'evidence/nonlinear/stokes_run/stage_03_beta_8/selected'
paths={'original':v/'stage_03_beta_8/original/geometry.json','raw':v/'stage_03_beta_8/selected/geometry.json','smoothed':v/'smoothed_geometry.json','stokes':n/'geometry.json'}
def read(p):return json.loads(p.read_text())
for name,p in paths.items():
 d=read(p);poly=shape(d.get('geometry',d));
 # Sample actual horizontal branch sections; preserve width, not centerline asymmetry.
 ys=np.linspace(4.00001,6.99999,601);r=[]
 for y in ys:
  sect=poly.intersection(LineString([(7,y),(13,y)]));r.append(sect.length/2)
 step=.065; x=np.arange(-.15,20.16,step); y=np.arange(-.15,7.16,step); z=np.arange(-2.2,2.21,step)
 X=x[:,None,None];Y=y[None,:,None];Z=z[None,None,:]
 main=np.maximum(np.sqrt((Y-2)**2+Z**2)-2,np.maximum(-X,X-20))
 # Extend base radius into the main cylinder to form a single joined fluid-volume surface.
 rad=np.interp(y,ys,r,left=r[0],right=.5)[None,:,None]
 side=np.maximum(np.sqrt((X-10)**2+Z**2)-rad,np.maximum(2-Y,Y-7))
 sdf=np.minimum(main,side).astype('float32')
 verts,faces,_,_=marching_cubes(sdf,level=0,spacing=(step,)*3);verts+=np.array([x[0],y[0],z[0]])
 # Remove cap faces at physical openings. Display fluid-volume skin, not a wall-thickness model.
 f=verts[faces];keep=~((np.max(np.abs(f[:,:,0]),axis=1)<.075)|(np.max(np.abs(f[:,:,0]-20),axis=1)<.075)|(np.max(np.abs(f[:,:,1]-7),axis=1)<.075));faces=faces[keep]
 np.asarray(verts,dtype='<f4').tofile(out/f'{name}-cyl-vertices.bin');np.asarray(faces,dtype='<u4').tofile(out/f'{name}-cyl-faces.bin')
 (out/f'{name}-geometry.json').write_text(json.dumps({'contour':list(poly.exterior.coords),'profile_y_mm':ys.tolist(),'radius_mm':r,'cyl_mesh_step_mm':step,'interpretation':'sampled equivalent circular width and union; fluid volume surface, open port caps; not CFD'}))
 print(name,verts.shape,faces.shape)
stages=[]
for i,folder in enumerate(['stage_01_beta_0','stage_02_beta_2','stage_03_beta_8']):
 for label in (['original','selected'] if i==0 else ['selected']):
  m=read(v/folder/label/'metrics.json')['results'];stages.append({'label':'Original' if label=='original' else f'β = {[0,2,8][i]}','power':m['power_w_per_m']*1e6,'side':m['side_fraction']*100,'image':f'evidence/pilot/local/{folder}/images/'+('original_physical_geometry.png' if label=='original' else 'optimized_physical_geometry.png'),'scope':'immersed checkpoint'})
(out/'stages.json').write_text(json.dumps(stages))

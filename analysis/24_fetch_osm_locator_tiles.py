"""Cache the OpenStreetMap standard tiles behind the Zimbabwe locator panel.

Tiles are fetched once at zoom 7 and reused, following the OSM tile usage policy
(identifying user agent, no bulk download). Map data © OpenStreetMap contributors,
ODbL; the attribution must accompany any figure that shows these tiles.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,math,time,urllib.request
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/public/osm_tiles_z7';OUT.mkdir(parents=True,exist_ok=True)
Z=7;WEST,EAST,SOUTH,NORTH=24.5,34.0,-23.5,-14.0
def tile(lon,lat):
    n=2**Z;r=math.radians(lat)
    return int((lon+180)/360*n),int((1-math.log(math.tan(r)+1/math.cos(r))/math.pi)/2*n)
x0,y0=tile(WEST,NORTH);x1,y1=tile(EAST,SOUTH)
records=[]
for x in range(x0,x1+1):
    for y in range(y0,y1+1):
        f=OUT/f'{Z}_{x}_{y}.png'
        if not f.exists():
            req=urllib.request.Request(f'https://tile.openstreetmap.org/{Z}/{x}/{y}.png',
                headers={'User-Agent':'MidlandsPM25-locator/1.0 (academic figure, single cached download)'})
            f.write_bytes(urllib.request.urlopen(req,timeout=60).read());time.sleep(.5)
        records.append(dict(file=f.name,z=Z,x=x,y=y,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
meta=dict(source='https://tile.openstreetmap.org/{z}/{x}/{y}.png',style='OpenStreetMap Standard',
    attribution='© OpenStreetMap contributors',licence='ODbL 1.0 (data); tiles CC BY-SA 2.0',
    retrieved_utc=datetime.now(timezone.utc).isoformat(),zoom=Z,x_range=[x0,x1],y_range=[y0,y1],tiles=records)
(OUT/'osm_tiles_metadata.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
print(len(records),'tiles',x0,x1,y0,y1)

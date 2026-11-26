import {buildPipe} from "./pipe.js";

const $=id=>document.getElementById(id);
const fetchJSON=async p=>{const r=await fetch(p);if(!r.ok)throw Error(p+': '+r.status);return r.json()};
const fetchBin=async(p,T)=>{const r=await fetch(p);if(!r.ok)throw Error(p+': '+r.status);return new T(await r.arrayBuffer())};
const names={raw:'Navier–Stokes optimized',stokes:'Stokes optimized',smoothed:'Smoothed (provisional)',original:'Original'};
const ISOMETRIC = {yaw:-.53,pitch:.43,zoom:1,panX:0,panY:0};
const camera = {...ISOMETRIC};
const datasets={},cylinders={};let result,stages,pipeData,panels=[],timer,limits={speed:1,pressure:1},lastStage=0,refreshToken=0,refreshing=true;
function cross(a,b){return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]}function sub(a,b){return a.map((v,i)=>v-b[i])}function norm(a){const d=Math.hypot(...a)||1;return a.map(v=>v/d)}
function rot(p){let [x,y,z]=p;let c=Math.cos(camera.yaw),s=Math.sin(camera.yaw);let xx=c*x+s*z,zz=-s*x+c*z;c=Math.cos(camera.pitch);s=Math.sin(camera.pitch);return [xx,c*y-s*zz,s*y+c*zz]}
function basis(b){let [a,c,d]=b;return [a*(2*a-1),c*(2*c-1),d*(2*d-1),4*a*c,4*c*d,4*d*a]}
function evalFE(d,i,b){const f=basis(b),vi=d.order[i],p=d.triangles;let u=0,v=0,pr=0;for(let k=0;k<6;k++){u+=f[k]*d.velocity_m_per_s[2*vi[k]];v+=f[k]*d.velocity_m_per_s[2*vi[k]+1]}for(let k=0;k<3;k++)pr+=b[k]*d.pressure_pa[p[3*i+k]];return [u,v,pr]}
function xy(d,i,b){let x=0,y=0;for(let k=0;k<3;k++){const p=d.triangles[3*i+k];x+=b[k]*d.pressure_points_mm[2*p];y+=b[k]*d.pressure_points_mm[2*p+1]}return [x,y]}
async function loadField(name){if(datasets[name])return datasets[name];let root=`assets/fields/${name}/`,meta=await fetchJSON(root+'metadata.json'),d={meta};await Promise.all(Object.entries(meta.arrays).map(async([k,a])=>d[k]=await fetchBin(root+a.file,a.dtype==='uint32'?Uint32Array:Float32Array)));const count=d.triangles.length/3;d.order=[];const verts=[],lines=[],edges=new Map();let maxS=0,maxP=0;
 for(let i=0;i<count;i++){const points=[0,1,2].map(k=>xy(d,i,[0,1,2].map(j=>j===k?1:0)));let target=[...points,...[[0,1],[1,2],[2,0]].map(([a,b])=>points[a].map((x,j)=>(x+points[b][j])/2))];let raw=Array.from(d.velocity_dofs.slice(i*6,i*6+6));d.order[i]=target.map(([x,y])=>raw.reduce((a,b)=>Math.hypot(d.velocity_points_mm[2*a]-x,d.velocity_points_mm[2*a+1]-y)<Math.hypot(d.velocity_points_mm[2*b]-x,d.velocity_points_mm[2*b+1]-y)?a:b));if(new Set(d.order[i]).size!==6)throw Error('P2 node order ambiguity');
 const add=(b,z,normal)=>{let [x,y]=xy(d,i,b),[u,v,p]=evalFE(d,i,b),s=Math.hypot(u,v);maxS=Math.max(maxS,s);maxP=Math.max(maxP,Math.abs(p));verts.push(x-10,y-3.5,z,...normal,s,p)};
 const N=3, B=(a,b)=>[1-(a+b)/N,a/N,b/N];for(let a=0;a<N;a++)for(let b=0;b<N-a;b++){let t=[B(a,b),B(a+1,b),B(a,b+1)];for(let z of [.3,-.3])t.forEach(q=>add(q,z,[0,0,z>0?1:-1]));if(a+b<N-1){t=[B(a+1,b),B(a+1,b+1),B(a,b+1)];for(let z of [.3,-.3])t.forEach(q=>add(q,z,[0,0,z>0?1:-1]));}}
 for(let [a,b] of [[0,1],[1,2],[2,0]]){let ia=d.triangles[i*3+a],ib=d.triangles[i*3+b],key=[ia,ib].sort((x,y)=>x-y).join(':');if(edges.has(key))edges.delete(key);else edges.set(key,[i,a,b]);for(let k of [a,b]){let [x,y]=points[k];lines.push(x-10,y-3.5,.304,0,0,1,0,0)}}
 }
 for(let [i,a,b] of edges.values()){let ba=[0,0,0],bb=[0,0,0];ba[a]=1;bb[b]=1;let A=xy(d,i,ba),B=xy(d,i,bb);let normal=norm([B[1]-A[1],A[0]-B[0],0]);for(let [p,z] of [[ba,.3],[bb,.3],[bb,-.3],[ba,.3],[bb,-.3],[ba,-.3]]){let [x,y]=xy(d,i,p),[u,v,pr]=evalFE(d,i,p);verts.push(x-10,y-3.5,z,...normal,Math.hypot(u,v),pr)}}
 d.vertices=new Float32Array(verts);d.lines=new Float32Array(lines);d.maxS=maxS;d.maxP=maxP;datasets[name]=d;return d}
function loadCylinder(name) {
    return cylinders[name] ??= {vertices: buildPipe(pipeData.profiles[name].sections_mm, pipeData.wall_mm)};
}

function preparePlane(field) {
    if (field.plane) return;
    const plane = [];
    for (let i = 0; i < field.vertices.length; i += 24) {
        if (field.vertices[i+5] !== 1 || field.vertices[i+2] < 0) continue;
        for (let k = i; k < i + 24; k += 8) {
            plane.push(field.vertices[k], field.vertices[k+1], 0, 0, 0, 1, field.vertices[k+6], field.vertices[k+7]);
        }
    }
    field.plane = new Float32Array(plane);
}

function sampleVectors(field) {
    if (field.arrows) return field.arrows;
    // Keep one real P2 sample near the center of each grid tile.
    const samples = new Map();
    const thirds = [1/3, 1/3, 1/3];
    for (let i = 0; i < field.triangles.length / 3; i++) {
        const [x,y] = xy(field, i, thirds);
        const ix = Math.floor(x / 0.8), iy = Math.floor(y / 0.65);
        const key = ix + ":" + iy;
        const distance = (x - (ix + 0.5)*0.8)**2 + (y - (iy + 0.5)*0.65)**2;
        if (!samples.has(key) || distance < samples.get(key).distance) samples.set(key,{x,y,i,distance});
    }
    const arrows = [];
    for (const {x,y,i} of samples.values()) {
        const [u,v] = evalFE(field,i,thirds);
        const magnitude = Math.hypot(u,v);
        if (magnitude < 1e-7) continue;
        const dx = u / magnitude, dy = v / magnitude;
        const length = 22 * magnitude;
        const a = [x - 10, y - 3.5], b = [a[0]+dx*length, a[1]+dy*length];
        const tip = 0.3 * length;
        const c = [b[0] - tip*(dx - 0.48*dy), b[1] - tip*(dy + 0.48*dx)];
        const d = [b[0] - tip*(dx + 0.48*dy), b[1] - tip*(dy - 0.48*dx)];
        for (const [px,py] of [a,b,b,c,b,d]) arrows.push(px,py,0.008,0,0,1,0,0);
    }
    return field.arrows = new Float32Array(arrows);
}

const VS=`attribute vec3 aPos;attribute vec3 aNormal;attribute vec2 aField;uniform vec4 uRot;uniform vec2 uScale;uniform vec2 uPan;varying vec2 field;varying vec3 normal;varying float sectionZ;vec3 rotate(vec3 p){vec3 q=vec3(uRot.x*p.x+uRot.y*p.z,p.y,-uRot.y*p.x+uRot.x*p.z);return vec3(q.x,uRot.z*q.y-uRot.w*q.z,uRot.w*q.y+uRot.z*q.z);}void main(){vec3 p=rotate(aPos);gl_Position=vec4(p.xy*uScale+uPan,-p.z/30.0,1.0);normal=rotate(aNormal);field=aField;sectionZ=aPos.z;}`;
const FS = `
precision mediump float;
varying vec2 field;
varying vec3 normal;
varying float sectionZ;
uniform int uMode;
uniform vec2 uLimits;
uniform float uSection;
void main() {
    if (uMode == 4 && sectionZ > uSection) discard;
    vec3 c = vec3(.49,.65,.68);
    if (uMode == 1) {
        float t = clamp(field.x/uLimits.x,0.0,1.0);
        c = mix(vec3(.91,.96,1.),vec3(.05,.23,.55),t);
    }
    if (uMode == 2) {
        float t = clamp(field.y/uLimits.y,-1.,1.);
        c = t < 0. ? mix(vec3(.96),vec3(.13,.34,.65),-t)
                   : mix(vec3(.96),vec3(.67,.18,.15),t);
    }
    if (uMode == 3) c = vec3(.035,.11,.17);
    if (uMode == 4) {
        vec3 n = normalize(normal);
        vec3 light = normalize(vec3(-.3,.65,1.));
        float diffuse = abs(dot(n,light));
        float shine = pow(abs(dot(n,normalize(vec3(-.15,.36,1.)))),16.);
        vec3 metal = field.x < .5 ? vec3(.44,.64,.72) : vec3(.16,.28,.32);
        if (field.x > 1.5) metal = vec3(.77,.83,.86);
        c = metal * (.40 + .60 * diffuse) + .30 * shine;
    }
    gl_FragColor = vec4(c,1.);
}`;

function init(canvas){let gl=canvas.getContext('webgl',{antialias:true});if(!gl)throw Error('WebGL unavailable');let prog=gl.createProgram();for(let [type,src] of [[gl.VERTEX_SHADER,VS],[gl.FRAGMENT_SHADER,FS]]){let sh=gl.createShader(type);gl.shaderSource(sh,src);gl.compileShader(sh);if(!gl.getShaderParameter(sh,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(sh));gl.attachShader(prog,sh)}gl.linkProgram(prog);if(!gl.getProgramParameter(prog,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(prog));gl.useProgram(prog);gl.enable(gl.DEPTH_TEST);let buf=gl.createBuffer(), line=gl.createBuffer(),p={canvas,gl,prog,buf,line,pipe:gl.createBuffer(),arrows:gl.createBuffer(),labels:canvas.parentElement.querySelectorAll('[data-port]'),name:null,data:null};let dragging=false,last=null,moved=false;canvas.addEventListener('pointerdown',e=>{dragging=true;last=[e.clientX,e.clientY];moved=false;canvas.setPointerCapture(e.pointerId)});canvas.addEventListener('pointermove',e=>{if(!dragging)return;let dx=e.clientX-last[0],dy=e.clientY-last[1];if(Math.abs(dx)+Math.abs(dy)>1)moved=true;if(e.shiftKey||$('pan').getAttribute('aria-pressed')==='true'){camera.panX+=dx/canvas.clientWidth*2;camera.panY-=dy/canvas.clientHeight*2}else{camera.yaw+=dx*.012;camera.pitch=Math.max(-1.3,Math.min(1.3,camera.pitch+dy*.012))}last=[e.clientX,e.clientY];drawAll()});canvas.addEventListener('pointerup',e=>{dragging=false;if(!moved)inspect(p,e)});canvas.addEventListener('pointercancel',()=>{dragging=false;last=null});canvas.addEventListener('wheel',e=>{e.preventDefault();camera.zoom=Math.max(.55,Math.min(4,camera.zoom*Math.exp(-e.deltaY*.001)));drawAll()},{passive:false});return p}
function upload(panel, buffer, data) {
    const gl = panel.gl;
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
    return data.length / 8;
}

function setData(p, name, pipe, fields, mode) {
    p.name = name;
    p.data = fields;
    p.pipeCount = pipe ? upload(p, p.pipe, pipe.vertices) : 0;
    p.count = p.lineCount = p.arrowCount = 0;
    if (!fields) return;
    const cutaway = mode === "cutaway";
    if (cutaway) preparePlane(fields);
    p.count = upload(p, p.buf, cutaway ? fields.plane : fields.vertices);
    if (cutaway) {
        if (!fields.centerLines) {
            fields.centerLines = fields.lines.slice();
            for (let i = 2; i < fields.centerLines.length; i += 8) fields.centerLines[i] = .006;
        }
        p.lineCount = upload(p, p.line, fields.centerLines);
    } else {
        p.lineCount = upload(p, p.line, fields.lines);
    }
    const vectors = sampleVectors(fields);
    if (cutaway) p.arrowCount = upload(p,p.arrows,vectors);
    else {
        if (!fields.extrudedArrows) {
            fields.extrudedArrows = vectors.slice();
            for (let i = 2; i < fields.extrudedArrows.length; i += 8) fields.extrudedArrows[i] += .3;
        }
        p.arrowCount = upload(p,p.arrows,fields.extrudedArrows);
    }
}

function scales(p) {
    const aspect = p.canvas.width / p.canvas.height;
    const x = Math.min(1/11, 1/(5.0*aspect)) * camera.zoom;
    return [x,x*aspect];
}

function draw(p) {
    const gl = p.gl, canvas = p.canvas;
    const dpr = Math.min(devicePixelRatio, 2);
    const width = Math.max(1,Math.round(canvas.clientWidth*dpr)), height = Math.max(1,Math.round(canvas.clientHeight*dpr));
    if (canvas.width !== width || canvas.height !== height) { canvas.width = width; canvas.height = height; }
    gl.viewport(0,0,width,height);
    gl.clearColor(.965,.98,.985,1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    const locate = name => gl.getUniformLocation(p.prog,name);
    const [sx,sy] = scales(p);
    gl.uniform4f(locate("uRot"),Math.cos(camera.yaw),Math.sin(camera.yaw),Math.cos(camera.pitch),Math.sin(camera.pitch));
    gl.uniform2f(locate("uScale"),sx,sy);
    gl.uniform2f(locate("uPan"),camera.panX,camera.panY);
    gl.uniform2f(locate("uLimits"),limits.speed,limits.pressure);
    const mode = $("mode").value, field = $("field").value, cutaway = mode === "cutaway";
    gl.uniform1f(locate("uSection"),cutaway ? Number($("section-depth").value) : 3);
    const render = (buffer,count,type,shader) => {
        if (!count) return;
        gl.uniform1i(locate("uMode"),shader);
        gl.bindBuffer(gl.ARRAY_BUFFER,buffer);
        for (const [name,size,offset] of [["aPos",3,0],["aNormal",3,12],["aField",2,24]]) {
            const attribute = gl.getAttribLocation(p.prog,name);
            gl.enableVertexAttribArray(attribute);
            gl.vertexAttribPointer(attribute,size,gl.FLOAT,false,32,offset);
        }
        gl.drawArrays(type,0,count);
    };
    render(p.pipe,p.pipeCount,gl.TRIANGLES,4);
    render(p.buf,p.count,gl.TRIANGLES,field==="speed" ? 1 : field==="pressure" ? 2 : 0);
    if (field==="mesh") render(p.line,p.lineCount,gl.LINES,3);
    if ($("vectors").checked) render(p.arrows,p.arrowCount,gl.LINES,3);
    for (const label of p.labels) {
        const [x,y] = rot(({inlet:[-10,-1.5,0],outlet:[10,-1.5,0],side:[0,3.5,0]})[label.dataset.port]);
        label.style.left = Math.min(86,Math.max(14,(x*sx+camera.panX+1)*50)) + "%";
        const offset = label.dataset.port === "side" ? -7 : 13;
        label.style.top = Math.min(94,Math.max(6,(1-y*sy-camera.panY)*50+offset)) + "%";
    }
}
function drawAll(){if(!refreshing)panels.forEach(draw)}
function inspect(p,e){if(refreshing||!p.data)return;let d=p.data,r=p.canvas.getBoundingClientRect(),x=(e.clientX-r.left)/r.width*2-1,y=1-(e.clientY-r.top)/r.height*2,[sx,sy]=scales(p);let best=null;for(let i=0;i<d.triangles.length/3;i++)for(let z of ($('mode').value==='cutaway'?[0]:[.3,-.3])){let q=[0,1,2].map(k=>{let [a,b]=xy(d,i,[0,1,2].map(j=>j===k?1:0));let v=rot([a-10,b-3.5,z]);return [v[0]*sx+camera.panX,v[1]*sy+camera.panY,v[2]]});let [a,b,c]=q,den=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1]);if(Math.abs(den)<1e-10)continue;let l=((b[1]-c[1])*(x-c[0])+(c[0]-b[0])*(y-c[1]))/den,m=((c[1]-a[1])*(x-c[0])+(a[0]-c[0])*(y-c[1]))/den,n=1-l-m;if(Math.min(l,m,n)<-1e-6)continue;let dep=l*a[2]+m*b[2]+n*c[2];if(!best||dep>best.dep)best={i,b:[l,m,n],dep}}if(!best){$('inspect').textContent='Click a visible planar face to inspect an element.';return}let [u,v,pr]=evalFE(d,best.i,best.b),[X,Y]=xy(d,best.i,best.b);$('inspect').textContent=`${names[p.name]} · element ${best.i} · (${X.toFixed(3)}, ${Y.toFixed(3)}) mm · |u| ${Math.hypot(u,v).toExponential(4)} m/s · u=(${u.toExponential(3)}, ${v.toExponential(3)}) m/s · p=${pr.toExponential(4)} Pa (P2/P1 evaluation)`}
function metrics(name){let a=result.nonlinear.fitted_wall.series[name==='raw'?'navier_stokes':name];if(!a)a=result.mesh.series['fitted-'+name];let x=a.at(-1),side=x.side_fraction*100,ref=result.headline['fitted-original'].side_percent,ok=side>=ref-1e-8;let originalNS=result.nonlinear.cross_evaluation.cases.find(x=>x.model==='navier_stokes'&&x.design==='original');let c=result.nonlinear.cross_evaluation.cases.find(x=>x.model==='navier_stokes'&&x.design===(name==='raw'?'navier_stokes_optimized':name==='stokes'?'stokes_optimized':name));return `<div class="metrics"><b>Finest fitted mesh:</b> ${(x.mechanical_power_W_per_m*1e6).toFixed(6)} µW/m · side ${side.toFixed(4)}% · area ${x.physical_fluid_area_mm2.toFixed(3)} mm²<br>Resistance ${x.hydraulic_resistance_pa_s_per_m2.toFixed(2)} Pa·s/m²<br>Physical reference ${ref.toFixed(4)}%: ${ok?'PASS':'FAIL'}${name==='smoothed'?' · mesh screen FAIL (provisional)':''}<br>${c?`Immersed NS floor ${(c.side_requirement*100).toFixed(4)}%: ${c.side_flow_percent.toFixed(4)}% — ${c.constraint_satisfied?'PASS':'FAIL'}; density bounds PASS<br>Immersed objective f = ${(c.power_w_per_m/originalNS.power_w_per_m).toFixed(6)}; scaled g = ${((c.side_requirement-c.side_flow_percent/100)/0.001).toFixed(5)}`:'Immersed constraints not applicable to postprocessed shape.'}</div>`}
async function refresh() {
    const token = ++refreshToken;
    const mode = $("mode").value, name = $("candidate").value;
    const displayField = mode !== "realistic";
    refreshing = true;
    $("inspect").textContent = "Loading selected comparison.";
    $("right-title").textContent = names[name];
    $("field").disabled = !displayField;
    $("vector-control").hidden = !displayField;
    $("section-control").hidden = mode !== "cutaway";
    $("mode-note").textContent = ({
        realistic: "Cylindrical engineering visualization from the planar contour. The display wall is 0.22 mm thick. The bore is 4 mm and the side outlet is 1 mm. The inner and outer lofts preserve open ports.",
        cutaway: "The section exposes the archived planar FEM solution at z = 0. Fields use the 0.25 mm fitted mesh and a shared scale. Arrows show sampled velocity with display-scaled length.",
        numerical: "Exact planar fitted geometry shown as a 0.6 mm extrusion. The mesh is the original P2/P1 analysis mesh at 0.25 mm. The extrusion depth and arrow lengths are for display."
    })[mode];
    const fields = displayField ? await Promise.all([loadField("original"),loadField(name)]) : [null,null];
    if (token !== refreshToken) return;
    if (displayField) limits = {
        speed: Math.max(...fields.map(d=>d.maxS), 1e-12),
        pressure: Math.max(...fields.map(d=>d.maxP), 1e-12)
    };
    ["original",name].forEach((design,i) => {
        setData(panels[i],design,mode === "numerical" ? null : loadCylinder(design),fields[i],mode);
        $(i===0 ? "left-metrics" : "right-metrics").innerHTML = metrics(design);
    });
    refreshing = false;
    $("viewer").dataset.ready = mode;
    $("viewer-fallback").hidden = true;
    $("inspect").textContent = displayField
        ? "Click a field element for coordinates, pressure and velocity. Drag to orbit. Use Shift or Pan to move both views."
        : "Drag to orbit. Use Shift or Pan to move both views. Select Cutaway or Numerical geometry to inspect the FEM fields.";
    legend();
    drawAll();
}
function legend(){let f=$('field').value;if($('mode').value==='realistic'||f==='geometry'||f==='mesh'){$('legend').innerHTML='Geometry / mesh display; no field scale.';return}let p=f==='pressure',v=p?limits.pressure:limits.speed;$('legend').innerHTML=`<div class="bar" style="background:${p?'linear-gradient(90deg,#2157a6,#f5f5f5,#ab2e26)':'linear-gradient(90deg,#e8f5ff,#0d3b8c)'}"></div>${p?(-v).toExponential(3)+' · 0 · ': '0 · '}${v.toExponential(3)} ${p?'Pa (signed pressure)':'m/s (speed)'} · shared range`}
function stage(i) {
    lastStage = i;
    $("stage").value = i;
    const checkpoint = stages[i], images = pipeData.history[i];
    $("stage-label").textContent = images.label;
    $("stage-image").src = images.geometry;
    $("stage-density").src = images.density;
    $("stage-metrics").textContent = "Immersed checkpoint · " + checkpoint.power.toFixed(6) +
        " µW/m · side flow " + checkpoint.side.toFixed(4) + "%";
}

function stopPlayback() {
    clearInterval(timer);
    timer = null;
    $("play").textContent = "Play";
}

async function main() {
    [result, stages, pipeData] = await Promise.all(
        ["data/results.json","assets/viewer/stages.json","assets/viewer/pipe-profiles.json"].map(fetchJSON));
    const h = result.headline;
    $("headline").innerHTML = [
        [h.raw_power_saving_percent.toFixed(2)+"%","Fitted pumping-power saving"],
        [h.adjoint_serial_speedup.toFixed(2)+"×","Serial adjoint speedup"],
        [h.kkt_stationarity.toExponential(3),"Projected stationarity"],
        [h["fitted-raw"].side_percent.toFixed(3)+"%","Optimized physical side flow"]
    ].map(([value,label]) => "<div class=\"stat\"><strong>"+value+"</strong>"+label+"</div>").join("");
    stage(0);
    $("stage").oninput = e => { stopPlayback(); stage(+e.target.value); };
    $("play").onclick = () => {
        if (timer) { stopPlayback(); return; }
        $("play").textContent = "Pause";
        timer = setInterval(() => stage((lastStage+1)%stages.length),1200);
    };
    $("rewind").onclick = () => { stopPlayback(); stage(0); };
    panels = [$("left"),$("right")].map(init);
    $("mode").onchange = $("candidate").onchange = () => refresh().catch(error);
    $("field").onchange = () => { legend(); drawAll(); };
    $("section-depth").oninput = $("vectors").onchange = drawAll;
    $("zoom-in").onclick = () => { camera.zoom = Math.min(4,camera.zoom * 1.25); drawAll(); };
    $("zoom-out").onclick = () => { camera.zoom = Math.max(.55,camera.zoom / 1.25); drawAll(); };
    $("pan").onclick = e => {
        const button = e.currentTarget;
        button.setAttribute("aria-pressed",button.getAttribute("aria-pressed") === "false");
    };
    for (const [button,angles] of Object.entries({front:[0,0],top:[0,Math.PI/2],isometric:[ISOMETRIC.yaw,ISOMETRIC.pitch]})) {
        $(button).onclick = () => { [camera.yaw,camera.pitch] = angles; drawAll(); };
    }
    $("reset").onclick = () => { Object.assign(camera,ISOMETRIC); drawAll(); };
    window.addEventListener("resize",drawAll);
    await refresh();
}
function error(e){console.error(e);refreshing=false;delete $('viewer').dataset.ready;$('viewer-fallback').hidden=false;$('inspect').textContent='Viewer unavailable: '+e.message+' — see archived figures below.'}main().catch(error);

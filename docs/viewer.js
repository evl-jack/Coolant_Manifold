const bg = [0.914, 0.941, 0.929];
const mean = (a, b) => (a + b) / 2;
const neutral = [0.96, 0.96, 0.92];
const low = [0.18, 0.36, 0.64];
const high = [0.66, 0.23, 0.12];
const speedHigh = [0.02, 0.35, 0.31];
const colorCSS = color => "rgb(" + color.map(value => Math.round(value * 255)).join(" ") + ")";
const depth = 0.8;
const rad = Math.PI / 180;
const defaults = {yaw: -18 * rad, pitch: 31 * rad, zoom: 1, x: 0, y: 0};

function program(gl, vertex, fragment) {
  const result = gl.createProgram();
  for (const [type, code] of [[gl.VERTEX_SHADER, vertex], [gl.FRAGMENT_SHADER, fragment]]) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, code);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) throw Error(gl.getShaderInfoLog(shader));
    gl.attachShader(result, shader);
    gl.deleteShader(shader);
  }
  gl.linkProgram(result);
  if (!gl.getProgramParameter(result, gl.LINK_STATUS)) throw Error(gl.getProgramInfoLog(result));
  return result;
}

const vertexStart = `
attribute vec3 point;
uniform vec3 right, up, forward;
uniform vec2 gain, pan;
void place() {
  vec3 p = point - vec3(10., 3.1, 0.);
  gl_Position = vec4(vec2(dot(p, right), dot(p, up)) * gain + pan,
                    -dot(p, forward) / 24., 1.);
}`;

function gpu(canvas) {
  const gl = canvas.getContext("webgl", {antialias: true});
  if (!gl) throw Error("WebGL is unavailable on this device");
  const shell = program(gl, vertexStart + `
    attribute vec3 rgb;
    varying vec3 tint;
    void main(){ place(); tint = rgb; }`,
    "precision mediump float; varying vec3 tint; void main(){ gl_FragColor=vec4(tint,1.); }");
  const field = program(gl, vertexStart + `
    attribute vec3 corner, uxv, uxm, uyv, uym, pv;
    varying vec3 b, xv, xm, yv, ym, p;
    void main(){ place(); b=corner; xv=uxv; xm=uxm; yv=uyv; ym=uym; p=pv; }`,
    `precision highp float;
    varying vec3 b, xv, xm, yv, ym, p;
    uniform int mode;
    uniform float limit;
    uniform vec3 bottomColor, middleColor, topColor;
    void main(){
      vec3 v=b*(2.*b-1.);
      vec3 m=4.*b*b.yzx;
      vec2 u=vec2(dot(v,xv)+dot(m,xm),dot(v,yv)+dot(m,ym));
      float s=mode == 1 ? length(u) : dot(b,p);
      float t=clamp(s/limit,-1.,1.);
      vec3 rgb=t<0. ? mix(middleColor,bottomColor,-t) : mix(middleColor,topColor,t);
      gl_FragColor=vec4(rgb,1.);
    }`);
  const buffers = new WeakMap();
  gl.enable(gl.DEPTH_TEST);
  gl.depthFunc(gl.LEQUAL);
  function draw(vertices, attribs, mode, camera, lines, range) {
    if (!vertices.length) return;
    const shader = mode ? field : shell;
    gl.useProgram(shader);
    const loc = name => gl.getUniformLocation(shader, name);
    gl.uniform3fv(loc("right"), camera.right);
    gl.uniform3fv(loc("up"), camera.up);
    gl.uniform3fv(loc("forward"), camera.forward);
    gl.uniform2fv(loc("gain"), camera.gain);
    gl.uniform2fv(loc("pan"), camera.pan);
    if (mode) {
      gl.uniform1i(loc("mode"), mode);
      gl.uniform1f(loc("limit"), range);
      gl.uniform3fv(loc("bottomColor"), low);
      gl.uniform3fv(loc("middleColor"), neutral);
      gl.uniform3fv(loc("topColor"), mode === 1 ? speedHigh : high);
    }
    if (!buffers.has(vertices)) {
      const buffer = gl.createBuffer();
      buffers.set(vertices, buffer);
      gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
      gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);
    }
    gl.bindBuffer(gl.ARRAY_BUFFER, buffers.get(vertices));
    const stride = 12 * attribs.length;
    const locations = attribs.map((name, i) => {
      const index = gl.getAttribLocation(shader, name);
      gl.enableVertexAttribArray(index);
      gl.vertexAttribPointer(index, 3, gl.FLOAT, false, stride, 12 * i);
      return index;
    });
    gl.drawArrays(lines ? gl.LINES : gl.TRIANGLES, 0, vertices.length / (3 * attribs.length));
    locations.forEach(index => gl.disableVertexAttribArray(index));
  }
  return {
    clear() {
      const ratio = Math.min(window.devicePixelRatio, 1.6);
      const w = Math.max(1, Math.round(canvas.clientWidth * ratio));
      const h = Math.max(1, Math.round(canvas.clientHeight * ratio));
      if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
      gl.viewport(0, 0, w, h);
      gl.clearColor(...bg, 1);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    },
    shell: (vertices, camera, lines = false) => draw(vertices, ["point", "rgb"], 0, camera, lines),
    field: (vertices, mode, camera, range) => draw(vertices, ["point", "corner", "uxv", "uxm", "uyv", "uym", "pv"], mode, camera, false, range)
  };
}

function meshBuilder() {
  const vertices = [];
  return {
    triangle(a, b, c, base) {
      const u = b.map((v, i) => v - a[i]), v = c.map((w, i) => w - a[i]);
      const normal = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
      const length = Math.hypot(...normal);
      if (length < 1e-13) return;
      const light = .63 + .37 * Math.abs((.22 * normal[0] + .48 * normal[1] + .85 * normal[2]) / length);
      const rgb = base.map(channel => channel * light);
      for (const point of [a,b,c]) vertices.push(...point, ...rgb);
    },
    quad(a, b, c, d, color) {
      this.triangle(a,b,c,color);
      this.triangle(a,c,d,color);
    },
    finish: () => new Float32Array(vertices)
  };
}

function extrude(data, color) {
  const mesh = meshBuilder();
  for (const face of data.caps) {
    for (const z of [-depth / 2, depth / 2]) mesh.triangle(...face.map(([x,y]) => [x,y,z]), color);
  }
  for (const ring of data.rings) {
    for (let i = 1; i < ring.length; i++) {
      const a = ring[i-1], b = ring[i];
      if ((a[0] === b[0] && (a[0] === 0 || a[0] === 20)) || (a[1] === 7 && b[1] === 7)) continue;
      mesh.quad([...a, -depth/2], [...b,-depth/2], [...b,depth/2], [...a,depth/2], color);
    }
  }
  return mesh.finish();
}

function cylinder(data, color) {
  const mesh = meshBuilder(), n = 64, r = Math.min(2, data.profile[0][3]), tau = 2*Math.PI;
  const phi = Array.from({length:n+1}, (_, i) => tau*i/n);
  const branch = phi.map(a => [r*Math.cos(a), r*Math.sin(a)]);
  const shoulder = ([x,z]) => [10+x, 2+Math.sqrt(Math.max(0, 4-z*z)), z];
  // Remove the opening from the upper pipe wall and use the same seam samples for the branch skirt.
  const edge = phi.slice(n/4,3*n/4+1).map((a,i) => [branch[i+n/4][1], Math.abs(r*Math.cos(a))]).reverse();
  const onPipe = (x, a) => [x, 2+2*Math.cos(a), 2*Math.sin(a)];
  function strip(a,b,left,right) {
    mesh.quad(onPipe(left[0],a),onPipe(right[0],a),onPipe(right[1],b),onPipe(left[1],b),color);
  }
  for (let i=1; i<edge.length; i++) {
    const [[za,xa],[zb,xb]] = edge.slice(i-1,i+1);
    const a = Math.asin(za/2), b = Math.asin(zb/2);
    strip(a,b,[0,0],[10-xa,10-xb]);
    strip(a,b,[10+xa,10+xb],[20,20]);
  }
  const end = Math.asin(r/2);
  for (let i=0; i<n; i++) {
    const a = end+(tau-2*end)*i/n, b = end+(tau-2*end)*(i+1)/n;
    strip(a,b,[0,0],[10,10]);
    strip(a,b,[10,10],[20,20]);
  }
  let previous = branch.map(shoulder);
  for (const [y,,,radius] of data.profile) {
    const next = phi.map(a => [10+radius*Math.cos(a),y,radius*Math.sin(a)]);
    for (let i=0; i<n; i++) mesh.quad(previous[i], previous[i+1],next[i+1],next[i],color);
    previous = next;
  }
  return mesh.finish();
}

const dot = (a,b) => a.reduce((sum,value,i) => sum+value*b[i],0);
export function interpolate(data, key, triangle, b) {
  const dofs = data.velocityDofs[triangle], f = data.fields[key], indices = data.triangles[triangle];
  const q = [...b.map(v => v*(2*v-1)), 4*b[0]*b[1], 4*b[1]*b[2], 4*b[2]*b[0]];
  return {ux:dot(q,dofs.map(i=>f.ux[i])), uy:dot(q,dofs.map(i=>f.uy[i])),
          pressure:dot(b,indices.map(i=>f.pressure[i]))};
}

function prep(data) {
  let speed = 0, pressure = 0;
  const buffer = {};
  for (const [key, f] of Object.entries(data.fields)) {
    const vertices = [];
    data.triangles.forEach((ids,t) => {
      const dofs = data.velocityDofs[t], ux = dofs.map(i=>f.ux[i]), uy = dofs.map(i=>f.uy[i]);
      const p = ids.map(i=>f.pressure[i]);
      pressure = Math.max(pressure,...p.map(Math.abs));
      // The quadratic Bernstein control polygon bounds the entire vector polynomial.
      for (let j=0;j<3;j++) {
        const k=(j+1)%3;
        speed=Math.max(speed,Math.hypot(ux[j],uy[j]),
          Math.hypot(2*ux[j+3]-mean(ux[j],ux[k]),2*uy[j+3]-mean(uy[j],uy[k])));
      }
      ids.forEach((i,j) => vertices.push(...data.points[i],0, ...[0,1,2].map(k=>Number(k===j)), ...ux,...uy,...p));
    });
    buffer[key]=new Float32Array(vertices);
  }
  const lines=[];
  for(const ids of data.triangles) {
    for(let j=0;j<3;j++) {
      for(const i of [ids[j],ids[(j+1)%3]]) lines.push(...data.points[i],0,.32,.41,.39);
    }
  }
  return {buffer, lines:new Float32Array(lines), speed:speed||1, pressure:pressure||1};
}

function axes(state, canvas) {
  const {yaw,pitch}=state, c=Math.cos(yaw), s=Math.sin(yaw), d=Math.cos(pitch), t=Math.sin(pitch);
  const scale=state.zoom*Math.min(canvas.clientWidth/23, canvas.clientHeight/13);
  return {right:[c,0,s], up:[s*t,d,-c*t], forward:[-s*d,t,c*d],
          gain:[2*scale/canvas.clientWidth,2*scale/canvas.clientHeight],pan:[state.x,state.y]};
}

function atPlane(event, canvas, camera) {
  const rect=canvas.getBoundingClientRect();
  const a=(2*(event.clientX-rect.left)/rect.width-1-camera.pan[0])/camera.gain[0];
  const b=(1-2*(event.clientY-rect.top)/rect.height-camera.pan[1])/camera.gain[1];
  const xyz=camera.right.map((v,i)=>a*v+b*camera.up[i]);
  const n=-xyz[2]/camera.forward[2];
  return [xyz[0]+n*camera.forward[0]+10,xyz[1]+n*camera.forward[1]+3.1];
}

function locate(data,x,y) {
  for (let i=0;i<data.triangles.length;i++) {
    const [[ax,ay],[bx,by],[cx,cy]]=data.triangles[i].map(j=>data.points[j]);
    const det=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
    const u=((x-ax)*(cy-ay)-(y-ay)*(cx-ax))/det, v=((bx-ax)*(y-ay)-(by-ay)*(x-ax))/det;
    if (u>=-1e-8&&v>=-1e-8&&u+v<=1+1e-8) return [i,[1-u-v,u,v]];
  }
  return null;
}

export async function initViewer(root, data, get) {
  const status = root.querySelector("#viewer-status"), state={...defaults};
  const panels=[...root.querySelectorAll(".view-panel")].map((panel,i)=>{
    const canvas=panel.querySelector("canvas");
    return {canvas, key:i ? "selected":"original", gpu:gpu(canvas),
            output:panel.querySelector(".inspect"), metric:panel.querySelector(".view-metrics")};
  });
  const mode=root.querySelector("#mode"), field=root.querySelector("#field"), mesh=root.querySelector("#mesh");
  const slider=root.querySelector("#checkpoint"), play=root.querySelector("#play"), legend=root.querySelector("#legend");
  const checkpoints=data.checkpoints, final=checkpoints.length-1;
  const geometry=new Map(), shell=new Map();
  let active=final, fields, prepared, fieldRequest, timer, loaded=false;
  slider.max=final;
  slider.value=final;
  function draw() {
    if (!loaded) return;
    const shading=mode.value==="planar"&&active===final&&prepared ? field.selectedIndex : 0;
    for (const panel of panels) {
      const camera=axes(state,panel.canvas);
      panel.gpu.clear();
      if (shading) {
        panel.gpu.field(prepared.buffer[panel.key],shading,camera,shading===1?prepared.speed:prepared.pressure);
        if(mesh.checked) panel.gpu.shell(prepared.lines,camera,true);
      } else {
        panel.gpu.shell(shell.get(mode.value+"/"+active+"/"+panel.key),camera);
      }
    }
    legend.hidden=!shading;
    if(shading) {
      const range=shading===1?prepared.speed:prepared.pressure;
      root.querySelector(".legend-bar").style.background="linear-gradient(90deg,"+
        (shading===1 ? [neutral,speedHigh] : [low,neutral,high]).map(colorCSS).join(",")+")";
      root.querySelector(".legend-labels").textContent=shading===1
        ? "Speed · 0 to "+range.toPrecision(3)+" m/s"
        : "Pressure · ±"+range.toPrecision(3)+" Pa · 0 at neutral midpoint";
    }
  }
  async function update() {
    const index=Number(slider.value), current=checkpoints[index], view=mode.value;
    field.disabled=view!=="planar"||index!==final;
    mesh.disabled=field.disabled||field.selectedIndex===0;
    status.textContent="Loading archived checkpoint…";
    if(!geometry.has(index)) geometry.set(index,get(current.geometry));
    const designs=await geometry.get(index);
    if(index!==Number(slider.value)||view!==mode.value) return;
    for(const [key,contour] of Object.entries(designs)) {
      const tag=view+"/"+index+"/"+key;
      const color=key==="original"?[.31,.46,.58]:[.15,.56,.44];
      if(!shell.has(tag)) shell.set(tag,(view==="planar"?extrude:cylinder)(contour,color));
    }
    active=index;
    loaded=true;
    for(const panel of panels) {
      const run=panel.key==="original"?current.reference:current.selected;
      const pass=run.side_fraction>=current.reference.side_fraction-1e-9;
      panel.metric.textContent="Immersed pilot · "+(1e6*run.power_w_per_m).toFixed(4)+" µW/m · Side flow "+
        (100*run.side_fraction).toFixed(3)+"% · "+(pass?"Pass":"Fail")+" (limit "+
        (100*current.reference.side_fraction).toFixed(3)+"%)";
      panel.output.textContent="Drag to orbit · Shift+drag to pan · Wheel to zoom";
    }
    const title="β = "+current.beta+" · archived stage "+current.stage+" of "+checkpoints.length;
    status.textContent=title+(view==="cylinder" ? " · Conceptual 3D interpretation. Planar results are shown below." :
      index!==final ? " · Planar extrusion. Fields are available for the final stage." : " · Planar extrusion, 0.8 mm display depth.");
    draw();
    if(!field.disabled&&field.selectedIndex) {
      status.textContent=title+" · Loading archived immersed FEM fields…";
      fieldRequest ??= get(data.fields).then(value=>{ fields=value; prepared=prep(value); })
        .catch(error=>{fieldRequest=null;throw error;});
      await fieldRequest;
      if(index===Number(slider.value)&&view===mode.value&&field.selectedIndex) {
        status.textContent=title+" · Immersed FEM plane includes the penalized region. Pressure uses the solver datum.";
        draw();
      }
    }
  }
  const refresh=()=>update().catch(error=>{status.textContent=error.message;console.error(error);});
  for(const control of [mode,field,mesh,slider]) control.addEventListener("input",refresh);
  function stop() {clearInterval(timer);timer=null;play.textContent="Play checkpoints";}
  play.addEventListener("click",()=>{
    if(timer) return stop();
    slider.value=0;
    refresh();
    play.textContent="Pause";
    timer=setInterval(()=>{
      if(Number(slider.value)===final) return stop();
      slider.stepUp();
      refresh();
    },950);
  });
  root.querySelector("#reset").addEventListener("click",()=>{Object.assign(state,defaults);draw();});
  root.querySelectorAll("[data-zoom]").forEach(button=>button.addEventListener("click",()=>{
    state.zoom=Math.min(12,Math.max(.45,state.zoom*Number(button.dataset.zoom)));draw();
  }));
  const observer=new ResizeObserver(draw);
  for(const {canvas,key,output} of panels) {
    observer.observe(canvas);
    let dragging, locked=false;
    function inspect(event) {
      if(!fields||field.disabled||!field.selectedIndex) return;
      const camera=axes(state,canvas);
      if(Math.abs(camera.forward[2])<1e-5) { output.textContent="Orbit to face the FEM plane.";return; }
      const [x,y]=atPlane(event,canvas,camera), result=locate(fields,x,y);
      if(!result) {output.textContent="Outside the immersed analysis mesh.";return;}
      const [i,b]=result, value=interpolate(fields,key,i,b);
      output.textContent="Element "+(i+1)+" · ("+x.toFixed(2)+", "+y.toFixed(2)+") mm · u = ("+
        value.ux.toPrecision(3)+", "+value.uy.toPrecision(3)+") m/s · p = "+value.pressure.toPrecision(3)+" Pa";
    }
    canvas.addEventListener("pointerdown",event=>{
      dragging={x:event.clientX,y:event.clientY,move:false,pan:event.shiftKey||event.button===1||event.button===2};
      canvas.setPointerCapture(event.pointerId);
    });
    canvas.addEventListener("pointermove",event=>{
      if(!dragging){if(!locked)inspect(event);return;}
      const dx=event.clientX-dragging.x, dy=event.clientY-dragging.y;
      dragging.move ||= Math.abs(dx)+Math.abs(dy)>2;
      if(dragging.pan){state.x+=2*dx/canvas.clientWidth;state.y-=2*dy/canvas.clientHeight;}
      else {state.yaw=(state.yaw+dx*.012)%(tauTurn);state.pitch=Math.max(-1.48,Math.min(1.48,state.pitch+dy*.012));}
      dragging.x=event.clientX;dragging.y=event.clientY;draw();
    });
    canvas.addEventListener("pointerup",event=>{if(dragging&&!dragging.move){locked=!locked;inspect(event);}dragging=null;});
    canvas.addEventListener("pointercancel",()=>{dragging=null;});
    canvas.addEventListener("contextmenu",event=>event.preventDefault());
    canvas.addEventListener("wheel",event=>{
      event.preventDefault();state.zoom=Math.min(12,Math.max(.45,state.zoom*Math.exp(-event.deltaY*.001)));draw();
    },{passive:false});
    canvas.addEventListener("keydown",event=>{
      const axesKeys={ArrowLeft:[-.08,0],ArrowRight:[.08,0],ArrowUp:[0,-.08],ArrowDown:[0,.08]};
      const delta=axesKeys[event.key];
      if(!delta)return;
      event.preventDefault();
      if(event.shiftKey){state.x+=delta[0];state.y-=delta[1];}
      else {state.yaw+=delta[0];state.pitch=Math.max(-1.48,Math.min(1.48,state.pitch+delta[1]));}
      draw();
    });
  }
  await update();
  root.dataset.ready="true";
}
const tauTurn=2*Math.PI;

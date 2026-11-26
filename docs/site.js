import {initViewer} from "./viewer.js";

const $ = selector => document.querySelector(selector);
const fmt = (n, places) => n == null ? "Pending" : n.toFixed(places);
async function get(path) {
  const response = await fetch(path);
  if (!response.ok) throw Error("Could not load " + path + " (" + response.status + ")");
  return response.json();
}

function table(container, columns, rows) {
  const element=document.createElement("table");
  const head=element.createTHead().insertRow();
  columns.forEach(label=>{
    const cell=document.createElement("th");
    cell.scope="col";
    cell.textContent=label;
    head.append(cell);
  });
  const body=element.createTBody();
  for(const row of rows) {
    const target=body.insertRow();
    for(const cell of row) target.insertCell().textContent=cell??"Pending";
  }
  container.append(element);
}

function chart(container, journal) {
  const series=journal.filter(point=>point.phase==="mma");
  if(!series.length) return;
  const svgNS="http://www.w3.org/2000/svg", svg=document.createElementNS(svgNS,"svg");
  svg.setAttribute("viewBox","0 0 640 256");
  svg.setAttribute("class","trace");
  svg.setAttribute("role","img");
  svg.setAttribute("aria-label","MMA normalized pumping power by archived evaluation, grouped by continuation stage");
  const lo=Math.min(1,...series.map(p=>p.normalized_power)), hi=Math.max(1,...series.map(p=>p.normalized_power));
  const x=i=>90+520*i/Math.max(1,series.length-1);
  const y=v=>205-155*(v-lo)/Math.max(1e-8,hi-lo);
  function add(tag, attributes, text) {
    const shape=document.createElementNS(svgNS,tag);
    for(const [key,value] of Object.entries(attributes)) shape.setAttribute(key,value);
    if(text) shape.textContent=text;
    svg.append(shape);
  }
  add("text",{x:22,y:27,fill:"#172c32","font-size":17},"Normalized power by MMA evaluation");
  add("path",{d:"M90 45V205H610",stroke:"#62746f",fill:"none"});
  for(let i=0;i<=4;i++) {
    const value=lo+(hi-lo)*i/4;
    add("text",{x:16,y:y(value)+4,fill:"#42585e","font-size":13},value.toFixed(5));
  }
  add("text",{x:89,y:231,fill:"#42585e","font-size":13},"First evaluation");
  add("text",{x:547,y:231,fill:"#42585e","font-size":13},"Final evaluation");
  const inks=["#206558","#8c572a","#344d90"];
  for(const stage of new Set(series.map(p=>p.stage))) {
    const coordinates=series.flatMap((p,i)=>p.stage===stage?[[x(i),y(p.normalized_power)]]:[]);
    add("polyline",{points:coordinates.map(p=>p.join(",")).join(" "),fill:"none",stroke:inks[stage-1],"stroke-width":2});
    add("text",{x:coordinates[0][0],y:45,fill:inks[stage-1],"font-size":12},"Stage "+stage);
  }
  container.append(svg);
}

function releaseSections(data) {
  for(const card of document.querySelectorAll("[data-evidence]")) {
    const record=data.pending[card.dataset.evidence];
    if(!record) {
      const pending=document.createElement("p");
      pending.className="pending";
      pending.textContent="[PENDING: "+card.querySelector("strong").textContent.toUpperCase()+"]";
      card.append(pending);
      table(card,["Result","Value","Units","Status"],[["Numerical release",null,null,"Awaiting source record"]]);
      continue;
    }
    table(card,record.columns,record.rows);
    const source=document.createElement("p");
    source.className="muted";
    source.textContent=[record.status,record.experiment,record.source].join(" · ");
    card.append(source);
    if(record.figure) {
      const img=document.createElement("img"), figure=document.createElement("figure");
      img.alt=record.figureAlt||card.querySelector("strong").textContent;
      img.loading="lazy";
      img.src="assets/"+record.figure;
      figure.append(img);
      card.append(figure);
    }
  }
}

async function report(root) {
  const data=await get("assets/results.json");
  const {fitted}=data;
  const start=fitted.Original, raw=fitted["Raw optimized"];
  $("#headline-power").textContent=fmt(100*(1-raw.power.value/start.power.value),2)+"%";
  $("#headline-start").textContent=fmt(start.flow.value,3)+"%";
  $("#headline-end").textContent=fmt(raw.flow.value,3)+"%";
  $("#headline-area").textContent=fmt(raw.area.value,3)+" "+raw.area.units;
  $("#headline-turn").textContent=fmt(raw.turn.value,2)+raw.turn.units;
  const keys=Object.keys(fitted);
  const definitions=[["Mechanical port power","power",4],["Side-outlet flow","flow",3],
    ["Effective hydraulic resistance","resistance",1],["Physical fluid area","area",3],
    ["Maximum wall turn","turn",2],["Element count","mesh",0]];
  table($("#fitted-table"),["Metric",...keys],definitions.map(([name,key,places])=>
    [name+" ("+start[key].units+")",...keys.map(design=>fmt(fitted[design][key].value,places))]));
  table($("#checkpoints"),["β","Unique evaluations","Pilot exit","P/P₀","Side flow (%)","Minimum (%)","Flow constraint"],
    data.checkpoints.map(({beta,status,evaluations,selected,reference})=>[beta,evaluations,status,fmt(selected.normalized_power,6),
      fmt(100*selected.side_fraction,4),fmt(100*reference.side_fraction,4),
      selected.side_fraction>=reference.side_fraction-1e-9?"Pass":"Fail"]));
  if(data.kkt) table($("#endpoint-kkt"),["Archived replay","Multiplier","Stationarity ∞-norm","Complementarity"],
    data.kkt.map(k=>[k.verified?"Passed":"Failed",k.multipliers[0],k.stationarity_inf.toExponential(3),
                      k.complementarity_inf.toExponential(3)]));
  $("#smoke-figure").hidden=!data.adjointSmoke;
  if(data.adjointSmoke) $("#smoke-figure img").src="assets/adjoint-smoke.png";
  chart($("#history"),data.history);
  releaseSections(data);
  await initViewer(root,data,get);
}

const viewer=$("#viewer");
if(viewer) report(viewer).catch(error=>{
  $("#viewer-status").textContent=error.message;
  console.error(error);
});

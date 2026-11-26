const assert = require("node:assert/strict");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const http = require("node:http");
const {spawn} = require("node:child_process");

const base = "/coolant-manifold/";
const root = path.resolve(process.argv[2] || path.join(__dirname, "dist"));
const chrome = process.env.PROJECT4_CHROME || (process.platform === "win32" ?
  "C:/Program Files/Google/Chrome/Application/chrome.exe" : "chromium");
const mime = {".html":"text/html",".css":"text/css",".js":"text/javascript",".json":"application/json",".png":"image/png"};
const loads = [], missing = [];
const listen = http.createServer((req,res) => {
  const pathname = new URL(req.url,"http://localhost").pathname;
  const file = path.resolve(root,pathname.slice(base.length)||"index.html");
  if (!pathname.startsWith(base)||path.relative(root,file).startsWith("..")||!fs.existsSync(file)||!fs.statSync(file).isFile()) {
    if(pathname!=="/favicon.ico") missing.push(pathname);
    res.writeHead(404).end();
    return;
  }
  loads.push(pathname);
  res.writeHead(200,{"Content-Type":mime[path.extname(file)] || "application/octet-stream"});
  fs.createReadStream(file).pipe(res);
});
const sleep = ms => new Promise(resolve=>setTimeout(resolve,ms));
async function until(check, description) {
  for(let i=0;i<120;i++){
    if(await check()) return;
    await sleep(100);
  }
  throw Error(description);
}

async function run() {
  await new Promise(resolve=>listen.listen(0,"127.0.0.1",resolve));
  const url="http://127.0.0.1:"+listen.address().port+base;
  const profile=fs.mkdtempSync(path.join(os.tmpdir(),"fixed-pipe-browser-"));
  const child=spawn(chrome,["--headless=new","--incognito","--no-first-run","--window-size=1440,1050",
    "--use-gl=angle","--use-angle=swiftshader","--enable-unsafe-swiftshader","--remote-debugging-port=0",
    "--user-data-dir="+profile,url],{stdio:["ignore","ignore","pipe"]});
  let logs="", socket;
  child.on("error",error=>{logs+=error.message;});
  child.stderr.on("data",chunk=>{logs+=chunk;});
  try {
    const portFile=path.join(profile,"DevToolsActivePort");
    await until(()=>fs.existsSync(portFile),"Chromium startup timeout").catch(error=>{throw Error(error.message+" "+logs);});
    const port=fs.readFileSync(portFile,"utf8").split("\n")[0];
    const tabs=await (await fetch("http://127.0.0.1:"+port+"/json/list")).json();
    socket=new WebSocket(tabs.find(tab=>tab.url===url||tab.type==="page").webSocketDebuggerUrl);
    await new Promise((resolve,reject)=>{socket.addEventListener("open",resolve);socket.addEventListener("error",reject);});
    let counter=0;
    const waiting=new Map(), errors=[];
    socket.addEventListener("message",event=>{
      const data=JSON.parse(event.data);
      if(data.method==="Runtime.exceptionThrown") errors.push(data.params.exceptionDetails.text);
      if(waiting.has(data.id)){waiting.get(data.id)(data);waiting.delete(data.id);}
    });
    const send=(method,params={})=>new Promise(resolve=>{
      const id=++counter; waiting.set(id,resolve); socket.send(JSON.stringify({id,method,params}));
    });
    const evaluate=async fn=>{
      const response=await send("Runtime.evaluate",{expression:"("+fn.toString()+")()",awaitPromise:true,returnByValue:true});
      if(response.error?.message.includes("context")) return null;
      assert(!response.error,JSON.stringify(response.error));
      assert(!response.result.exceptionDetails,response.result.exceptionDetails?.exception?.description);
      return response.result.result.value;
    };
    await send("Runtime.enable");
    await until(()=>evaluate(()=>document.querySelector("#viewer")?.dataset.ready==="true"),"Viewer did not initialize")
      .catch(async error=>{throw Error(error.message+" "+await evaluate(()=>document.querySelector("#viewer-status")?.textContent));});
    assert(!loads.includes(base+"assets/fields.json"),"Fields loaded eagerly");
    const desktop=await evaluate(()=>[document.querySelector("#headline-power").textContent,
      document.querySelectorAll("math").length,document.querySelectorAll(".pending").length,
      document.querySelectorAll(".view-panel canvas").length]);
    assert.equal(desktop[0],"4.83%");
    assert(desktop[1]>=6);
    assert(desktop[2]<=7);
    assert.equal(desktop[3],2);
    const features=await evaluate(async ()=>{
      const $=selector=>document.querySelector(selector), wait=ms=>new Promise(resolve=>setTimeout(resolve,ms));
      const select=(id,value)=>{const control=$(id);control.value=value;control.dispatchEvent(new Event("input"));};
      const sample=()=>{
        const image=document.createElement("canvas");
        image.width=image.height=72;
        const ctx=image.getContext("2d");
        ctx.drawImage($("canvas"),0,0,72,72);
        const pixels=ctx.getImageData(0,0,72,72).data;
        return new Set(Array.from({length:72*72},(_,i)=>pixels[i*4]+256*pixels[i*4+1])).size;
      };
      // Capture in the same event loop as the WebGL redraw.
      $("#reset").click();
      const planar=sample();
      select("#mode","cylinder");await wait(150);$("#reset").click();
      const cylinder=sample(), conceptual=$("#viewer-status").textContent.includes("Conceptual");
      select("#mode","planar");await wait(100);
      select("#field","Velocity magnitude");
      for(let i=0;i<40&&$("#legend").hidden;i++)await wait(100);
      $("#reset").click();
      const velocity=sample(), velocityUnit=$(".legend-labels").textContent;
      select("#field","Pressure");await wait(100);$("#reset").click();
      const pressure=sample(), pressureUnit=$(".legend-labels").textContent;
      const canvas=$("canvas"), r=canvas.getBoundingClientRect();
      canvas.dispatchEvent(new PointerEvent("pointermove",{clientX:r.left+r.width/2,clientY:r.top+r.height*.53}));
      const inspection=$(".inspect").textContent;
      $("#play").click();await wait(140);const playback=$("#viewer-status").textContent.includes("β = 0");
      $("#play").click();
      select("#checkpoint","2");await wait(120);
      return {planar,cylinder,conceptual,velocity,pressure,velocityUnit,pressureUnit,inspection,playback};
    });
    assert(features.planar>2 && features.cylinder>10 && features.conceptual,JSON.stringify(features));
    assert(features.velocity>80 && features.pressure>80,JSON.stringify(features));
    assert(features.velocityUnit.includes("m/s") && features.pressureUnit.includes("Pa"),JSON.stringify(features));
    assert(features.inspection.includes("Element") && features.playback,JSON.stringify(features));
    assert(loads.includes(base+"assets/fields.json"));
    await send("Emulation.setDeviceMetricsOverride",{width:390,height:844,deviceScaleFactor:1,mobile:true});
    await sleep(250);
    const phone=await evaluate(()=>[document.documentElement.scrollWidth,
      getComputedStyle(document.querySelector(".pair")).gridTemplateColumns.split(" ").length]);
    assert(phone[0]<=390&&phone[1]===1,JSON.stringify(phone));
    for(const page of ["optimization.html","handoff.html"]) {
      const response=await fetch(url+page);
      assert(response.ok);
      assert((await response.text()).includes("</main>"));
    }
    assert.deepEqual(missing,[]);
    assert.deepEqual(errors,[]);
    console.log("PASS private Chromium profile, repository subpath, WebGL geometry and fields, interpolation inspection, playback, equations, phone layout and local links");
  } finally {
    socket?.close();
    child.kill();
    listen.close();
  }
}
run().catch(error=>{console.error(error);listen.close();process.exitCode=1;});

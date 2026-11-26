import {readFileSync} from "node:fs";
import {buildPipe} from "./assets/pipe.js";
const {profiles, wall_mm} = JSON.parse(readFileSync(new URL("./assets/viewer/pipe-profiles.json", import.meta.url), "utf8"));
for (const [name, profile] of Object.entries(profiles)) {
    const {ports, ...audit} = buildPipe(profile.sections_mm, wall_mm, true);
    const radii = ports.map((ring, index) => Math.max(...ring.map(([x,y,z])=>index===2?Math.hypot(x-10,z):Math.hypot(y-2,z))));
    console.log(name, audit, radii);
    if (audit.openEdges || !audit.connected || radii.some((r,i)=>Math.abs(r-(i===2?.5:2))>1e-8)) process.exitCode = 1;
}

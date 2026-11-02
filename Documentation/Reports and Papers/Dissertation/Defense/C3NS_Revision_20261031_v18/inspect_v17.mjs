import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "file:///C:/Users/Ben%20Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";

const input = "../C3NS_Revision_20261031_v17/output/Bolen_PhD_Defense_RoseGreen_20261031_v17.pptx";
const deck = await PresentationFile.importPptx(await FileBlob.load(input));
const result = await deck.inspect({
  select: "$.slides[0,1,2,3,4,5,6,7,8,9,10,43,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69]",
  include: ["id", "index", "title", "elements.id", "elements.name", "elements.kind", "elements.text", "elements.position", "notes.text"],
  fileName: "build/v17-target-slides.ndjson",
});
await fs.writeFile("build/v17-target-slides.stdout.ndjson", result.ndjson, "utf8");
const montage = await deck.export({ format: "webp", montage: true, scale: 0.35 });
await fs.writeFile("build/v17-montage.webp", new Uint8Array(await montage.arrayBuffer()));
console.log(`slides=${deck.slides.items.length} records=${result.recordCount}`);

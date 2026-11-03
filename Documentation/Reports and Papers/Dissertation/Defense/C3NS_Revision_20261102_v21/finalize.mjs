import fs from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {execFileSync} from 'node:child_process';
import {finalizePresentation} from 'file:///C:/Users/Ben%20Bolen/.codex/plugins/cache/openai-primary-runtime/presentations/26.1017.22136/skills/presentations/container_tools/artifact_tool_utils.mjs';

process.env.RUNTIME_NODE_MODULES='C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES=process.env.RUNTIME_NODE_MODULES;

const W='D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Defense/C3NS_Revision_20261102_v21';
const S='C:/Users/Ben Bolen/.codex/plugins/cache/openai-primary-runtime/presentations/26.1017.22136/skills/presentations';
const reference='D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Defense/C3NS_Revision_20261102_v20/output/Bolen_PhD_Defense_RoseGreen_20261102_v20.pptx';
const python='C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe';
execFileSync(python,['-X','utf8',W+'/verify_delivery.py','--target',W+'/build/candidate.pptx','--skip-pdf'],{stdio:'inherit'});
const structure=JSON.parse(await fs.readFile(W+'/build/structure.json','utf8'));
const order=JSON.parse(await fs.readFile(W+'/build/slide_order.json','utf8'));
const owners=[15,25,30,52,101,102].map(n=>order.indexOf(n)+1).sort((a,b)=>a-b);
const result=await finalizePresentation({
  workspaceDir:W,
  candidatePath:W+'/build/candidate.pptx',
  finalPath:W+'/output/Bolen_PhD_Defense_RoseGreen_20261102_v21.pptx',
  explicitTotalSlideCount:structure.total,
  requiredNativeChartOwnerSlides:owners,
  requiredEmbeddedWorkbookChartOwnerSlides:owners,
  pythonExecutable:'C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe',
  integrityValidatorPath:S+'/container_tools/inspect_presentation_package_integrity.py',
  layoutValidatorPath:S+'/container_tools/inspect_presentation_layout_geometry.py',
  layoutArgs:['--expected-slide-size-emu','9144000,5143500','--validate-bullet-geometry','--validate-heading-fit'],
  fontPolicy:{
    basis:'reference',
    families:['Calibri','Arial'],
    referencePath:reference,
    referenceSha256:createHash('sha256').update(await fs.readFile(reference)).digest('hex')
  },
  materializeLiteralChartWorkbooks:false,
  verifyArtifactToolImport:true,
  receiptPath:W+'/build/validation_delivery.json'
});
execFileSync(python,['-X','utf8',W+'/verify_delivery.py','--target',result.finalPath,'--skip-pdf'],{stdio:'inherit'});
console.log(JSON.stringify({path:result.finalPath,sha256:result.finalSha256,integrity:result.packageIntegrity.passed,import:result.firstPartyImport.passed}));


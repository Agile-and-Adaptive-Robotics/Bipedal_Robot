import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {finalizePresentation} from 'file:///C:/Users/Ben%20Bolen/.codex/plugins/cache/openai-primary-runtime/presentations/26.1017.22136/skills/presentations/container_tools/artifact_tool_utils.mjs';
const W=process.cwd(),B=path.join(W,'build/plot_rotation_180');
process.env.RUNTIME_NODE_MODULES='C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES=process.env.RUNTIME_NODE_MODULES;
const S='C:/Users/Ben Bolen/.codex/plugins/cache/openai-primary-runtime/presentations/26.1017.22136/skills/presentations';
const reference=path.resolve(W,'../C3NS_Revision_20261102_v20/output/Bolen_PhD_Defense_RoseGreen_20261102_v20.pptx');
await fs.mkdir(path.join(B,'validated'),{recursive:true});
const result=await finalizePresentation({
 workspaceDir:W,candidatePath:path.join(B,'candidate.pptx'),finalPath:path.join(B,'validated/current.pptx'),
 explicitTotalSlideCount:145,requiredNativeChartOwnerSlides:[15,25,32,56,104,105],requiredEmbeddedWorkbookChartOwnerSlides:[15,25,32,56,104,105],
 pythonExecutable:'C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe',
 integrityValidatorPath:S+'/container_tools/inspect_presentation_package_integrity.py',
 layoutValidatorPath:S+'/container_tools/inspect_presentation_layout_geometry.py',
 layoutArgs:['--expected-slide-size-emu','9144000,5143500','--validate-bullet-geometry','--validate-heading-fit'],
 fontPolicy:{basis:'reference',families:['Calibri','Arial'],referencePath:reference,referenceSha256:createHash('sha256').update(await fs.readFile(reference)).digest('hex')},
 materializeLiteralChartWorkbooks:false,verifyArtifactToolImport:true,receiptPath:path.join(B,'validation.json')
});
console.log(JSON.stringify({sha256:result.finalSha256,path:result.finalPath}));

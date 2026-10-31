import fs from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {finalizePresentation} from 'file:///C:/Users/Ben%20Bolen/.codex/plugins/cache/openai-primary-runtime/presentations/26.1017.22136/skills/presentations/container_tools/artifact_tool_utils.mjs';
const W='D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Defense/C3NS_Revision_20261031_v17';
const S='C:/Users/Ben Bolen/.codex/plugins/cache/openai-primary-runtime/presentations/26.1017.22136/skills/presentations';
const reference='D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Defense/C3NS_Revision_20261030_v16/output/Bolen_PhD_Defense_RoseGreen_20261030_v16.pptx';
const result=await finalizePresentation({
  workspaceDir:W,candidatePath:W+'/build/candidate.pptx',finalPath:W+'/output/Bolen_PhD_Defense_RoseGreen_20261031_v17.pptx',
  explicitTotalSlideCount:123,requiredNativeChartOwnerSlides:[16,26,31,44,97,98],requiredEmbeddedWorkbookChartOwnerSlides:[16,26,31,44,97,98],
  pythonExecutable:'C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe',
  integrityValidatorPath:S+'/container_tools/inspect_presentation_package_integrity.py',
  layoutValidatorPath:S+'/container_tools/inspect_presentation_layout_geometry.py',
  layoutArgs:['--expected-slide-size-emu','9144000,5143500','--validate-bullet-geometry','--validate-heading-fit'],
  fontPolicy:{basis:'reference',families:['Calibri','Arial'],referencePath:reference,referenceSha256:createHash('sha256').update(await fs.readFile(reference)).digest('hex')},
  materializeLiteralChartWorkbooks:false,verifyArtifactToolImport:true,receiptPath:W+'/build/validation_delivery5.json'
});
console.log(JSON.stringify({path:result.finalPath,sha256:result.finalSha256,integrity:result.packageIntegrity.passed,import:result.firstPartyImport.passed}));




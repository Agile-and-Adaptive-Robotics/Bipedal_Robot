import path from 'node:path';
import {createRequire} from 'node:module';
const require=createRequire('C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/package.json');
const sharp=require('sharp');
const root=path.resolve('assets/neural/slide_panels');
for (const name of ['NEURAL_RIGHT_HIERARCHY','CONTACT_PF_AND_FOOT_LOADING','RF_MOTOR_FORCE_CROP','RF_LENGTH_SENSORY_CROP']) {
 await sharp(path.join(root,name+'.svg'),{density:180}).png().toFile(path.join(root,name+'.png'));
 console.log(name);
}

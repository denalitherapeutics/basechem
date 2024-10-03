import {checkTask} from "../../main/static/js/tasks.js"
import { showLoadingModal } from "../../mni_common/static/mni_common/js/loading.js"
import PropCalc from "../../main/static/js/propcalc.js"
import LigandAlign, {parseAndShowAlignment} from "../../main/static/js/align.js"
import EspMap from "../../main/static/js/esp.js"
import Dock from "../../main/static/js/dock.js"
import Torsion from "../../main/static/js/torsion.js"
import Mmp from "../../main/static/js/mmp.js"
import {updateSavedGems} from "../../main/static/js/saveGems.js"
import {addComps, addCompsEventListeners} from "../../main/static/js/addComps.js"
import {enableRowReorder} from "../../main/static/js/reorderCos.js"
import {copyDataId} from "../../mni_common/static/mni_common/js/select2.js"
import {toggleDisableChildren, toggleElement} from "../../mni_common/static/mni_common/js/toggle.js"
import AjaxForms from "../../mni_common/static/mni_common/js/ajaxForms.js"
import {copyToClipboard} from "../../mni_common/static/mni_common/js/copyToClipboard.js"
import {mmpKetcherHandler} from "../../main/static/js/mmp_ketcher.js"


window.checkTask = checkTask;
window.showLoadingModal = showLoadingModal;
window.parseAndShowAlignment = parseAndShowAlignment;
window.toggleDisableChildren = toggleDisableChildren;
window.toggleElement = toggleElement;
window.copyDataId = copyDataId;
window.updateSavedGems = updateSavedGems;
window.addComps = addComps;
window.copyToClipboard = copyToClipboard;
window.enableRowReorder = enableRowReorder;
window.mmpKetcherHandler = mmpKetcherHandler;
window.PropCalc = PropCalc;
window.LigandAlign = LigandAlign;
window.Dock = Dock;
window.EspMap = EspMap;
window.Torsion = Torsion;
window.Mmp = Mmp;
document.addEventListener('DOMContentLoaded', (event) => {
  AjaxForms();
  addCompsEventListeners();
});

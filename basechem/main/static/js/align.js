import {PRIMARY_OPACITY, CARBON_COLORS, STYLE, download, removeModel, addModel, fitViewerToModels, colorTags, styleReceptor, toggleTorsionAlerts} from './viewer.js'
import {saveViewerGems} from "./saveGems.js"
import {handleNewCompsAdded} from "./addComps.js"
import {copyDataId} from "../../../mni_common/static/mni_common/js/select2.js"
import {showConfOptions, addSaveGemsModal, manageNewAnalysisForm, collectAnalysisTask, collectRefsTask,} from "./tasks.js"

/***************************************************
 *             Constants/Global Vars               *
 ***************************************************/
// Keeps track of how many times each color is used
const colorFrequency = CARBON_COLORS.reduce(function(obj, x) {
    obj[x] = 0;
    return obj;
}, {});
// Keeps track of which color is used for each model
const colorMap = {}
const models = {}
/* parsedData is an object with all the results from from the backend alignment. It takes the form:
 * {compounds: {comp_id: {moltext, r_mmff_rel_energy, r_bc_rmsd_to_lsalign}}, 
 * references: {series_id: moltext},
 * receptors: {series_id: {"moltext": moltext, "residues": int_array}}}
 */
var parsedData = {"compounds": {}, "references": {}, "receptors": {}}
const SLAB_NEAR = -10
const SLAB_FAR = 10


/**
 * This function handles all of the front-end interactions for the ligand aligner
 * after the backend alignment process is complete
 */
 export default function(){
    const groupName = JSON.parse(document.getElementById('group-name').textContent);
    if(groupName){
        // Disable the new analysis form while tasks are running
        manageNewAnalysisForm()
        // Create 3DmolJS viewer to display compounds
        const viewer = $3Dmol.createViewer('align-viewer')
        // Wait for the alignrefs task to finish and add the results to the dropdown when ready
        collectRefsTask(parsedData, viewer, 'align', collectAlignRefs);        
        // Add event listeners to the download buttons
        download(models, '.align-download');
        // Add event listener to the fit-viewer button
        fitViewerToModels(viewer, models, '#fit-viewer')
        // Add event listener to the toggle-torsion-alerts button
        toggleTorsionAlerts(viewer, models, '#toggle-torsion-alerts')
        // Add event listener to the save-viewer-gems button
        saveViewerGems(models, '#save-viewer-gems')
        // Add event listener to the toggle-receptor button
        toggleReceptor(viewer);
        // Wait for each align task to finish and add the results to the frontend when ready
        collectAlignTasks(parsedData, viewer)    
    }  
}

/**
 * For each `select-conformer` widget in the table, wait for the associated DjangoQ task to complete
 * and display the results.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the displayed 3dmol viewer
 */
function collectAlignTasks(parsedData, viewer){
    document.querySelectorAll(".select-conformer").forEach(element => {
        const taskName = element.getAttribute("data-task-name")
        collectAnalysisTask(taskName, parsedData, viewer, collectSuccessfulAlignTask)
    })
    handleNewCompsAdded("align-table", parsedData, viewer, newCompAddedHelper)
}

/**
 * Called in `collectAnalysisTask` when the task returns successfully - this function uses
 * the result from the backend to update parsedData, the viewer, and the DOM
 * @param {Number} coPk the PK of the CompoundOccurrence whose task returned successfully
 * @param {Object} ajaxData the data returned from the backend, including task results and HTML strings to add to the DOM
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the displayed 3dmol viewer
 */
function collectSuccessfulAlignTask(coPk, ajaxData, parsedData, viewer){
    parsedData.compounds[`co-${coPk}`] = ajaxData.taskResult
    const selectConformerElement = document.querySelector(`.select-conformer#co-${coPk}`)
    selectConformerElement.insertAdjacentHTML('beforeend', ajaxData.confOptions)
    selectConformerElement.parentElement.style.display = ""
    addSaveGemsModal(coPk, ajaxData.saveGemsModal)
    selectConformers(coPk, parsedData, viewer)
}

/**
 * Used in handleNewCompsAdded, this function collects task results for new compounds.
 * @param {Object} trElement a table row element for a newly added compound
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the displayed 3dmol viewer
 */
function newCompAddedHelper(trElement, parsedData, viewer){
    const taskName = trElement.querySelector(".select-conformer").getAttribute("data-task-name")
    collectAnalysisTask(taskName, parsedData, viewer, collectSuccessfulAlignTask)
    $('.select-conformer').select2({
        placeholder: "Conformers",
        width: "100%",
        templateResult: copyDataId,
        templateSelection: copyDataId,
    })
}

/**
 * Waits for the alignrefs DjangoQ task to finish. Once it is complete, populates the reference dropdown
 * with results. If the task is not successful, the reference dropdown will remain disabled.
 * @param {Object} taskData a dictionary with the result of the alignrefs task
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data.
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
 async function collectAlignRefs(taskData, parsedData, viewer){
    // Add taskData to parsedData so that models can be added to the viewer
    parsedData["references"] = taskData.taskResult.references
    parsedData["receptors"] = taskData.taskResult.receptors
    // Add options to the display-reference select element
    const selectElement = document.getElementById('display-reference')
    selectElement.insertAdjacentHTML('beforeend', taskData.refOptions)
    selectElement.disabled=false;
    $('#display-reference').trigger('change')
    // Add selected references to the viewer
    const m = renderReferences(parsedData, viewer)
    if(m){
        viewer.zoomTo({model:m})
        viewer.render()
    }
    // Add event listeners to the selectReferences dropdown
    await selectReferences(parsedData, viewer);
    $('#display-reference').trigger('change')
 }

/***************************************************
 *                 Event Listeners                 *
 ***************************************************/

/**
 * Adds event listeners to a select-conformer dropdown for a particular CompoundOccurrence.
 * When a new conformer is selected, adds a model to the viewer with uniquely colored carbons
 * (if possible). When a conformer is de-selected, removes the model from the viewer and
 * frees up the color of the model's carbons to be used by the next model
 * @param {Number} coPk the PK of a CompoundOccurrence object
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the 3dmol viewer that is displayed 
 */
function selectConformers(coPk, parsedData, viewer){
    const querySelector = `.select-conformer#co-${coPk}`
    $(querySelector).on('select2:select', function (e) {
        var confId = e.params.data.id;
        const moltext = parsedData.compounds[this.id][confId].moltext
        const torsionAlerts = parsedData.compounds[this.id][confId].torsion_alerts
        addModelWrapper(viewer, moltext, confId, {torsionAlerts:torsionAlerts})
        viewer.render()
        // Resize the table header to fit the new table content
        $('#align-table').dataTable().fnAdjustColumnSizing();
    });
    $(querySelector).on('select2:unselect', function (e) {
        var confId = e.params.data.id;
        removeModelWrapper(viewer, confId)
        viewer.render()
        // Resize the table header to fit the new table content
        $('#align-table').dataTable().fnAdjustColumnSizing();
    });
    // Ensure that all tags are colored when the select element's options are updated.
    $(querySelector).on("change", function (e){
        colorTags(colorMap)
    })
}


/**
 * Adds an event listener to the select-references dropdown. When a new reference is selected,
 * a model is added to the viewer with uniquely colored carbons (if possible). When a pose
 * is de-selected, the model is removed from the viewer.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
function selectReferences(parsedData, viewer){
    $('#display-reference').on('select2:select', function (e) {
        var data = e.params.data;
        const moltext = parsedData.references[data.id]
        addModelWrapper(viewer, moltext, data.id, {})
        viewer.render()
    });
    $('#display-reference').on('select2:unselect', function (e) {
        var data = e.params.data;
        removeModelWrapper(viewer, data.id)
        viewer.render()
    });
    // Ensure that all tags are colored when the select element's options are updated.
    $('#display-reference').on("change", function (e){
        colorTags(colorMap)
    })
}

/**
 * Show/hide the receptor in the viewer
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
function toggleReceptor(viewer){
    $('#toggle-receptor').on("click", function (e) {
        if("receptor" in models){
            removeModelWrapper(viewer, 'receptor')
        }else{
            const residues = addReceptor(viewer)
            styleReceptorWrapper(viewer, residues)
        }
        viewer.setSlab(SLAB_NEAR,SLAB_FAR)
        viewer.render()
    })
}


/***************************************************
 *                     Helpers                     *
 ***************************************************/

/**
 * Adds models for the selected references to the viewer.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @returns the model added to the viewer
 */
 function renderReferences(parsedData, viewer){
    const selected = $('#display-reference').select2('data');
    var m = null
    for(var i=0; i<selected.length; i++){
        const moltext = parsedData.references[selected[i].id]
        m = addModelWrapper(viewer, moltext, selected[i].id, {})
    }
    return m
}

/**
 * Adds a receptor model to the viewer
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @returns an integer array of residues in the binding pocket for the selected receptor
 */
 function addReceptor(viewer){
    const receptorSeries = getReceptorSeries()
    const receptorMoltext = parsedData.receptors[receptorSeries].moltext
    addModelWrapper(viewer, receptorMoltext, "receptor", {format: "pdb"})
    viewer.render()
    return parsedData.receptors[receptorSeries].residues
}

 /**
  * Removes the receptor from the viewer (if one is displayed)
  * @param {Object} viewer the 3dmol viewer that is displayed
  */
  function removeReceptor(viewer){
    if("receptor" in models){
        removeModelWrapper(viewer, 'receptor')
    }
 }

/**
 * A wrapper around `removeModel` that passes the `colorFrequency`, `colorMap`,
 * and `models` objects in addition to the `viewer` and `id`
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {String} id the basechem identifier of this model (ex s-3 for series or 121-0 for a pose)
 */
function removeModelWrapper(viewer, id){
    removeModel(models, colorMap, colorFrequency, viewer, id)
    toggleReceptorButton(viewer)
}


/**
 * A wrapper around `addModel` that passes the `colorFrequency`, `colorMap`,
 * and `models` objects in addition to the required fields
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {String} moltext the moltext of the model to add
 * @param {String} id the basechem identifier of this model (ex s-3 for series or 121-0 for a pose)
 * @param {String} format "sdf" or "pdb", what is the format of the moltext
 * @param {Float} opacity the opacity of the model
 * @param {String} style "cartoon" or "stick", the atom style of the model
 * @returns the model object that was added
 */
function addModelWrapper(viewer, moltext, id, {format="sdf", opacity=PRIMARY_OPACITY, style=STYLE, torsionAlerts={}}= {}){
    const m = addModel(models, colorMap, colorFrequency, viewer, moltext, id,{format:format, opacity:opacity, style:style, torsionAlerts:torsionAlerts})
    toggleReceptorButton(viewer)
    return m
}

/**
 * A wrapper around `styleReceptor` that passes `models`, `SLAB_NEAR`, and `SLAB_FAR`
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {Array} residues an array of binding pocket residues (integers)
 */
 function styleReceptorWrapper(viewer, residues){
    styleReceptor(models, viewer, residues, SLAB_NEAR, SLAB_FAR)
}

/**
 * Given a model ID, return the ID of the most related series
 * @param {String} mId the basechem ID of a a 3Dmoljs model (ex: "s-12", "co-52-s-12")
 * @returns a string, the basechem ID of the series that corresponds to the given model (ex: "s-12")
 */
 function getSeriesForModel(mId){
    return mId.match(/s-\d+/)[0];
}

/**
 * Looks through all models currently in the viewer and returns the basechem ID of the
 * series model whose receptor can be shown in the viewer. A series ID is only returned
 * if all models in the viewer have the same series (either ARE the series or were aligned to the series)
 * @returns a string, the series ID (ex: "s-12") or null if no receptor can be shown
 */
 function getReceptorSeries(){
    // Collect all series relevant to the current models 
    const seriesSet = new Set()
    for(const mId in models){
        if(mId != "receptor"){
            seriesSet.add(getSeriesForModel(mId))
        }
    }
    // Return the series ID only if it's the only relevant series
    if(seriesSet.size == 1){
        const series = seriesSet.values().next().value
        if(series in parsedData.receptors){
            return series
        }
    }
    return null
 }

 /**
  * Disables/enables the "ToggleReceptor" based on the models currently in the viewer.
  * If all models have the same relevant series and that series has a receptor, the button
  * is enabled. If not, it is disabled.
  * @param {Object} viewer the 3dmol viewer that is displayed
  */
 function toggleReceptorButton(viewer){
    const receptorSeries = getReceptorSeries()
    const receptorBtn = document.getElementById("toggle-receptor")
    if(receptorSeries){
        receptorBtn.disabled = false
    }else{
        receptorBtn.disabled = true
        removeReceptor(viewer)
    }
 }


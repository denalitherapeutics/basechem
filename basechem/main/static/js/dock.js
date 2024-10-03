import { PRIMARY_OPACITY, STYLE, CARBON_COLORS, fitViewerToModels, download, removeModel, addModel, colorTags, styleReceptor, toggleTorsionAlerts, toggleToklatAnnotations} from './viewer.js'
import {saveViewerGems} from "./saveGems.js"
import {handleNewCompsAdded} from "./addComps.js"
import {copyDataId} from "../../../mni_common/static/mni_common/js/select2.js"
import {addSaveGemsModal, manageNewAnalysisForm, collectAnalysisTask, collectRefsTask} from "./tasks.js"

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
const SLAB_NEAR = -10
const SLAB_FAR = 10
/* parsedData is an object with all the results from the backend docking. It takes the form:
 * {compounds: {comp_id: {"moltext": moltext, "dockingScore": number}},
 *  references: {series_id: moltext},
 *  receptors: {series_id: {"moltext": moltext, "residues": int_array}}}
 */
var parsedData = {"compounds": {}, "references": {}, "receptors": {}}

/**
 * This function handles all of the front-end interactions for docking
 * after the backend docking process is complete
 */
 export default function(){
    const groupName = JSON.parse(document.getElementById("group-name").textContent);
    if(groupName){
        // Disable the new analysis form while tasks are running
        manageNewAnalysisForm()
        // Create 3DmolJS viewer to display compounds
        const viewer = $3Dmol.createViewer('dock-viewer');
        // Wait for dockrefs task to finish and add results to the dropdown when ready
        collectRefsTask(parsedData, viewer, 'dock', collectDockRefs);
        // Add event listeners to the download buttons
        download(models, '.dock-download')
        // Add event listener to the fit-viewer button
        fitViewerToModels(viewer, models, '#fit-viewer', SLAB_NEAR, SLAB_FAR)
        // Add event listener to the toggle-torsion-alerts button
        toggleTorsionAlerts(viewer, models, '#toggle-torsion-alerts')
        // Add event listener to the toggle-interactions button
        toggleToklatAnnotations(viewer, models, '#toggle-toklat-annotations')
        // Add event listener to the save-viewer-gems button
        saveViewerGems(models, '#save-viewer-gems')
        // Wait for each dock task to finish and add the results to the frontend when ready
        collectDockTasks(parsedData, viewer)    
    }
}

/**
 * For each `select-pose` widget in the table, wait for the associated DjangoQ task to complete
 * and display the results.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the displayed 3dmol viewer
 */
function collectDockTasks(parsedData, viewer){
    document.querySelectorAll(".select-pose").forEach(element => {
        const taskName = element.getAttribute("data-task-name")
        collectAnalysisTask(taskName, parsedData, viewer, collectSuccessfulDockTask) 
    })
    handleNewCompsAdded("dock-table", parsedData, viewer, newCompAddedHelper)
}

/**
 * Called in `collectAnalysisTask` when the task returns successfully - this function uses
 * the result from the backend to update parsedData, the viewer, and the DOM
 * @param {Number} coPk the PK of the CompoundOccurrence whose task returned successfully
 * @param {Object} ajaxData the data returned from the backend, including task results and HTML strings to add to the DOM
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the displayed 3dmol viewer 
 */
function collectSuccessfulDockTask(coPk, ajaxData, parsedData, viewer){
    parsedData.compounds[`co-${coPk}`] = ajaxData.taskResult
    const selectPoseElement = document.querySelector(`.select-pose#co-${coPk}`)
    selectPoseElement.insertAdjacentHTML('beforeend', ajaxData.confOptions);
    selectPoseElement.parentElement.style.display = ""
    addSaveGemsModal(coPk, ajaxData.saveGemsModal)
    selectPoses(coPk, parsedData, viewer)
}

function newCompAddedHelper(trElement, parsedData, viewer){
    // Collect the task results
    const taskName = trElement.querySelector(".select-pose").getAttribute("data-task-name")
    collectAnalysisTask(taskName, parsedData, viewer, collectSuccessfulDockTask)
    // Disable the row if it doesn't match the currently displayed receptor
    const receptorId = $('#display-receptor').select2('data')[0].id
    toggleAvailableCompounds(viewer, receptorId)
    $(".select-pose").select2({
        placeholder: "Poses",
        width: "100%",
        templateResult: copyDataId,
        templateSelection: copyDataId,
      })
}

/**
 * Waits for the dockrefs DjangoQ task to finish. Once it is complete, populates the reference dropdown
 * with results. If the task is not successful, the reference dropdown will remain disabled.
 * @param {Object} taskData a dictionary with the result of the dockrefs task
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data.
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
async function collectDockRefs(taskData, parsedData, viewer){
    // Add taskData to parsedData so that models can be added to the viewer
    parsedData["references"] = taskData.taskResult.references
    parsedData["receptors"] = taskData.taskResult.receptors
    // Add options to references select element
    var selectElement = document.getElementById('display-reference')
    selectElement.insertAdjacentHTML('beforeend', taskData.refOptions)
    selectElement.disabled=false;
    $('#display-reference').trigger('change')
    // Add options to receptor select element
    selectElement = document.getElementById('display-receptor')
    selectElement.insertAdjacentHTML('beforeend', taskData.recOptions)
    selectElement.disabled=false;
    $('#display-receptor').trigger('change')
    // Display receptor
    const residues = renderReceptor(parsedData, viewer)
    const receptorId = $('#display-receptor').select2('data')[0].id
    styleReceptorWrapper(viewer, residues)
    toggleAvailableCompounds(viewer, receptorId)
    // Display the matching reference
    const m = matchReferenceToReceptor(viewer, parsedData, receptorId)
    viewer.zoomTo({model: m})
    viewer.setSlab(SLAB_NEAR,SLAB_FAR)
    viewer.render()
    // Add event listeners to "Select References" dropdown
    selectReferences(parsedData, viewer)
    // Add event listeners to "Select Receptor" dropdown
    await selectReceptor(parsedData, viewer)
    $('#display-reference').trigger('change')
}

/***************************************************
 *                 Event Listeners                 *
 ***************************************************/

/**
 * Adds event listeners to the select-pose dropdown for a particular CompoundOccurrence.
 * When a new pose is selected, adds a model to the viewer with uniquely colored carbons
 * (if possible). When a pose is de-selected, removes the model from the viewer and frees
 * up the color of the model's carbons to be used by the next model
 * @param {Number} coPk the PK of a CompoundOccurrence object
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the 3dmol viewer that is displayed 
 */
 function selectPoses(coPk, parsedData, viewer){
    const querySelector = `.select-pose#co-${coPk}`
    $(querySelector).on('select2:select', function (e) {
        var confId = e.params.data.id;
        const moltext = parsedData.compounds[this.id][confId].moltext
        const torsionAlerts = parsedData.compounds[this.id][confId].torsionAlerts
        const toklatAnnotations = parsedData.compounds[this.id][confId].toklatAnnotations
        addModelWrapper(viewer, moltext, confId, {torsionAlerts:torsionAlerts, toklatAnnotations:toklatAnnotations})
        viewer.setSlab(SLAB_NEAR,SLAB_FAR)
        viewer.render()
        // Resize the table header to fit the new table content
        $('#dock-table').dataTable().fnAdjustColumnSizing();
    });
    $(querySelector).on('select2:unselect', function (e) {
        var confId = e.params.data.id;
        removeModelWrapper(viewer, confId)
        viewer.setSlab(SLAB_NEAR,SLAB_FAR)
        viewer.render()
        // Resize the table header to fit the new table content
        $('#dock-table').dataTable().fnAdjustColumnSizing();
    });
    // Ensure that all tags are colored when the select element's options are updated.
    $(querySelector).on("change", function (e){
        colorTags(colorMap)
    })
}

/**
 * Adds an event listener to the select-receptor dropdown. When a receptor is selected,
 * the model is added to the viewer and the previous receptor is removed. If the default
 * series was used for docking, then selecting a receptor also disables all compounds
 * in the table that were not docked using the selected receptor.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
 function selectReceptor(parsedData, viewer){
    //render selected receptor
    $('#display-receptor').on('select2:select', function (e) {
        const residues = renderReceptor(parsedData, viewer)
        styleReceptorWrapper(viewer, residues)
        // if default, disable any compounds w/ a different series
        const receptorId = $('#display-receptor').select2('data')[0].id
        toggleAvailableCompounds(viewer, receptorId)
        matchReferenceToReceptor(viewer, parsedData, receptorId)
        viewer.setSlab(SLAB_NEAR,SLAB_FAR)
        viewer.render()
    });
    $('#display-receptor').on('select2:unselect', function (e) {
        removeModelWrapper(viewer, 'receptor')
        toggleAvailableCompounds(viewer, null)
        viewer.setSlab(SLAB_NEAR,SLAB_FAR)
        viewer.render()
    });
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
        addModelWrapper(viewer, moltext, data.id)
        viewer.setSlab(SLAB_NEAR,SLAB_FAR)
        viewer.render()
    });
    $('#display-reference').on('select2:unselect', function (e) {
        var data = e.params.data;
        removeModelWrapper(viewer, data.id)
        viewer.setSlab(SLAB_NEAR,SLAB_FAR)
        viewer.render()
    });
    // Ensure that all tags are colored when the select element's options are updated.
    $('#display-reference').on("change", function (e){
        colorTags(colorMap)
    })
}


/***************************************************
 *                     Helpers                     *
 ***************************************************/

/**
 * Displays the reference compound with the same id as the currently displayed receptor.
 * This function is called when a new receptor is chosen from the dropdown, so that the reference
 * automatically changes with it.
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {String} receptorId the id (ex. s-3) of series whose receptor is displayed
 * @returns the model added to the viewer
 */
function matchReferenceToReceptor(viewer, parsedData, receptorId){
    // remove currently selected references from viewer
    const references = $('#display-reference').select2('data')
    for(var i=0; i< references.length; i++){
        if(references[i].id in models){
            removeModelWrapper(viewer,references[i].id)
        }
    }
    // update the `display-reference` select2 widget to select the new reference
    $('#display-reference').val([receptorId]).trigger("change")
    return renderReferences(parsedData, viewer)
}

/**
 * Disables/Enables compounds in the table based on the selected receptor. If the default
 * series was used for each compound, only enable compounds that were docked with the 
 * selected receptor. 
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {String} receptorId the id (ex. s-3) of series whose receptor is displayed
 */
function toggleAvailableCompounds(viewer, receptorId){
    const groupNameContainer = document.getElementById("group-name");
    if(groupNameContainer){
        const rows = document.querySelectorAll(`.compound-row`)
        const isDefault = groupNameContainer.innerHTML.includes("default")
        for(var i=0; i < rows.length; i++){
            if(isDefault & rows[i].getAttribute('series') !== receptorId){
                // disable row
                rows[i].style.opacity=0.4;
                toggleDisableChildren(rows[i], true, false)
                // If there is a select-poses dropdown, remove selected poses from the viewer
                // and clear the select widget
                const selectPoseWrapper = document.querySelector(`.select-pose-wrapper#${rows[i].id}`)
                if(selectPoseWrapper.style.display != "none"){
                    const poses = $(`.select-pose#${rows[i].id}`).select2('data')
                    for(var j=0; j< poses.length; j++){
                        removeModelWrapper(viewer,poses[j].id)
                    }
                    $(`.select-pose#${rows[i].id}`).val([]).trigger("change")
                }
            }else{
                // enable row
                rows[i].style.opacity=1;
                toggleDisableChildren(rows[i], false, false)
            }
        }
    }
}

/**
 * Adds the model for the selected receptor to the viewer.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @returns an integer array of residues in the binding pocket for the selected receptor
 */
 function renderReceptor(parsedData, viewer){
    const selected = $('#display-receptor').select2('data');
    const receptorId = selected[0].id
    const moltext = parsedData.receptors[receptorId].moltext
    const residues = parsedData.receptors[receptorId].residues
    addModelWrapper(viewer, moltext, "receptor", {format:'pdb'})
    return residues
}



/**
 * Adds models for the selected references to the viewer.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @returns the model added to the viewer
 */
 function renderReferences(parsedData, viewer){
    const selected = $('#display-reference').select2('data');
    var m = {}
    for(var i=0; i<selected.length; i++){
        const moltext = parsedData.references[selected[i].id]
        m = addModelWrapper(viewer, moltext, selected[i].id)
    }
    return m
}

/**
 * A wrapper around `removeModel` that passes the `colorFrequency`, `colorMap`,
 * and `models` objects in addition to the `viewer` and `id`
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {String} id the basechem identifier of this model (ex s-3 for series or 121-0 for a pose)
 */
 function removeModelWrapper(viewer, id){
    removeModel(models, colorMap, colorFrequency, viewer, id)
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
 * @param {Object} torsionAlerts a dictionary with torsion alert data for this model
 * @param {Object} toklatAnnotations a dictionary with toklat annotation data for this model
 * @returns the model object that was added
 */
function addModelWrapper(viewer, moltext, id, {format="sdf", opacity=PRIMARY_OPACITY, style=STYLE, torsionAlerts={}, toklatAnnotations={"cylinders": [], "spheres": []}} = {}){
    const m = addModel(models, colorMap, colorFrequency, viewer, moltext, id, {format: format, opacity: opacity, style:style, torsionAlerts:torsionAlerts, toklatAnnotations:toklatAnnotations})
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
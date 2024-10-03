import {PRIMARY_OPACITY, CARBON_COLORS, STYLE, download, removeModel, addModel, colorTags, fitViewerToModels} from './viewer.js'
import {handleNewCompsAdded} from "./addComps.js"
import {manageNewAnalysisForm, collectAnalysisTask, collectRefsTask} from './tasks.js';

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
const models = {"single": {}, "0": {}, "1":{}, "2": {}, "3":{}}
const surfaces = {"single": {}, "0": {}, "1":{}, "2": {}, "3":{}}
const SLAB_NEAR = -10
const SLAB_FAR = 10
var showReceptors = false
/* parsedData is an object with all the results from the backend esp map generation. It takes the form:
 * {compounds: {comp_id: {pqr, related_series}},
 * references: {series_id: {pqr, dx}},
 * receptors: {series_id: {pqr, residues}}}
 */
var parsedData = {"compounds": {}, "references": {}, "receptors": {}}


/**
 * This function handles all of the front-end interactions after the backend ESP process is complete
 */
 export default function(){
    const groupName = JSON.parse(document.getElementById("group-name").textContent);
    if(groupName){
        // Disable the new analysis form while tasks are running
        manageNewAnalysisForm()
        // Create 3DmolJS viewers to display compounds
        const viewers = {"single": $3Dmol.createViewer('esp-viewer')}
        const backgroundColors = ["#ededed","#ffffff","#ffffff","#ededed"]
        for(let i=0; i<4;i++){
            viewers[i]= $3Dmol.createViewer(`esp-viewer-${i}`, {backgroundColor:backgroundColors[i]})
        }
        // Wait for esprefs task to finish and add results to the dropdown when ready
        collectRefsTask(parsedData, viewers, 'esp', collectEspRefs);
        // Add event listener to the "Toggle Viewer" button
        toggleViewer();
        // Add event listeners to the download buttons
        download(models["single"], '.esp-download');
        // Add event listener to the surface color max input
        toggleSurfaceColors(parsedData, viewers)
        // Add event listeners to the fit-viewer button
        fitViewerToModels(viewers["single"], models["single"], '#fit-viewer', SLAB_NEAR, SLAB_FAR)
        fitViewerToModels(viewers["0"], models["0"], '#fit-viewer', SLAB_NEAR, SLAB_FAR)
        fitViewerToModels(viewers["1"], models["1"], '#fit-viewer', SLAB_NEAR, SLAB_FAR)
        fitViewerToModels(viewers["2"], models["2"], '#fit-viewer', SLAB_NEAR, SLAB_FAR)
        fitViewerToModels(viewers["3"], models["3"], '#fit-viewer', SLAB_NEAR, SLAB_FAR)
        // Wait for each ESP task to finish and add the results to the frontend when ready
        collectEspTasks(parsedData, viewers)
    }
       
}

/**
 * For each `toggle-compound` widget in the table, wait for the associated DjangoQ task to complete
 * and display the results.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 */
function collectEspTasks(parsedData, viewers){
    document.querySelectorAll(".toggle-compound").forEach(element => {
        const taskName = element.getAttribute("data-task-name")
        collectAnalysisTask(taskName, parsedData, viewers, collectSuccessfulEspTask)
    })
    handleNewCompsAdded("esp-table", parsedData, viewers, newCompAddedHelper)
}

/**
 * Used in handleNewCompsAdded, this function collects task results for new compounds.
 * @param {Object} trElement a table row element for a newly added compound
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 */
function newCompAddedHelper(trElement, parsedData, viewers){
    const taskName = trElement.querySelector(".toggle-compound").getAttribute("data-task-name")
    collectAnalysisTask(taskName, parsedData, viewers, collectSuccessfulEspTask)
}

/**
 * Called in `collectAnalysisTask` when the task returns successfully - this function uses
 * the result from the backend to update parsedData, the viewer, and the DOM.
 * @param {Number} coPk the PK of a CompoundOccurrence object
 * @param {Object} ajaxData the data returned from the backend, including task results and HTML strings to add to the DOM
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 */
function collectSuccessfulEspTask(coPk, ajaxData, parsedData, viewers){
    // Add task result to parseData
    parsedData.compounds[`co-${coPk}`] = ajaxData.taskResult
    // Show the toggle button and add its event listener
    toggleCompound(coPk, parsedData, viewers)
    toggleOpacity(`co-${coPk}`, parsedData, viewers)
    // Resize the table header to fit the new table content
    $('#esp-table').dataTable().fnAdjustColumnSizing();
    const loadingModalShowing = document.getElementById("loading-modal").classList.contains("show")
    const showHelpModal = JSON.parse(document.getElementById("show-help-modal").textContent);
    if(loadingModalShowing && showHelpModal){
        $('#esp-help-modal').modal('show');
    }
}

/**
 * Waits for the esprefs DjangoQ task to finish. Once it is complete, populates the reference dropdown
 * with results. If the task is not successful, the reference dropdown will remain disabled.
 * @param {Object} taskData a dictionary with the result of the esprefs task
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data.
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
async function collectEspRefs(taskData, parsedData, viewers){
    // Add taskData to parsedData so that models can be added to the viewer
    parsedData["references"] = taskData.taskResult.references
    parsedData["receptors"] = taskData.taskResult.receptors
    // Add options to display-reference select element
    const selectElement = document.getElementById('display-reference')
    selectElement.insertAdjacentHTML('beforeend', taskData.refOptions)
    selectElement.disabled=false;
    $('#display-reference').trigger('change')
    // Add selected references to the viewer
    renderReferences(parsedData, viewers)
    // Add event listeners to the selectReferences dropdown
    await selectReferences(parsedData, viewers);
    // Add event listener to the reference opacity input
    toggleOpacity('reference-opacity', parsedData, viewers)
    $('#display-reference').trigger('change')
    // Add event listener to the "Toggle Receptors" button
    toggleReceptors(parsedData, viewers)
}

/***************************************************
 *                 Event Listeners                 *
 ***************************************************/

/**
 * Adds an event listeners to the toggle-compound button for a particular CompoundOccurrence.
 * When the button is clicked, this adds a model to the viewer with uniquely colored carbons (if possible).
 * When the button is clicked again, this removes the model from the viewer and frees up the color of
 * the model's carbons to be used by the next model
 * @param {Number} coPk the pk of the CompoundOccurrence object
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 */
function toggleCompound(coPk, parsedData, viewers){
    const toggleCompoundBtn = document.querySelector(`.toggle-compound#co-${coPk}`)
    toggleCompoundBtn.parentElement.style.display = ""
    const opacityInput = document.querySelector(`.toggle-opacity#co-${coPk}`)
    opacityInput.style.display = ""
    $(`.toggle-compound#co-${coPk}`).on('click', function (e) {
        const coData = parsedData.compounds[this.id]
        if(this.id in models['single']){
            removeModelWrapper(viewers, this.id)
            // Remove background color from button
            const tag = document.querySelector(`button.toggle-compound#${this.id}`)
            tag.style.backgroundColor = "";
            this.innerText = "Show"
        }else{
            const added = addModelWrapper(parsedData, viewers,coData.pqr, coData.dx, this.id)
            if(added){
                this.innerText="Hide"
            }
        }
    });
}

/**
 * Changes the opacity of a surface based on changes to the user input
 * @param {String} elementId the ID of a toggle-opacity input
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 */
function toggleOpacity(elementId, parsedData, viewers){
    $(`.toggle-opacity#${elementId}`).on('change', function (e) {
        if(this.id in surfaces['single']){
            removeModelWrapper(viewers, this.id)
            const coData = parsedData.compounds[this.id]
            addModelWrapper(parsedData, viewers, coData.pqr, coData.dx, this.id)
        }else if(this.id == "reference-opacity"){
            // Update opacity of all references
            for(let surfaceId in surfaces['single']){
                if(surfaceId.includes("s-")){
                    removeModelWrapper(viewers, surfaceId)
                    addModelWrapper(parsedData, viewers, parsedData.references[surfaceId], "", surfaceId)
                }
            }
        }
    });
}


/**
 * Changes the coloring boundaries for all surfaces based on the user input
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 */
function toggleSurfaceColors(parsedData, viewers){
    $("#toggle-surface-color-max").on('change', function(){
        for(const surfaceId in surfaces["single"]){
            let data = {}
            if(surfaceId.includes("s-")){
                data = parsedData.references[surfaceId]
            }else{
                data = parsedData.compounds[surfaceId]
            }
            removeModelWrapper(viewers, surfaceId)
            addModelWrapper(parsedData, viewers, data.pqr, data.dx, surfaceId)
        }
    })
    
}

/**
 * Toggles which viewer (single or quad) is currently visible
 */
function toggleViewer(){
    $('#toggle-viewer').on("click", function(e){
        const single = document.getElementById("esp-viewer")
        const quad = document.getElementById("quad-esp-viewer-container")
        const toggleSingleOn = single.style.visibility == "hidden"
        single.style.visibility = toggleSingleOn ? "" : "hidden"
        quad.style.visibility = toggleSingleOn ? "hidden" : ""
    })
}

/**
 * If `showReceptors` is true, changes it to false and removes all receptors. If
 * `showReceptors` is false, changes it to true and adds all receptors
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 */
function toggleReceptors(parsedData, viewers){
    $('#toggle-receptors').on("click", function(e){
        showReceptors = !showReceptors
        if(showReceptors){
            for(const viewerId in models){
                addReceptor(parsedData, viewers, viewerId)
            }
        }else{
            for(const viewerId in models){
                removeReceptor(viewers, viewerId)
            }
        }
    })
}


/**
 * Adds an event listener to the select-references dropdown. When a new reference is selected,
 * a model is added to the viewer with uniquely colored carbons (if possible). When a pose
 * is de-selected, the model is removed from the viewer.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 */
function selectReferences(parsedData, viewers){
    $('#display-reference').on('select2:select', function (e) {
        var data = e.params.data;
        const refData = parsedData.references[data.id]
        const added = addModelWrapper(parsedData, viewers, refData.pqr, refData.dx, data.id)
        if(!added){
            // Model was not added (too many models showing), so de-select the selected value
            const vals = $('#display-reference').select2('data');
            var new_vals = []
            for(var i=0; i<vals.length;i++){
                if(vals[i].id != data.id){
                    new_vals.push(vals[i].id)
                }
            }
            $('#display-reference').val(new_vals);
            $('#display-reference').trigger('change');
            // Update the tag colors
            colorTags(colorMap)
        }
    });
    $('#display-reference').on('select2:unselect', function (e) {
        var data = e.params.data;
        removeModelWrapper(viewers, data.id)
    });
    // Ensure that all tags are colored when the select element's options are updated.
    $('#display-reference').on("change", function (e){
        colorTags(colorMap)
    })
}

/***************************************************
 *                    Surfaces                     *
 ***************************************************/


/**
 * Given a 3dmol object, return the gradient object that should be used to color its ESP surface
 * @param {Object} model the 3dmol model whose surface gradient is being constructed
 * @returns a 3dmol red to blue gradient object, normalized to the given model
 */
 function surfGradient(model){
    const state = model.getInternalState()
    var min = 0
    var max = 0
    for(var i=0; i<state.atoms.length; i++){
        const charge = state.atoms[i].properties.partialCharge
        if(charge < min){
            min = charge
        }else if(charge > max){
            max = charge
        }
    }
    return new $3Dmol.Gradient.RWB(min,max)
}

/**
 * Adds a surface to the viewer
 * @param {Object} surfaces a dictionary mapping model ids to their 3dmol surface
 * @param {Object} viewer the 3dmol viewer to add this surface to
 * @param {String} id the basechem identifier of the compound whose surface is being added
 * @param {Object} model the 3dmol model whose surface this is
 * @returns a promise to the surface that was added to the viewer
 */
function addSurface(surfaces, viewer, id, model, dxData){
    // Determine opacity
    let opacity = document.querySelector('.toggle-opacity#reference-opacity').value
    if(document.querySelector(`.toggle-opacity#${id}`)){
        opacity = document.querySelector(`.toggle-opacity#${id}`).value
    }
    const max = Number(document.querySelector('#toggle-surface-color-max').value)
    const min = -max
    // Create new surface
    const styleSpec = {opacity:opacity, voldata: dxData, volformat: 'dx', volscheme: {gradient:'rwb', min:min, max:max, mid:0}} 
    const surfPromise = viewer.addSurface($3Dmol.SurfaceType.VDW, styleSpec, {model: model.getID()})
    
    // Remove existing surface if it exists
    if(surfaces[id]){viewer.removeSurface(surfaces[id])}
    // Add surface promise to dictionary
    surfaces[id] = surfPromise;
    return surfPromise;
}

/**
 * Removes a surface from the viewer
 * @param {Object} surfaces a dictionary mapping model ids to their 3dmol surface
 * @param {Object} viewer the 3dmol viewer to remove the surface from
 * @param {String} id the basechem identifier of the compound whose surface is being removed
 */
function removeSurface(surfaces, viewer, id){
    // Remove model from viewer and delete from models dict
    viewer.removeSurface(surfaces[id].surfid)
    delete surfaces[id];
}

/***************************************************
 *                     Helpers                     *
 ***************************************************/

/**
 * Adds models for the selected references to the viewer.
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 */
 function renderReferences(parsedData, viewers){
    const selected = $('#display-reference').select2('data');
    for(var i=0; i<selected.length; i++){
        const refData = parsedData.references[selected[i].id]
        addModelWrapper(parsedData, viewers, refData.pqr, refData.dx, selected[i].id)
    }
}

/**
 * A wrapper around `removeModel` that passes the `colorFrequency`, `colorMap`,
 * and `models` objects in addition to the `viewer` and `id`
 * @param {Object} viewers a dictionary of 3dmol viewers
 * @param {String} id the basechem identifier of this model (ex s-3 for series or 121-0 for a pose)
 */
function removeModelWrapper(viewers, id){
    let quadViewerId = null
    for(const viewerId in models){
        if(viewerId != 'single' & id in models[viewerId]){
            quadViewerId = viewerId
        }
    }
    if(quadViewerId){
        // Remove single
        removeModel(models['single'], colorMap, colorFrequency, viewers.single, id)
        removeSurface(surfaces['single'], viewers.single, id)
        viewers["single"].render()
        // Remove quad
        removeModel(models[quadViewerId], colorMap, colorFrequency, viewers[quadViewerId], id)
        removeSurface(surfaces[quadViewerId], viewers[quadViewerId], id)
        removeReceptor(viewers,quadViewerId)
        if("receptor" in models[quadViewerId]){
            // Remove receptor if it exists
            removeModel(models[quadViewerId], colorMap, colorFrequency, viewers[quadViewerId], "receptor")
        }
        viewers[quadViewerId].render()
    }
}


/**
 * A wrapper around `addModel` that passes the `colorFrequency`, `colorMap`,
 * and `models` objects in addition to the required fields
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 * @param {String} pqr the pqr data of the model to add
 * @param {String} dx the dx data that should be used to color the model's surface
 * @param {String} id the basechem identifier of this model (ex s-3 for series or 121-0 for a pose)
 * @param {String} format "sdf" or "pdb" or "pqr", what is the format of the moltext
 * @param {Float} opacity the opacity of the model
 * @param {String} style "cartoon" or "stick", the atom style of the model
 * @returns true if the model was added successfully, false if it failed
 */
function addModelWrapper(parsedData, viewers, pqr, dx, id, {format="pqr", opacity=PRIMARY_OPACITY, style=STYLE} = {}){
    let quadViewerId = null
    for(const viewerId in models){
        if(viewerId != 'single' & Object.keys(models[viewerId]).length == 0){
            quadViewerId = viewerId
            break
        }
    }
    if(quadViewerId){
        // Add single
        const singleM = addModel(models['single'], colorMap, colorFrequency, viewers.single, pqr, id, {format:format, opacity:opacity, style:style})
        addSurface(surfaces['single'], viewers.single, id, singleM, dx)
        // Add quad
        const quadM = addModel(models[quadViewerId], colorMap, colorFrequency, viewers[quadViewerId], pqr, id, {format:format, opacity: opacity, style:style})
        addSurface(surfaces[quadViewerId], viewers[quadViewerId], id, quadM, dx)
        if(showReceptors){
            addReceptor(parsedData, viewers, quadViewerId)
        }
        viewers[quadViewerId].zoomTo({model:quadM})
        viewers[quadViewerId].setSlab(SLAB_NEAR,SLAB_FAR)
        viewers[quadViewerId].render()
        if(Object.keys(models.single).length == 1){
            // When adding the first model to the viewer, zoom and set slab
            viewers.single.zoomTo({model:singleM})
            viewers.single.setSlab(SLAB_NEAR,SLAB_FAR)
        }
        viewers.single.render()
        return true
    }else{
        alert("You can only view 4 ESP maps at a time. Hide an existing map before showing this one")
        return false
    }
}

/**
 * Given an array of atoms, create a 3DMol gradient object that goes from red (negative)
 * to blue (positive), using the partial charge to determine the maximum and minimum values
 * @param {Array} atoms a list of 3DMol atoms
 */
 function createRWBGradient(atoms){
    let minCharge = 0
    let maxCharge = 0
    for(const atom of atoms){
        const charge = atom.properties.partialCharge
        if(charge < minCharge){
            minCharge = charge
        }
        if(charge > maxCharge){
            maxCharge = charge
        }
    }
    return new $3Dmol.Gradient.RWB(minCharge, maxCharge)
}

/**
 * Given a viewer ID, this function determines which series is most relevant for the viewer. The most relevant
 * series is the one that appears the most in the viewer (either the series itself, or a compound belonging to that series)
 * AND also has a valid receptor
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {String} viewerId "0" - "3" or "single", the viewer whose series is being picked
 * @returns a string, the series for this viewer
 */
 function chooseSeriesForViewer(parsedData, viewerId){
    const modelsInViewer = Object.keys(models[viewerId])
    if(modelsInViewer.length  == 0){
        return ""
    }
    const seriesCount = {}
    // Count how many times each series appears
    for(let mIndex = 0; mIndex < modelsInViewer.length; mIndex++){
        const cId = modelsInViewer[mIndex]
        var series = cId
        // If the model isn't itself a series, get its related series
        if(!cId.includes("s-")){
            series = parsedData.compounds[cId].related_series
        }
        // If a series was found and it has a receptor, count it
        if(series && (series in parsedData.receptors)){
            if(!(series in seriesCount)){
                seriesCount[series] = 0
            }
            seriesCount[series] += 1
        }
    }
    // Return the series that appears the most
    const sortedSeries = Object.entries(seriesCount).sort((x, y) => y[1] - x[1])
    if(sortedSeries.length > 0){
        return sortedSeries[0][0]
    }
    return ""
}


/***************************************************
 *                    Receptors                    *
 ***************************************************/

/**
 * Adds a receptor to the viewer with the given `viewerId`
 * @param {Object} parsedData a dictionary of the form {"compounds": {}, "receptors": {}, "references": {}} with molecule data
 * @param {Object} viewers a dictionary of 3dmol viewers
 * @param {String} viewerId "0" - "3" or "single", the viewer that should have a receptor added
 * @returns the seriesId of the added receptor ("" if no receptor added)
 */
function addReceptor(parsedData, viewers, viewerId){
    // Pick which series to use
    const seriesId = chooseSeriesForViewer(parsedData, viewerId)
    if(!seriesId){
        // No series: don't add a receptor
        return ""
    }
    // Add the receptor
    const pqr = parsedData.receptors[seriesId].pqr
    const m = addModel(models[viewerId], colorMap, colorFrequency, viewers[viewerId], pqr, "receptor",{format:"pqr"})
    // Create styles for the pocket
    const pocket = {resi:parsedData.receptors[seriesId].residues}
    const surfaceHs = {elem: ["Hsh","Hlp","Hh","Hp", "Hl", "Hs"]}
    const pocketStyle = {sphere:{hidden:false, opacity:0.7}}
    const hiddenStyle = {stick:{hidden:true}, cartoon:{hidden:true}}

    // Apply styles
    m.setStyle({}, hiddenStyle);
    m.setStyle(pocket, pocketStyle);
    m.setStyle(surfaceHs, hiddenStyle)
    // Apply gradient based on partial charge
    const gradient = createRWBGradient(m.selectedAtoms(pocket))
    m.setColorByFunction(pocket, (atom) => gradient.valueToHex(atom.properties.partialCharge))

    viewers[viewerId].render()
    return seriesId
}

/**
 * Removes the receptor from the viewer with the given "viewerId"
 * @param {Object} viewers a dictionary of 3dmol viewers
 * @param {String} viewerId "0" - "3" or "single", the viewer that should have a receptor removed
 */
function removeReceptor(viewers, viewerId){
    if("receptor" in models[viewerId]){
        removeModel(models[viewerId], colorMap, colorFrequency, viewers[viewerId], "receptor")
        viewers[viewerId].render()
    }
}

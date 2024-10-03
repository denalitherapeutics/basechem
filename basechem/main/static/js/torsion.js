import { PRIMARY_OPACITY, STYLE, CARBON_COLORS, pickCarbonColor, download, removeModel, addModel, ATOM_COLORS, fitViewerToModels, getTorsionConfId} from './viewer.js'
import {saveViewerGems} from "./saveGems.js"
import {handleNewCompsAdded} from "./addComps.js"
import {addSaveGemsModal, collectAnalysisTask} from "./tasks.js"

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
const SELECTED_COLOR = 0xFC03DF
const ACTIVE_POINT_SIZE = 10
const INACTIVE_POINT_SIZE = 3
/* parsedData is an object with all the results from the backend torsion scan. It takes the form:
 * {co-id: {
 *     "delta_energy": Number,
 *     "initial_dihedral": Number,
 *     "torsions":{conf-id:{
 *         "moltext": String,
 *         "rel_energy": Number,
 *         "dihedral": Number,
 *      }}
 * }}
 */
var parsedData = {}

/**
 * This function handles all of the front-end interactions for torsion scans
 * after the backend torsion scan process is complete
 */
export default function(){
    const groupName = JSON.parse(document.getElementById('group-name').textContent);
    if(groupName){
        // TORSION RESULTS
        // Create 3DmolJS viewer to display compounds
        const viewer = $3Dmol.createViewer('torsion-viewer', {backgroundColor:"#ededed"});
        // Wait for each torsion task to finish and add the results to the frontend when ready
        collectTorsionTasks(parsedData, viewer)  
        // Setup initial blank chart
        setupChart(parsedData, viewer)
        // Add event listener to the "Show Lowest Energy" buttons
        showLowestEnergy(parsedData,viewer)
        // Add event listener to the "Match Dihedral" button
        showMatchDihedral(parsedData,viewer)
        // Add event listeners to the download buttons
        download(models, '.torsion-download')
        // Add event listener to the save-viewer-gems button
        saveViewerGems(models, '#save-viewer-gems')
        // Add event listener to the fit-viewer button
        fitViewerToModels(viewer, models, '#fit-viewer')
    }else{
        // TORSION SUBMISSION
        const viewer = $3Dmol.createViewer('torsion-viewer');
        // Add selected pioneer to the viewer
        addPioneerModel(viewer)
        // Add event listener to the choose pioneer dropdown
        choosePioneer(viewer)
        // Add event listener to the fit-viewer button
        fitViewerToModels(viewer, models, '#fit-viewer')
    }
}

/***************************************************
 *                Torsion Submission               *
 ***************************************************/

/**
 * Adds an event listener to the `choose-pioneer` dropdown so that the chosen pioneer
 * is the one displayed in the viewer
 * @param {Object} viewer the 3dmol viewer that is displayed 
 */
function choosePioneer(viewer){
    $('#choose-pioneer').on('change', function(){addPioneerModel(viewer)});
}

/**
 * Adds the pioneer chosen in the `choose-pioneer` dropdown to the viewer
 * @param {Object} viewer the 3dmol viewer that is displayed 
 */
function addPioneerModel(viewer){
    // Remove any existing models from the viewer
    for(var comp in models){
        removeModelWrapper(viewer, comp)
    }
    // Add the current model to the viewer
    const dropdown = $('#choose-pioneer')[0]
    const moltext = dropdown.options[dropdown.selectedIndex].getAttribute("moltext")
    const m = addModelWrapper(viewer, moltext, dropdown.value)
    // Make it clickable
    viewer.setClickable({}, true, clickAtom)
    // Render
    viewer.zoomTo({model: m})
    viewer.setSlab(SLAB_NEAR,SLAB_FAR)
    viewer.render()
}

/**
 * This function runs each time an atom is clicked in the viewer. If not currently selected, the atom
 * will be selected and highlight SELECTED_COLOR. If already selected, the atom will be deselected
 * and will change back to it's original color. This function also enables/disables the submit button
 * depending on how many atoms are selected.
 * @param {Object} atom the 3dmoljs atom object that was clicked
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
function clickAtom(atom, viewer){
    const m = viewer.getModel(atom.m)
    // Make sure there are not already 4 atoms selected
    var numSelected = m.selectedAtoms({selected:true}).length
    if(!atom.selected && numSelected >=4){
        alert("You can only select 4 atoms for a torsion scan! De-select one of your atoms before selecting another.")
        return
    }
    atom.selected = !atom.selected

    // Update color of clicked atom
    m.setColorByFunction({serial:atom.serial},function(atom){
        if(atom.selected){
            return SELECTED_COLOR;
        }else{
            // Return to original color
            const atomColors = ATOM_COLORS
            atomColors['C'] = colorMap[$('#choose-pioneer')[0].value];
            return atomColors[atom.elem];
        }
    })
    viewer.render()

    // Enable/disable the submit button depending on how many atoms are selected
    const selectedAtoms = m.selectedAtoms({selected:true})
    if(selectedAtoms.length == 4){
        // Update the `selected_atoms` input with the currently selected atom numbers
        var atomNums = []
        for(var i=0; i<selectedAtoms.length; i++){
            const atom = selectedAtoms[i]
            atomNums.push(`${atom.serial}`)
        }
        $("input[name='selected_atoms']")[0].value = atomNums.join(",")
        // Enable submit
        $("#torsion-btn")[0].disabled = false
    }else{
        $("#torsion-btn")[0].disabled = true
    }
}


/***************************************************
 *                 Torsion Results                 *
 ***************************************************/

/**
 * For each `toggle-compound` button in the table, wait for the associated DjangoQ task to complete
 * and display the results.
 * @param {Object} parsedData output from the backend torsion scan
 * @param {Object} viewer the displayed 3dmol viewer
 */
function collectTorsionTasks(parsedData, viewer){
    document.querySelectorAll(".toggle-compound").forEach(element => {
        const taskName = element.getAttribute("data-task-name")
        collectAnalysisTask(taskName, parsedData, viewer, collectSuccessfulTorsionTask)
    })
    handleNewCompsAdded("torsion-table", parsedData, viewer, newCompAddedHelper)
}

/**
 * Used in handleNewCompsAdded, this function collects task results for new compounds.
 * @param {Object} trElement a table row element for a newly added compound
 * @param {Object} parsedData output from the backend torsion scan
 * @param {Object} viewer the displayed 3dmol viewer
 */
function newCompAddedHelper(trElement, parsedData, viewer){
    const taskName = trElement.querySelector(".toggle-compound").getAttribute("data-task-name")
    collectAnalysisTask(taskName, parsedData, viewer, collectSuccessfulTorsionTask)
}

/**
 * Called in `collectAnalysisTask` when the task returns successfully - this function uses
 * the result from the backend to update parsedData, the viewer, and the DOM.
 * @param {Number} coPk the PK of a CompoundOccurrence object
 * @param {Object} ajaxData the data returned from the backend, including task results and HTML strings to add to the DOM
 * @param {Object} parsedData output from the backend torsion scan
 * @param {Object} viewer the displayed 3dmol viewer
 */
function collectSuccessfulTorsionTask(coPk, ajaxData, parsedData, viewer ){
    parsedData[`co-${coPk}`] = ajaxData.taskResult
    toggleCompound(coPk, parsedData, viewer)
    addSaveGemsModal(coPk, ajaxData.saveGemsModal)
    // Add delta energy to the table
    const deltaEnergy = parseFloat(parsedData[`co-${coPk}`]["delta_energy"]).toFixed(3);
    document.querySelector(`.delta-energy#co-${coPk}`).innerText = deltaEnergy
    // Resize the table header to fit the new table content
    $('#torsion-table').dataTable().fnAdjustColumnSizing();
    const loadingModalShowing = document.getElementById("loading-modal").classList.contains("show")
    const showHelpModal = JSON.parse(document.getElementById("show-help-modal").textContent);
    if(loadingModalShowing && showHelpModal){
        $('#torsion-help-modal').modal('show');
    }
}

/**
 * Initializes the chart and defines click behavior
 * @param {Object} parsedData output from the backend torsion scan
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
function setupChart(parsedData, viewer){
    // Initialize chart
    const chart = new Chart(
        document.getElementById('dihedral-energy-chart'),
        {
            type: 'scatter',
            // add initial dataset so that the chart renders
            data: {datasets: [{data: [], label:"initial"}]},
            options: {
                events: ['click'],
                interaction: {
                    mode: 'point'
                },
                plugins:{
                    tooltip: {
                        enabled: false
                    },
                    legend: {
                        onClick: null,
                        labels: {
                            filter: function(legendItem, chartData) {
                                const label = chartData.datasets[legendItem.datasetIndex].label
                                if(label.includes("initial")){
                                    return false
                                }
                                return true;
                            }
                        }
                    },
                    // Annotations allow us to add lines to specify the initial dihedrals
                    annotation: {annotations: {}}
                },
                onClick(evt, points, chart) {
                    for(var i=0; i<points.length; i++){
                        const point = points[i]
                        const dataset = chart.data.datasets[point.datasetIndex]
                        togglePoint(parsedData, viewer, chart, dataset, point.index)   
                    }
                }
            }
        }
    )
    chart.update()
}  

/**
 * Updates the chart to either add or remove the data for `coId` (removes if coId is already
 * in the chart, adds if not)
 * @param {Object} chart a chartJs object
 * @param {Object} parsedData output from the backend torsion scan
 * @param {Object} coId the compound occurrence basechem ID of the CO being added to the chart
 */
function updateChart(chart, parsedData, coId) {
    if(coIdDisplayedInGraph(coId, chart)){
      // Remove dataset
      for(var i=0; i<chart.data.datasets.length; i++){
        if(chart.data.datasets[i].id == coId){
          chart.data.datasets.splice(i, 1);
          break;
        }
      }
      // Remove vertical line annotation for initial dihedral
      delete chart.options.plugins.annotation.annotations[coId] 
    }else{
      // Construct data
      const data = []
      for(var dihedral in parsedData[coId].torsions){
        const torsion = parsedData[coId].torsions[dihedral]
        data.push({
          x: parseInt(torsion.dihedral),
          y: torsion.rel_energy
        })
      }
      // Sort data in order of dihedral so the lines connect the correct points
      data.sort((a, b) => (a.x > b.x) ? 1 : -1)
      // Add dataset
      const carbonColor = pickCarbonColor(colorMap,colorFrequency, coId, true)
      const dataset = {
        label: getChartLabel(coId),
        id: coId,
        data: data,
        borderColor: carbonColor,
        backgroundColor: carbonColor,
        pointBackgroundColor: carbonColor,
        showLine: true,
        pointRadius: Array(data.length).fill(INACTIVE_POINT_SIZE),
        pointHoverRadius: Array(data.length).fill(ACTIVE_POINT_SIZE),
      }
      chart.data.datasets.push(dataset)
      // Add a vertical line annotation for the initial dihedral
      if(parsedData[coId]["initial_dihedral"]){
        chart.options.plugins.annotation.annotations[coId] = {
            type: 'line',
            xMin: parseFloat(parsedData[coId]["initial_dihedral"]),
            xMax: parseFloat(parsedData[coId]["initial_dihedral"]),
            borderColor: carbonColor,
            borderWidth: 2,
          }
      }
    }
    chart.update();
    toggleButtons(chart)
  }

/**
 * Enables/Disables the "Lowest Energy" and "Match Dihedral" buttons depending on the state
 * of what's clicked. The "Lowest Energy" button requires that at least one CO is displayed
 * in the graph. The "Match Dihedral" button requires that at least one CO is displayed
 * in the graph ant that exactly 1 pose is clicked and shown in the viewer
 * @param {Object} chart a chartJs object 
 */
function toggleButtons(chart){
    // Enable/Disable the buttons depending on what's in the viewer
    const numDatasets = chart.data.datasets.length -1 // Subtract the initial, blank dataset
    if(numDatasets == 0){
        $('#lowest-energy-btn')[0].disabled = true
        $('#match-dihedral-btn')[0].disabled = true
    }else if(numDatasets > 0){
        $('#lowest-energy-btn')[0].disabled = false
        $('#match-dihedral-btn')[0].disabled = !(Object.keys(models).length == 1)
    }
}

/***************************************************
 *                 Event Listeners                 *
 ***************************************************/

/**
 * Displays the toggle-compound button for a particular CompoundOccurrence.
 * Toggling a CO that is already displayed removes it from the chart, removes its models
 * from the viewer, and frees up the color that was being used. Toggling a CO that is not
 * yet displayed adds it to the chart and picks a color to use 
 * @param {Number} coPk the pk of the CompoundOccurrence object
 * @param {Object} parsedData output from the backend torsion scan
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
function toggleCompound(coPk, parsedData, viewer){
    const toggleBtn = document.querySelector(`.toggle-compound#co-${coPk}`)
    toggleBtn.parentElement.style.display = ""
    const chart = Chart.getChart('dihedral-energy-chart');
    $(`.toggle-compound#co-${coPk}`).on('click', function (e) {
        const coId = this.id
        const tag = document.querySelector(`button.toggle-compound#${coId}`)
        updateChart(chart, parsedData, this.id)
        if(coIdDisplayedInGraph(this.id, chart)){
            // Add background color to button
            this.innerText="Hide"
            tag.style.backgroundColor = pickCarbonColor(colorMap,colorFrequency, coId, true);
        }else{
            // Remove background color from button
            tag.style.backgroundColor = "";
            colorFrequency[colorMap[coId]] -= 1
            delete colorMap[coId]
            this.innerText = "Show"
            // Remove all models for this CO
            for(var modelId in models){
                if(modelId.includes(coId)){
                    removeModelWrapper(viewer, modelId)  
                }
            }
            viewer.render()
        }
    });
}

/**
 * Adds an event listener to the "Lowest Energy" button. When clicked this button will remove
 * all models from the viewer and display the lowest energy pose for each CO that is currently
 * displayed in the chart. This also updates the chart to reflect which poses are currently
 * being displayed in the viewer.
 * @param {Object} parsedData output from the backend torsion scan
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
function showLowestEnergy(parsedData, viewer){
    const chart = Chart.getChart('dihedral-energy-chart');
    $('#lowest-energy-btn').on('click', function(){
        for(var datasetIndex=0; datasetIndex<chart.data.datasets.length;datasetIndex++){
            const dataset = chart.data.datasets[datasetIndex]
            for(var pointIndex=0; pointIndex< dataset.data.length; pointIndex++){
                const dihedral = dataset.data[pointIndex].x
                const relEnergy = dataset.data[pointIndex].y
                const confId = getTorsionConfId(dataset, dihedral)
                if(relEnergy == 0){
                    // Select this point and add it to the viewer
                    togglePoint(parsedData, viewer, chart, dataset, pointIndex, "select")
                }else if(models[confId]){
                    // Deselect this point and remove it from the viewer
                    togglePoint(parsedData, viewer, chart, dataset, pointIndex, "deselect")
                }
            }
        }
    })
}

/**
 * Adds an event listener to the "Match Dihedral" button. When clicked this button will find
 * the dihedral of the currently displayed model and display the poses of the same
 * dihedral for all other COs that are currently displayed. This also updates the chart
 * to reflect which poses are currently displayed
 * @param {Object} parsedData output from the backend torsion scan
 * @param {Object} viewer the 3dmol viewer that is displayed
 */
function showMatchDihedral(parsedData,viewer){
    const chart = Chart.getChart('dihedral-energy-chart');
    $('#match-dihedral-btn').on('click', function () {
        // Check that there is only one pose currently selected in the viewer
        // (the HTML enforces this, this is just a safety check)
        if(Object.keys(models).length != 1){
            return
        }
        // Get basechem identifier of selected model and get the CO id and dihedral from it
        const modelId = Object.keys(models)[0]
        // matchDihedral pulls the dihedral (ex: -180) out the modelId (ex: co-25-neg180)
        var matchDihedral = modelId.match(/co-\d+-((neg)?\d+)/)[1]
        if(matchDihedral.includes("neg")){
            matchDihedral = `-${matchDihedral.substring(3)}`
        }
        matchDihedral = Number(matchDihedral)

        // For each dataset, select the point with the same dihedral as `matchDihedral`
        for(var datasetIndex=0; datasetIndex<chart.data.datasets.length;datasetIndex++){
            const dataset = chart.data.datasets[datasetIndex]
            for(var pointIndex=0; pointIndex< dataset.data.length; pointIndex++){
                const dihedral = dataset.data[pointIndex].x
                if(dihedral == matchDihedral){
                    togglePoint(parsedData, viewer, chart, dataset, pointIndex, "select")
                }
            }
        }
    })
}



/***************************************************
 *                      Utils                      *
 ***************************************************/


/**
 * Given a compound occurrence ID, return the label for the chart legend
 * @param {String} coId the compound occurrence ID
 */
function getChartLabel(coId){
    return($(`#${coId}-name`)[0].innerText)
}

/**
 * @param {String} coId a compound occurrence basechem id (ex: "co-12")
 * @param {Object} chart a chartJs object
 * @returns a boolean, is this compound occurrence currently displayed in the chart
 */
function coIdDisplayedInGraph(coId, chart){
    for(let i=0; i<chart.data.datasets.length; i++){
        if(chart.data.datasets[i].id == coId){
            return true
        }
    }
    return false
}

/**
 * Handles the action of a user clicking a point on the graph. If the point is already clicked,
 * it deselects it and removes its model from the viewer. If the point is
 * not yet clicked, it selects it and adds its model to the viewer.
 * @param {Object} parsedData output from the backend torsion scan
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {Object} chart a chartJs object
 * @param {Object} dataset the chartJs dataset whose point is being clicked
 * @param {Number} pointIndex the index of the point being clicked
 * @param {String} force "select" if you want this point to be selected (even if it's already
 *  selected); "deselect" if you want this point to be deselected (even if it's already
 *  deselected); null if you want the default behavior (selects if currently deselected, deselects
 *  if currently selected)
 */
function togglePoint(parsedData, viewer, chart, dataset, pointIndex, force=null){
    const dihedral = dataset.data[pointIndex].x
    var confId = getTorsionConfId(dataset, dihedral)
    const select = (dataset.pointRadius[pointIndex] == INACTIVE_POINT_SIZE && force != "deselect")|| force == "select"
    if(select){
        // Select point
        dataset.pointRadius[pointIndex] = ACTIVE_POINT_SIZE
        dataset.pointHoverRadius[pointIndex] = INACTIVE_POINT_SIZE
        // Add to viewer
        if(!models[confId]){
            const moltext = parsedData[dataset.id].torsions[confId].moltext
            const m = addModelWrapper(viewer, moltext, confId)
            viewer.zoomTo({model: m})
        }
    }else{
        // Deselect point
        dataset.pointRadius[pointIndex] = INACTIVE_POINT_SIZE
        dataset.pointHoverRadius[pointIndex] = ACTIVE_POINT_SIZE
        // Remove from viewer
        if(models[confId]){
            removeModelWrapper(viewer, confId)
        }
    }
    viewer.render()
    toggleButtons(chart)
    chart.update()
    // After clicking, this point, click somewhere in the background of the chart so that
    // the point sizes update
    let canvas = document.querySelector("#dihedral-energy-chart");
    let rect = canvas.getBoundingClientRect();
    canvas.dispatchEvent(new MouseEvent("click",{
        clientX: rect.left + 50,
        clientY: rect.top + 50
    }))
}

/**
 * A wrapper around `removeModel` that passes the `colorFrequency`, `colorMap`,
 * and `models` objects in addition to the `viewer` and `id`
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {String} id the basechem identifier of this model (ex s-3 for series or 121-0 for a pose)
 */
function removeModelWrapper(viewer, id){
    removeModel(models, colorMap, colorFrequency, viewer, id, false)
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
function addModelWrapper(viewer, moltext, id, {format="sdf", opacity=PRIMARY_OPACITY, style=STYLE} = {}){
    const m = addModel(models, colorMap, colorFrequency, viewer, moltext, id, {format:format, opacity:opacity, style:style})
    m.addAtomSpecs(["selected"])
    return m
}

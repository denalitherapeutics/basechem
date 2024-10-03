import { collectAnalysisTask, manageNewAnalysisForm } from './tasks.js'
import { handleNewCompsAdded, addCompOnClick } from "./addComps.js"

/**
 * This function handles all of the front-end interactions for propcalc
 */
export default function () {
  // Toggle between grid and list view
  addToggleViewEventListeners()
  // Disable the new analysis form while tasks are running
  manageNewAnalysisForm()
  // Wait for each propcalc task to finish and add the results to the DOM when ready
  collectPropCalcTasks()
}

/**
 * Adds event listeners to the show-table-btn and show-grid-btn
 */
function addToggleViewEventListeners() {
  function toggleView() {
    const showTableBtn = document.getElementById("show-table-btn")
    const showGridBtn = document.getElementById("show-grid-btn")
    const tableDisplay = showGridBtn.style.display == "" ? 'none' : ''
    const gridDisplay = showTableBtn.style.display == "" ? 'none' : ''
    document.getElementById("propcalc-grid-container").style.display = gridDisplay
    showTableBtn.style.display = gridDisplay
    document.getElementById("grid-view-helptext").style.display = gridDisplay
    document.getElementById("propcalc-table-container").style.display = tableDisplay
    showGridBtn.style.display = tableDisplay
  }
  $("#show-grid-btn").on("click", toggleView)
  $('#show-table-btn').on("click", toggleView)
}

/**
 * For each `task-status` in the table, wait for the associated DjangoQ task to complete
 * and add the results to the table/grid
 */
function collectPropCalcTasks() {
  document.querySelectorAll(".task-status").forEach(element => {
    const taskName = element.getAttribute("data-task-name")
    collectAnalysisTask(taskName, {}, null, collectSuccessfulPropCalcTask)
  })
  handleNewCompsAdded("propcalc-table", {}, null, newCompAddedHelper)
}

/**
 * Called in `collectAnalysisTask` when the task returns successfully - this function uses
 * the result from the backend to update parsedData, the viewer, and the DOM. This function additionally
 * accepts `parsedData` and `viewer` as arguments because they are relevant to all other basechem analyses,
 * but they are not used here.
 * @param {Number} coPk the PK of the CompoundOccurrence whose task returned successfully
 * @param {Object} ajaxData the data returned from the backend, including task results and HTML strings to add to the DOM
 */
function collectSuccessfulPropCalcTask(coPk, ajaxData, _parsedData, _viewer) {
  // Replace old row
  const dataTable = $('#propcalc-table').DataTable()
  const oldRow = dataTable.row(`#co-${coPk}`)
  const orderNum = oldRow.data()[0]
  oldRow.remove()
  const newRow = dataTable.row.add($(ajaxData.tableRow))
  dataTable.cell(newRow.index(), 0).data(orderNum)
  dataTable.draw()
  $(document.getElementById(`add-comp-btn-${coPk}`)).on("click", addCompOnClick)
  // Replace old grid item
  const grid = document.getElementById("propcalc-grid-container")
  const oldGridItem = grid.querySelector(`.propcalc-grid-item#co-${coPk}`)
  oldGridItem.insertAdjacentHTML("afterend", ajaxData.gridItem)
  oldGridItem.remove()
  $('.hover-tooltip').tooltip({trigger:"hover"})
  // Update tooltip text with the latest training date (if inductive props exist)
  if(ajaxData.taskResult["latest_lm_data_date"]){
    const hlmTooltipElement = document.getElementById("hlm-tooltip")
    if(hlmTooltipElement){
      const hlmTooltip = hlmTooltipElement.getAttribute("data-original-title").replace("?",ajaxData.taskResult.latest_lm_data_date)
      hlmTooltipElement.setAttribute("data-original-title", hlmTooltip)
    }
    const rlmTooltipElement = document.getElementById("rlm-tooltip")
    if(rlmTooltipElement){
      const rlmTooltip = rlmTooltipElement.getAttribute("data-original-title").replace("?",ajaxData.taskResult.latest_lm_data_date)
      rlmTooltipElement.setAttribute("data-original-title", rlmTooltip)
    }  
  }
  if(ajaxData.taskResult["latest_logd_data_date"]){
    const logdTooltipElement = document.getElementById("logd-tooltip")
    if(logdTooltipElement){
      const logdTooltip = logdTooltipElement.getAttribute("data-original-title").replace("?",ajaxData.taskResult.latest_logd_data_date)
      logdTooltipElement.setAttribute("data-original-title", logdTooltip)
    }
  }
}

/**
 * Used in handleNewCompsAdded, this function collects task results for new compounds.
 * This function additionally accepts `parsedData` and `viewer` as arguments because they
 * are relevant to all other basechem analyses, but they are not used here.
 * @param {Object} trElement a table row element for a newly added compound
 */
function newCompAddedHelper(trElement, _parsedData, _viewer) {
  const taskName = trElement.querySelector(".task-status").getAttribute("data-task-name")
    collectAnalysisTask(taskName, {}, null, collectSuccessfulPropCalcTask)
}
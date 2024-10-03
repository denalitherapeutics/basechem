import { collectAnalysisTask, manageNewAnalysisForm } from './tasks.js'

/**
 * This function handles all of the front-end interactions for MMP
 */
export default function () {
  // Disable the new analysis form while tasks are running
  manageNewAnalysisForm()
  // Wait for each propcalc task to finish and add the results to the DOM when ready
  collectMMPTasks()
}


/**
 * For each `task-status` in the table, wait for the associated DjangoQ task to complete
 * and add the results to the grid
 */
function collectMMPTasks() {
  document.querySelectorAll(".task-status").forEach(element => {
    const taskName = element.getAttribute("data-task-name")
    collectAnalysisTask(taskName, {}, null, collectSuccessfulMMPTask)
  })
}

/**
 * Called in `collectAnalysisTask` when the task returns successfully - this function uses
 * the result from the backend to update parsedData, the viewer, and the DOM. This function additionally
 * accepts `parsedData` and `viewer` as arguments because they are relevant to all other basechem analyses,
 * but they are not used here.
 * @param {Number} _coPk the PK of the CompoundOccurrence whose task returned successfully
 * @param {Object} ajaxData the data returned from the backend, including task results and HTML strings to add to the DOM
 */
function collectSuccessfulMMPTask(_coPk, ajaxData, _parsedData, _viewer) {
  // Replace old row
  const dnGrid = document.getElementById("dns")
  dnGrid.insertAdjacentHTML("beforeEnd", ajaxData.dnGrid)
  const ideaGrid = document.getElementById("ideas")
  ideaGrid.insertAdjacentHTML("beforeEnd", ajaxData.ideaGrid)
  $('.hover-tooltip').tooltip({trigger:"hover"})
}
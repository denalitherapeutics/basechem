import { hideLoadingModal } from "../../../mni_common/static/mni_common/js/loading.js"
import { postAjaxForm } from "../../../mni_common/static/mni_common/js/ajaxForms.js"
import { toggleDisableChildren } from "../../../mni_common/static/mni_common/js/toggle.js";

/**
 * Checks to see if an async task is complete. If it is, refreshes the page
 * @param {Number} millisecs The amount of time to wait between ajax calls
 * @param {String} identifier The group name or task ID of the async task(s)
 * @param {String} task_type The type of task being checked, "group" or "task"
 */
export function checkTask(millisecs, identifier, taskType) {
  if (identifier) {
    window.setInterval(function () {
      $.ajax({
        url: '/ajax/check-group/',
        data: { "identifier": identifier, "task_type": taskType },
        success: function (data) {
          if ("result" in data) {
            location.reload()
          }
        }
      })
    }, millisecs);
  };
};

/**
 * This is a general-purpose function for ALIGN, DOCK, ESP, and TORSION. The provided
 * taskName must have a data-task-name attribute. This function waits for the DjangoQ
 * task specified by data-task-name to complete and updates the DOM based on the results.
 * @param {Object} taskName a select or button HTML element that should be enabled upon
 *  receiving task results. Must have a `data-task-name` attribute with a DjangoQ task name.
 * @param {Object} parsedData an object with all results for the given analysis. The results from
 *  this task are added to parsedData.
 * @param {Object} viewer the displayed 3dmol viewer
 * @param {Object} successMethod an analysis-specific helper function that should accept coPk, parsedData, and viewer
 *  as arguments. This function should do more specific actions such as add event listeners and display results.
 */
export function collectAnalysisTask(taskName, parsedData, viewer, successMethod) {
  const coPk = taskName.split("_")[1]
  const collectionPk = JSON.parse(document.getElementById('collection-pk').textContent);
  const groupName = JSON.parse(document.getElementById('group-name').textContent);
  var startNextInterval = true
  const interval = window.setInterval(function () {
    // If the previous ajax call is still running, don't start a new one
    if (!startNextInterval) { return }
    startNextInterval = false
    $.ajax({
      url: '/ajax/collect-task/',
      data: {
        "collection_pk": collectionPk,
        "co_pk": coPk,
        "group_name": groupName,
        "task_name": taskName,
      },
      success: function (data) {
        if (!(data.taskStatus == "in progress")) {
          // Stop making ajax calls if the task is no longer inprogress
          clearInterval(interval)
        }
        if (data.taskStatus == "complete") {
          document.querySelector(`.task-status#co-${coPk}`).style.display = "none";
          successMethod(coPk, data, parsedData, viewer);
          hideLoadingModal()
        } else if (data.taskStatus == "dropped") {
          const status = "Dropped"
          document.querySelector(`.task-status#co-${coPk}`).innerText = status
          document.querySelectorAll(`.task-status-update#co-${coPk}`).forEach(element => {
            element.innerText = status
          })
        } else if (data.taskStatus == "error") {
          const status = "Error"
          document.querySelector(`.task-status#co-${coPk}`).innerText = status
          document.querySelectorAll(`.task-status-update#co-${coPk}`).forEach(element => {
            element.innerText = status
          })
        }
        startNextInterval = true
      }
    })
  }, 1000)
}

/**
* This is a general-purpose function for ALIGN, DOCK, and ESP. It waits for the references
* (ex:"alignrefs") DjangoQ task to finish and updates the DOM based on the results.
* @param {Object} parsedData an object with all results for the given analysis. The results from
*  this task are added to parsedData.
* @param {Object} viewer the displayed 3dmol viewer
* @param {String} analysis the analysis being run (ex: "align")
* @param {Object} collectRefsFunc an analysis-specific helper function that should accept taskData, parsedData, and viewer
*  as arguments. This function should do more specific actions such as add event listeners and display results.
*/
export function collectRefsTask(parsedData, viewer, analysis, collectRefsFunc) {
  const groupName = JSON.parse(document.getElementById('group-name').textContent);
  const collectionPk = groupName.split("_")[1]
  var startNextInterval = true
  const interval = window.setInterval(function () {
    if (!startNextInterval) { return }
    startNextInterval = false
    $.ajax({
      url: '/ajax/collect-task/',
      data: {
        "task_name": `${analysis}refs_${collectionPk}`,
        "group_name": groupName,
        "collection_pk": collectionPk,
      },
      success: async function (taskData) {
        if ("taskResult" in taskData) {
          clearInterval(interval)
          await collectRefsFunc(taskData, parsedData, viewer)
        }
        startNextInterval = true;
      }
    })
  }, 1000)
}

/**
 * Given a string of HTML for a save-gems-modal, add this modal to the DOM and enable the relevant
 * toggle-save-gems-modal button.
 * @param {Number} coPk the PK of the CompoundOccurrence whose save-gems-modal is being added
 * @param {String} modalHtml a string of HTML encoding a modal
 */
export function addSaveGemsModal(coPk, modalHtml) {
  // Append the saveGemsModal (modal opened by backpack button) to the DOM
  document.getElementById(`save-gems-modals`).insertAdjacentHTML('beforeend', modalHtml);
  $(`#save-gems-form-${coPk}`).on('submit', function (event) {
    event.preventDefault();
    postAjaxForm.call(this);
  });
  // Enable backpack button
  const backpackButton = document.getElementById(`toggle-save-gems-modal-${coPk}`)
  backpackButton.disabled = false
}


/**
* This function disables the `new-analysis-form` and enables it once all running tasks have returned.
*/
export function manageNewAnalysisForm() {
  // Disable new analysis form
  const submitForm = document.getElementById('new-analysis-form')
  if(submitForm){
    toggleDisableChildren(submitForm, true)
  }
  // When all tasks have returned, enable new analysis form
  const interval = window.setInterval(function () {
    var allTasksReturned = true;
    document.querySelectorAll(".task-status").forEach(element => {
      // If a task status is still showing "Job running", the not all tasks have returned
      if ((element.style.display != "none") && (element.innerText == "Job running")) {
        allTasksReturned = false;
      }
    })
    if (allTasksReturned) {
      clearInterval(interval)
      if(submitForm){
        toggleDisableChildren(submitForm, false)
      }
      // If the loading modal is still showing (because all tasks failed), hide it
      hideLoadingModal()
    }
  }, 5000)
}